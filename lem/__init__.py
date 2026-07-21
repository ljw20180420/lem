import bioframe as bf
import cooler
import cooltools
import gsd.hoomd
import hoomd
import matplotlib as mpl
import matplotlib.pyplot as plt
import numpy as np
import pandas as pd
from matplotlib.ticker import EngFormatter

import polychrom_hoomd.build as build
import polychrom_hoomd.log as log
from polychrom_hoomd import forces, render
from polykit.generators.initial_conformations import grow_cubic

mpl.use("agg")


def get_E1(cfg: dict) -> None:
    clr = cooler.Cooler(cfg["hic"])
    bins = clr.bins().fetch((cfg["chrom"], cfg["start"], cfg["end"]))

    genome = bf.load_fasta(cfg["genome"])
    gc_cov = bf.frac_gc(bins[["chrom", "start", "end"]], genome)

    cis_eigs = cooltools.eigs_cis(
        clr,
        gc_cov,
        view_df=pd.DataFrame({
            "chrom": [cfg["chrom"]],
            "start": [cfg["start"]],
            "end": [cfg["end"]],
            "name": [cfg["chrom"]],
        }),
        n_eigs=1,
    )
    df_E1 = cis_eigs[1][["chrom", "start", "end", "E1"]].assign(
        E1=lambda df: df["E1"].interpolate()
    )

    (cfg["data_dir"] / "output").mkdir(exist_ok=True, parents=True)
    df_E1.to_csv(cfg["data_dir"] / "output" / "E1.csv", index=False)

    bp_formatter = EngFormatter("b")
    ax = df_E1.plot(x="start", y="E1")
    ax.xaxis.set_major_formatter(bp_formatter)
    ax.set_xticks(
        np.arange(cfg["start"], cfg["end"], 100000),
    )
    ax.tick_params(axis="x", labelrotation=45)
    fig = ax.get_figure()
    fig.tight_layout()
    fig.savefig(cfg["data_dir"] / "output" / "E1.pdf")
    plt.close(fig)


def align_E1(cfg: dict) -> None:
    df_E1 = pd.read_csv(cfg["data_dir"] / "output" / "E1.csv", header=0).assign(
        AB=lambda df: pd.cut(
            df["E1"], bins=[-float("inf"), 0, float("inf")], labels=["B", "A"]
        )
    )

    df_bin = bf.binnify(
        chromsizes=pd.Series(
            data=[cfg["end"] - cfg["start"]],
            index=[cfg["chrom"]],
        ),
        binsize=cfg["bin"],
    ).assign(
        start=lambda df: df["start"] + cfg["start"],
        end=lambda df: df["end"] + cfg["start"],
    )

    df_AB = bf.closest(df_bin, df_E1)[["chrom", "start", "end", "AB_"]].rename(
        columns={"AB_": "AB"}
    )
    df_AB.to_csv(cfg["data_dir"] / "output" / "AB.csv", index=False)


class Frame:
    def __init__(self, cfg: dict) -> None:
        self.cfg = cfg

    def __call__(self) -> gsd.hoomd.Frame:
        frame = self.set_box()
        frame = self.populate_particles(frame)
        frame = self.populate_bonds(frame)
        frame = self.populate_angles(frame)

        return frame

    def set_box(self) -> gsd.hoomd.Frame:
        assert (self.cfg["end"] - self.cfg["start"]) % self.cfg["bin"] == 0, (
            "bin does not divide locus size"
        )
        # Initialize simulation with the appropriate box size
        number_of_monomers = (self.cfg["end"] - self.cfg["start"]) // self.cfg["bin"]
        number_of_monomers *= self.cfg["n_copies"]
        box_length = (number_of_monomers / self.cfg["density"]) ** (1 / 3.0)
        frame = gsd.hoomd.Frame()
        box = [box_length] * 3 + [0] * 3
        frame.configuration.box = np.asarray(box, dtype=np.float32)

        return frame

    def populate_particles(self, frame: gsd.hoomd.Frame) -> gsd.hoomd.Frame:
        # Get monomer types of the locus as A(0)/B(1) compartment.
        frame.particles.typeid, frame.particles.types = pd.factorize(
            pd.read_csv(self.cfg["data_dir"] / "output" / "AB.csv", header=0)["AB"],
            sort=True,
        )
        # Replicate the locus to optimize statistical sampling and mitigate potential finite-size effects.
        frame.particles.typeid = np.tile(frame.particles.typeid, self.cfg["n_copies"])
        frame.particles.N = len(frame.particles.typeid)
        # Check whether A -> 0 and B -> 1.
        assert frame.particles.types[0] == "A", "A to 0 and B to 1"
        frame.particles.types = frame.particles.types.to_list()
        frame.particles.diameter = np.ones(frame.particles.N, dtype=np.float32)

        # Build random, dense initial conformations.
        box_length = frame.configuration.box.max().item()
        frame.particles.position = grow_cubic(
            N=frame.particles.N, boxSize=int(box_length - 1)
        ).astype(np.float32)
        # Centralize
        frame.particles.position -= frame.particles.position.mean(axis=0, keepdims=True)

        return frame

    def populate_bonds(self, frame: gsd.hoomd.Frame) -> gsd.hoomd.Frame:
        frame.bonds.types = list(self.cfg["force"]["Bonded forces"].keys())

        one_rep_size = frame.particles.N // self.cfg["n_copies"]
        start_ids = np.add.outer(
            np.arange(0, frame.particles.N, one_rep_size),
            np.arange(one_rep_size - 1),
        ).flatten()
        frame.bonds.N = len(start_ids)
        frame.bonds.group = np.add.outer(start_ids, np.arange(2))

        typeid = frame.bonds.types.index("Backbone")
        frame.bonds.typeid = np.full(frame.bonds.N, typeid, dtype=np.uint32)

        return frame

    def populate_angles(self, frame: gsd.hoomd.Frame) -> gsd.hoomd.Frame:
        frame.angles.types = list(self.cfg["force"]["Angular forces"].keys())

        one_rep_size = frame.particles.N // self.cfg["n_copies"]
        start_ids = np.add.outer(
            np.arange(0, frame.particles.N, one_rep_size),
            np.arange(one_rep_size - 2),
        ).flatten()
        frame.angles.N = len(start_ids)
        frame.angles.group = np.add.outer(start_ids, np.arange(3))

        typeid = frame.angles.types.index("Curvature")
        frame.angles.typeid = np.full(frame.angles.N, typeid, dtype=np.uint32)

        return frame


def init_trajectory(cfg: dict) -> None:
    frame = Frame(cfg)()
    with gsd.hoomd.open(
        name=cfg["data_dir"] / "output" / "trajectory.gsd", mode="w"
    ) as fd:
        fd.append(frame)


def append_trajectory(cfg: dict, frame: gsd.hoomd.Frame) -> None:
    with gsd.hoomd.open(
        name=cfg["data_dir"] / "output" / "trajectory.gsd", mode="a"
    ) as fd:
        fd.append(frame)


def draw_trajectory(cfg: dict) -> None:
    # Visualize starting conformation using the Fresnel backend (A compartments in blue, B in red)
    (cfg["data_dir"] / "output" / "frames").mkdir(exist_ok=True, parents=True)
    with gsd.hoomd.open(
        name=cfg["data_dir"] / "output" / "trajectory.gsd", mode="r"
    ) as fd:
        for i, frame in enumerate(fd):
            render.fresnel(frame, show="compartments", cmap="coolwarm").static(
                pathtrace=False,
                png_output_file=cfg["data_dir"]
                / "output"
                / "frames"
                / f"frame-{i}.png",
            )


class Integrator:
    def __init__(self, cfg: dict) -> None:
        self.cfg

    def __call__(self) -> hoomd.md.Integrator:
        # Setup neighbor list
        nl = hoomd.md.nlist.Cell(buffer=0.4)

        # Set chromosome excluded volume
        repulsion_forces = forces.get_repulsion_forces(nl, **self.cfg["force"])

        # Set bonded/angular potentials
        bonded_forces = forces.get_bonded_forces(**self.cfg["force"])
        angular_forces = forces.get_angular_forces(**self.cfg["force"])

        # Set attractive/DPD forces
        dpd_forces = forces.get_dpd_forces(nl, **self.cfg["force"])
        attraction_forces = forces.get_attraction_forces(nl, **self.cfg["force"])

        # Define full force_field
        dpd_force_field = (
            repulsion_forces
            + bonded_forces
            + angular_forces
            + attraction_forces
            + dpd_forces
        )

        # Setup integrator methods
        nve = hoomd.md.methods.ConstantVolume(filter=hoomd.filter.All())
        dpd_integrator = hoomd.md.Integrator(
            dt=5e-3, methods=[nve], forces=dpd_force_field
        )

        return dpd_integrator


def simulate(cfg: dict, steps: int) -> None:
    assert steps % cfg["period"] == 0, "period does not divide steps"
    if cfg["device"] == "gpu":
        device = hoomd.device.GPU(notice_level=3)
    else:
        assert cfg["device"] == "cpu", "only support devices cpu and gpu"
        device = hoomd.device.CPU(notice_level=3)
    simulation = hoomd.Simulation(device=device, seed=cfg["seed"])
    trajectory = (cfg["data_dir"] / "output" / "trajectory.gsd").as_posix()
    simulation.create_state_from_gsd(filename=trajectory)
    dpd_integrator = Integrator(cfg)()
    simulation.operations.integrator = dpd_integrator
    simulation.operations.writers.append(
        hoomd.write.GSD(
            trigger=hoomd.trigger.Periodic(period=cfg["period"]),
            filename=trajectory,
            logger=log.get_logger(simulation),
        )
    )
    simulation.run(steps)
