import bioframe as bf
import cooler
import cooltools
import cupy as cp
import gsd.hoomd
import hoomd
import matplotlib as mpl
import matplotlib.pyplot as plt
import numpy as np
import pandas as pd
from matplotlib.ticker import EngFormatter

import polychrom_hoomd.forces
import polychrom_hoomd.log
import polychrom_hoomd.render
from fundenberg.onestate_extruder import LEFTranslocatorDirectional, compute_LEF_pos
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
        LEF_num = int(
            frame.particles.N // (self.cfg["1d"]["LEF_separation"] / self.cfg["bin"])
        )
        frame.bonds.N = len(start_ids) + LEF_num
        lef1 = np.random.randint(low=0, high=frame.particles.N - 1, size=LEF_num)
        lef2 = lef1 + 1
        frame.bonds.group = np.concatenate(
            (
                np.add.outer(start_ids, np.arange(2)),
                np.stack((lef1, lef2), axis=1),
            ),
            axis=0,
        )

        backbone_typeid = frame.bonds.types.index("Backbone")
        LEF_dummy_typeid = frame.bonds.types.index("LEF_dummy")
        frame.bonds.typeid = np.concatenate(
            (
                np.full(frame.bonds.N - LEF_num, backbone_typeid, dtype=np.uint32),
                np.full(LEF_num, LEF_dummy_typeid, dtype=np.uint32),
            ),
            axis=0,
        )

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
    if cfg["device"] == "gpu":
        device = hoomd.device.GPU(notice_level=3)
    else:
        assert cfg["device"] == "cpu", "only support devices cpu and gpu"
        device = hoomd.device.CPU(notice_level=3)
    simulation = hoomd.Simulation(device=device, seed=cfg["seed"])
    simulation.create_state_from_snapshot(snapshot=frame)
    hoomd.write.GSD.write(
        state=simulation.state,
        filename=cfg["data_dir"] / "output" / "trajectory.gsd",
        mode="wb",  # Use 'wb' to write/overwrite a clean file
    )


class Integrator:
    def __init__(self, cfg: dict) -> None:
        self.cfg = cfg

    def __call__(self) -> hoomd.md.Integrator:
        # Setup neighbor list
        nl = hoomd.md.nlist.Cell(buffer=0.4)

        # Set chromosome excluded volume
        repulsion_forces = polychrom_hoomd.forces.get_repulsion_forces(
            nl, **self.cfg["force"]
        )

        # Set bonded/angular potentials
        bonded_forces = polychrom_hoomd.forces.get_bonded_forces(**self.cfg["force"])
        angular_forces = polychrom_hoomd.forces.get_angular_forces(**self.cfg["force"])

        # Set attractive/DPD forces
        dpd_forces = polychrom_hoomd.forces.get_dpd_forces(nl, **self.cfg["force"])
        attraction_forces = polychrom_hoomd.forces.get_attraction_forces(
            nl, **self.cfg["force"]
        )

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
            dt=self.cfg["tau_3d"], methods=[nve], forces=dpd_force_field
        )

        return dpd_integrator


def get_logger(simulation: hoomd.Simulation) -> hoomd.logging.Logger:
    logger = hoomd.logging.Logger(categories=["scalar", "string"])
    logger.add(simulation, quantities=["timestep", "tps"])
    status = polychrom_hoomd.log.Status(simulation)
    logger[("Status", "etr")] = (status, "etr", "string")
    thermo = hoomd.md.compute.ThermodynamicQuantities(filter=hoomd.filter.All())
    simulation.operations.computes.append(thermo)
    logger.add(thermo, quantities=["kinetic_temperature"])

    return logger


def load_simulation(cfg: dict) -> hoomd.Simulation:
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
            logger=get_logger(simulation),
        )
    )

    return simulation


def warmup(cfg: dict) -> None:
    assert cfg["warmup"] % cfg["period"] == 0, "period does not divide warmup"
    simulation = load_simulation(cfg)
    simulation.run(steps=cfg["warmup"])


def loop_extrusion_1d(
    cfg: dict, simulation: hoomd.Simulation, steps: int
) -> np.ndarray:
    trajectory_1d = compute_LEF_pos(
        extrusion_engine=LEFTranslocatorDirectional,
        n_tot=simulation.state.N_particles,
        trajectory_length=steps,
        dummy_steps=cfg["1d"]["warmup"],
        LEF_lifetime=cfg["1d"]["LEF_lifetime"] / cfg["1d"]["tau_1d"],
        LEF_separation=cfg["1d"]["LEF_separation"],
        kb_per_site=cfg["bin"],
    )

    return trajectory_1d


def update_topology(
    simulation, bond_list, xp: cp.__class__ | np.__class__, thermalize=False
):
    """Update topology on either GPU or CPU, based on availability"""

    LEF_typeid = simulation.state.bond_types.index("LEF")
    LEF_dummy_typeid = simulation.state.bond_types.index("LEF_dummy")

    if len(bond_list) > 0:
        # Discard contiguous loops
        bond_array = xp.array(bond_list, dtype=xp.int32)
        type_array = xp.full(len(bond_array), LEF_typeid, dtype=xp.int32)

        redundant_bonds = xp.less(bond_array[:, 1] - bond_array[:, 0], 1)
        n_prune = int(xp.count_nonzero(redundant_bonds))

        ids = xp.random.randint(
            low=0, high=simulation.state.N_particles - 1, size=n_prune, dtype=xp.int32
        )

        bond_array[redundant_bonds] = xp.stack((ids, ids + 1), axis=1)
        type_array[redundant_bonds] = LEF_dummy_typeid

    else:
        bond_array = xp.empty(0, dtype=xp.int32)
        type_array = xp.empty(0, dtype=xp.int32)

    _update_topology_local(
        simulation, bond_array, type_array, LEF_typeid, LEF_dummy_typeid, xp
    )

    if thermalize:
        simulation.state.thermalize_particle_momenta(filter=hoomd.filter.All(), kT=1.0)


def _update_topology_local(system, bond_array, type_array, type_id, dummy_id, xp):
    """Update topology locally on the GPU"""
    device = "gpu" if xp.__name__ == "cupy" else "cpu"
    with getattr(system.state, f"{device}_local_snapshot") as local_snap:
        bond_ids = xp.asarray(local_snap.bonds.typeid)

        is_bound = xp.equal(bond_ids, type_id)
        is_unbound = xp.equal(bond_ids, dummy_id)

        is_LEF = xp.logical_or(is_bound, is_unbound)

        if bond_array.shape[0] == type_array.shape[0] == xp.count_nonzero(is_LEF):
            local_snap.bonds.group[is_LEF] = bond_array.astype(xp.uint32)
            local_snap.bonds.typeid[is_LEF] = type_array.astype(xp.uint32)

        else:
            raise RuntimeError("Unable to dynamically resize bond arrays on the GPU")


def loop_extrusion_3d(cfg: dict, steps: int) -> None:
    n_3d_to_1d = cfg["1d"]["tau_1d"] / cfg["tau_3d"]
    assert n_3d_to_1d % cfg["period"] == 0, "period does not divide n_3d_to_1d"
    simulation = load_simulation(cfg)
    trajectory_1d = loop_extrusion_1d(cfg, simulation, steps)
    xp = np if cfg["device"] == "cpu" else cp
    for bond_list in trajectory_1d:
        update_topology(simulation, bond_list, xp, thermalize=False)
        simulation.run(n_3d_to_1d)


def draw_trajectory(cfg: dict) -> None:
    # Visualize starting conformation using the Fresnel backend (A compartments in blue, B in red)
    (cfg["data_dir"] / "output" / "frames").mkdir(exist_ok=True, parents=True)
    with gsd.hoomd.open(
        name=cfg["data_dir"] / "output" / "trajectory.gsd", mode="r"
    ) as fd:
        for i, frame in enumerate(fd):
            polychrom_hoomd.render.fresnel(
                frame, show="compartments", cmap="coolwarm"
            ).static(
                pathtrace=False,
                png_output_file=cfg["data_dir"]
                / "output"
                / "frames"
                / f"frame-{i}.png",
            )
