import bioframe as bf
import cooler
import cooltools
import gsd.hoomd
import matplotlib as mpl
import matplotlib.pyplot as plt
import numpy as np
import pandas as pd
from matplotlib.ticker import EngFormatter

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


def set_snapshot_box(cfg: dict) -> gsd.hoomd.Frame:
    assert (cfg["end"] - cfg["start"]) % cfg["bin"] == 0, (
        "bin does not divide locus size"
    )
    # Initialize simulation with the appropriate box size
    number_of_monomers = (cfg["end"] - cfg["start"]) // cfg["bin"]
    number_of_monomers *= cfg["n_copies"]
    box_length = (number_of_monomers / cfg["density"]) ** (1 / 3.0)
    snapshot = gsd.hoomd.Frame()
    box = [box_length] * 3 + [0] * 3
    snapshot.configuration.box = np.asarray(box, dtype=np.float32)

    return snapshot


def populate_particles_to_snapshot(
    cfg: dict, snapshot: gsd.hoomd.Frame
) -> gsd.hoomd.Frame:
    # Get monomer types of the locus as A(0)/B(1) compartment.
    snapshot.particles.typeid, snapshot.particles.types = pd.factorize(
        pd.read_csv(cfg["data_dir"] / "output" / "AB.csv", header=0)["AB"], sort=True
    )
    # Replicate the locus to optimize statistical sampling and mitigate potential finite-size effects.
    snapshot.particles.typeid = np.tile(snapshot.particles.typeid, cfg["n_copies"])
    snapshot.particles.N = len(snapshot.particles.typeid)
    # Check whether A -> 0 and B -> 1.
    assert snapshot.particles.types[0] == "A", "A to 0 and B to 1"
    snapshot.particles.types = snapshot.particles.types.to_list()

    # Build random, dense initial conformations.
    box_length = snapshot.configuration.box.max().item()
    snapshot.particles.position = grow_cubic(
        N=snapshot.particles.N, boxSize=int(box_length - 1)
    ).astype(np.float32)
    # Centralize
    snapshot.particles.position -= snapshot.particles.position.mean(
        axis=0, keepdims=True
    )

    return snapshot


def populate_bonds_to_snapshot(cfg: dict, snapshot: gsd.hoomd.Frame) -> gsd.hoomd.Frame:
    snapshot.bonds.types = list(cfg["force"]["Bonded forces"].keys())

    one_rep_size = snapshot.particles.N // cfg["n_copies"]
    start_ids = np.add.outer(
        np.arange(0, snapshot.particles.N, one_rep_size), np.arange(one_rep_size - 1)
    ).flatten()
    snapshot.bonds.group = np.add.outer(start_ids, np.arange(2))

    return snapshot


def populate_angles_to_snapshot(
    cfg: dict, snapshot: gsd.hoomd.Frame
) -> gsd.hoomd.Frame:
    snapshot.angles.types = list(cfg["force"]["Angular forces"].keys())

    one_rep_size = snapshot.particles.N // cfg["n_copies"]
    start_ids = np.add.outer(
        np.arange(0, snapshot.particles.N, one_rep_size), np.arange(one_rep_size - 2)
    ).flatten()
    snapshot.angles.group = np.add.outer(start_ids, np.arange(3))

    return snapshot
