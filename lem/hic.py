import os

import bioframe as bf
import cooler
import cooltools
import matplotlib.pyplot as plt
import pandas as pd
from Bio import Seq
from matplotlib.axes import Axes
from matplotlib.colors import LinearSegmentedColormap
from matplotlib.figure import Figure
from matplotlib.ticker import EngFormatter
from mpl_toolkits.axes_grid1 import make_axes_locatable
from numpy.typing import ArrayLike


def _white2red() -> LinearSegmentedColormap:
    return LinearSegmentedColormap(
        name="white2red",
        segmentdata={
            "red": [(0.0, 1.0, 1.0), (1.0, 1.0, 1.0)],
            "green": [(0.0, 1.0, 1.0), (1.0, 0.0, 0.0)],
            "blue": [(0.0, 1.0, 1.0), (1.0, 0.0, 0.0)],
        },
    )


def _heatmap(mat: ArrayLike, chrom: str, start: int, end: int) -> tuple[Figure, Axes]:
    fig, ax = plt.subplots(figsize=(7, 6))
    im = ax.matshow(mat, extent=[start, end, start, end], vmin=0, cmap=_white2red())
    fig.colorbar(im, fraction=0.046, pad=0.04)
    ax.xaxis.set_major_formatter(EngFormatter("b"))
    ax.yaxis.set_major_formatter(EngFormatter("b"))
    ax.xaxis.set_label_position("top")
    ax.xaxis.set_tick_params(rotation=45)
    ax.set_xlabel(chrom)
    ax.set_ylabel(chrom)
    fig.tight_layout()

    return fig, ax


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


def heatmap(cfg: dict) -> None:
    clr = cooler.Cooler(cfg["hic"])
    mat = clr.matrix(balance=False).fetch((cfg["chrom"], cfg["start"], cfg["end"]))
    fig, ax = _heatmap(mat, cfg["chrom"], cfg["start"], cfg["end"])

    divider = make_axes_locatable(ax)
    top_ax = divider.append_axes("top", size="10%", pad=0.1, sharex=ax)
    df_E1 = pd.read_csv(cfg["data_dir"] / "output" / "E1.csv", header=0)
    df_E1.melt(id_vars=["chrom", "E1"])
    breakpoint()
    # Add AB compartment
    top_ax.plot()

    fig.savefig(cfg["data_dir"] / "output" / "heatmap.png")
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


def bed2seq(
    bed: os.PathLike,
    genome: os.PathLike = "/home/ljw/.local/share/genomes/GRCm38/GRCm38.fa",
) -> None:
    df = pd.read_csv(bed)
    for chrom, start, end, strand in zip(
        df["chrom"], df["start"], df["end"], df["strand"]
    ):
        seq = genome.get(chrom).ff.fetch(chrom, start, end).upper()
        if strand == "-":
            seq = str(Seq.Seq(seq).reverse_complement())

        print(seq)


def find_CBS(cfg: dict):
    df = pd.read_csv(
        cfg["data_dir"] / "GSE235386" / "GSM7501570_esc.bed.gz",
        sep="\t",
        names=["chrom", "start", "end"],
    )
    range = pd.DataFrame({
        "chrom": [cfg["chrom"]],
        "start": [cfg["start"]],
        "end": [cfg["end"]],
    })
    df_CBS = bf.overlap(df, range, how="inner")[["chrom", "start", "end"]]
    pCBS = pd.read_csv(cfg["data_dir"] / "pCBS.csv", header=0)
    df_pCBS = bf.overlap(df_CBS, pCBS, how="inner")
    breakpoint()


def align_CBS(cfg: dict) -> None:
    df_AB = pd.read_csv(cfg["data_dir"] / "output" / "AB.csv", header=0)
    for name, file in cfg["CBS"].items():
        df_CBS = pd.read_csv(cfg["data_dir"] / file, header=0)
        df_strand = bf.closest(df_CBS, df_AB, return_index=True)[
            ["strand", "index_"]
        ].set_index(keys="index_")
        df_AB[name] = df_strand
        df_AB[name] = df_AB[name].fillna(".")

    df_AB.to_csv(cfg["data_dir"] / "output" / "AB_CBS.csv", index=False)
