import bioframe as bf
import cooler
import cooltools
import matplotlib.pyplot as plt
import numpy as np
import pandas as pd
from matplotlib.ticker import EngFormatter


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
