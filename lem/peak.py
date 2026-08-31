import os
import subprocess

import pandas as pd


def crossmap(cfg: dict) -> None:
    (cfg["data_dir"] / "result" / "peak").mkdir(exist_ok=True, parents=True)
    subprocess.run(
        args=[
            "CrossMap",
            "bigwig",
            os.fspath(cfg["data_dir"] / "data" / "mm9ToMm10.over.chain"),
            os.fspath(cfg["data_dir"] / "data" / "ESC.CTCF.merged.sort.bam_RPKM.bw"),
            os.fspath(
                cfg["data_dir"]
                / "result"
                / "peak"
                / "ESC.CTCF.merged.sort.bam_RPKM.mm10",
            ),
        ],
        check=True,
    )

    subprocess.run(
        args=[
            "CrossMap",
            "bed",
            os.fspath(cfg["data_dir"] / "data" / "mm39ToMm10.over.chain"),
            os.fspath(cfg["data_dir"] / "data" / "tCBS_ltj.bed"),
            os.fspath(cfg["data_dir"] / "result" / "peak" / "tCBS_ltj_mm10.bed"),
        ],
        check=True,
    )

    df = pd.read_csv(
        cfg["data_dir"] / "result" / "peak" / "tCBS_ltj_mm10.bed",
        sep="\t",
        names=["chrom", "start", "end", "name", "score", "strand"],
    )
    df["name"] = df["name"].map(
        lambda ele: (
            ele
            .replace("mmPcdhα", "Pcdha")
            .replace("mmPcdhβ", "Pcdhb")
            .replace("mmPcdhγ", "Pcdhg")
        )
    )
    df.to_csv(
        cfg["data_dir"] / "result" / "peak" / "tCBS_ltj_mm10.bed",
        sep="\t",
        index=False,
        header=False,
    )


def call_peak(
    cfg: dict,
) -> None:
    (cfg["data_dir"] / "result" / "peak").mkdir(exist_ok=True, parents=True)
    subprocess.run(
        args=[
            "bigWigToBedGraph",
            os.fspath(
                cfg["data_dir"]
                / "result"
                / "peak"
                / "ESC.CTCF.merged.sort.bam_RPKM.mm10.bw"
            ),
            os.fspath(
                cfg["data_dir"]
                / "result"
                / "peak"
                / "ESC.CTCF.merged.sort.bam_RPKM.mm10.bdg"
            ),
        ],
        check=True,
    )

    subprocess.run(
        args=[
            "macs3",
            "bdgpeakcall",
            "-i",
            os.fspath(
                cfg["data_dir"]
                / "result"
                / "peak"
                / "ESC.CTCF.merged.sort.bam_RPKM.mm10.bdg"
            ),
            "-o",
            os.fspath(
                cfg["data_dir"]
                / "result"
                / "peak"
                / "ESC.CTCF.merged.sort.bam_RPKM.mm10.bed"
            ),
        ],
        check=True,
    )

    df_peak = pd.read_csv(
        cfg["data_dir"] / "result" / "peak" / "ESC.CTCF.merged.sort.bam_RPKM.mm10.bed",
        sep="\t",
        names=[
            "chrom",
            "chromStart",
            "chromEnd",
            "name",
            "score",
            "strand",
            "signalValue",
            "pValue",
            "qValue",
            "peak",
        ],
        skiprows=1,
    )

    df_peak.to_csv(
        cfg["data_dir"]
        / "result"
        / "peak"
        / "ESC.CTCF.merged.sort.bam_RPKM.mm10.narrowPeak",
        sep="\t",
        index=False,
        header=False,
    )
