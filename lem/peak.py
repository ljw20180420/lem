import os
import subprocess


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
                / "ESC.CTCF.merged.sort.bam_RPKM.mm10.narrowPeak"
            ),
        ],
        check=True,
    )
