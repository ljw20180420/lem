import io
import os
import subprocess

import numpy as np
import pandas as pd
import py2bit
from Bio import motifs
from Bio.Seq import Seq


def get_pCBS(shift_file: os.PathLike, cpcdh_file: os.PathLike) -> pd.DataFrame:
    df_shift = pd.read_csv(shift_file, header=0)
    df_cpcdh = pd.read_csv(cpcdh_file, header=0)
    df = df_shift.merge(
        right=df_cpcdh, how="inner", left_on="gene", right_on="name", validate="1:1"
    )
    df = df.astype({
        "CDS_start": float,
        "CDS_end": float,
    }).astype({
        "CDS_start": int,
        "CDS_end": int,
    })
    df = df.assign(
        end=lambda df: df["CDS_start"] + df["shift"] + 2,
        start=lambda df: df["end"] - 27,
        score=0,
    )[["chrom", "start", "end", "name", "score", "strand"]]

    return df


def get_motif(seq_file: os.PathLike, matrix_file: os.PathLike, meme_file: os.PathLike):
    df = pd.read_csv(seq_file, header=0)
    m = motifs.create([Seq(seq) for seq in df["seq"]])

    matrix = np.array([m.counts[base] for base in ["A", "C", "G", "T"]])
    np.savetxt(matrix_file, matrix, fmt="%d")
    subprocess.run(
        args=[f"matrix2meme <{os.fspath(matrix_file)} >{os.fspath(meme_file)}"],
        shell=True,
        check=False,
    )


def correct_tCBS(cfg) -> None:
    df = pd.read_csv(
        "tCBS.bed", sep="\t", names=["chrom", "start", "end", "name", "score", "strand"]
    )

    seqs = []
    with py2bit.open("/home/ljw/sdb1/ucsc/hubs/myHub/lmm10/lmm10.2bit") as tb:
        for chrom, start, end, strand in zip(
            df["chrom"], df["start"], df["end"], df["strand"]
        ):
            seq = tb.sequence(chrom, start, end)
            if strand == "-":
                seq = str(Seq(seq).reverse_complement())
            seqs.append(seq)

    pd.DataFrame({"gene": df["name"], "seq": seqs}).to_csv(
        cfg["data_dir"] / "result" / "tCBS.csv", header=["gene", "seq"], index=False
    )

    get_motif(
        seq_file=cfg["data_dir"] / "result" / "tCBS.csv",
        matrix_file=cfg["data_dir"] / "result" / "tCBS.txt",
        meme_file=cfg["data_dir"] / "result" / "tCBS.meme",
    )

    for cluster, prefix in zip(["alpha", "beta", "gamma"], ["Pcdha", "Pcdhb", "Pcdhg"]):
        pd.read_csv(cfg["data_dir"] / "result" / "tCBS.csv", header=0).query(
            "gene.str.startswith(@prefix)"
        ).to_csv(cfg["data_dir"] / "result" / f"tCBS.{cluster}.csv", index=False)

        get_motif(
            seq_file=cfg["data_dir"] / "result" / f"tCBS.{cluster}.csv",
            matrix_file=cfg["data_dir"] / "result" / f"tCBS.{cluster}.txt",
            meme_file=cfg["data_dir"] / "result" / f"tCBS.{cluster}.meme",
        )


def update_tCBS_score() -> None:
    df = pd.read_csv(
        "tCBS.bed", sep="\t", names=["chrom", "start", "end", "name", "score", "strand"]
    )
    df.to_csv("tCBS.bed.bak", sep="\t", index=False, header=False)

    scores = []
    for chrom, start, end, name, score, strand in df.itertuples(index=False):
        if name.startswith("Pcdha"):
            cmd = "ta"
        elif name.startswith("Pcdhb"):
            cmd = "tb"
        elif name.startswith("Pcdhg"):
            cmd = "tg"
        else:
            raise ValueError("Unknown name prefix")

        proc = subprocess.run(
            args=["./call_cbs.sh", cmd, chrom, str(start), str(end)],
            check=True,
            capture_output=True,
            text=True,
        )

        df_fimo = pd.read_csv(io.StringIO(proc.stdout), sep="\t", header=0)
        max_row = df_fimo.loc[df_fimo["score"].idxmax()]
        if max_row["strand"] != strand:
            raise ValueError("Inconsistent strand")
        scores.append(max_row["score"])

    df["score"] = scores
    df.to_csv("tCBS.bed", sep="\t", index=False, header=False)
