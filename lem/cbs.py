import io
import os
import subprocess

import numpy as np
import pandas as pd
import py2bit
from Bio import motifs
from Bio.Seq import Seq


def get_pCBS(cfg: dict) -> None:
    df_shift = pd.read_csv("pCBS_shift.txt", header=0)
    df_cpcdh = pd.read_csv(cfg["data_dir"] / "result" / "cpcdh.csv", header=0)
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
    (cfg["data_dir"] / "result" / "CBS").mkdir(exist_ok=True, parents=True)
    df.to_csv(
        cfg["data_dir"] / "result" / "CBS" / "pCBS.bed",
        sep="\t",
        header=False,
        index=False,
    )


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


def get_eCBS(cfg: dict) -> None:
    (cfg["data_dir"] / "result" / "CBS").mkdir(exist_ok=True, parents=True)

    get_motif(
        seq_file="human_eCBS.csv",
        matrix_file=cfg["data_dir"] / "result" / "CBS" / "eCBS.txt",
        meme_file=cfg["data_dir"] / "result" / "CBS" / "eCBS.meme",
    )

    df_alpha = (
        pd
        .read_csv(cfg["data_dir"] / "result" / "cpcdh.csv", header=0)
        .query(r"name.str.match(r'^Pcdha\d')")
        .reset_index(drop=True)[["chrom", "start", "end", "name"]]
    )
    starts = []
    ends = []
    scores = []
    strands = []
    for chrom, start, end, name in df_alpha.itertuples(index=False):
        end = (start + end) // 2
        proc = subprocess.run(
            args=["./call_cbs.sh", "ecbs", chrom, str(start), str(end)],
            capture_output=True,
            text=True,
            check=True,
        )
        df_fimo = (
            pd
            .read_csv(io.StringIO(proc.stdout), sep="\t", header=0)
            .query("strand == '+'")
            .reset_index(drop=True)
        )
        max_row = df_fimo.loc[df_fimo["score"].idxmax()]
        starts.append(max_row["start"])
        ends.append(max_row["stop"] + 1)
        scores.append(max_row["score"])
        strands.append(max_row["strand"])

    pd.DataFrame({
        "chrom": df_alpha["chrom"],
        "start": starts,
        "end": ends,
        "name": df_alpha["name"],
        "score": scores,
        "strand": strands,
    }).to_csv(
        cfg["data_dir"] / "result" / "CBS" / "eCBS.bed",
        sep="\t",
        index=False,
        header=False,
    )


def correct_tCBS(cfg: dict) -> None:
    df = pd.read_csv(
        "tCBS.bed", sep="\t", names=["chrom", "start", "end", "name", "score", "strand"]
    )

    seqs = []
    with py2bit.open(cfg["2bit"]) as tb:
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

    (cfg["data_dir"] / "result" / "CBS").mkdir(exist_ok=True, parents=True)

    get_motif(
        seq_file=cfg["data_dir"] / "result" / "CBS" / "tCBS.csv",
        matrix_file=cfg["data_dir"] / "result" / "CBS" / "tCBS.txt",
        meme_file=cfg["data_dir"] / "result" / "CBS" / "tCBS.meme",
    )

    for cluster, prefix in zip(["alpha", "beta", "gamma"], ["Pcdha", "Pcdhb", "Pcdhg"]):
        pd.read_csv(cfg["data_dir"] / "result" / "tCBS.csv", header=0).query(
            "gene.str.startswith(@prefix)"
        ).to_csv(cfg["data_dir"] / "result" / f"tCBS.{cluster}.csv", index=False)

        get_motif(
            seq_file=cfg["data_dir"] / "result" / "CBS" / f"tCBS.{cluster}.csv",
            matrix_file=cfg["data_dir"] / "result" / "CBS" / f"tCBS.{cluster}.txt",
            meme_file=cfg["data_dir"] / "result" / "CBS" / f"tCBS.{cluster}.meme",
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


def bed2seq(
    cfg: dict,
    bed: os.PathLike,
) -> None:
    df = pd.read_csv(
        bed,
        sep="\t",
        names=["chrom", "start", "end", "name", "score", "strand"],
    )

    with py2bit.open(cfg["2bit"]) as tb:
        for chrom, start, end, name, strand in zip(
            df["chrom"], df["start"], df["end"], df["name"], df["strand"]
        ):
            seq = tb.sequence(chrom, start, end)
            if strand == "-":
                seq = str(Seq(seq).reverse_complement())
            print(seq, name)
