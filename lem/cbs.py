import io
import os
import pathlib
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


def correct_motif(site_file: os.PathLike, two_bit_file: os.PathLike) -> None:
    df = pd.read_csv(
        site_file, sep="\t", names=["chrom", "start", "end", "name", "score", "strand"]
    )

    seqs = []
    with py2bit.open(two_bit_file) as tb:
        for chrom, start, end, strand in zip(
            df["chrom"], df["start"], df["end"], df["strand"]
        ):
            seq = tb.sequence(chrom, start, end)
            if strand == "-":
                seq = str(Seq(seq).reverse_complement())
            seqs.append(seq)

    pd.DataFrame({"gene": df["name"], "seq": seqs}).to_csv(
        "correct.csv", header=["gene", "seq"], index=False
    )
    get_motif(
        seq_file="correct.csv",
        matrix_file=pathlib.Path("correct.txt"),
        meme_file=pathlib.Path("correct.meme"),
    )

    chroms = []
    starts = []
    ends = []
    scores = []
    strands = []

    chroms_changed = []
    starts_changed = []
    ends_changed = []
    names_changed = []
    scores_changed = []
    strands_changed = []
    for chrom, start, end, name, score, strand in df.itertuples(index=False):
        peak_start = start - 500
        peak_end = end + 500
        proc = subprocess.run(
            args=[
                "./call_cbs.sh",
                "correct.meme",
                chrom,
                f"{peak_start}",
                f"{peak_end}",
            ],
            check=True,
            capture_output=True,
            text=True,
        )

        try:
            df_fimo = (
                pd
                .read_csv(io.StringIO(proc.stdout), sep="\t")
                .query("strand == @strand")
                .reset_index(drop=True)
            )
        except pd.errors.EmptyDataError as e:
            print(name)
            raise e
        if len(df_fimo) == 0:
            print(name)
            raise Exception("No correct strand found")

        max_row = df_fimo.loc[df_fimo["score"].idxmax()]

        chroms.append(max_row["sequence_name"])
        starts.append(max_row["start"])
        ends.append(max_row["stop"] + 1)
        scores.append(max_row["score"])
        strands.append(max_row["strand"])

        if (
            max_row["sequence_name"] != chrom
            or max_row["start"] != start
            or max_row["stop"] + 1 != end
            or max_row["strand"] != strand
        ):
            chroms_changed.append(chrom)
            starts_changed.append(start)
            ends_changed.append(end)
            names_changed.append(name)
            scores_changed.append(score)
            strands_changed.append(strand)

    pd.DataFrame({
        "chrom": chroms,
        "start": starts,
        "end": ends,
        "name": df["name"],
        "score": scores,
        "strand": strands,
    }).to_csv(site_file, sep="\t", header=False, index=False)

    pd.DataFrame({
        "chrom": chroms_changed,
        "start": starts_changed,
        "end": ends_changed,
        "name": names_changed,
        "score": scores_changed,
        "strand": strands_changed,
    }).to_csv(f"{os.fspath(site_file)}.changed", sep="\t", header=False, index=False)
