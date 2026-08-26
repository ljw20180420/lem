import os
import subprocess

import numpy as np
import pandas as pd
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
    )[["chrom", "start", "end", "name", "score", "strand"]]

    return df


def get_motif(
    cfg: dict, seq_file: os.PathLike, matrix_file: os.PathLike, meme_file: os.PathLike
):
    df = pd.read_csv(seq_file, header=0)
    m = motifs.create([Seq(seq) for seq in df["seq"]])

    matrix = np.array([m.counts[base] for base in ["A", "C", "G", "T"]])
    np.savetxt(matrix_file, matrix, fmt="%d")
    subprocess.run(
        args=[f"matrix2meme <{os.fspath(matrix_file)} >{os.fspath(meme_file)}"],
        shell=True,
    )
