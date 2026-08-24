import os

import pandas as pd


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
        end=lambda df: df["CDS_start"] + df["shift"] + 1,
        start=lambda df: df["end"] - 42,
    )[["chrom", "start", "end", "name", "score", "strand"]]

    return df
