import os

import pandas as pd


def get_cpcdh_exon(gtffile: os.PathLike) -> pd.DataFrame:
    df = pd.read_csv(
        gtffile,
        sep="\t",
        names=[
            "chrom",
            "source",
            "feature",
            "start",
            "end",
            "score",
            "strand",
            "frame",
            "attributes",
        ],
    )
    df = df.query(
        "chrom == 'chr18' and attributes.str.contains(r'Pcdh[abg][abc]?[0-9]{1,2}') and feature=='exon'"
    ).reset_index(drop=True)
    attributes = df["attributes"].str.split(expand=True)
    df = df.assign(
        name=attributes[9].str.strip('";'),
        transcript_id=attributes[3].str.strip('";'),
        exon_number=attributes[5].str.strip('";').astype(int),
        total_exon_number=lambda df: df.groupby("transcript_id")[
            "exon_number"
        ].transform(max),
    ).drop(
        columns=[
            "feature",
            "source",
            "score",
            "frame",
            "attributes",
        ]
    )

    df = (
        df
        .query(
            "exon_number == 1 and (total_exon_number == 4 or name.str.contains(r'^Pcdhb')) or name == 'Pcdha1' or name == 'Pcdhga1'"
        )
        .reset_index(drop=True)
        .assign(score=".")
    )
    for exon in ["Pcdha1", "Pcdhga1"]:
        for exon_number in [2, 3, 4]:
            prefix = "ace" if exon == "Pcdha1" else "gce"
            df.loc[
                (df["exon_number"] == exon_number) & (df["name"] == exon), "name"
            ] = f"{prefix}{exon_number - 1}"

    df = df[["chrom", "start", "end", "name", "score", "strand"]].sort_values(
        by=["start", "end"], ignore_index=True
    )

    return df
