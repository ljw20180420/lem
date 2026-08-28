import pandas as pd


def get_cpcdh_exon(cfg: dict) -> None:
    df = pd.read_csv(
        cfg["data_dir"] / "data" / "mm10.refGene.gtf.gz",
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
        "chrom == 'chr18' and attributes.str.contains(r'Pcdh[abg][abc]?[0-9]{1,2}') and (feature=='exon' or feature=='CDS')"
    ).reset_index(drop=True)
    attributes = df["attributes"].str.split(expand=True)
    df = df.assign(
        start=lambda df: df["start"] - 1,
        name=attributes[9].str.strip('";'),
        transcript_id=attributes[3].str.strip('";'),
        exon_number=attributes[5].str.strip('";').astype(int),
        total_exon_number=lambda df: df.groupby("transcript_id")[
            "exon_number"
        ].transform(max),
    ).drop(
        columns=[
            "source",
            "score",
            "frame",
            "attributes",
        ]
    )

    df = (
        df
        .query("transcript_id.str.startswith('NM')")
        .query(
            "exon_number == 1 and (total_exon_number == 4 or name.str.contains(r'^Pcdhb')) or name == 'Pcdha1' or name == 'Pcdhga1'"
        )
        .reset_index(drop=True)
    )
    for exon in ["Pcdha1", "Pcdhga1"]:
        for exon_number in [2, 3, 4]:
            prefix = "ace" if exon == "Pcdha1" else "gce"
            df.loc[
                (df["exon_number"] == exon_number) & (df["name"] == exon), "name"
            ] = f"{prefix}{exon_number - 1}"

    df = df.pivot_table(
        values=["start", "end"],
        index=["chrom", "strand", "name", "transcript_id"],
        columns="feature",
    )
    df.columns = df.columns.to_flat_index().map(lambda tp: f"{tp[1]}_{tp[0]}")

    df = (
        df
        .reset_index()
        .rename(columns={"exon_start": "start", "exon_end": "end"})
        .assign(score=".")[
            ["chrom", "start", "end", "name", "score", "strand", "CDS_start", "CDS_end"]
        ]
        .assign(
            CDS_start=lambda df: df["CDS_start"].fillna("."),
            CDS_end=lambda df: df["CDS_end"].fillna("."),
        )
        .sort_values(by=["start", "end"], ignore_index=True)
        .astype({
            "start": int,
            "end": int,
        })
    )

    (cfg["data_dir"] / "result").mkdir(exist_ok=True, parents=True)
    df.to_csv(cfg["data_dir"] / "result" / "cpcdh.csv", index=False)
