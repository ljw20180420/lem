import os
from pathlib import Path


def cpcdh() -> dict:
    data_dir = Path("/home/ljw/sdc1/cpcdh")

    return {
        "data_dir": data_dir,
        "genome": "/home/ljw/.local/share/genomes/GRCm38/GRCm38.fa",
        "hic": os.fspath(data_dir / "wt.mcool::resolutions/10000"),
        "chrom": "chr18",
        "start": 36900000,
        "end": 37900000,
        "bin": 2500,
        "density": 0.2,
        "n_copies": 10,
    }
