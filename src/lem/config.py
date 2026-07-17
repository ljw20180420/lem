import os
from pathlib import Path


def force_fudenberg2026() -> dict:
    return {
        "Angular forces": {"Curvature": {"Stiffness": 2, "Type": "Harmonic"}},
        "Bonded forces": {
            "Backbone": {
                "Rest length": 1.0,
                "Type": "Harmonic",
                "Wiggle distance": 0.1,
            },
            "LEF": {"Rest length": 0.5, "Type": "Harmonic", "Wiggle distance": 0.2},
            "LEF_dummy": {"Rest length": 0, "Type": "Harmonic", "Wiggle distance": 0},
        },
        "External forces": {"Confinement": {}},
        "Non-bonded forces": {
            "Repulsion": {
                "Cutoff": 1.0,
                "Matrix": {"A": {"A": 5.0, "B": 5.0}, "B": {"A": 5.0, "B": 5.0}},
                "Type": "Polychrom",
            },
            "Attraction": {
                "Cutoff": 1.5,
                "Matrix": {"A": {"A": 0, "B": 0}, "B": {"A": 0, "B": 0.05}},
                "Type": "Polychrom",
            },
        },
    }


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
        "n_copies": 100,
        "seed": 63036,
        "force": force_fudenberg2026(),
    }
