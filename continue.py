import json

import polychrom_hoomd.extrude as extrude
import polychrom_hoomd.render as render
from onestate_extruder import LEFTranslocatorDirectional, compute_LEF_pos

# Run
system.run(1e5)

# Visualize new conformation
snapshot_relaxed = system.state.get_snapshot()
render.fresnel(snapshot_relaxed, show="compartments", cmap="coolwarm").static(
    pathtrace=True
)

# Parse extrusion parameters
with open("onestate_extrusion_dict.json", "r") as dict_file:
    extrusion_dict = json.load(dict_file)

# Lattice and polymer time units (in seconds)
tau_1d = 20  # [s], yields a baseline extrusion rate of 2*2.5/20 = 0.25 kb/s
tau_3d = 0.005  # [s], from fitting polymer kinetics to experimental MSDs

# Rescale baseline velocity to 1 kb/s
velocity_multiplier = 4

tau_1d_rescaled = tau_1d / velocity_multiplier
n_3d_to_1d = tau_1d_rescaled / tau_3d

# Convert residence times from seconds to lattice units
extrusion_dict["LEF_lifetime"] /= tau_1d_rescaled

# Precompute 1D loop extrusion trajectory
lef_trajectory = compute_LEF_pos(
    LEFTranslocatorDirectional, len(chromatin_types), **extrusion_dict
)

# Run 3D simulation
for lef_positions in lef_trajectory:
    extrude.update_topology(system, lef_positions, thermalize=False)
    system.run(n_3d_to_1d)

# Visualize new conformation, highlighting extruder positions
snapshot_final = system.state.get_snapshot()
render.fresnel(
    snapshot_final, show="compartments", cmap="coolwarm", rescale_backbone_bonds=0.05
).static(pathtrace=True)
