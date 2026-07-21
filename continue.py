import json

import hoomd
import numpy as np

import polychrom_hoomd.build as build
import polychrom_hoomd.extrude as extrude
import polychrom_hoomd.forces as forces
import polychrom_hoomd.log as log
import polychrom_hoomd.render as render
from onestate_extruder import LEFTranslocatorDirectional, compute_LEF_pos

# Generate RNG seed
rng_seed = np.random.randint([0, 2**16])

# Initialise HooMD on the CPU or GPU, based on availability
hoomd_device = build.get_hoomd_device()
# Initialize empty simulation object
system = hoomd.Simulation(device=hoomd_device, seed=rng_seed)

# Setup neighbor list
nl = hoomd.md.nlist.Cell(buffer=0.4)


# Setup HooMD simulation object
system.create_state_from_snapshot(snapshot)

# Set chromosome excluded volume
repulsion_forces = forces.get_repulsion_forces(nl, **force_dict)

# Set bonded/angular potentials
bonded_forces = forces.get_bonded_forces(**force_dict)
angular_forces = forces.get_angular_forces(**force_dict)

# Set attractive/DPD forces
dpd_forces = forces.get_dpd_forces(nl, **force_dict)
attraction_forces = forces.get_attraction_forces(nl, **force_dict)

# Define full force_field
dpd_force_field = (
    repulsion_forces + bonded_forces + angular_forces + attraction_forces + dpd_forces
)

# Setup integrator methods
nve = hoomd.md.methods.NVE(filter=hoomd.filter.All())
dpd_integrator = hoomd.md.Integrator(dt=5e-3, methods=[nve], forces=dpd_force_field)

# Set up logs and integrator objects
logger = log.get_logger(system)

system.operations.integrator = dpd_integrator
system.operations.writers.append(log.table_formatter(logger, period=1000))

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
