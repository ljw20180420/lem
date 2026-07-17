import fresnel as fl
import numpy as np
from matplotlib.cm import get_cmap
from matplotlib.colors import Normalize

import polychrom_hoomd.utils as utils


def fresnel(
    snap,
    cmap="viridis",
    rescale_backbone_bonds=1.0,
    show=None,
    color_array=None,
    **kwargs,
):
    """
    Wrapper around polykit.renderers.backends for HooMD rendering using the Fresnel library
    """

    bonds = snap.bonds.group.copy()
    positions = utils.unwrap_coordinates(snap)

    bond_mask = np.ones(snap.particles.N, dtype=bool)
    polymer_mask = np.ones(snap.particles.N, dtype=bool)

    LEF_typeid = snap.bonds.types.index("LEF") if "LEF" in snap.bonds.types else -1

    polymer_mask[bonds] = False
    bond_mask[bonds[snap.bonds.typeid == LEF_typeid]] = False

    bond_mask[polymer_mask] = False

    radii = snap.particles.diameter.copy() * 0.5
    radii[bond_mask] *= rescale_backbone_bonds

    if isinstance(color_array, np.ndarray):
        colorscale = color_array

    else:
        colorscale = np.zeros(snap.particles.N)

        if show == "chromosomes":
            chrom_bounds = utils.get_chrom_bounds(snap)

            for i, bounds in enumerate(chrom_bounds):
                colorscale[bounds[0] : bounds[1] + 1] = i + 1

        elif show == "loops":
            loop_bounds = bonds[snap.bonds.typeid == LEF_typeid]

            for i, bounds in enumerate(loop_bounds):
                colorscale[bounds[0] : bounds[1] + 1] = i + 1

        elif show == "compartments":
            colorscale = snap.particles.typeid.copy()

        elif show == "strains":
            strains = np.diff(positions[bonds], axis=1)
            colorscale = np.linalg.norm(strains, axis=-1).flatten()

        else:
            colorscale = np.arange(snap.particles.N)

    colors = get_cmap(cmap)(Normalize()(colorscale))

    return _fresnel(positions, bonds, colors, radii, **kwargs)


def _fresnel(
    self,
    positions,
    bonds,
    colors,
    radii,
    intensity=0.0,
    metal=0.7,
    specular=0.8,
    spec_trans=0.1,
    roughness=0.1,
    outline=0.04,
):
    """
    Render individual polymer/particle configurations using the Fresnel backend library

    Parameters
    ----------
    positions : Nx3 float array
        List of 3D positions of the monomers to be displayed
    bonds : Mx2 int array
        List of pairwise inter-monomer bonds to be displayed
    colors : Nx3 or Nx4 or Mx3 or Mx4 float array
        List of RGB colors to be assigned to each monomer or bond
    radii : Mx1 float array
        List of bond/particle radii
    intensity : float
        Intensity of extra light for gamma correction
    roughness : float
        Roughness of the rendering material. Nominally in the range [0.1,1]
    metal : float
        Set to 0 for dielectric materials, or 1 for metals. Intermediate values interpolate between the 2
    specular : float
        Controls the strength of specular highlights. Nominally in the range [0.1,1]
    spec_trans : float
        Controls the amount of specular light transmission. In the range [0,1]
    outline : float
        Width of the outline material

    Returns
    -------
        IPython image object suitable for embedding in Jupyter notebooks
    """

    scene = fl.Scene()

    geometry = fl.geometry.Cylinder(scene, N=bonds.shape[0], outline_width=outline)

    geometry.points[:] = positions[bonds]
    geometry.radius[:] = radii[bonds].min(axis=1)

    corrected_colors = fl.color.linear(colors)

    if corrected_colors.shape[0] == positions.shape[0]:
        geometry.color[:] = corrected_colors[bonds]
    elif corrected_colors.shape[0] == bonds.shape[0]:
        geometry.color[:] = corrected_colors[:, None, :]
    else:
        raise ValueError("Color array does not match particle or bond dimensions")

    geometry.material = fl.material.Material(
        color=fl.color.linear([0.25, 0.25, 0.25]),
        roughness=roughness,
        metal=metal,
        specular=specular,
        spec_trans=spec_trans,
        primitive_color_mix=1.0,
        solid=0.0,
    )
    geometry.outline_material = fl.material.Material(
        color=fl.color.linear([0.25, 0.25, 0.25]),
        roughness=roughness,
        metal=metal,
        specular=specular,
        spec_trans=spec_trans,
        primitive_color_mix=0.0,
        solid=0.0,
    )

    polymer_mask = np.ones(positions.shape[0], dtype=bool)
    polymer_mask[bonds] = False

    num_unbound_atoms = np.count_nonzero(polymer_mask)

    if num_unbound_atoms > 0:
        geometry2 = fl.geometry.Sphere(
            scene, N=num_unbound_atoms, outline_width=outline
        )

        geometry2.radius[:] = radii[polymer_mask]
        geometry2.position[:] = positions[polymer_mask]

        geometry2.color[:] = corrected_colors[polymer_mask]

        geometry2.material = geometry.material
        geometry2.outline_material = geometry.outline_material

    scene.camera = fl.camera.Orthographic.fit(scene, view="isometric", margin=0)
    scene.lights.append(
        fl.light.Light(direction=[0, 0, 1], color=[intensity] * 3, theta=np.pi)
    )

    return scene
