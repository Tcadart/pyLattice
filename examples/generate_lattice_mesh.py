"""
Example script to generate a lattice mesh
"""

from pyLattice.lattice import Lattice
from pyLattice.plotting_lattice import LatticePlotting

name_file = "hybrid_cell"

lattice_object = Lattice.from_json(name_file)

lattice_object.generate_mesh_lattice_Gmsh(cutMeshAtBoundary=True, volume_computation=True)

lattice_object.timing.summary()

vizualizer = LatticePlotting()
vizualizer.visualize_lattice_3D(lattice_object, beam_color_type="radii")
