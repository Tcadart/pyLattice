"""
Example of lattice structure design with all entries from the package.
"""

from lattice import Lattice
from plotting_lattice import LatticePlotting

name_file = "full_entry_lattice"

lattice_object = Lattice.from_json(name_file, verbose=1)

vizualizer = LatticePlotting()
vizualizer.visualize_lattice_3D(lattice_object, beam_color_type="radii")
