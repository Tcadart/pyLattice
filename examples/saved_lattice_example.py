"""
This example shows how to load a saved lattice and visualize it.
"""
from lattice import Lattice
from plotting_lattice import LatticePlotting

name_lattice = "Kelvin_helmet"

lattice = Lattice.open_pickle_lattice(name_lattice)

vizualizer = LatticePlotting()
vizualizer.visualize_lattice_3D(lattice, beam_color_type="radii")
vizualizer.show()
