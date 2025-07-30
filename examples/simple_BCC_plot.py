"""
Simple example of how to plot a simple BCC lattice using Matplotlib.
"""

from lattice import Lattice
from plotting_lattice import LatticePlotting


name_file = "simple_BCC"

lattice_object = Lattice.from_json(name_file)

vizualizer = LatticePlotting()
vizualizer.visualize_lattice_3D(lattice_object, beam_color_type="radii")
