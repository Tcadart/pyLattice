"""
Simple example of how to plot a simple BCC lattice using Matplotlib.
"""

from pyLattice.lattice import Lattice
from pyLattice.plotting_lattice import LatticePlotting


name_file = "design/simple_BCC"

lattice_object = Lattice(name_file)

vizualizer = LatticePlotting()
vizualizer.visualize_lattice(lattice_object, beam_color_type="radii")
