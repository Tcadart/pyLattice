"""
Simple example of simulation of a beam in flexion using pyLatticeSim.
"""

from pyLattice.lattice import Lattice
from pyLattice.plotting_lattice import LatticePlotting
from pyLatticeSim.utils_simulation import solve_FEM_FenicsX


name_file = "simulation_beam_flexion"

lattice_object = Lattice.from_json(name_file)

solve_FEM_FenicsX(lattice_object)

vizualizer = LatticePlotting()
vizualizer.visualize_lattice(lattice_object, beam_color_type="radii", enable_boundary_conditions=True,
                             deformedForm=True)
