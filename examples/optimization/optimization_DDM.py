"""
Examples of a simple optimization case.
"""
from pyLattice.plotting_lattice import LatticePlotting
from pyLatticeOpti.lattice_opti import LatticeOpti
from pyLattice.utils import save_JSON_to_Grasshopper

path = "optimization/"
name_file = "optimization_DDM_surrogate"

lattice_object = LatticeOpti(path + name_file, verbose=1, convergence_plotting = True)

lattice_object.solve_DDM()
lattice_object.optimize_lattice()

# lattice_object.reset_penalized_beams()
# save_JSON_to_Grasshopper(lattice_object, name_file + "_optimized")

# Visualization optimized lattice
vizualizer = LatticePlotting()
vizualizer.visualize_lattice(lattice_object, beam_color_type="radii", enable_boundary_conditions=True,
                             deformedForm = True)
