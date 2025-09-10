"""
Examples of a simple optimization case.
"""
import numpy as np

from pyLattice.plotting_lattice import LatticePlotting
from pyLatticeOpti.lattice_opti import LatticeOpti
from pyLattice.utils import save_JSON_to_Grasshopper
from pyLatticeSim.lattice_sim import LatticeSim
from pyLatticeSim.utils_simulation import solve_FEM_FenicsX

path = "optimization/"
name_file = "optimization_DDM_surrogate"

lattice_Sim_object = LatticeSim(path + name_file, enable_domain_decomposition_solver = False)

x_fem = solve_FEM_FenicsX(lattice_Sim_object)[0]
print("FEM displacement", x_fem)

# vizualizer = LatticePlotting()
# vizualizer.visualize_lattice(lattice_Sim_object, beam_color_type="radii",
#                              enable_boundary_conditions=True,
#                              deformedForm=True)

# lattice_object = LatticeOpti(path + name_file, verbose=1, convergence_plotting = True)
lattice_object = LatticeSim(path + name_file, enable_domain_decomposition_solver = True)

x_DDM = lattice_object.solve_DDM()[0]

print("DDM displacement", x_DDM)
relative_error = np.linalg.norm(x_fem - x_DDM) / np.linalg.norm(x_fem)
print("Relative error between FEM and DDM", relative_error)

# lattice_object.reset_penalized_beams()
# save_JSON_to_Grasshopper(lattice_object, name_file + "_optimized")

# Visualization optimized lattice
# vizualizer = LatticePlotting()
# vizualizer.visualize_lattice(lattice_object, beam_color_type="radii", enable_boundary_conditions=True,
#                              deformedForm = True)
