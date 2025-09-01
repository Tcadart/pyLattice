"""
Example of a domain decomposition simulation using pyLatticeSim.
"""
import numpy as np

from pyLatticeSim.lattice_sim import LatticeSim
from pyLattice.plotting_lattice import LatticePlotting
from pyLatticeSim.utils_schur import load_schur_complement_dataset
from pyLatticeSim.utils_simulation import solve_FEM_FenicsX


name_file = "simulation/3PointBending"

solver_DDM = LatticeSim(name_file, verbose=1, enable_domain_decomposition_solver=True)

enable_precondioner = True
number_iteration_max = 1500
solver_DDM.define_parameters(enable_precondioner, number_iteration_max)

xsol = solver_DDM.solve_DDM()[0]


# Visualization with matplotlib
vizualizer = LatticePlotting()
vizualizer.visualize_lattice(solver_DDM, beam_color_type="radii",
                             enable_boundary_conditions=True,
                             deformedForm=True)

lattice_Sim_object = LatticeSim(name_file)

sol, simulation_lattice = solve_FEM_FenicsX(lattice_Sim_object)

np.set_printoptions(precision=6, suppress=False, formatter={'float_kind': '{: .6e}'.format})
for xi, si in zip(np.ravel(xsol), np.ravel(sol)):
    print(f"{xi: .6e}   {si: .6e}")

error_rel = np.linalg.norm(xsol - sol, ord=2) / np.linalg.norm(sol, ord=2)
print("Relative error between DDM and direct solver:", error_rel)


# Visualization with matplotlib
vizualizer = LatticePlotting()
vizualizer.visualize_lattice(lattice_Sim_object, beam_color_type="radii",
                             enable_boundary_conditions=True,
                             deformedForm=True)
