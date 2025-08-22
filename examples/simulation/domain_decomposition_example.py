"""
Example of a domain decomposition simulation using pyLatticeSim.
"""

from pyLatticeSim.lattice_sim import LatticeSim
from pyLatticeSim.domain_decomposition_solver import DomainDecompositionSolver
from pyLattice.plotting_lattice import LatticePlotting


name_file = "3PointBending"

solver_DDM = DomainDecompositionSolver(name_file, verbose=1)

enable_precondioner = True
number_iteration_max = 100
solver_DDM.define_parameters(enable_precondioner, number_iteration_max)

solver_DDM.solve_DDM()


# Visualization with matplotlib
vizualizer = LatticePlotting()
vizualizer.visualize_lattice(solver_DDM, beam_color_type="radii",
                             enable_boundary_conditions=True,
                             deformedForm=True)