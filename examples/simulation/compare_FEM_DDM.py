"""
Example script to compare FEM simulation with DDM simulation for a beam under flexion.
"""
import numpy as np

from pyLatticeSim.lattice_sim import LatticeSim
from pyLattice.plotting_lattice import LatticePlotting
from pyLatticeSim.utils_simulation import solve_FEM_FenicsX

path = "simulation/"
name_file = "simulation_beam_flexion"

lattice_Sim_object = LatticeSim(path + name_file)

sol_FEM = solve_FEM_FenicsX(lattice_Sim_object)[0]

# Visualization with matplotlib
vizualizer = LatticePlotting()
vizualizer.visualize_lattice(lattice_Sim_object, beam_color_type="radii",
                             enable_boundary_conditions=True,
                             deformedForm=True)

lattice_object = LatticeSim(path + name_file, enable_domain_decomposition_solver = True)

sol_DDM = lattice_object.solve_DDM()[0]

relative_error = np.linalg.norm(sol_FEM - sol_DDM) / np.linalg.norm(sol_FEM)
print("Relative error between FEM and DDM", relative_error)

vizualizer = LatticePlotting()
vizualizer.visualize_lattice(lattice_object, beam_color_type="radii",
                             enable_boundary_conditions=True,
                             deformedForm=True)