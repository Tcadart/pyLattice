"""
Simple example of simulation of a beam in flexion using pyLatticeSim.
"""

from pyLattice.lattice import Lattice
from pyLattice.plotting_lattice import LatticePlotting
from pyLatticeSim.utils_simulation import solve_FEM_FenicsX
from pyLatticeSim.export_simulation_results import exportSimulationResults


name_file = "simulation_beam_flexion"

lattice_object = Lattice.from_json(name_file)

sol, simulation_lattice = solve_FEM_FenicsX(lattice_object)

# Visualization with matplotlib
vizualizer = LatticePlotting()
vizualizer.visualize_lattice(lattice_object, beam_color_type="radii",
                             enable_boundary_conditions=True,
                             deformedForm=True)

# Export the results to Paraview
export_results = exportSimulationResults(simulation_lattice, name_file)
export_results.export_displacement_rotation()
export_results.write_function()
