"""
Simple homogenization example of a hybrid cell
"""

from pyLattice.lattice import Lattice
from pyLatticeSim.utils import create_homogenization_figure
from pyLatticeSim.utils_simulation import get_homogenized_properties
from pyLatticeSim.export_simulation_results import exportSimulationResults


name_file = "hybrid_cell_homogenization"

lattice_object = Lattice.from_json(name_file)

mat_S_orthotropic, homogenization_analysis = get_homogenized_properties(lattice_object)

create_homogenization_figure(mat_S_orthotropic, save=True)

# Export simulations to Paraview
exportData = exportSimulationResults(homogenization_analysis, name_file)
exportData.export_data_homogenization()