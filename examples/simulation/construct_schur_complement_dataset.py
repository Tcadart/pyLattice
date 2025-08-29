"""
Construct a dataset of Schur complements for parametried cells.
"""

from pyLatticeSim.lattice_sim import LatticeSim
from pyLatticeSim.utils_simulation import get_schur_complement


name_file = "simulation/hybrid_cell_simulation"

lattice_object = LatticeSim(name_file)

schur_complement = get_schur_complement(lattice_object)

print("Schur complement matrix:\n", schur_complement)

lattice_object.reset_cell_with_new_radii([0.01, 0.02, 0.03])

schur_complement_updated = get_schur_complement(lattice_object)

print("Updated Schur complement matrix:\n", schur_complement_updated)
