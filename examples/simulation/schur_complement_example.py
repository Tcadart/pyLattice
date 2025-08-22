"""
Simple schur complement example of an hybrid cell
"""

from pyLatticeSim.lattice_sim import LatticeSim
from pyLatticeSim.utils_simulation import get_schur_complement


name_file = "hybrid_cell_simulation"

lattice_object = LatticeSim(name_file)

schur_complement = get_schur_complement(lattice_object)

print("Schur complement matrix:\n", schur_complement)
