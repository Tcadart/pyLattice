"""
Simple schur complement example of an hybrid cell
"""

from pyLattice.lattice import Lattice
from pyLatticeSim.utils_simulation import get_schur_complement


name_file = "hybrid_cell_simulation"

lattice_object = Lattice.from_json(name_file)

schur_complement = get_schur_complement(lattice_object)

print("Schur complement matrix:\n", schur_complement)
