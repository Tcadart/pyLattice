"""
Construct a dataset of Schur complements for parametried cells.
"""
from itertools import product

import numpy as np

from pyLatticeSim.lattice_sim import LatticeSim
from pyLatticeSim.utils_schur import get_schur_complement, save_schur_complement_npz


name_file = "simulation/hybrid_cell_simulation"

lattice_object = LatticeSim(name_file)

start_radius = 0.01
end_radius = 0.11
step_radius = 0.01
radius_range = np.round(np.arange(start_radius, end_radius, step_radius), 3)

radius_values_batch = []
schur_matrix_batch = []

for i, radius_combinations in enumerate(product(radius_range, repeat=len(lattice_object.geom_types)), start=1):
    if sum(radius_combinations) <= 0.003:
        continue
    print(f"Combination {i}: {radius_combinations}")
    lattice_object.reset_cell_with_new_radii(list(radius_combinations))
    schur_complement = get_schur_complement(lattice_object)

    radius_values_batch.append(list(radius_combinations))
    schur_matrix_batch.append(schur_complement)
    if sum(radius_combinations) == 0.1:
        print(schur_complement)

save_schur_complement_npz(lattice_object, radius_values_batch, schur_matrix_batch)
