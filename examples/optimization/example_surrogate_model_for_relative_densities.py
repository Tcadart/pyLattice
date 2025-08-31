"""
Example file for generating a surrogate model to predict relative denities based on cell properties.
"""
from pyLattice.lattice import Lattice
from pyLatticeOpti.surrogate_model_relative_densities import compute_relative_densities_dataset, plot_3D_iso_surface

name_cell = "hybrid_cell"
name_file = "design/" + name_cell

hybrid_cell = Lattice(name_file, verbose=-1)

compute_relative_densities_dataset(hybrid_cell)

plot_3D_iso_surface(hybrid_cell)





