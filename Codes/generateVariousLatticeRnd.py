import random
import os
import numpy as np
import matplotlib.pyplot as plt
from mpl_toolkits.mplot3d import Axes3D
from Lattice import Lattice
from saveLattice import save_lattice


# def generate_lattices(num_lattices, output_dir):
#     if not os.path.exists(output_dir):
#         os.makedirs(output_dir)
#
#     lattice_types = [0, 1, 4, 5, 11]  # Example lattice types
#     radius_values = [0.05, 0.1, 0.15]  # Different radii
#     grad_rules = ['constant', 'linear', 'parabolic', 'sinusoide', 'exponential']
#     num_cells_options = [2, 3, 4]  # Different lattice sizes
#     multimat_options = [0, 1, -1]
#     grad_directions = [[1, 0, 0], [0, 1, 0], [0, 0, 1]]
#     grad_parameters = [[1.0, 0.0, 2.0], [0.5, 0.0, 1.0]]
#
#     for i in range(num_lattices):
#         lattice_params = {
#             'cell_size_x': 1,
#             'cell_size_y': 1,
#             'cell_size_z': 1,
#             'num_cells_x': random.choice(num_cells_options),
#             'num_cells_y': random.choice(num_cells_options),
#             'num_cells_z': random.choice(num_cells_options),
#             'Lattice_Type': random.choice(lattice_types),
#             'Radius': random.choice(radius_values),
#             'gradRadiusProperty': [random.choice(grad_rules), random.choice(grad_directions),
#                                    random.choice(grad_parameters)],
#             'gradDimProperty': [random.choice(grad_rules), random.choice(grad_directions),
#                                 random.choice(grad_parameters)],
#             'gradMatProperty': [random.choice(multimat_options), random.randint(1, 3)],
#             'simMethod': random.choice([0, 1]),
#             'uncertaintyNode': random.choice([0, 1])
#         }
#
#         lattice = Lattice(**lattice_params)
#         filename = f"{output_dir}/lattice_{i + 1}.png"
#         save_lattice(lattice, filename)
#         print(f"Lattice {i + 1} saved as '{filename}'")

def generate_lattices(num_lattices, output_dir):
    if not os.path.exists(output_dir):
        os.makedirs(output_dir)

    lattice_types = [0, 1, 4, 5, 11]
    radius_values = [0.05, 0.1, 0.15]
    grad_rules = ['constant', 'linear', 'parabolic', 'sinusoide', 'exponential']
    num_cells_options = [2, 3, 4]
    multimat_options = [0, 1, -1]
    grad_directions = [[1, 0, 0], [0, 1, 0], [0, 0, 1]]
    grad_parameters = [[1.0, 0.0, 2.0], [0.5, 0.0, 1.0]]
    hybridLatticeDataOptions = [[0.1, 0.2, 0.3], [0.05, 0.15, 0.25]]

    for i in range(num_lattices):
        # Randomly decide if this should be a simple or hybrid lattice
        if random.choice([True, False]):  # 50% chance for each
            lattice_params = {
                'cell_size_x': 1,
                'cell_size_y': 1,
                'cell_size_z': 1,
                'num_cells_x': random.choice(num_cells_options),
                'num_cells_y': random.choice(num_cells_options),
                'num_cells_z': random.choice(num_cells_options),
                'Lattice_Type': random.choice(lattice_types),
                'Radius': random.choice(radius_values),
                'gradRadiusProperty': [random.choice(grad_rules), random.choice(grad_directions), random.choice(grad_parameters)],
                'gradDimProperty': [random.choice(grad_rules), random.choice(grad_directions), random.choice(grad_parameters)],
                'gradMatProperty': [random.choice(multimat_options), random.randint(1, 3)],
                'simMethod': random.choice([0, 1]),
                'uncertaintyNode': random.choice([0, 1])
            }
            lattice = Lattice(**lattice_params)
        else:
            # Generate a hybrid lattice
            lattice_params = {
                'cell_size_x': 1,
                'cell_size_y': 1,
                'cell_size_z': 1,
                'simMethod': random.choice([0, 1]),
                'uncertaintyNode': random.choice([0, 1]),
                'hybridLatticeData': random.choice(hybridLatticeDataOptions)
            }
            lattice = Lattice.hybridgeometry(**lattice_params)

        filename = f"{output_dir}/lattice_{i + 1}.png"
        save_lattice(lattice, filename)
        print(f"Lattice {i + 1} saved as '{filename}'")