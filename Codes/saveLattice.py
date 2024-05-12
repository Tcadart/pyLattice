import matplotlib.pyplot as plt
from mpl_toolkits.mplot3d import Axes3D
import numpy as np
import os
from Lattice import *

def save_lattice(lattice, filename):
    fig = plt.figure()
    ax = fig.add_subplot(111, projection='3d')
    for beam in lattice.beams_obj:
        p1 = beam.point1
        p2 = beam.point2
        ax.plot([p1.x, p2.x], [p1.y, p2.y], [p1.z, p2.z], 'b-')
    ax.set_xlabel('X')
    ax.set_ylabel('Y')
    ax.set_zlabel('Z')
    plt.savefig(filename)
    plt.close(fig)


def generate_and_save_lattices(num_samples, lattice_params, path):
    if not os.path.exists(path):
        os.makedirs(path)

    for i in range(num_samples):
        lattice = Lattice(**lattice_params)
        save_lattice(lattice, os.path.join(path, f'lattice_{i + 1}.png'))
        print(f"Lattice {i + 1} saved at '{path}'")