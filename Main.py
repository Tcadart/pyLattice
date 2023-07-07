from Lattice import *
import matplotlib.pyplot as plt


def visu_Cellule(ax):
    ax.set_title("Cellule choisie")
    lattice.cells[0].visualize_3d(ax)


def visu_lattice(ax):
    ax.set_title("Lattice généré")
    lattice.visualize_3d(ax)


def display_lattice_points_beams():
    fig = plt.figure()
    ax1 = fig.add_subplot(1, 2, 1, projection='3d')
    visu_Cellule(ax1)
    ax2 = fig.add_subplot(1, 2, 2, projection='3d')
    visu_lattice(ax2)
    plt.show()


lattice = Lattice(1, 1, 1, 2, 2, 2)

cells = lattice.generate_custom_lattice(7)
# cells = lattice.generate_random_lattice(20, 40, 0.01)
display_lattice_points_beams()
# cells = lattice.generate_random_lattice(50, 80, 0.01)
# print(lattice.affichage_points_console())
# # print(lattice.affichage_points_console()[1])
# # print(lattice.affichage_points_console()[1][1])
# print("b", lattice.affichage_beams_console())
