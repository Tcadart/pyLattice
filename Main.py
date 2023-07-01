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

print("Entrez les paramètres du Lattice :")

# num_cells_x = int(input("Nombre de cellules (en X) : "))
# num_cells_y = int(input("Nombre de cellules (en Y) : "))
# num_cells_z = int(input("Nombre de cellules (en Z) : "))
#
# lattice = Lattice(1, 1, 1, num_cells_x, num_cells_y, num_cells_z)
lattice = Lattice(1, 1, 1, 2, 2, 2)

print("Choisissez le mode : 1 pour rnd, 2 pour choisi")
choix = int(input())
cells = []

while choix != 999:
    if choix == 1:
        # min_points = int(input("Nombre minimum de points : "))
        # max_points = int(input("Nombre maximum de points : "))
        # cells = lattice.generate_random_lattice(min_points, max_points, 0.01)
        cells = lattice.generate_random_lattice(10, 15, 0.01)

    elif choix == 2:
        latticeChoix = int(input("Entrez le nombre du lattice choisi : "))
        cells = lattice.generate_custom_lattice(latticeChoix)
    lattice.cells = cells
    # lattice.remove_cell(3)
    display_lattice_points_beams()

print("Programme terminé.")
