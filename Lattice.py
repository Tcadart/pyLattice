import copy
import matplotlib.pyplot as plt
from mpl_toolkits.mplot3d import Axes3D
from Cellule import *
from Point import *
from Beam import *


class Lattice:
    def __init__(self, cell_size_x, cell_size_y, cell_size_z, num_cells_x, num_cells_y, num_cells_z):
        self.cell_size_x = cell_size_x
        self.cell_size_y = cell_size_y
        self.cell_size_z = cell_size_z
        self.num_cells_x = num_cells_x
        self.num_cells_y = num_cells_y
        self.num_cells_z = num_cells_z
        self.cells = []

        self.size_x = self.cell_size_x * self.num_cells_x
        self.size_y = self.cell_size_y * self.num_cells_y
        self.size_z = self.cell_size_z * self.num_cells_z

    def generate_random_cells(self, Tmin, Tmax, padding):
        cell = Cellule(self.cell_size_x, self.cell_size_y, self.cell_size_z, 0, 0, 0)
        _, _, _, _ = cell.generate_points(Tmin, Tmax, padding)
        cell.merge_points()
        cell.generate_beams()
        cell.remove_unused_points()
        return cell

    def generate_custom_cells(self, lattice_choix):
        cell = Cellule(self.cell_size_x, self.cell_size_y, self.cell_size_z, 0, 0, 0)
        BCC = cell.Lattice_geometry(lattice_choix, 1.0)
        # cell.generate_beams_from_given_point_list(BCC, 0)
        cell.generate_beams_from_given_point_list(BCC)
        self.cells.append(cell)
        return self.cells

    def generate_random_lattice(self, Tmin, Tmax, padding):
        self.cells = []
        random_cell = self.generate_random_cells(Tmin, Tmax, padding)
        for i in range(self.num_cells_x):
            for j in range(self.num_cells_y):
                for k in range(self.num_cells_z):
                    cell = copy.deepcopy(random_cell)
                    cell.translate((i, j, k))
                    self.cells.append(cell)
        return self.cells

    def generate_random_lattice(self, Tmin, Tmax, padding):
        self.cells = []
        random_cell = self.generate_random_cells(Tmin, Tmax, padding)  # Step 1: Generate a random cell

        for i in range(self.num_cells_x):
            for j in range(self.num_cells_y):
                for k in range(self.num_cells_z):
                    new_cell = copy.deepcopy(random_cell)  # Step 2: Copy the random cell
                    # Step 2: Translate the copied cell
                    new_points = new_cell.points
                    new_cell.points = new_points
                    new_cell.merge_points()
                    new_beams = new_cell.beams
                    new_cell.beams = new_beams
                    new_cell.remove_unused_points()  # Step 3: Remove unused points
                    new_cell.translate((i, j, k))
                    self.cells.append(new_cell)

        return self.cells

    def generate_custom_lattice(self, lattice_choix):
        self.cells = []
        custom_cell = self.generate_custom_cells(lattice_choix)[0]  # Generate a single custom cell
        BCC = custom_cell.Lattice_geometry(lattice_choix, 1.0)
        for i in range(self.num_cells_x):
            for j in range(self.num_cells_y):
                for k in range(self.num_cells_z):
                    new_points, new_beams = custom_cell.generate_beams_from_given_point_list(BCC)
                    new_cell = Cellule(custom_cell.cell_size_x, custom_cell.cell_size_y, custom_cell.cell_size_z,
                                       custom_cell.x, custom_cell.y, custom_cell.z)
                    new_cell.points = new_points
                    new_cell.beams = new_beams
                    new_cell.translate((i, j, k))
                    self.cells.append(new_cell)
        return self.cells

    def remove_cell(self, index):
        if 0 <= index < len(self.cells):
            del self.cells[index]
        else:
            raise IndexError("Invalid cell index.")

    def count_unique_points(self):
        list_points_lattice = []
        for cell in self.cells:
            for point in cell.points:
                list_points_lattice.append(point)
        point_counts = {}
        for point in list_points_lattice:
            if point in point_counts:
                point_counts[point] += 1
            else:
                point_counts[point] = 1
        unique_points = [point for point, count in point_counts.items() if count > 1]
        return unique_points

    def display_all_points(self):
        unique_points = self.count_unique_points()
        for index, point in enumerate(unique_points):
            print(f"[{index}, {point.x}, {point.y}, {point.z}],")

    def display_all_beams(self):
        updated_beams = set()
        unique_points = self.count_unique_points()
        for index, cell in enumerate(self.cells):
            for beam in cell.beams:
                point1_index = self.get_unique_point_index(beam.point1, unique_points)
                point2_index = self.get_unique_point_index(beam.point2, unique_points)
                if point1_index is not None and point2_index is not None:
                    updated_beams.add(tuple(sorted([point1_index, point2_index])))

        for index, beam in enumerate(updated_beams):
            print(f" [{index}, {beam[0]}, {beam[1]}],")

    def get_unique_point_index(self, point, unique_points):
        for index, unique_point in enumerate(unique_points):
            if point.x == unique_point.x and point.y == unique_point.y and point.z == unique_point.z:
                return index
        return None

    def visualize_3d(self, ax):
        self.display_all_points()
        self.display_all_beams()
        for index, cell in enumerate(self.cells):
            cell.display_point(ax, 'r', 'pink', 'black')
            cell.display_beams(ax, 'b', 'r')
            ax.set_xlabel('X')
            ax.set_ylabel('Y')
            ax.set_zlabel('Z')
        ax.set_xlim3d(0, self.size_x)
        ax.set_ylim3d(0, self.size_y)
        ax.set_zlim3d(0, self.size_z)
