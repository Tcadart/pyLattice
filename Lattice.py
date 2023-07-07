import copy

from Cellule import *
import time
import numpy as np

class Lattice:
    def __init__(self, cell_size_x, cell_size_y, cell_size_z, num_cells_x, num_cells_y, num_cells_z):
        self.cell_size_x = cell_size_x
        self.cell_size_y = cell_size_y
        self.cell_size_z = cell_size_z
        self.num_cells_x = num_cells_x
        self.num_cells_y = num_cells_y
        self.num_cells_z = num_cells_z
        self.cells = []
        self._beams = []
        self.size_x = self.cell_size_x * self.num_cells_x
        self.size_y = self.cell_size_y * self.num_cells_y
        self.size_z = self.cell_size_z * self.num_cells_z

    # def generate_random_cells(self, Tmin, Tmax, padding):
    #     cell = Cellule(self.cell_size_x, self.cell_size_y, self.cell_size_z, 0, 0, 0)
    #     cell.generate_points(Tmin, Tmax, padding)
    #     # cell.points = cell.generate_points(Tmin, Tmax, padding)
    #     cell.merge_points()
    #     cell.generate_beams()
    #     # cell.beams = cell.generate_beams()
    #     cell.remove_unused_points()
    #     cell.remove_unused_points()
    #     self.cells.append(cell)
    #     return self.cells
    def generate_random_cells(self, Tmin, Tmax, padding):
        cell = Cellule(self.cell_size_x, self.cell_size_y, self.cell_size_z, 0, 0, 0)
        cell.generate_points(Tmin, Tmax, padding)
        # _, _, _, _ = cell.generate_points(Tmin, Tmax, padding)
        time.sleep(1)
        cell.merge_points()
        cell.generate_beams()
        cell.remove_unused_points()
        return cell

    def generate_custom_cells(self, lattice_choix):
        cell = Cellule(self.cell_size_x, self.cell_size_y, self.cell_size_z, 0, 0, 0)
        BCC = cell.Lattice_geometry(lattice_choix, 1.0)
        cell.generate_beams_from_given_point_list(BCC)
        self.cells.append(cell)
        return self.cells

    @property
    def beams(self):
        return self._beams

    def get_all_beams(self):
        all_beams = []
        for cell in self.cells:
            all_beams.extend(cell.beams)
        return all_beams

    def generate_random_lattice(self, Tmin, Tmax, padding):
        self.cells = []
        random_cell = self.generate_random_cells(Tmin, Tmax, padding)

        for i in range(self.num_cells_x):
            for j in range(self.num_cells_y):
                for k in range(self.num_cells_z):
                    new_cell = copy.deepcopy(random_cell)
                    # Step 2: Translate the copied cell
                    new_points = new_cell.points
                    new_cell.points = new_points
                    new_cell.merge_points()
                    new_beams = new_cell.beams
                    new_cell.beams = new_beams
                    new_cell.remove_unused_points()
                    new_cell.translate((i, j, k))
                    self.cells.append(new_cell)

        return self.cells

    # def generate_random_lattice(self, Tmin, Tmax, padding):
    #     self.cells = []
    #     random_cell = self.generate_random_cells(Tmin, Tmax, padding)
    #     # random_cell = self.generate_random_cells(Tmin, Tmax, padding)[0]
    #
    #     for i in range(self.num_cells_x):
    #         for j in range(self.num_cells_y):
    #             for k in range(self.num_cells_z):
    #                 new_cell = copy.deepcopy(random_cell)
    #                 new_points = new_cell.points
    #                 new_cell.points = new_points
    #                 new_cell.merge_points()
    #                 new_beams = new_cell.beams
    #                 new_cell.beams = new_beams
    #                 # new_cell.remove_unused_points()
    #                 new_cell.translate((i, j, k))
    #                 self.cells.append(new_cell)
    #
    #     return self.cells

    # def generate_random_lattice(self, Tmin, Tmax, padding):
    #     self.cells = []
    #     random_cell = self.generate_random_cells(Tmin, Tmax, padding)[0]
    #     for i in range(self.num_cells_x):
    #         for j in range(self.num_cells_y):
    #             for k in range(self.num_cells_z):
    #                 new_points = random_cell.generate_points(Tmin, Tmax, padding)
    #                 new_beams = random_cell.generate_beams()
    #                 new_cell = Cellule(random_cell.cell_size_x, random_cell.cell_size_y, random_cell.cell_size_z,
    #                                    random_cell.x, random_cell.y, random_cell.z)
    #                 new_cell.points = new_points
    #                 new_cell.beams = new_beams
    #                 new_cell.remove_unused_points()
    #                 new_cell.translate((i, j, k))
    #                 self.cells.append(new_cell)
    #     return self.cells

    def generate_custom_lattice(self, lattice_choix):
        self.cells = []
        custom_cell = self.generate_custom_cells(lattice_choix)[0]  # Generate a single custom cell
        BCC = custom_cell.Lattice_geometry(lattice_choix, 0.15)
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

    def affichage_points_console(self):
        unique_points = self.count_unique_points()
        node_data = []
        for index, point in enumerate(unique_points):
            node_data.append([index, point.x, point.y, point.z])
        return node_data

    # def affichage_points_console(self):
    #     unique_points = self.count_unique_points()
    #     for index, point in enumerate(unique_points):
    #         print(f"[{index}, {point.x}, {point.y}, {point.z}],")

    def count_unique_beams(self):
        updated_beams = set()
        unique_points = self.count_unique_points()
        for index, cell in enumerate(self.cells):
            for beam in cell.beams:
                point1_index = self.get_unique_point_index(beam.point1, unique_points)
                point2_index = self.get_unique_point_index(beam.point2, unique_points)
                if point1_index is not None and point2_index is not None:
                    unique_beam = tuple(sorted([point1_index, point2_index]))
                    updated_beams.add(unique_beam)
        return list(updated_beams)

    # def affichage_beams_console(self):
    #     unique_beams = self.count_unique_beams()
    #     for index, beam in enumerate(unique_beams):
    #         print(f"[{index}, {beam[0]}, {beam[1]}],")

    def affichage_beams_console(self):
        unique_beams = self.count_unique_beams()
        Beam_data = []
        for index, beam in enumerate(unique_beams):
            Beam_data.append([index, beam[0], beam[1]])
        return Beam_data

    def get_unique_point_index(self, point, unique_points):
        for index, unique_point in enumerate(unique_points):
            if point.x == unique_point.x and point.y == unique_point.y and point.z == unique_point.z:
                return index
        return None

    def visualize_3d(self, ax):
        # self.affichage_points_console()
        # print("h")
        # self.affichage_beams_console()
        print("p", self.affichage_points_console())
        print("b", self.affichage_beams_console())
        for index, cell in enumerate(self.cells):
            cell.display_point(ax, 'r', 'pink', 'black')
            cell.display_beams(ax, 'b', 'r')
            ax.set_xlabel('X')
            ax.set_ylabel('Y')
            ax.set_zlabel('Z')
        ax.set_xlim3d(0, self.size_x)
        ax.set_ylim3d(0, self.size_y)
        ax.set_zlim3d(0, self.size_z)


    def Getangle(liaison, Lattice_geom):
        angle = []
        angle_deg = []
        for j in range(len(liaison)):
            u = [Lattice_geom[liaison[0]][3] - Lattice_geom[liaison[0]][0],
                 Lattice_geom[liaison[0]][4] - Lattice_geom[liaison[0]][1],
                 Lattice_geom[liaison[0]][5] - Lattice_geom[liaison[0]][2]]
            v = [Lattice_geom[liaison[j]][3] - Lattice_geom[liaison[j]][0],
                 Lattice_geom[liaison[j]][4] - Lattice_geom[liaison[j]][1],
                 Lattice_geom[liaison[j]][5] - Lattice_geom[liaison[j]][2]]
            if np.dot(u, v) < (np.linalg.norm(u) * np.linalg.norm(v)):
                angle_rad = acos(np.dot(u, v) / (np.linalg.norm(u) * np.linalg.norm(v)))
                angle_deg.append(degrees(angle_rad))
            else:
                angle_deg.append(0)
        angle_deg = np.array(angle_deg)
        if np.all(angle_deg == 0):
            angle.append(0)
        else:
            non_zero_angle = [x for x in angle_deg if x >= 0.01]
            angle.append(min(non_zero_angle))
        angle = np.array(angle)
        return angle
