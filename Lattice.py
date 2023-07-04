from Cellule import *


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

    def generate_custom_cells(self, lattice_choix):
        cell = Cellule(self.cell_size_x, self.cell_size_y, self.cell_size_z, 0, 0, 0)
        BCC = cell.Lattice_geometry(lattice_choix, 1.0)
        # cell.generate_beams_from_given_point_list(BCC, 0)
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
        thomas = []
        for index, point in enumerate(unique_points):
            thomas.append([index, float(point.x), float(point.y), float(point.z)])
        return thomas

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

    def affichage_beams_console(self):
        unique_beams = self.count_unique_beams()
        thomasB = []
        for index, beam in enumerate(unique_beams):
            # thomasB.append(beam)
            thomasB.append([index, beam[0], beam[1]])
            # print(f"[{index}, {beam[0]}, {beam[1]}],")
        return thomasB

    def get_unique_point_index(self, point, unique_points):
        for index, unique_point in enumerate(unique_points):
            if point.x == unique_point.x and point.y == unique_point.y and point.z == unique_point.z:
                return index
        return None
