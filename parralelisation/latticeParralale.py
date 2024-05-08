import copy
from Codes.Cellule import *
import multiprocessing


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

    def generate_cells_parallel(self, min_points, max_points, step_size, num_processes):
        def generate_cells_single_process(start_x, end_x):
            cells_partial = []
            for i in range(start_x, end_x):
                for j in range(self.num_cells_y):
                    for k in range(self.num_cells_z):
                        cell = Cellule(self.cell_size_x, self.cell_size_y, self.cell_size_z, i, j, k)
                        cell.generate_points(min_points, max_points, step_size)
                        cell.generate_beams()
                        cells_partial.append(cell)
            return cells_partial

        self.cells = []
        process_pool = multiprocessing.Pool(processes=num_processes)

        chunk_size_x = self.num_cells_x // num_processes
        start_indices = [i * chunk_size_x for i in range(num_processes)]
        end_indices = [(i + 1) * chunk_size_x for i in range(num_processes - 1)]
        end_indices.append(self.num_cells_x)

        results = process_pool.starmap(generate_cells_single_process, zip(start_indices, end_indices))

        for cells_partial in results:
            self.cells.extend(cells_partial)

        process_pool.close()
        process_pool.join()

        return self.cells

    def generate_random_cells(self, Tmin, Tmax, padding):
        cell = Cellule(self.cell_size_x, self.cell_size_y, self.cell_size_z, 0, 0, 0)
        _, _, _, _ = cell.generate_points(Tmin, Tmax, padding)
        cell.merge_points()
        cell.generate_beams()
        cell.remove_unused_points()
        return cell

    def generate_custom_cells(self, lattice_choix):
        cell = Cellule(self.cell_size_x, self.cell_size_y, self.cell_size_z, 0, 0, 0)
        BCC = cell.Lattice_geometry(lattice_choix)
        # cell.generate_beams_from_given_point_list(BCC, 0)
        cell.generate_beams_from_given_point_list(BCC)
        self.cells.append(cell)
        return self.cells

    def generate_lattice(self, Tmin, Tmax, padding):
        self.cells = []
        random_cell = self.generate_random_cells(Tmin, Tmax, padding)
        for i in range(self.num_cells_x):
            for j in range(self.num_cells_y):
                for k in range(self.num_cells_z):
                    cell = copy.deepcopy(random_cell)  # Utiliser une copie de la cellule générée aléatoirement
                    cell.translate((i, j, k))
                    self.cells.append(cell)

    def remove_cell(self, index):
        if 0 <= index < len(self.cells):
            del self.cells[index]
        else:
            raise IndexError("Invalid cell index.")
    def visualize_3d(self, ax):
        for index, cell in enumerate(self.cells):
            # cell.affichage_points()
            # cell.display_point(ax)
            # cell.display_beams(ax)
            cell.display_point(ax, 'r', 'pink','black')
            cell.display_beams(ax, 'b', 'r')
            ax.set_xlabel('X')
            ax.set_ylabel('Y')
            ax.set_zlabel('Z')
        ax.set_xlim3d(0, self.size_x)
        ax.set_ylim3d(0, self.size_y)
        ax.set_zlim3d(0, self.size_z)
