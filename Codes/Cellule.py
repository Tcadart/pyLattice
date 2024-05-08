from Point import *
from Beam import *
import math
import random


class Cellule:
    def __init__(self, cell_size_x, cell_size_y, cell_size_z, x, y, z):
        """
        Initialize a Cellule (cell) with its dimensions and position.

        :param cell_size_x: x-dimension of the cell
        :param cell_size_y: y-dimension of the cell
        :param cell_size_z: z-dimension of the cell
        :param x: x-coordinate of the cell's position
        :param y: y-coordinate of the cell's position
        :param z: z-coordinate of the cell's position
        """
        self.cell_size_x = cell_size_x
        self.cell_size_y = cell_size_y
        self.cell_size_z = cell_size_z
        self.x = x
        self.y = y
        self.z = z
        self.nodes = []
        self.beams = []
        self._centerCell = self.calculateCenterCell()
    
    @property
    def centerCell(self):
        return self._centerCell

    # def nbPointsTBCI(self, Tmin, Tmax):
    #     """
    #     Generate a random distribution of points within the cell.

    #     :param Tmin: Minimum number of points
    #     :param Tmax: Maximum number of points
    #     :return: Tuple containing total points, boundary points, corner points, and interior points
    #     """
    #     pointtotal = random.randint(Tmin, Tmax)
    #     pointinterieur = random.randint(0, pointtotal)
    #     proportion_coins = random.uniform(0, 1)
    #     pointBordOuCoins = pointtotal - pointinterieur
    #     pointcoin = int(pointBordOuCoins * proportion_coins)
    #     pointbord = pointBordOuCoins - pointcoin
    #     if pointcoin > 8:
    #         pointbord += pointcoin - 8
    #         pointcoin = 8
    #     elif pointcoin < 0:
    #         pointcoin = 0
    #         pointbord += pointcoin
    #     return pointtotal, pointbord, pointcoin, pointinterieur

    def translate(self, translation):
        """
        Translate the cell and its points by a given translation vector.

        :param translation: Tuple (dx, dy, dz) specifying the translation along x, y, and z axes
        """
        self.x += translation[0]
        self.y += translation[1]
        self.z += translation[2]

        for point in self.nodes:
            point.x += translation[0]
            point.y += translation[1]
            point.z += translation[2]

    # def generate_nodes(self, Tmin, Tmax, padding):
    #     self.nodes = []
    #     pointtotal, pointbord, pointcoin, pointinterieur = self.nbPointsTBCI(Tmin, Tmax)
    #     nodes_interior = []
    #     nodes_border = []
    #     nodes_corner = []
    #     for _ in range(pointinterieur):
    #         x = round(random.uniform(0 + padding, self.cellSizeX - padding), 2)
    #         y = round(random.uniform(0 + padding, self.cellSizeY - padding), 2)
    #         z = round(random.uniform(0 + padding, self.cellSizeZ - padding), 2)
    #         point = Point(x, y, z)
    #         nodes_interior.append(point)
    #     for _ in range(pointcoin):
    #         x = random.choice([0, self.cellSizeX])
    #         y = random.choice([0, self.cellSizeY])
    #         z = random.choice([0, self.cellSizeZ])
    #         point = Point(x, y, z)
    #         nodes_corner.append(point)
    #     for _ in range(pointbord):
    #         if random.random() < 0.5:
    #             x = random.choice([0, self.cellSizeX])
    #             y = random.uniform(0, self.cellSizeY)
    #             z = random.uniform(0, self.cellSizeZ)
    #         else:
    #             x = random.uniform(0, self.cellSizeX)
    #             y = random.choice([0, self.cellSizeY])
    #             z = random.uniform(0, self.cellSizeZ)
    #         point = Point(round(x, 2), round(y, 2), round(z, 2))
    #         nodes_border.append(point)
    #     self.nodes = nodes_interior + nodes_border + nodes_corner
    #     # return self.nodes, nodes_border, nodes_corner, nodes_interior
    #     return self.nodes

    # def generate_nodes(self, Tmin, Tmax, padding):
    #     """
    #     Generate nodes within the cell.

    #     :param Tmin: Minimum number of points
    #     :param Tmax: Maximum number of points
    #     :param padding: Minimum distance from the cell boundary
    #     :return: List of generated Node objects
    #     """
    #     pointtotal, pointbord, pointcoin, pointinterieur = self.nbPointsTBCI(Tmin, Tmax)
    #     nodes_interior = []
    #     nodes_border = []
    #     nodes_corner = []
    #     min_distance = 0.02

    #     def check_min_distance(x, y, z, existing_nodes):
    #         for point in existing_nodes:
    #             distance = ((point.x - x) ** 2 + (point.y - y) ** 2 + (point.z - z) ** 2) ** 0.5
    #             if distance < min_distance:
    #                 return False
    #         return True

    #     for _ in range(pointinterieur):
    #         x, y, z = None, None, None
    #         while True:
    #             x = round(random.uniform(0 + padding, self.cellSizeX - padding), 2)
    #             y = round(random.uniform(0 + padding, self.cellSizeY - padding), 2)
    #             z = round(random.uniform(0 + padding, self.cellSizeZ - padding), 2)
    #             if check_min_distance(x, y, z, nodes_interior):
    #                 break
    #         point = Point(x, y, z)
    #         nodes_interior.append(point)

    #     for _ in range(pointcoin):
    #         x, y, z = None, None, None
    #         while True:
    #             x = random.choice([0, self.cellSizeX])
    #             y = random.choice([0, self.cellSizeY])
    #             z = random.choice([0, self.cellSizeZ])
    #             if check_min_distance(x, y, z, nodes_corner):
    #                 break
    #         point = Point(x, y, z)
    #         nodes_corner.append(point)

    #     for _ in range(pointbord):
    #         x, y, z = None, None, None
    #         while True:
    #             if random.random() < 0.5:
    #                 x = random.choice([0, self.cellSizeX])
    #                 y = random.uniform(0, self.cellSizeY)
    #                 z = random.uniform(0, self.cellSizeZ)
    #             else:
    #                 x = random.uniform(0, self.cellSizeX)
    #                 y = random.choice([0, self.cellSizeY])
    #                 z = random.uniform(0, self.cellSizeZ)
    #             if check_min_distance(x, y, z, nodes_border):
    #                 break
    #         point = Point(round(x, 2), round(y, 2), round(z, 2))
    #         nodes_border.append(point)

    #     self.nodes = nodes_interior + nodes_border + nodes_corner
    #     return self.nodes

    def merge_nodes(self):
        """
        Merge nodes on the same plane to simplify the structure.
        """
        z_values = [point.z for point in self.nodes]
        z_min = min(z_values)
        z_max = max(z_values)
        merged_nodes = set()
        for point in self.nodes:
            x, y, z = point.x, point.y, point.z
            if z == z_min or (x in [0, self.cell_size_x] and y in [0, self.cell_size_y]) or (
                    z == z_max and x != 0 and x != self.cell_size_x and y != 0 and y != self.cell_size_y):
                merged_nodes.add(Point(x, y, z_max))
            elif z == z_max or (x in [0, self.cell_size_x] and y in [0, self.cell_size_y]) or (
                    z == z_min and x != 0 and x != self.cell_size_x and y != 0 and y != self.cell_size_y):
                merged_nodes.add(Point(x, y, z_min))
            else:
                merged_nodes.add(point)
        self.nodes = list(merged_nodes)

    # def generate_beams(self):
    #     """
    #     Generate beams between nodes.

    #     :return: List of generated Beam objects
    #     """
    #     num_nodes = len(self.nodes)
    #     min_beams = num_nodes // 2
    #     max_beams = num_nodes
    #     num_beams = random.randint(min_beams, max_beams)
    #     used_pairs = set()
    #     used_beams = set()
    #     while num_beams > 0:
    #         index1 = random.randint(0, num_nodes - 1)
    #         index2 = random.randint(0, num_nodes - 1)
    #         if index1 != index2 and (index1, index2) not in used_pairs and (index2, index1) not in used_pairs:
    #             point1 = self.nodes[index1]
    #             point2 = self.nodes[index2]
    #             beam = Beam(point1, point2)
    #             if beam not in used_beams:
    #                 self.beams.append(beam)
    #                 used_pairs.add((index1, index2))
    #                 used_beams.add(beam)
    #                 num_beams -= 1
    #     return self.beams

    def remove_unused_nodes(self):
        """
        Remove nodes not used in any beams.

        :return: List of remaining Node objects
        """
        used_nodes = set()
        for beam in self.beams:
            used_nodes.add(beam.point1)
            used_nodes.add(beam.point2)
        self.nodes = list(used_nodes)
        return self.nodes

    def Lattice_geometry(self, Lattice):
        BCC = [(0.0, 0.0, 0.0, 0.5, 0.5, 0.5),
               (0.5, 0.5, 0.5, 1.0, 1.0, 1.0),
               (0.5, 0.5, 0.5, 1.0, 1.0, 0.0),
               (0.5, 0.5, 0.5, 0.0, 0.0, 1.0),
               (0.5, 0.5, 0.5, 0.0, 1.0, 0.0),
               (0.5, 0.5, 0.5, 0.0, 1.0, 1.0),
               (1.0, 0.0, 1.0, 0.5, 0.5, 0.5),
               (0.5, 0.5, 0.5, 1.0, 0.0, 0.0)]
        Octet = [(0.0, 0.0, 0.0, 0.5, 0.0, 0.5),
                 (1.0, 0.0, 1.0, 0.5, 0.0, 0.5),
                 (0.0, 0.0, 1.0, 0.5, 0.0, 0.5),
                 (1.0, 0.0, 0.0, 0.5, 0.0, 0.5),
                 (0.0, 0.0, 0.0, 0.0, 0.5, 0.5),
                 (0.0, 1.0, 1.0, 0.0, 0.5, 0.5),
                 (0.0, 0.0, 1.0, 0.0, 0.5, 0.5),
                 (0.0, 1.0, 0.0, 0.0, 0.5, 0.5),
                 (0.0, 0.0, 0.0, 0.5, 0.5, 0.0),
                 (1.0, 1.0, 0.0, 0.5, 0.5, 0.0),
                 (1.0, 0.0, 0.0, 0.5, 0.5, 0.0),
                 (0.0, 1.0, 0.0, 0.5, 0.5, 0.0),
                 (0.0, 0.0, 1.0, 0.5, 0.5, 1.0),
                 (1.0, 1.0, 1.0, 0.5, 0.5, 1.0),
                 (1.0, 0.0, 1.0, 0.5, 0.5, 1.0),
                 (0.0, 1.0, 1.0, 0.5, 0.5, 1.0),
                 (1.0, 0.5, 0.5, 1.0, 1.0, 1.0),
                 (1.0, 0.0, 0.0, 1.0, 0.5, 0.5),
                 (1.0, 0.5, 0.5, 1.0, 1.0, 0.0),
                 (1.0, 0.0, 1.0, 1.0, 0.5, 0.5),
                 (0.5, 1.0, 0.5, 1.0, 1.0, 1.0),
                 (0.0, 1.0, 0.0, 0.5, 1.0, 0.5),
                 (0.5, 1.0, 0.5, 1.0, 1.0, 0.0),
                 (0.0, 1.0, 1.0, 0.5, 1.0, 0.5),
                 (0.5, 0.0, 0.5, 0.5, 0.5, 0.0),
                 (0.5, 0.0, 0.5, 0.0, 0.5, 0.5),
                 (0.5, 0.0, 0.5, 1.0, 0.5, 0.5),
                 (0.5, 0.0, 0.5, 0.5, 0.5, 1.0),
                 (0.5, 1.0, 0.5, 0.5, 0.5, 0.0),
                 (0.5, 1.0, 0.5, 0.0, 0.5, 0.5),
                 (0.5, 1.0, 0.5, 1.0, 0.5, 0.5),
                 (0.5, 1.0, 0.5, 0.5, 0.5, 1.0),
                 (1.0, 0.5, 0.5, 0.5, 0.5, 0.0),
                 (0.5, 0.5, 0.0, 0.0, 0.5, 0.5),
                 (0.5, 0.5, 1.0, 0.0, 0.5, 0.5),
                 (0.5, 0.5, 1.0, 1.0, 0.5, 0.5)]
        OctetExt = [(0.0, 0.0, 0.0, 0.5, 0.0, 0.5),
                    (1.0, 0.0, 1.0, 0.5, 0.0, 0.5),
                    (0.0, 0.0, 1.0, 0.5, 0.0, 0.5),
                    (1.0, 0.0, 0.0, 0.5, 0.0, 0.5),
                    (0.0, 0.0, 0.0, 0.0, 0.5, 0.5),
                    (0.0, 1.0, 1.0, 0.0, 0.5, 0.5),
                    (0.0, 0.0, 1.0, 0.0, 0.5, 0.5),
                    (0.0, 1.0, 0.0, 0.0, 0.5, 0.5),
                    (0.0, 0.0, 0.0, 0.5, 0.5, 0.0),
                    (1.0, 1.0, 0.0, 0.5, 0.5, 0.0),
                    (1.0, 0.0, 0.0, 0.5, 0.5, 0.0),
                    (0.0, 1.0, 0.0, 0.5, 0.5, 0.0),
                    (0.0, 0.0, 1.0, 0.5, 0.5, 1.0),
                    (1.0, 1.0, 1.0, 0.5, 0.5, 1.0),
                    (1.0, 0.0, 1.0, 0.5, 0.5, 1.0),
                    (0.0, 1.0, 1.0, 0.5, 0.5, 1.0),
                    (1.0, 0.5, 0.5, 1.0, 1.0, 1.0),
                    (1.0, 0.0, 0.0, 1.0, 0.5, 0.5),
                    (1.0, 0.5, 0.5, 1.0, 1.0, 0.0),
                    (1.0, 0.0, 1.0, 1.0, 0.5, 0.5),
                    (0.5, 1.0, 0.5, 1.0, 1.0, 1.0),
                    (0.0, 1.0, 0.0, 0.5, 1.0, 0.5),
                    (0.5, 1.0, 0.5, 1.0, 1.0, 0.0),
                    (0.0, 1.0, 1.0, 0.5, 1.0, 0.5)]
        OctetInt = [(0.5, 0.0, 0.5, 0.5, 0.5, 0.0),
                    (0.5, 0.0, 0.5, 0.0, 0.5, 0.5),
                    (0.5, 0.0, 0.5, 1.0, 0.5, 0.5),
                    (0.5, 0.0, 0.5, 0.5, 0.5, 1.0),
                    (0.5, 1.0, 0.5, 0.5, 0.5, 0.0),
                    (0.5, 1.0, 0.5, 0.0, 0.5, 0.5),
                    (0.5, 1.0, 0.5, 1.0, 0.5, 0.5),
                    (0.5, 1.0, 0.5, 0.5, 0.5, 1.0),
                    (1.0, 0.5, 0.5, 0.5, 0.5, 0.0),
                    (0.5, 0.5, 0.0, 0.0, 0.5, 0.5),
                    (0.5, 0.5, 1.0, 0.0, 0.5, 0.5),
                    (0.5, 0.5, 1.0, 1.0, 0.5, 0.5)]
        BCCZ = [(0.5, 0.5, 0.5, 1.0, 1.0, 1.0),
                (0.0, 0.0, 0.0, 0.5, 0.5, 0.5),
                (0.5, 0.5, 0.5, 1.0, 1.0, 0.0),
                (0.0, 0.0, 1.0, 0.5, 0.5, 0.5),
                (0.5, 0.5, 0.5, 0.0, 1.0, 0.0),
                (1.0, 0.0, 1.0, 0.5, 0.5, 0.5),
                (0.5, 0.5, 0.5, 0.0, 1.0, 1.0),
                (1.0, 0.0, 0.0, 0.5, 0.5, 0.5),
                (0.5, 0.5, 0.0, 0.5, 0.5, 0.5),
                (0.5, 0.5, 0.5, 0.5, 0.5, 1.0)]
        Cubic = [(0.0, 0.0, 0.0, 0.0, 0.0, 1.0),
                 (1.0, 0.0, 0.0, 1.0, 0.0, 1.0),
                 (0.0, 1.0, 0.0, 0.0, 1.0, 1.0),
                 (1.0, 1.0, 0.0, 1.0, 1.0, 1.0),
                 (0.0, 0.0, 0.0, 1.0, 0.0, 0.0),
                 (0.0, 0.0, 0.0, 0.0, 1.0, 0.0),
                 (1.0, 1.0, 0.0, 0.0, 1.0, 0.0),
                 (1.0, 1.0, 0.0, 1.0, 0.0, 0.0),
                 (0.0, 0.0, 1.0, 1.0, 0.0, 1.0),
                 (0.0, 0.0, 1.0, 0.0, 1.0, 1.0),
                 (1.0, 1.0, 1.0, 0.0, 1.0, 1.0),
                 (1.0, 1.0, 1.0, 1.0, 0.0, 1.0)]
        OctahedronZ = [(0.5, 0.0, 0.5, 0.5, 0.5, 0.0),
                       (0.5, 0.0, 0.5, 0.0, 0.5, 0.5),
                       (0.5, 0.0, 0.5, 1.0, 0.5, 0.5),
                       (0.5, 0.0, 0.5, 0.5, 0.5, 1.0),
                       (0.5, 1.0, 0.5, 0.5, 0.5, 0.0),
                       (0.5, 1.0, 0.5, 0.0, 0.5, 0.5),
                       (0.5, 1.0, 0.5, 1.0, 0.5, 0.5),
                       (0.5, 1.0, 0.5, 0.5, 0.5, 1.0),
                       (1.0, 0.5, 0.5, 0.5, 0.5, 0.0),
                       (0.5, 0.5, 0.0, 0.0, 0.5, 0.5),
                       (0.5, 0.5, 1.0, 0.0, 0.5, 0.5),
                       (0.5, 0.5, 1.0, 1.0, 0.5, 0.5),
                       (0.5, 0.5, 0.0, 0.5, 0.5, 1.0)]
        OctahedronZcross = [(0.5, 0.0, 0.5, 0.5, 0.5, 0.0),
                            (0.5, 0.0, 0.5, 0.0, 0.5, 0.5),
                            (0.5, 0.0, 0.5, 1.0, 0.5, 0.5),
                            (0.5, 0.0, 0.5, 0.5, 0.5, 1.0),
                            (0.5, 1.0, 0.5, 0.5, 0.5, 0.0),
                            (0.5, 1.0, 0.5, 0.0, 0.5, 0.5),
                            (0.5, 1.0, 0.5, 1.0, 0.5, 0.5),
                            (0.5, 1.0, 0.5, 0.5, 0.5, 1.0),
                            (1.0, 0.5, 0.5, 0.5, 0.5, 0.0),
                            (0.5, 0.5, 0.0, 0.0, 0.5, 0.5),
                            (0.5, 0.5, 1.0, 0.0, 0.5, 0.5),
                            (0.5, 0.5, 1.0, 1.0, 0.5, 0.5),
                            (0.5, 0.5, 0.0, 0.5, 0.5, 0.5),
                            (0.5, 0.5, 1.0, 0.5, 0.5, 0.5),
                            (0.5, 0.0, 0.5, 0.5, 0.5, 0.5),
                            (0.0, 0.5, 0.5, 0.5, 0.5, 0.5),
                            (1.0, 0.5, 0.5, 0.5, 0.5, 0.5),
                            (0.5, 1.0, 0.5, 0.5, 0.5, 0.5)]
        Kelvin = [(0.5, 0.25, 0, 0.25, 0.5, 0),
                  (0.5, 0.25, 0, 0.75, 0.5, 0),
                  (0.5, 0.75, 0, 0.25, 0.5, 0),
                  (0.5, 0.75, 0, 0.75, 0.5, 0),
                  (0.5, 0.25, 1, 0.25, 0.5, 1),
                  (0.5, 0.25, 1, 0.75, 0.5, 1),
                  (0.5, 0.75, 1, 0.25, 0.5, 1),
                  (0.5, 0.75, 1, 0.75, 0.5, 1),
                  (0.5, 0, 0.25, 0.25, 0, 0.5),
                  (0.5, 0, 0.25, 0.75, 0, 0.5),
                  (0.5, 0, 0.75, 0.25, 0, 0.5),
                  (0.5, 0, 0.75, 0.75, 0, 0.5),
                  (0.5, 1, 0.25, 0.25, 1, 0.5),
                  (0.5, 1, 0.25, 0.75, 1, 0.5),
                  (0.5, 1, 0.75, 0.25, 1, 0.5),
                  (0.5, 1, 0.75, 0.75, 1, 0.5),
                  (0, 0.5, 0.25, 0, 0.25, 0.5),
                  (0, 0.5, 0.25, 0, 0.75, 0.5),
                  (0, 0.5, 0.75, 0, 0.25, 0.5),
                  (0, 0.5, 0.75, 0, 0.75, 0.5),
                  (1, 0.5, 0.25, 1, 0.25, 0.5),
                  (1, 0.5, 0.25, 1, 0.75, 0.5),
                  (1, 0.5, 0.75, 1, 0.25, 0.5),
                  (1, 0.5, 0.75, 1, 0.75, 0.5),
                  (0.5, 0.25, 0, 0.5, 0, 0.25),
                  (0.25, 0.5, 0, 0, 0.5, 0.25),
                  (0.75, 0.5, 0, 1, 0.5, 0.25),
                  (0.5, 0.75, 0, 0.5, 1, 0.25),
                  (0.25, 0, 0.5, 0, 0.25, 0.5),
                  (0.75, 0, 0.5, 1, 0.25, 0.5),
                  (0.75, 1, 0.5, 1, 0.75, 0.5),
                  (0.25, 1, 0.5, 0, 0.75, 0.5),
                  (0.5, 0, 0.75, 0.5, 0.25, 1),
                  (0, 0.5, 0.75, 0.25, 0.5, 1),
                  (0.5, 1, 0.75, 0.5, 0.75, 1),
                  (1, 0.5, 0.75, 0.75, 0.5, 1)]
        CubicV2 = [(0.5, 0.0, 0.5, 0.5, 0.5, 0.5),
                   (0.0, 0.5, 0.5, 0.5, 0.5, 0.5),
                   (0.5, 1.0, 0.5, 0.5, 0.5, 0.5),
                   (1.0, 0.5, 0.5, 0.5, 0.5, 0.5),
                   (0.5, 0.5, 0.0, 0.5, 0.5, 0.5),
                   (0.5, 0.5, 1.0, 0.5, 0.5, 0.5)]
        CubicV3 = [(0.0, 0.0, 0.0, 0.5, 0.0, 0.0),
                   (0.5, 0.0, 0.0, 1.0, 0.0, 0.0),
                   (0.0, 1.0, 0.0, 0.5, 1.0, 0.0),
                   (0.5, 1.0, 0.0, 1.0, 1.0, 0.0),
                   (0.5, 0.0, 0.0, 0.5, 1.0, 0.0),
                   (0.0, 0.0, 1.0, 0.5, 0.0, 1.0),
                   (0.5, 0.0, 1.0, 1.0, 0.0, 1.0),
                   (0.0, 1.0, 1.0, 0.5, 1.0, 1.0),
                   (0.5, 1.0, 1.0, 1.0, 1.0, 1.0),
                   (0.5, 0.0, 1.0, 0.5, 1.0, 1.0),
                   (0.5, 0.0, 0.0, 0.5, 0.0, 1.0),
                   (0.5, 1.0, 0.0, 0.5, 1.0, 1.0)]
        CubicV4 = [(0.5, 0.0, 0.0, 0.5, 0.5, 0.0),
                   (0.0, 0.5, 0.0, 0.5, 0.5, 0.0),
                   (0.5, 1.0, 0.0, 0.5, 0.5, 0.0),
                   (1.0, 0.5, 0.0, 0.5, 0.5, 0.0),
                   (0.5, 0.0, 1.0, 0.5, 0.5, 1.0),
                   (0.0, 0.5, 1.0, 0.5, 0.5, 1.0),
                   (0.5, 1.0, 1.0, 0.5, 0.5, 1.0),
                   (1.0, 0.5, 1.0, 0.5, 0.5, 1.0),
                   (0.5, 0.5, 0.0, 0.5, 0.5, 1.0)]
        Newlattice = [
            (0.0, 0.0, 0.0, 0.25, 0.25, 0.25),
            (0.25, 0.25, 0.25, 0.5, 0.0, 0.0),
            (0.25, 0.25, 0.25, 0.0, 0.0, 0.5),
            (0.25, 0.25, 0.25, 0.0, 0.5, 0.0),
            (1.0, 0.0, 0.0, 0.75, 0.25, 0.25),
            (0.75, 0.25, 0.25, 0.5, 0.0, 0.0),
            (0.75, 0.25, 0.25, 1.0, 0.5, 0.0),
            (0.75, 0.25, 0.25, 1.0, 0.0, 0.5),
            (0.0, 1.0, 0.0, 0.25, 0.75, 0.25),
            (0.25, 0.75, 0.25, 0.0, 0.5, 0.0),
            (0.25, 0.75, 0.25, 0.5, 1.0, 0.0),
            (0.25, 0.75, 0.25, 0.0, 1.0, 0.5),
            (1.0, 1.0, 0.0, 0.75, 0.75, 0.25),
            (0.75, 0.75, 0.25, 1.0, 0.5, 0.0),
            (0.75, 0.75, 0.25, 0.5, 1.0, 0.0),
            (0.75, 0.75, 0.25, 1.0, 1.0, 0.5),
            (0.0, 0.0, 1.0, 0.25, 0.25, 0.75),
            (0.25, 0.25, 0.75, 0.0, 0.5, 1.0),
            (0.25, 0.25, 0.75, 0.5, 0.0, 1.0),
            (0.25, 0.25, 0.75, 0.0, 0.0, 0.5),
            (1.0, 0.0, 1.0, 0.75, 0.25, 0.75),
            (0.75, 0.25, 0.75, 1.0, 0.0, 0.5),
            (0.75, 0.25, 0.75, 0.5, 0.0, 1.0),
            (0.75, 0.25, 0.75, 1.0, 0.5, 1.0),
            (0.0, 1.0, 1.0, 0.25, 0.75, 0.75),
            (0.25, 0.75, 0.75, 0.0, 1.0, 0.5),
            (0.25, 0.75, 0.75, 0.5, 1.0, 1.0),
            (0.25, 0.75, 0.75, 0.0, 0.5, 1.0),
            (1.0, 1.0, 1.0, 0.75, 0.75, 0.75),
            (0.75, 0.75, 0.75, 1.0, 1.0, 0.5),
            (0.75, 0.75, 0.75, 0.5, 1.0, 1.0),
            (0.75, 0.75, 0.75, 1.0, 0.5, 1.0)
        ]
        Diamond = [
            (0.0, 0.0, 0.0, 0.25, 0.25, 0.25),
            (0.25, 0.25, 0.25, 0.5, 0.5, 0.0),
            (0.25, 0.25, 0.25, 0.0, 0.5, 0.5),
            (0.25, 0.25, 0.25, 0.5, 0.0, 0.5),
            (1.0, 0.0, 0.0, 0.75, 0.25, 0.25),
            (0.75, 0.25, 0.25, 0.5, 0.5, 0.0),
            (0.75, 0.25, 0.25, 1.0, 0.5, 0.5),
            (0.75, 0.25, 0.25, 0.5, 0.0, 0.5),
            (1.0, 1.0, 0.0, 0.75, 0.75, 0.25),
            (0.75, 0.75, 0.25, 0.5, 0.5, 0.0),
            (0.75, 0.75, 0.25, 1.0, 0.5, 0.5),
            (0.75, 0.75, 0.25, 0.5, 1.0, 0.5),
            (0.0, 1.0, 0.0, 0.25, 0.75, 0.25),
            (0.25, 0.75, 0.25, 0.5, 0.5, 0.0),
            (0.25, 0.75, 0.25, 0.0, 0.5, 0.5),
            (0.25, 0.75, 0.25, 0.5, 1.0, 0.5),
            (0.0, 0.0, 1.0, 0.25, 0.25, 0.75),
            (0.25, 0.25, 0.75, 0.5, 0.5, 1.0),
            (0.25, 0.25, 0.75, 0.0, 0.5, 0.5),
            (0.25, 0.25, 0.75, 0.5, 0.0, 0.5),
            (1.0, 0.0, 1.0, 0.75, 0.25, 0.75),
            (0.75, 0.25, 0.75, 0.5, 0.5, 1.0),
            (0.75, 0.25, 0.75, 1.0, 0.5, 0.5),
            (0.75, 0.25, 0.75, 0.5, 0.0, 0.5),
            (1.0, 1.0, 1.0, 0.75, 0.75, 0.75),
            (0.75, 0.75, 0.75, 0.5, 0.5, 1.0),
            (0.75, 0.75, 0.75, 1.0, 0.5, 0.5),
            (0.75, 0.75, 0.75, 0.5, 1.0, 0.5),
            (0.0, 1.0, 1.0, 0.25, 0.75, 0.75),
            (0.25, 0.75, 0.75, 0.5, 0.5, 1.0),
            (0.25, 0.75, 0.75, 0.0, 0.5, 0.5),
            (0.25, 0.75, 0.75, 0.5, 1.0, 0.5)
        ]
        angleGeom = 20 #Angle en degres
        hGeom = 0.35
        valGeom = hGeom-math.tan(angleGeom*math.pi/180)/2
        Auxetic = [(0.5, 0.0, 0.0, 0.5, 0.0, hGeom),
                   (0.5, 0.0, 1.0, 0.5, 0.0, 1-hGeom),
                   (0.0, 0.0, valGeom, 0.0, 0.0, 1-valGeom),
                   (1.0, 0.0, valGeom, 1.0, 0.0, 1-valGeom),
                   (0.0, 0.0, valGeom, 0.5, 0.0, hGeom),
                   (0.0, 0.0, 1-valGeom, 0.5, 0.0, 1-hGeom),
                   (1.0, 0.0, 1-valGeom, 0.5, 0.0, 1-hGeom),
                   (1.0, 0.0, valGeom, 0.5, 0.0, hGeom),
                   (0.5, 1.0, 0.0, 0.5, 1.0, hGeom),
                   (0.5, 1.0, 1.0, 0.5, 1.0, 1-hGeom),
                   (0.0, 1.0, valGeom, 0.0, 1.0, 1-valGeom),
                   (1.0, 1.0, valGeom, 1.0, 1.0, 1-valGeom),
                   (0.0, 1.0, valGeom, 0.5, 1.0, hGeom),
                   (0.0, 1.0, 1-valGeom, 0.5, 1.0, 1-hGeom),
                   (1.0, 1.0, 1-valGeom, 0.5, 1.0, 1-hGeom),
                   (1.0, 1.0, valGeom, 0.5, 1.0, hGeom),
                   (1.0, 0.0, valGeom, 1.0, 0.5, hGeom),
                   (1.0, 1.0, valGeom, 1.0, 0.5, hGeom),
                   (1.0, 0.5, 0.0, 1.0, 0.5, hGeom),
                   (1.0, 0.5, 1-hGeom, 1.0, 1.0, 1-valGeom),
                   (1.0, 0.5, 1-hGeom, 1.0, 0.0, 1-valGeom),
                   (1.0, 0.5, 1-hGeom, 1.0, 0.5, 1.0),
                   (0.0, 0.0, valGeom, 0.0, 0.5, hGeom),
                   (0.0, 1.0, valGeom, 0.0, 0.5, hGeom),
                   (0.0, 0.5, 0.0, 0.0, 0.5, hGeom),
                   (0.0, 0.5, 1-hGeom, 0.0, 1.0, 1-valGeom),
                   (0.0, 0.5, 1-hGeom, 0.0, 0.0, 1-valGeom),
                   (0.0, 0.5, 1-hGeom, 0.0, 0.5, 1.0)]
        Hichem = [(0.0, 0.0, 0.0, 0.5, 0.5, 0.5),
                   (0.5, 0.5, 0.5, 1.0, 1.0, 1.0),
                   (0.5, 0.5, 0.5, 1.0, 1.0, 0.0),
                   (0.5, 0.5, 0.5, 0.0, 0.0, 1.0),
                   (0.5, 0.5, 0.5, 0.0, 1.0, 0.0),
                   (0.5, 0.5, 0.5, 0.0, 1.0, 1.0),
                   (1.0, 0.0, 1.0, 0.5, 0.5, 0.5),
                   (0.5, 0.5, 0.5, 1.0, 0.0, 0.0),
                  (0.0, 0.0, 0.0, 0.5, 0.0, 0.5),
                  (0.5, 0.0, 0.0, 0.5, 0.0, 0.5),
                  (1.0, 0.0, 0.0, 0.5, 0.0, 0.5),
                  (1.0, 0.0, 0.5, 0.5, 0.0, 0.5),
                  (1.0, 0.0, 1.0, 0.5, 0.0, 0.5),
                  (0.5, 0.0, 1.0, 0.5, 0.0, 0.5),
                  (0.0, 0.0, 1.0, 0.5, 0.0, 0.5),
                  (0.0, 0.0, 0.5, 0.5, 0.0, 0.5),
                  (0.0, 1.0, 0.0, 0.5, 1.0, 0.5),
                  (0.5, 1.0, 0.0, 0.5, 1.0, 0.5),
                  (1.0, 1.0, 0.0, 0.5, 1.0, 0.5),
                  (1.0, 1.0, 0.5, 0.5, 1.0, 0.5),
                  (1.0, 1.0, 1.0, 0.5, 1.0, 0.5),
                  (0.5, 1.0, 1.0, 0.5, 1.0, 0.5),
                  (0.0, 1.0, 1.0, 0.5, 1.0, 0.5),
                  (0.0, 1.0, 0.5, 0.5, 1.0, 0.5),
                  (1.0, 0.0, 0.0, 1.0, 0.5, 0.5),
                  (1.0, 0.5, 0.0, 1.0, 0.5, 0.5),
                  (1.0, 1.0, 0.0, 1.0, 0.5, 0.5),
                  (1.0, 1.0, 0.5, 1.0, 0.5, 0.5),
                  (1.0, 1.0, 1.0, 1.0, 0.5, 0.5),
                  (1.0, 0.5, 1.0, 1.0, 0.5, 0.5),
                  (1.0, 0.0, 1.0, 1.0, 0.5, 0.5),
                  (1.0, 0.0, 0.5, 1.0, 0.5, 0.5),
                  (0.0, 0.0, 0.0, 0.0, 0.5, 0.5),
                  (0.0, 0.5, 0.0, 0.0, 0.5, 0.5),
                  (0.0, 1.0, 0.0, 0.0, 0.5, 0.5),
                  (0.0, 1.0, 0.5, 0.0, 0.5, 0.5),
                  (0.0, 1.0, 1.0, 0.0, 0.5, 0.5),
                  (0.0, 0.5, 1.0, 0.0, 0.5, 0.5),
                  (0.0, 0.0, 1.0, 0.0, 0.5, 0.5),
                  (0.0, 0.0, 0.5, 0.0, 0.5, 0.5),
                  (0.0, 0.0, 0.0, 0.5, 0.5, 0.0),
                  (0.5, 0.0, 0.0, 0.5, 0.5, 0.0),
                  (1.0, 0.0, 0.0, 0.5, 0.5, 0.0),
                  (1.0, 0.5, 0.0, 0.5, 0.5, 0.0),
                  (1.0, 1.0, 0.0, 0.5, 0.5, 0.0),
                  (0.5, 1.0, 0.0, 0.5, 0.5, 0.0),
                  (0.0, 1.0, 0.0, 0.5, 0.5, 0.0),
                  (0.0, 0.5, 0.0, 0.5, 0.5, 0.0),
                  (0.0, 0.0, 1.0, 0.5, 0.5, 1.0),
                  (0.5, 0.0, 1.0, 0.5, 0.5, 1.0),
                  (1.0, 0.0, 1.0, 0.5, 0.5, 1.0),
                  (1.0, 0.5, 1.0, 0.5, 0.5, 1.0),
                  (1.0, 1.0, 1.0, 0.5, 0.5, 1.0),
                  (0.5, 1.0, 1.0, 0.5, 0.5, 1.0),
                  (0.0, 1.0, 1.0, 0.5, 0.5, 1.0),
                  (0.0, 0.5, 1.0, 0.5, 0.5, 1.0)
                  ]
        Hybrid1 = [
            (0.25, 0.25, 0.25, 0.5, 0.0, 0.0),
            (0.25, 0.25, 0.25, 0.0, 0.0, 0.5),
            (0.25, 0.25, 0.25, 0.0, 0.5, 0.0),
            (0.75, 0.25, 0.25, 0.5, 0.0, 0.0),
            (0.75, 0.25, 0.25, 1.0, 0.5, 0.0),
            (0.75, 0.25, 0.25, 1.0, 0.0, 0.5),
            (0.25, 0.75, 0.25, 0.0, 0.5, 0.0),
            (0.25, 0.75, 0.25, 0.5, 1.0, 0.0),
            (0.25, 0.75, 0.25, 0.0, 1.0, 0.5),
            (0.75, 0.75, 0.25, 1.0, 0.5, 0.0),
            (0.75, 0.75, 0.25, 0.5, 1.0, 0.0),
            (0.75, 0.75, 0.25, 1.0, 1.0, 0.5),
            (0.25, 0.25, 0.75, 0.0, 0.5, 1.0),
            (0.25, 0.25, 0.75, 0.5, 0.0, 1.0),
            (0.25, 0.25, 0.75, 0.0, 0.0, 0.5),
            (0.75, 0.25, 0.75, 1.0, 0.0, 0.5),
            (0.75, 0.25, 0.75, 0.5, 0.0, 1.0),
            (0.75, 0.25, 0.75, 1.0, 0.5, 1.0),
            (0.25, 0.75, 0.75, 0.0, 1.0, 0.5),
            (0.25, 0.75, 0.75, 0.5, 1.0, 1.0),
            (0.25, 0.75, 0.75, 0.0, 0.5, 1.0),
            (0.75, 0.75, 0.75, 1.0, 1.0, 0.5),
            (0.75, 0.75, 0.75, 0.5, 1.0, 1.0),
            (0.75, 0.75, 0.75, 1.0, 0.5, 1.0)
        ]
        Hybrid2 = [
            (0.5, 0.0, 0.0, 0.5, 0.5, 0.5),
            (1.0, 0.0, 0.5, 0.5, 0.5, 0.5),
            (0.5, 0.0, 1.0, 0.5, 0.5, 0.5),
            (0.0, 0.0, 0.5, 0.5, 0.5, 0.5),
            (0.5, 1.0, 0.0, 0.5, 0.5, 0.5),
            (1.0, 1.0, 0.5, 0.5, 0.5, 0.5),
            (0.5, 1.0, 1.0, 0.5, 0.5, 0.5),
            (0.0, 1.0, 0.5, 0.5, 0.5, 0.5),
            (0.0, 0.5, 0.0, 0.5, 0.5, 0.5),
            (0.0, 0.5, 1.0, 0.5, 0.5, 0.5),
            (1.0, 0.5, 0.0, 0.5, 0.5, 0.5),
            (1.0, 0.5, 1.0, 0.5, 0.5, 0.5)
        ]
        if (Lattice == -1):
            Lattice = random.randint(0, 10)
        if (Lattice == 0):
            return BCC
        if (Lattice == 1):
            return Octet
        if (Lattice == 2):
            return OctetExt
        if (Lattice == 3):
            return OctetInt
        if (Lattice == 4):
            return BCCZ
        if (Lattice == 5):
            return Cubic
        if (Lattice == 6):
            return OctahedronZ
        if (Lattice == 7):
            return OctahedronZcross
        if (Lattice == 8):
            return Kelvin
        if (Lattice == 9):
            return CubicV2
        if (Lattice == 10):
            return CubicV3
        if (Lattice == 11):
            return CubicV4
        if (Lattice == 12):
            return Newlattice
        if (Lattice == 13):
            return Diamond
        if (Lattice == 14):
            return Auxetic
        if (Lattice == 15):
            return Hichem
        if (Lattice == 16):
            return Hybrid1
        if (Lattice == 17):
            return Hybrid2
        


    def generate_beams_from_given_point_list(self, latticeType, Radius, gradRadius, gradDim, gradMat, posCell):
        """
        Generate beams and nodes using a given point list and lattice parameters.

        :param latticeType: Type of lattice geometry
        :param Radius: Beam radius
        :param gradRadius: Gradient of beam radius
        :param gradDim: Gradient of dimensions
        :param gradMat: Gradient of material
        :param posCell: Position of the cell
        :return: List of generated Node objects and Beam objects
        """
        self.beams = []
        self.nodes = []
        for line in self.Lattice_geometry(latticeType):
            x1, y1, z1, x2, y2, z2 = map(float, line)
            point1 = Point((x1)*self.cell_size_x+self.x, (y1)*self.cell_size_y+self.y, (z1)*self.cell_size_z+self.z)
            point2 = Point((x2)*self.cell_size_x+self.x, (y2)*self.cell_size_y+self.y, (z2)*self.cell_size_z+self.z)
            self.nodes.append(point1)
            self.nodes.append(point2)
            beam = Beam(point1, point2, Radius, self.cell_size_x, self.cell_size_y, self.cell_size_z, gradRadius, gradMat,posCell,0)
            self.beams.append(beam)
        return self.nodes, self.beams

    def random_coordinate(self,coord, mu, sigma):
        mod_coord = []
        for pos in coord:
            # mod_coord.append(pos + random.uniform(0, 0.5) - 0.25)
            mod_coord.append(pos + random.gauss(mu, sigma))
        return mod_coord

    def add_point(self, point):
        x, y, z = map(float, point)
        point_obj = Point((x) * self.cell_size_x + self.x, (y) * self.cell_size_y + self.y,
                      (z) * self.cell_size_z + self.z)
        self.nodes.append(point_obj)

    def generate_beams_random(self, Radius, gradRadius, gradDim, gradMat, posCell):
        self.beams = []
        self.nodes = []
        # Corner node
        corner_node = random.randint(0,1)
        if corner_node == 1: # Corner nodes
            map_corner = [(0.0, 0.0, 0.0),
                   (1.0, 1.0, 1.0),
                   (1.0, 1.0, 0.0),
                   (0.0, 0.0, 1.0),
                   (0.0, 1.0, 0.0),
                   (0.0, 1.0, 1.0),
                   (1.0, 0.0, 1.0),
                   (1.0, 0.0, 0.0)]
            for idx, point_corner in enumerate(map_corner):
                self.add_point(point_corner)
        # Edge
        map_edge = [(0.0, 0.0), (0.0, 1.0), (1.0, 0.0), (1.0, 1.0)]
        for i in range(3): # 3 direction of edge node
            Edge_node = random.randint(0,1)
            if Edge_node == 1:
                point_mod = 0.5 + random.uniform(0, 0.5) - 0.25
                for idx, point_edge in enumerate(map_edge):
                    point = list(point_edge)
                    point.insert(i, point_mod)
                    self.add_point(point)
        # Face
        map_face = [(0.25, 0.25), (0.75, 0.25), (0.5, 0.5), (0.25, 0.75), (0.75, 0.75)]
        for i in range(3):  # 3 direction of face node
            Face_node = [random.randint(0,1) for _ in range(5)]
            for idx, point_face in enumerate(map_face):
                if Face_node[idx] == 1:
                    point = list(self.random_coordinate(map_face[idx], 0, 0.25))
                    point_sym = list(point)
                    point.insert(i, 0.0)
                    point_sym.insert(i, 1.0)
                    self.add_point(point)
                    self.add_point(point_sym)
        # Interior
        map_interior = [(0.25, 0.25, 0.25), (0.5, 0.25, 0.25), (0.75, 0.25, 0.25),
                        (0.25, 0.5, 0.25), (0.5, 0.5, 0.25), (0.75, 0.5, 0.25),
                        (0.25, 0.75, 0.25), (0.5, 0.75, 0.25), (0.75, 0.75, 0.25),

                        (0.25, 0.25, 0.5), (0.5, 0.25, 0.5), (0.75, 0.25, 0.5),
                        (0.25, 0.5, 0.5), (0.5, 0.5, 0.5), (0.75, 0.5, 0.5),
                        (0.25, 0.75, 0.5), (0.5, 0.75, 0.5), (0.75, 0.75, 0.5),

                        (0.25, 0.25, 0.75), (0.5, 0.25, 0.75), (0.75, 0.25, 0.75),
                        (0.25, 0.5, 0.75), (0.5, 0.5, 0.75), (0.75, 0.5, 0.75),
                        (0.25, 0.75, 0.75), (0.5, 0.75, 0.75), (0.75, 0.75, 0.75)]
        interior_node = [random.randint(0, 1) for _ in range(len(map_interior))]
        for idx, point_interior in enumerate(map_interior):
            if interior_node[idx] == 1:
                point = list(self.random_coordinate(map_interior[idx], 0, 0.1))
                self.add_point(point)

        self.remove_nodes_outside_unit_cube()
        print(len(self.nodes))

        self.remove_nodes_too_close_advanced()
        print(len(self.nodes))



        # Generate beam at least 2 beams per node
        for beam in range(1):
            beam = Beam(self.nodes[random.randint(0, len(self.nodes)-1)], self.nodes[random.randint(0, len(self.nodes)-1)], Radius,
                        self.cell_size_x, self.cell_size_y, self.cell_size_z, gradRadius, gradMat,
                        posCell, 0)
            if self.beam_already_exist(beam):
                self.beams.append(beam)
        return self.nodes, self.beams

    def beam_already_exist(self, beam_test):
        # print(beam_test.point1, beam_test.point2)
        for beam in self.beams:
            if beam.point1 == beam_test.point1:
                if beam.point2 == beam_test.point2:
                    # print(beam.point1,beam.point2)
                    # print(False)
                    return False
            if beam.point1 == beam_test.point2:
                if beam.point2 == beam_test.point1:
                    # print(beam.point1,beam.point2)
                    # print(False)
                    return False
        # print(True)
        return True

    def remove_nodes_outside_unit_cube(self):
        """
        Remove nodes that are outside of the unit cube defined by the cell dimensions.
        """
        inside_nodes = []
        for point in self.nodes:
            if 0 <= point.x <= self.cell_size_x and 0 <= point.y <= self.cell_size_y and 0 <= point.z <= self.cell_size_z:
                inside_nodes.append(point)
        self.nodes = inside_nodes

    def is_node_in_corner(self, node):
        """
        Check if a node is in a corner of the unit cube.

        :param node: The node to check.
        :return: True if the node is in a corner, False otherwise.
        """
        return ((node.x in [0, self.cell_size_x]) and
                (node.y in [0, self.cell_size_y]) and
                (node.z in [0, self.cell_size_z]))

    def remove_nodes_too_close_advanced(self, min_distance=0.1):
        """
        Remove nodes that are too close together with a refined decision process.
        """
        to_remove = set()


        def is_surface_not_corner(node):
            return self.is_node_on_surface(node) and not self.is_node_in_corner(node)

        for i, node_i in enumerate(self.nodes):
            if i in to_remove:
                continue

            for j, node_j in enumerate(self.nodes[i + 1:], start=i + 1):
                if j in to_remove:
                    continue

                # Calculate Euclidean distance
                distance = math.sqrt((node_i.x - node_j.x) ** 2 +
                                     (node_i.y - node_j.y) ** 2 +
                                     (node_i.z - node_j.z) ** 2)

                if distance < min_distance:

                    if self.is_node_in_corner(node_i):
                        to_remove.add(j)
                        continue
                    elif self.is_node_in_corner(node_j):
                        to_remove.add(i)
                        continue


                    if self.is_node_on_surface(node_i) and not self.is_node_on_surface(node_j):
                        to_remove.add(i)
                    elif not self.is_node_on_surface(node_i) and self.is_node_on_surface(node_j):
                        to_remove.add(j)
                    else:

                        to_remove.add(j)  # Choix arbitraire pour supprimer
                        if is_surface_not_corner(node_j):
                            symmetric_node = self.find_symmetric_node(node_j)
                            if symmetric_node and not self.is_node_in_corner(symmetric_node):
                                symmetric_index = self.nodes.index(symmetric_node)
                                to_remove.add(symmetric_index)

        self.nodes = [node for index, node in enumerate(self.nodes) if index not in to_remove]

    def is_node_on_surface(self, node):
        """
        Check if a node is on the surface of the unit cube.

        :param node: The node to check.
        :return: True if the node is on the surface, False otherwise.
        """
        return (node.x in [0, self.cell_size_x] or
                node.y in [0, self.cell_size_y] or
                node.z in [0, self.cell_size_z])

    def find_symmetric_node(self, node):
        """
        Find a symmetric node to a given node, based on the cube's faces. If a node is on the face X=0,
        its symmetric is on the face X=cellSizeX, and analogously for Y and Z axes.

        :param node: The node for which to find a symmetric counterpart.
        :return: The symmetric node, if found.
        """

        symmetric_x = self.cell_size_x - node.x if node.x in [0, self.cell_size_x] else node.x
        symmetric_y = self.cell_size_y - node.y if node.y in [0, self.cell_size_y] else node.y
        symmetric_z = self.cell_size_z - node.z if node.z in [0, self.cell_size_z] else node.z


        for n in self.nodes:
            if n.x == symmetric_x and n.y == symmetric_y and n.z == symmetric_z:
                return n
        return None

    def calculate_connections(self):
        """
        Calculate the number of connections for each node.

        :return: Dictionary containing the number of connections for each node index
        """
        connections = {}
        for i, point in enumerate(self.nodes):
            connections[i] = 0
        for beam in self.beams:
            start_point_index = self.nodes.index(beam.point1)
            end_point_index = self.nodes.index(beam.point2)
            connections[start_point_index] += 1
            connections[end_point_index] += 1
        return connections

    def remove_isolated_beams(self):
        """
        Remove isolated beams (beams connected to only one node).
        """
        point_counts = {}

        for beam in self.beams:
            point1 = beam.point1
            point2 = beam.point2

            if point1 not in point_counts:
                point_counts[point1] = 0
            if point2 not in point_counts:
                point_counts[point2] = 0

            point_counts[point1] += 1
            point_counts[point2] += 1

        self.beams = [beam for beam in self.beams if point_counts[beam.point1] > 1 and point_counts[beam.point2] > 1]

    def num_connections(self):
        """
        Calculate the total number of connections in the cell.

        :return: Total number of connections
        """
        total_connections = 0
        for beam in self.beams:
            total_connections += beam.num_connections()
        return total_connections

    def display_point(self, ax, color1=None, color2=None, color3=None):
        """
        Display nodes in the 3D plot.

        :param ax: Matplotlib 3D axis object
        :param color1: Color for corner points
        :param color2: Color for corner-edge points
        :param color3: Color for edge points
        """
        for point in self.nodes:
            x, y, z = point.x, point.y, point.z
            color = color1
            if z == 0 or z == self.cell_size_z:
                if (x == 0 or x == self.cell_size_x) and (y == 0 or y == self.cell_size_y):
                    color = color2
                else:
                    color = 'y'
            elif (x == 0 or x == self.cell_size_x) or (y == 0 or y == self.cell_size_y):
                color = color3
            ax.scatter(x, y, z, c=color)

    def display_point_simple(self, ax, color):
        """
        Display nodes in the 3D plot using a simple color.

        :param ax: Matplotlib 3D axis object
        :param color: Color for points
        """
        for point in self.nodes:
            x, y, z = point.x, point.y, point.z
            ax.scatter(x, y, z, c=color)
            
    def display_beams(self, ax):
        """
        Display beams in the 3D plot.

        :param ax: Matplotlib 3D axis object
        """
        color = ['blue','green','black','yellow','orange']
        for index, beam in enumerate(self.beams):
            point1 = beam.point1
            point2 = beam.point2
            ax.plot([point1.x, point2.x], [point1.y, point2.y], [point1.z, point2.z], color=color[beam.material])

    def visualize_3d(self, ax):
        """
        Visualize the cell in a 3D plot.

        :param ax: Matplotlib 3D axis object
        """
        self.display_point(ax, 'r', 'pink', 'black')
        self.display_beams(ax, 'b', 'r')
    
    def calculateCenterCell(self):
        return [self.x+self.cell_size_x/2,self.y+self.cell_size_y/2,self.z+self.cell_size_z/2]
