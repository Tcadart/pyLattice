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
    #         x = round(random.uniform(0 + padding, self.cell_size_x - padding), 2)
    #         y = round(random.uniform(0 + padding, self.cell_size_y - padding), 2)
    #         z = round(random.uniform(0 + padding, self.cell_size_z - padding), 2)
    #         point = Point(x, y, z)
    #         nodes_interior.append(point)
    #     for _ in range(pointcoin):
    #         x = random.choice([0, self.cell_size_x])
    #         y = random.choice([0, self.cell_size_y])
    #         z = random.choice([0, self.cell_size_z])
    #         point = Point(x, y, z)
    #         nodes_corner.append(point)
    #     for _ in range(pointbord):
    #         if random.random() < 0.5:
    #             x = random.choice([0, self.cell_size_x])
    #             y = random.uniform(0, self.cell_size_y)
    #             z = random.uniform(0, self.cell_size_z)
    #         else:
    #             x = random.uniform(0, self.cell_size_x)
    #             y = random.choice([0, self.cell_size_y])
    #             z = random.uniform(0, self.cell_size_z)
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
    #             x = round(random.uniform(0 + padding, self.cell_size_x - padding), 2)
    #             y = round(random.uniform(0 + padding, self.cell_size_y - padding), 2)
    #             z = round(random.uniform(0 + padding, self.cell_size_z - padding), 2)
    #             if check_min_distance(x, y, z, nodes_interior):
    #                 break
    #         point = Point(x, y, z)
    #         nodes_interior.append(point)

    #     for _ in range(pointcoin):
    #         x, y, z = None, None, None
    #         while True:
    #             x = random.choice([0, self.cell_size_x])
    #             y = random.choice([0, self.cell_size_y])
    #             z = random.choice([0, self.cell_size_z])
    #             if check_min_distance(x, y, z, nodes_corner):
    #                 break
    #         point = Point(x, y, z)
    #         nodes_corner.append(point)

    #     for _ in range(pointbord):
    #         x, y, z = None, None, None
    #         while True:
    #             if random.random() < 0.5:
    #                 x = random.choice([0, self.cell_size_x])
    #                 y = random.uniform(0, self.cell_size_y)
    #                 z = random.uniform(0, self.cell_size_z)
    #             else:
    #                 x = random.uniform(0, self.cell_size_x)
    #                 y = random.choice([0, self.cell_size_y])
    #                 z = random.uniform(0, self.cell_size_z)
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
            (0.75, 0.75, 0.75, 1.0, 0.5, 1.0),
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
                  (0.0, 0.5, 1.0, 0.5, 0.5, 1.0),
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
