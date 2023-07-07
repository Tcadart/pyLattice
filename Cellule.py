import matplotlib.pyplot as plt
from mpl_toolkits.mplot3d import Axes3D
from Point import *
from Beam import *


class Cellule:
    def __init__(self, cell_size_x, cell_size_y, cell_size_z, x, y, z):
        self.cell_size_x = cell_size_x
        self.cell_size_y = cell_size_y
        self.cell_size_z = cell_size_z
        self.x = x
        self.y = y
        self.z = z
        self.points = []
        self.beams = []

    def nbPointsTBCI(self, Tmin, Tmax):
        pointtotal = random.randint(Tmin, Tmax)
        pointinterieur = random.randint(0, pointtotal)
        proportion_coins = random.uniform(0, 1)
        pointBordOuCoins = pointtotal - pointinterieur
        pointcoin = int(pointBordOuCoins * proportion_coins)
        pointbord = pointBordOuCoins - pointcoin
        if pointcoin > 8:
            pointbord += pointcoin - 8
            pointcoin = 8
        elif pointcoin < 0:
            pointcoin = 0
            pointbord += pointcoin
        return pointtotal, pointbord, pointcoin, pointinterieur

    def translate(self, translation):
        self.x += translation[0]
        self.y += translation[1]
        self.z += translation[2]

        for point in self.points:
            point.x += translation[0]
            point.y += translation[1]
            point.z += translation[2]

    # def generate_points(self, Tmin, Tmax, padding):
    #     self.points = []
    #     pointtotal, pointbord, pointcoin, pointinterieur = self.nbPointsTBCI(Tmin, Tmax)
    #     points_interior = []
    #     points_border = []
    #     points_corner = []
    #     for _ in range(pointinterieur):
    #         x = round(random.uniform(0 + padding, self.cell_size_x - padding), 2)
    #         y = round(random.uniform(0 + padding, self.cell_size_y - padding), 2)
    #         z = round(random.uniform(0 + padding, self.cell_size_z - padding), 2)
    #         point = Point(x, y, z)
    #         points_interior.append(point)
    #     for _ in range(pointcoin):
    #         x = random.choice([0, self.cell_size_x])
    #         y = random.choice([0, self.cell_size_y])
    #         z = random.choice([0, self.cell_size_z])
    #         point = Point(x, y, z)
    #         points_corner.append(point)
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
    #         points_border.append(point)
    #     self.points = points_interior + points_border + points_corner
    #     # return self.points, points_border, points_corner, points_interior
    #     return self.points

    def generate_points(self, Tmin, Tmax, padding):
        pointtotal, pointbord, pointcoin, pointinterieur = self.nbPointsTBCI(Tmin, Tmax)
        points_interior = []
        points_border = []
        points_corner = []
        min_distance = 0.02

        def check_min_distance(x, y, z, existing_points):
            for point in existing_points:
                distance = ((point.x - x) ** 2 + (point.y - y) ** 2 + (point.z - z) ** 2) ** 0.5
                if distance < min_distance:
                    return False
            return True

        for _ in range(pointinterieur):
            x, y, z = None, None, None
            while True:
                x = round(random.uniform(0 + padding, self.cell_size_x - padding), 2)
                y = round(random.uniform(0 + padding, self.cell_size_y - padding), 2)
                z = round(random.uniform(0 + padding, self.cell_size_z - padding), 2)
                if check_min_distance(x, y, z, points_interior):
                    break
            point = Point(x, y, z)
            points_interior.append(point)

        for _ in range(pointcoin):
            x, y, z = None, None, None
            while True:
                x = random.choice([0, self.cell_size_x])
                y = random.choice([0, self.cell_size_y])
                z = random.choice([0, self.cell_size_z])
                if check_min_distance(x, y, z, points_corner):
                    break
            point = Point(x, y, z)
            points_corner.append(point)

        for _ in range(pointbord):
            x, y, z = None, None, None
            while True:
                if random.random() < 0.5:
                    x = random.choice([0, self.cell_size_x])
                    y = random.uniform(0, self.cell_size_y)
                    z = random.uniform(0, self.cell_size_z)
                else:
                    x = random.uniform(0, self.cell_size_x)
                    y = random.choice([0, self.cell_size_y])
                    z = random.uniform(0, self.cell_size_z)
                if check_min_distance(x, y, z, points_border):
                    break
            point = Point(round(x, 2), round(y, 2), round(z, 2))
            points_border.append(point)

        self.points = points_interior + points_border + points_corner
        return self.points

    def merge_points(self):
        z_values = [point.z for point in self.points]
        z_min = min(z_values)
        z_max = max(z_values)
        merged_points = set()
        for point in self.points:
            x, y, z = point.x, point.y, point.z
            if z == z_min or (x in [0, self.cell_size_x] and y in [0, self.cell_size_y]) or (
                    z == z_max and x != 0 and x != self.cell_size_x and y != 0 and y != self.cell_size_y):
                merged_points.add(Point(x, y, z_max))
            elif z == z_max or (x in [0, self.cell_size_x] and y in [0, self.cell_size_y]) or (
                    z == z_min and x != 0 and x != self.cell_size_x and y != 0 and y != self.cell_size_y):
                merged_points.add(Point(x, y, z_min))
            else:
                merged_points.add(point)
        self.points = list(merged_points)

    def generate_beams(self):
        num_points = len(self.points)
        min_beams = num_points // 2
        max_beams = num_points
        num_beams = random.randint(min_beams, max_beams)
        used_pairs = set()
        used_beams = set()
        while num_beams > 0:
            index1 = random.randint(0, num_points - 1)
            index2 = random.randint(0, num_points - 1)
            if index1 != index2 and (index1, index2) not in used_pairs and (index2, index1) not in used_pairs:
                point1 = self.points[index1]
                point2 = self.points[index2]
                beam = Beam(point1, point2)
                if beam not in used_beams:
                    self.beams.append(beam)
                    used_pairs.add((index1, index2))
                    used_beams.add(beam)
                    num_beams -= 1
        return self.beams

    def remove_unused_points(self):
        used_points = set()
        for beam in self.beams:
            used_points.add(beam.point1)
            used_points.add(beam.point2)
        self.points = list(used_points)
        return self.points

    def Lattice_geometry(self, Lattice, Radius_geom):
        BCC = [(0.0, 0.0, 0.0, 0.5, 0.5, 0.5, Radius_geom),
               (0.5, 0.5, 0.5, 1.0, 1.0, 1.0, Radius_geom),
               (0.5, 0.5, 0.5, 1.0, 1.0, 0.0, Radius_geom),
               (0.5, 0.5, 0.5, 0.0, 0.0, 1.0, Radius_geom),
               (0.5, 0.5, 0.5, 0.0, 1.0, 0.0, Radius_geom),
               (0.5, 0.5, 0.5, 0.0, 1.0, 1.0, Radius_geom),
               (1.0, 0.0, 1.0, 0.5, 0.5, 0.5, Radius_geom),
               (0.5, 0.5, 0.5, 1.0, 0.0, 0.0, Radius_geom)]
        Octet = [(0.0, 0.0, 0.0, 0.5, 0.0, 0.5, Radius_geom),
                 (1.0, 0.0, 1.0, 0.5, 0.0, 0.5, Radius_geom),
                 (0.0, 0.0, 1.0, 0.5, 0.0, 0.5, Radius_geom),
                 (1.0, 0.0, 0.0, 0.5, 0.0, 0.5, Radius_geom),
                 (0.0, 0.0, 0.0, 0.0, 0.5, 0.5, Radius_geom),
                 (0.0, 1.0, 1.0, 0.0, 0.5, 0.5, Radius_geom),
                 (0.0, 0.0, 1.0, 0.0, 0.5, 0.5, Radius_geom),
                 (0.0, 1.0, 0.0, 0.0, 0.5, 0.5, Radius_geom),
                 (0.0, 0.0, 0.0, 0.5, 0.5, 0.0, Radius_geom),
                 (1.0, 1.0, 0.0, 0.5, 0.5, 0.0, Radius_geom),
                 (1.0, 0.0, 0.0, 0.5, 0.5, 0.0, Radius_geom),
                 (0.0, 1.0, 0.0, 0.5, 0.5, 0.0, Radius_geom),
                 (0.0, 0.0, 1.0, 0.5, 0.5, 1.0, Radius_geom),
                 (1.0, 1.0, 1.0, 0.5, 0.5, 1.0, Radius_geom),
                 (1.0, 0.0, 1.0, 0.5, 0.5, 1.0, Radius_geom),
                 (0.0, 1.0, 1.0, 0.5, 0.5, 1.0, Radius_geom),
                 (1.0, 0.5, 0.5, 1.0, 1.0, 1.0, Radius_geom),
                 (1.0, 0.0, 0.0, 1.0, 0.5, 0.5, Radius_geom),
                 (1.0, 0.5, 0.5, 1.0, 1.0, 0.0, Radius_geom),
                 (1.0, 0.0, 1.0, 1.0, 0.5, 0.5, Radius_geom),
                 (0.5, 1.0, 0.5, 1.0, 1.0, 1.0, Radius_geom),
                 (0.0, 1.0, 0.0, 0.5, 1.0, 0.5, Radius_geom),
                 (0.5, 1.0, 0.5, 1.0, 1.0, 0.0, Radius_geom),
                 (0.0, 1.0, 1.0, 0.5, 1.0, 0.5, Radius_geom),
                 (0.5, 0.0, 0.5, 0.5, 0.5, 0.0, Radius_geom),
                 (0.5, 0.0, 0.5, 0.0, 0.5, 0.5, Radius_geom),
                 (0.5, 0.0, 0.5, 1.0, 0.5, 0.5, Radius_geom),
                 (0.5, 0.0, 0.5, 0.5, 0.5, 1.0, Radius_geom),
                 (0.5, 1.0, 0.5, 0.5, 0.5, 0.0, Radius_geom),
                 (0.5, 1.0, 0.5, 0.0, 0.5, 0.5, Radius_geom),
                 (0.5, 1.0, 0.5, 1.0, 0.5, 0.5, Radius_geom),
                 (0.5, 1.0, 0.5, 0.5, 0.5, 1.0, Radius_geom),
                 (1.0, 0.5, 0.5, 0.5, 0.5, 0.0, Radius_geom),
                 (0.5, 0.5, 0.0, 0.0, 0.5, 0.5, Radius_geom),
                 (0.5, 0.5, 1.0, 0.0, 0.5, 0.5, Radius_geom),
                 (0.5, 0.5, 1.0, 1.0, 0.5, 0.5, Radius_geom)]
        OctetExt = [(0.0, 0.0, 0.0, 0.5, 0.0, 0.5, Radius_geom),
                    (1.0, 0.0, 1.0, 0.5, 0.0, 0.5, Radius_geom),
                    (0.0, 0.0, 1.0, 0.5, 0.0, 0.5, Radius_geom),
                    (1.0, 0.0, 0.0, 0.5, 0.0, 0.5, Radius_geom),
                    (0.0, 0.0, 0.0, 0.0, 0.5, 0.5, Radius_geom),
                    (0.0, 1.0, 1.0, 0.0, 0.5, 0.5, Radius_geom),
                    (0.0, 0.0, 1.0, 0.0, 0.5, 0.5, Radius_geom),
                    (0.0, 1.0, 0.0, 0.0, 0.5, 0.5, Radius_geom),
                    (0.0, 0.0, 0.0, 0.5, 0.5, 0.0, Radius_geom),
                    (1.0, 1.0, 0.0, 0.5, 0.5, 0.0, Radius_geom),
                    (1.0, 0.0, 0.0, 0.5, 0.5, 0.0, Radius_geom),
                    (0.0, 1.0, 0.0, 0.5, 0.5, 0.0, Radius_geom),
                    (0.0, 0.0, 1.0, 0.5, 0.5, 1.0, Radius_geom),
                    (1.0, 1.0, 1.0, 0.5, 0.5, 1.0, Radius_geom),
                    (1.0, 0.0, 1.0, 0.5, 0.5, 1.0, Radius_geom),
                    (0.0, 1.0, 1.0, 0.5, 0.5, 1.0, Radius_geom),
                    (1.0, 0.5, 0.5, 1.0, 1.0, 1.0, Radius_geom),
                    (1.0, 0.0, 0.0, 1.0, 0.5, 0.5, Radius_geom),
                    (1.0, 0.5, 0.5, 1.0, 1.0, 0.0, Radius_geom),
                    (1.0, 0.0, 1.0, 1.0, 0.5, 0.5, Radius_geom),
                    (0.5, 1.0, 0.5, 1.0, 1.0, 1.0, Radius_geom),
                    (0.0, 1.0, 0.0, 0.5, 1.0, 0.5, Radius_geom),
                    (0.5, 1.0, 0.5, 1.0, 1.0, 0.0, Radius_geom),
                    (0.0, 1.0, 1.0, 0.5, 1.0, 0.5, Radius_geom)]
        OctetInt = [(0.5, 0.0, 0.5, 0.5, 0.5, 0.0, Radius_geom),
                    (0.5, 0.0, 0.5, 0.0, 0.5, 0.5, Radius_geom),
                    (0.5, 0.0, 0.5, 1.0, 0.5, 0.5, Radius_geom),
                    (0.5, 0.0, 0.5, 0.5, 0.5, 1.0, Radius_geom),
                    (0.5, 1.0, 0.5, 0.5, 0.5, 0.0, Radius_geom),
                    (0.5, 1.0, 0.5, 0.0, 0.5, 0.5, Radius_geom),
                    (0.5, 1.0, 0.5, 1.0, 0.5, 0.5, Radius_geom),
                    (0.5, 1.0, 0.5, 0.5, 0.5, 1.0, Radius_geom),
                    (1.0, 0.5, 0.5, 0.5, 0.5, 0.0, Radius_geom),
                    (0.5, 0.5, 0.0, 0.0, 0.5, 0.5, Radius_geom),
                    (0.5, 0.5, 1.0, 0.0, 0.5, 0.5, Radius_geom),
                    (0.5, 0.5, 1.0, 1.0, 0.5, 0.5, Radius_geom)]
        BCCZ = [(0.5, 0.5, 0.5, 1.0, 1.0, 1.0, Radius_geom),
                (0.0, 0.0, 0.0, 0.5, 0.5, 0.5, Radius_geom),
                (0.5, 0.5, 0.5, 1.0, 1.0, 0.0, Radius_geom),
                (0.0, 0.0, 1.0, 0.5, 0.5, 0.5, Radius_geom),
                (0.5, 0.5, 0.5, 0.0, 1.0, 0.0, Radius_geom),
                (1.0, 0.0, 1.0, 0.5, 0.5, 0.5, Radius_geom),
                (0.5, 0.5, 0.5, 0.0, 1.0, 1.0, Radius_geom),
                (1.0, 0.0, 0.0, 0.5, 0.5, 0.5, Radius_geom),
                (0.5, 0.5, 0.0, 0.5, 0.5, 0.5, Radius_geom),
                (0.5, 0.5, 0.5, 0.5, 0.5, 1.0, Radius_geom)]
        Cubic = [(0.0, 0.0, 0.0, 0.0, 0.0, 1.0, Radius_geom),
                 (1.0, 0.0, 0.0, 1.0, 0.0, 1.0, Radius_geom),
                 (0.0, 1.0, 0.0, 0.0, 1.0, 1.0, Radius_geom),
                 (1.0, 1.0, 0.0, 1.0, 1.0, 1.0, Radius_geom),
                 (0.0, 0.0, 0.0, 1.0, 0.0, 0.0, Radius_geom),
                 (0.0, 0.0, 0.0, 0.0, 1.0, 0.0, Radius_geom),
                 (1.0, 1.0, 0.0, 0.0, 1.0, 0.0, Radius_geom),
                 (1.0, 1.0, 0.0, 1.0, 0.0, 0.0, Radius_geom),
                 (0.0, 0.0, 1.0, 1.0, 0.0, 1.0, Radius_geom),
                 (0.0, 0.0, 1.0, 0.0, 1.0, 1.0, Radius_geom),
                 (1.0, 1.0, 1.0, 0.0, 1.0, 1.0, Radius_geom),
                 (1.0, 1.0, 1.0, 1.0, 0.0, 1.0, Radius_geom)]
        OctahedronZ = [(0.5, 0.0, 0.5, 0.5, 0.5, 0.0, Radius_geom),
                       (0.5, 0.0, 0.5, 0.0, 0.5, 0.5, Radius_geom),
                       (0.5, 0.0, 0.5, 1.0, 0.5, 0.5, Radius_geom),
                       (0.5, 0.0, 0.5, 0.5, 0.5, 1.0, Radius_geom),
                       (0.5, 1.0, 0.5, 0.5, 0.5, 0.0, Radius_geom),
                       (0.5, 1.0, 0.5, 0.0, 0.5, 0.5, Radius_geom),
                       (0.5, 1.0, 0.5, 1.0, 0.5, 0.5, Radius_geom),
                       (0.5, 1.0, 0.5, 0.5, 0.5, 1.0, Radius_geom),
                       (1.0, 0.5, 0.5, 0.5, 0.5, 0.0, Radius_geom),
                       (0.5, 0.5, 0.0, 0.0, 0.5, 0.5, Radius_geom),
                       (0.5, 0.5, 1.0, 0.0, 0.5, 0.5, Radius_geom),
                       (0.5, 0.5, 1.0, 1.0, 0.5, 0.5, Radius_geom),
                       (0.5, 0.5, 0.0, 0.5, 0.5, 1.0, Radius_geom)]
        OctahedronZcross = [(0.5, 0.0, 0.5, 0.5, 0.5, 0.0, Radius_geom),
                            (0.5, 0.0, 0.5, 0.0, 0.5, 0.5, Radius_geom),
                            (0.5, 0.0, 0.5, 1.0, 0.5, 0.5, Radius_geom),
                            (0.5, 0.0, 0.5, 0.5, 0.5, 1.0, Radius_geom),
                            (0.5, 1.0, 0.5, 0.5, 0.5, 0.0, Radius_geom),
                            (0.5, 1.0, 0.5, 0.0, 0.5, 0.5, Radius_geom),
                            (0.5, 1.0, 0.5, 1.0, 0.5, 0.5, Radius_geom),
                            (0.5, 1.0, 0.5, 0.5, 0.5, 1.0, Radius_geom),
                            (1.0, 0.5, 0.5, 0.5, 0.5, 0.0, Radius_geom),
                            (0.5, 0.5, 0.0, 0.0, 0.5, 0.5, Radius_geom),
                            (0.5, 0.5, 1.0, 0.0, 0.5, 0.5, Radius_geom),
                            (0.5, 0.5, 1.0, 1.0, 0.5, 0.5, Radius_geom),
                            (0.5, 0.5, 0.0, 0.5, 0.5, 0.5, Radius_geom),
                            (0.5, 0.5, 1.0, 0.5, 0.5, 0.5, Radius_geom),
                            (0.5, 0.0, 0.5, 0.5, 0.5, 0.5, Radius_geom),
                            (0.0, 0.5, 0.5, 0.5, 0.5, 0.5, Radius_geom),
                            (1.0, 0.5, 0.5, 0.5, 0.5, 0.5, Radius_geom),
                            (0.5, 1.0, 0.5, 0.5, 0.5, 0.5, Radius_geom)]
        Kelvin = [(0.5, 0.25, 0, 0.25, 0.5, 0, Radius_geom),
                  (0.5, 0.25, 0, 0.75, 0.5, 0, Radius_geom),
                  (0.5, 0.75, 0, 0.25, 0.5, 0, Radius_geom),
                  (0.5, 0.75, 0, 0.75, 0.5, 0, Radius_geom),
                  (0.5, 0.25, 1, 0.25, 0.5, 1, Radius_geom),
                  (0.5, 0.25, 1, 0.75, 0.5, 1, Radius_geom),
                  (0.5, 0.75, 1, 0.25, 0.5, 1, Radius_geom),
                  (0.5, 0.75, 1, 0.75, 0.5, 1, Radius_geom),
                  (0.5, 0, 0.25, 0.25, 0, 0.5, Radius_geom),
                  (0.5, 0, 0.25, 0.75, 0, 0.5, Radius_geom),
                  (0.5, 0, 0.75, 0.25, 0, 0.5, Radius_geom),
                  (0.5, 0, 0.75, 0.75, 0, 0.5, Radius_geom),
                  (0.5, 1, 0.25, 0.25, 1, 0.5, Radius_geom),
                  (0.5, 1, 0.25, 0.75, 1, 0.5, Radius_geom),
                  (0.5, 1, 0.75, 0.25, 1, 0.5, Radius_geom),
                  (0.5, 1, 0.75, 0.75, 1, 0.5, Radius_geom),
                  (0, 0.5, 0.25, 0, 0.25, 0.5, Radius_geom),
                  (0, 0.5, 0.25, 0, 0.75, 0.5, Radius_geom),
                  (0, 0.5, 0.75, 0, 0.25, 0.5, Radius_geom),
                  (0, 0.5, 0.75, 0, 0.75, 0.5, Radius_geom),
                  (1, 0.5, 0.25, 1, 0.25, 0.5, Radius_geom),
                  (1, 0.5, 0.25, 1, 0.75, 0.5, Radius_geom),
                  (1, 0.5, 0.75, 1, 0.25, 0.5, Radius_geom),

                  (1, 0.5, 0.75, 1, 0.75, 0.5, Radius_geom),
                  (0.5, 0.25, 0, 0.5, 0, 0.25, Radius_geom),
                  (0.25, 0.5, 0, 0, 0.5, 0.25, Radius_geom),
                  (0.75, 0.5, 0, 1, 0.5, 0.25, Radius_geom),
                  (0.5, 0.75, 0, 0.5, 1, 0.25, Radius_geom),
                  (0.25, 0, 0.5, 0, 0.25, 0.5, Radius_geom),
                  (0.75, 0, 0.5, 1, 0.25, 0.5, Radius_geom),
                  (0.75, 1, 0.5, 1, 0.75, 0.5, Radius_geom),
                  (0.25, 1, 0.5, 0, 0.75, 0.5, Radius_geom),
                  (0.5, 0, 0.75, 0.5, 0.25, 1, Radius_geom),
                  (0, 0.5, 0.75, 0.25, 0.5, 1, Radius_geom),
                  (0.5, 1, 0.75, 0.5, 0.75, 1, Radius_geom),
                  (1, 0.5, 0.75, 0.75, 0.5, 1, Radius_geom)]
        CubicV2 = [(0.5, 0.0, 0.5, 0.5, 0.5, 0.5, Radius_geom),
                   (0.0, 0.5, 0.5, 0.5, 0.5, 0.5, Radius_geom),
                   (0.5, 1.0, 0.5, 0.5, 0.5, 0.5, Radius_geom),
                   (1.0, 0.5, 0.5, 0.5, 0.5, 0.5, Radius_geom),
                   (0.5, 0.5, 0.0, 0.5, 0.5, 0.5, Radius_geom),
                   (0.5, 0.5, 1.0, 0.5, 0.5, 0.5, Radius_geom)]
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

    def generate_beams_from_given_point_list(self, BCC):
        self.beams = []
        self.points = []
        for line in BCC:
            x1, y1, z1, x2, y2, z2, radius = map(float, line)
            # x1, y1, z1, x2, y2, z2 = map(float, line)
            point1 = Point(x1, y1, z1)
            point2 = Point(x2, y2, z2)
            self.points.append(point1)
            self.points.append(point2)
            beam = Beam(point1, point2, radius=radius)
            # beam = Beam(point1, point2)
            self.beams.append(beam)
        return self.points, self.beams

    def calculate_connections(self):
        connections = {}
        for i, point in enumerate(self.points):
            connections[i] = 0
        for beam in self.beams:
            start_point_index = self.points.index(beam.point1)
            end_point_index = self.points.index(beam.point2)
            connections[start_point_index] += 1
            connections[end_point_index] += 1
        return connections

    def remove_isolated_beams(self):
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
        total_connections = 0
        for beam in self.beams:
            total_connections += beam.num_connections()
        return total_connections

    def display_point(self, ax, color1=None, color2=None, color3=None):
        for point in self.points:
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

    def display_beams(self, ax, line_color, text_color):
        for index, beam in enumerate(self.beams):
            point1 = beam.point1
            point2 = beam.point2

            ax.plot([point1.x, point2.x], [point1.y, point2.y], [point1.z, point2.z], color=line_color)

    def visualize_3d(self, ax):
        self.display_point(ax, 'r', 'pink', 'black')
        self.display_beams(ax, 'b', 'r')
