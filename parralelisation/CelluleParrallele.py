from Codes.Point import *
from Codes.Beam import *
import random
import concurrent.futures


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
        # Translater la cellule selon la translation spécifiée
        self.x += translation[0]
        self.y += translation[1]
        self.z += translation[2]

        for point in self.points:
            point.x += translation[0]
            point.y += translation[1]
            point.z += translation[2]

    def generate_points(self, Tmin, Tmax, padding):
        def generate_points_for_cell():
            pointtotal, pointbord, pointcoin, pointinterieur = self.nbPointsTBCI(Tmin, Tmax)
            points_interior = []
            points_border = []
            points_corner = []
            for _ in range(pointinterieur):
                x = round(random.uniform(0 + padding, self.cell_size_x - padding), 2)
                y = round(random.uniform(0 + padding, self.cell_size_y - padding), 2)
                z = round(random.uniform(0 + padding, self.cell_size_z - padding), 2)
                # point = Point(x, y, z)
                point = Point(round(x, 2), round(y, 2), round(z, 2))
                points_interior.append(point)
            for _ in range(pointcoin):
                x = round(random.choice([0, self.cell_size_x]), 2)
                y = round(random.choice([0, self.cell_size_y]), 2)
                z = round(random.choice([0, self.cell_size_z]), 2)
                # point = Point(x, y, z)
                point = Point(round(x, 2), round(y, 2), round(z, 2))
                points_corner.append(point)
            for _ in range(pointbord):
                if random.random() < 0.5:
                    x = round(random.choice([0, self.cell_size_x]), 2)
                    y = round(random.uniform(0, self.cell_size_y), 2)
                    z = round(random.uniform(0, self.cell_size_z), 2)
                else:
                    x = round(random.uniform(0, self.cell_size_x), 2)
                    y = round(random.choice([0, self.cell_size_y]), 2)
                    z = round(random.uniform(0, self.cell_size_z), 2)
                point = Point(round(x, 2), round(y, 2), round(z, 2))
                points_border.append(point)
            points = points_interior + points_border + points_corner
            return pointtotal, points_border, points_corner, points_interior, points

        with concurrent.futures.ThreadPoolExecutor() as executor:
            result = executor.submit(generate_points_for_cell)

        pointtotal, points_border, points_corner, points_interior, self.points = result.result()
        return pointtotal, points_border, points_corner, points_interior

    def merge_points(self):
        z_values = [point.z for point in self.points]
        z_min = min(z_values)
        z_max = max(z_values)
        merged_points = set()
        for point in self.points:
            x, y, z = point.x, point.y, point.z
            if z == z_min or (x in [0, self.cell_size_x] and y in [0, self.cell_size_y]) or (
                    z == z_max and x != 0 and x != self.cell_size_x and y != 0 and y != self.cell_size_y):
                merged_points.add(Point(round(x, 2), round(y, 2), round(z_max, 2)))
            elif z == z_max or (x in [0, self.cell_size_x] and y in [0, self.cell_size_y]) or (
                    z == z_min and x != 0 and x != self.cell_size_x and y != 0 and y != self.cell_size_y):
                merged_points.add(Point(round(x, 2), round(y, 2), round(z_min, 2)))
            else:
                merged_points.add(point)
        self.points = list(merged_points)  # Convert set back to a list

    def affichage_points(self):
        for index, point in enumerate(self.points):
            print(f"Point {index}: {point.x}, {point.y}, {point.z}")

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
            # Extraire les coordonnées des points et les paramètres du faisceau
            x1, y1, z1, x2, y2, z2, radius = map(float, line)
            # x1, y1, z1, x2, y2, z2 = map(float, line)
            # Créer les objets Point
            point1 = Point(x1, y1, z1)
            point2 = Point(x2, y2, z2)
            self.points.append(point1)
            self.points.append(point2)
            # Créer l'objet Beam avec les points et les paramètres
            beam = Beam(point1, point2, radius=radius)
            # beam = Beam(point1, point2)
            # Ajouter le faisceau à la liste des beams
            self.beams.append(beam)
        return self.points, self.beams

    def generate_beams(self):
        self.beams = []
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
        return self.points, self.beams

    def remove_unused_points(self):
        used_points = set()
        for beam in self.beams:
            used_points.add(beam.point1)
            used_points.add(beam.point2)
        self.points = list(used_points)
        return self.points

    def display_pointsAfterRemove(self):
        for index, point in enumerate(self.points):
            print(f"Point {index}: {point.x}, {point.y}, {point.z}")

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

    def num_points(self):
        return len(self.points)

    def num_beams(self):
        return len(self.beams)

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

    def display_beams(self, ax, color1=None, color2=None):
        connections = self.calculate_connections()
        for beam in self.beams:
            point1, point2 = beam.point1, beam.point2
            x = [point1.x, point2.x]
            y = [point1.y, point2.y]
            z = [point1.z, point2.z]
            color = color1
            if connections[self.points.index(point1)] == 1 and connections[self.points.index(point2)] == 1:
                color = color2
            ax.plot(x, y, z, c=color)

    def visualize_3d(self, ax):
        self.display_point(ax, 'r', 'pink', 'black')
        self.display_beams(ax, 'b', 'r')
