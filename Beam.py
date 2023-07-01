import random
import math


class Beam:
    def __init__(self, point1, point2, length=None, material=None, radius=None):
        self.point1 = point1
        self.point2 = point2
        self.length = length if length is not None else self.get_length(point1, point2)
        self.material = material if material is not None else random.choice([1, 2, 3])
        self.radius = radius if radius is not None else random.uniform(0.4, 1.0)

    def get_length(self, point1, point2):
        x1, y1, z1 = point1.x, point1.y, point1.z
        x2, y2, z2 = point2.x, point2.y, point2.z
        length = math.sqrt((x2 - x1) ** 2 + (y2 - y1) ** 2 + (z2 - z1) ** 2)
        return length
