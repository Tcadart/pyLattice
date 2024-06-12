import numpy as np
class Point:
    def __init__(self, x, y, z):
        self.x = x
        self.y = y
        self.z = z

    def __eq__(self, other):
        return isinstance(other, Point) and self.x == other.x and self.y == other.y and self.z == other.z

    def __hash__(self):
        return hash((self.x, self.y, self.z))

    def __sub__(self, other):
        return np.array([self.x - other.x, self.y - other.y, self.z - other.z])

    def __repr__(self):
        return f"Point({self.x}, {self.y}, {self.z})"

    def movePoint(self, x, y, z):
        self.x = x
        self.y = y
        self.z = z

