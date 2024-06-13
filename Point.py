import numpy as np


class Point:
    def __init__(self, x, y, z):
        self.x = x
        self.y = y
        self.z = z
        self.index = None

    def __eq__(self, other):
        return isinstance(other, Point) and self.x == other.x and self.y == other.y and self.z == other.z

    def __hash__(self):
        return hash((self.x, self.y, self.z))

    def __sub__(self, other):
        return np.array([self.x - other.x, self.y - other.y, self.z - other.z])

    def __repr__(self):
        return f"Point({self.x}, {self.y}, {self.z}, Index:{self.index})"

    def movePoint(self, xNew, yNew, zNew):
        """
        Move point at x, y, z
        """
        self.x = xNew
        self.y = yNew
        self.z = zNew

    def getPos(self):
        """
        Return list of node position
        """
        return self.x, self.y, self.z

    def setIndex(self, index):
        """
        Set node index
        """
        self.index = index

    def getData(self):
        """
        Return data structure to export lattice
        """
        return [self.index, self.x, self.y, self.z]
