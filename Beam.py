import math


class Beam:
    """
    Class Beam represente beam element in structures
    """

    def __init__(self, point1, point2, Radius: float, Material: int, Type: int):
        """
        Beam object represent a beam by 2 points a radius and a beam type

        Parameters:
        ------------
        point1: Point object
        point2: Point object
        Radius: float
        Type: int
            Center beam = 0; Modified beam = 1
        """
        self.point1 = point1
        self.point2 = point2
        self.radius = Radius
        self.material = Material

        # type = 0 => normal
        # type = 1 => beam mod
        # type = 2 => beam on boundary
        self.type = Type
        self.index = None
        self.angle1 = None
        self.angle2 = None

    def __repr__(self):
        return f"Beam({self.point1}, {self.point2}, Radius:{self.radius}, Type:{self.type}, Index:{self.index})"

    def __eq__(self, other):
        return isinstance(other, Beam) and self.point1 == other.point1 and self.point2 == other.point2

    def __hash__(self):
        return hash((self.point1, self.point2))

    def getLength(self):
        """
        Calculate the length of the beam.

        :return: Length of the beam
        """
        x1, y1, z1 = self.point1.x, self.point1.y, self.point1.z
        x2, y2, z2 = self.point2.x, self.point2.y, self.point2.z
        length = round(math.sqrt((x2 - x1) ** 2 + (y2 - y1) ** 2 + (z2 - z1) ** 2), 4)
        return length

    def getVolume(self):
        """
        Calculate the volume of the beam.

        :return: Volume of the beam
        """
        return math.pi * (self.radius ** 2) * self.length

    def changeBeamType(self, newType: int):
        """
        Change beam type

        Parameters:
        newtype: integer
            beam type wanted to assign
        """
        self.type = newType

    def getPointOnBeamFromDistance(self, distance, pointIndex):
        """
        Get the position [x,y,z] of a point on the beam at a distance of the point pointIndex

        Parameters:
        -----------
        distance: float
            Distance between pointIndex and new point
        pointIndex: int (1 or 2)
            Index of the beam point
        """
        beam_length = self.getLength()
        if pointIndex == 1:
            start_point = self.point1
            end_point = self.point2
        elif pointIndex == 2:
            start_point = self.point2
            end_point = self.point1
        else:
            raise ValueError("Point must be 1 or 2.")

        direction_ratios = [
            (end_point.x - start_point.x) / beam_length,
            (end_point.y - start_point.y) / beam_length,
            (end_point.z - start_point.z) / beam_length,
        ]

        factors = [dr * distance for dr in direction_ratios]

        point_mod = [
            start_point.x + factors[0],
            start_point.y + factors[1],
            start_point.z + factors[2]
        ]

        return point_mod

    def isPointOnBeam(self, node):
        """
        Find if input node is on the beam

        Return
        -------
        boolean: True => Point on line
        """
        vector1 = (self.point2.x - self.point1.x, self.point2.y - self.point1.y, self.point2.z - self.point1.z)
        vector2 = (node.x - self.point1.x, node.y - self.point1.y, node.z - self.point1.z)

        if (node.x == self.point1.x and node.y == self.point1.y and node.z == self.point1.z) or (
                node.x == self.point2.x and node.y == self.point2.y and node.z == self.point2.z):
            return False
        cross_product = (
            vector1[1] * vector2[2] - vector1[2] * vector2[1],
            vector1[2] * vector2[0] - vector1[0] * vector2[2],
            vector1[0] * vector2[1] - vector1[1] * vector2[0]
        )

        if cross_product == (0, 0, 0):
            dot_product = (vector2[0] * vector1[0] + vector2[1] * vector1[1] + vector2[2] * vector1[2])
            vector1_length_squared = (vector1[0] ** 2 + vector1[1] ** 2 + vector1[2] ** 2)
            return 0 <= dot_product <= vector1_length_squared
        else:
            return False

    def setIndex(self, index):
        """
        Set beam index
        """
        self.index = index

    def getData(self):
        """
        Return data structure to export lattice
        """
        return [self.index, self.point1.index, self.point2.index, self.type]

    def setAngle(self, AngleData):
        """
        Set angle data to beam
        """
        self.angle1 = AngleData[0:2]
        self.angle2 = AngleData[2:4]
