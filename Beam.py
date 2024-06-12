import math


class Beam:
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



    def __repr__(self):
        return f"Beam({self.point1}, {self.point2}, Radius:{self.radius}, Type:{self.type})"

    def __eq__(self, other):
        return isinstance(other, Beam) and self.point1 == other.point1 and self.point2 == other.point2


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