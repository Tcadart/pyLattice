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


    def changeBeamType(self, newType):
        self.type = newType


    def findPointMod(self, lengthMod):
        beamLength = self.getLength()
        # Direction ratio components
        DR = [(self.point2.x - self.point1.x) / beamLength,
              (self.point2.y - self.point1.y) / beamLength,
              (self.point2.z - self.point1.z) / beamLength]
        # Scaling the direction ratio components by lengthMod
        factor = [dr * lengthMod for dr in DR]
        # Calculating the new point position
        pointMod = [self.point1.x + factor[0], self.point1.y + factor[1], self.point1.z + factor[2]]
        return pointMod