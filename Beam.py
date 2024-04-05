import math


class Beam:
    def __init__(self, point1, point2, Radius,cell_size_x, cell_size_y, cell_size_z, gradRadius, gradMat, posCell,Type):
        """
        Initialize a Beam object.

        :param point1: First end point of the beam (Node object)
        :param point2: Second end point of the beam (Node object)
        :param Radius: Initial beam radius
        :param cell_size_x: Cell size in the x-direction
        :param cell_size_y: Cell size in the y-direction
        :param cell_size_z: Cell size in the z-direction
        :param gradRadius: Gradient of beam radius
        :param gradMat: Gradient of material
        :param posCell: Position of the cell
        :param Type: Center beam = 0; Modified beam = 1
        """
        self.point1 = point1
        self.point2 = point2
        self.Radius = Radius
        self.cell_size_x = cell_size_x
        self.cell_size_y = cell_size_y
        self.cell_size_z = cell_size_z
        self.gradRadius = gradRadius
        self.gradMat = gradMat
        self.posCell = posCell
        self.type = Type
        self.radius = self.setBeamRadius()
        self.length = self.get_length()
        self.volume = self.get_volume()
        self.material = self.setBeamMaterial()
        self.nodes = []

    def add_node(self, node):
        """
        Add a node to the beam's list of connected nodes.

        :param node: Node object to be added
        """
        self.nodes.append(node)

    def get_length(self):
        """
        Calculate the length of the beam.

        :return: Length of the beam
        """
        x1, y1, z1 = self.point1.x, self.point1.y, self.point1.z
        x2, y2, z2 = self.point2.x, self.point2.y, self.point2.z
        length = round(math.sqrt((x2 - x1) ** 2 + (y2 - y1) ** 2 + (z2 - z1) ** 2), 4)
        return length

    def get_volume(self):
        """
        Calculate the volume of the beam.

        :return: Volume of the beam
        """
        return math.pi * (self.radius ** 2) * self.length

    def setBeamRadius(self):
        """
        Calculate and set the beam radius.

        :return: Calculated beam radius
        """
        return self.Radius * self.gradRadius[self.posCell[0]][0] * self.gradRadius[self.posCell[1]][1] * self.gradRadius[self.posCell[2]][2]

    def setBeamMaterial(self):
        """
        Set the material of the beam based on the gradient and position.

        :return: Material index of the beam
        """
        Mat = self.gradMat[self.posCell[2]][self.posCell[1]][self.posCell[0]]
        return Mat
