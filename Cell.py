from Point import *
from Beam import *
import math
import random
from Geometry_Lattice import Lattice_geometry


class Cell(object):
    """
    Define Cell data for lattice structure
    """

    def __init__(self, posCell, initialCellSize, startCellPos, latticeType,
                 Radius, gradRadius, gradDim, gradMat):
        """
        Initialize a Cell with its dimensions and position

        Parameters:
        -----------
        posCell: list
            Position of the cell in the lattice
        initialCellSize: list
            Initial size of the cell
        startCellPos: list
            Position of the start of the cell
        latticeType: int
            Type of lattice geometry
        Radius: float
            Base radius of the beam
        gradRadius: list
            Gradient of the radius
        gradDim: list
            Gradient of the dimensions
        gradMat: list
            Gradient of the material
        """
        self.centerPoint = None
        self._beamMaterial = None
        self._beamRadius = None
        self.cellSize = None
        self.posCell = posCell
        self.coordinateCell = startCellPos
        self.beams = []
        self.index = None
        self.latticeType = latticeType
        self.hybridRadius = None

        self.getBeamMaterial(gradMat)
        self.getBeamRadius(gradRadius, Radius)
        self.getCellSize(initialCellSize, gradDim)
        self.generateBeamsInCell(latticeType, startCellPos)
        self.getCellCenter(initialCellSize)

    def generateBeamsInCell(self, latticeType, startCellPos, beamType=0):
        """
        Generate beams and nodes using a given lattice type and parameters.

        Parameters:
        -----------
        latticeType: int
            Type of lattice geometry
        startCellPos: list
            Position of the start of the cell
        beamType: int
            Type of beam
        """
        for line in Lattice_geometry(latticeType):
            x1, y1, z1, x2, y2, z2 = map(float, line)
            point1 = Point(x1 * self.cellSize[0] + startCellPos[0], y1 * self.cellSize[1] + startCellPos[1],
                           z1 * self.cellSize[2] + startCellPos[2])
            point2 = Point(x2 * self.cellSize[0] + startCellPos[0], y2 * self.cellSize[1] + startCellPos[1],
                           z2 * self.cellSize[2] + startCellPos[2])
            beam = Beam(point1, point2, self._beamRadius, self._beamMaterial, beamType)
            self.beams.append(beam)

    def getBeamMaterial(self, gradMat):
        """
        Get the material of the beam based on the gradient and position.

        Parameters:
        -----------
        gradMat: list
            Gradient of the material

        Returns:
        ---------
        materialType: int
            Material index of the beam
        """
        self._beamMaterial = gradMat[self.posCell[2]][self.posCell[1]][self.posCell[0]]

    def getBeamRadius(self, gradRadius, BaseRadius):
        """
        Calculate and return the beam radius

        Parameters:
        -----------
        gradRadius: list
            Gradient of the radius
        BaseRadius: float
            Base radius of the beam

        Returns:
        ---------
        actualBeamRadius: float
            Calculated beam radius
        """
        self._beamRadius = BaseRadius * gradRadius[self.posCell[0]][0] * gradRadius[self.posCell[1]][1] * \
                           gradRadius[self.posCell[2]][2]

    def getCellSize(self, initialCellSize, gradDim):
        """
        Calculate and return the cell size

        Parameters:
        -----------
        initialCellSize: 3-array
            Dimension of the initial cell without modification
        gradDim:

        Returns:
        ---------
        cellSize : float
            Calculated beam radius
        """
        self.cellSize = [initial_size * gradDim[pos][i] for i, (initial_size, pos) in
                         enumerate(zip(initialCellSize, self.posCell))]

    def getCellCenter(self, initialCellSize):
        """
        Calculate the center point of the cell
        """
        self.centerPoint = [initialCellSize[i] + self.cellSize[i] / 2 for i in range(3)]

    def getAllPoints(self):
        """
        Determine list of points in cell
        """
        pointList = []
        for beam in self.beams:
            for point in [beam.point1, beam.point2]:
                if point not in pointList:
                    pointList.append(point)
        return pointList

    def removeBeam(self, beamToDelete):
        """
        Removes a beam from the lattice

        Parameters:
        ------------
        beamToDelete: beam Object
            Beam to remove
        """
        try:
            self.beams.remove(beamToDelete)
        except ValueError:
            print("Beam not found in the list")

    def addBeam(self, beamToAdd):
        """
        Adding beam to cell
        """
        self.beams.append(beamToAdd)

    def setIndex(self, index):
        """
        Set cell index
        """
        self.index = index

    def random_coordinate(self, coord, mu, sigma):
        """
        Randomize coordinate of a node with a gaussian noise of parameter mu and sigma
        """
        mod_coord = []
        for pos in coord:
            # mod_coord.append(pos + random.uniform(0, 0.5) - 0.25)
            mod_coord.append(pos + random.gauss(mu, sigma))
        return mod_coord

    def add_point(self, point):
        x, y, z = map(float, point)
        point_obj = Point((x) * self.cellSizeX + self.x, (y) * self.cellSizeY + self.y,
                          (z) * self.cellSizeZ + self.z)
        self.nodes.append(point_obj)

    def generate_beams_random(self, Radius, gradRadius, gradDim, gradMat, posCell):
        self.beams = []
        self.nodes = []
        # Corner node
        corner_node = random.randint(0, 1)
        if corner_node == 1:  # Corner nodes
            map_corner = [(0.0, 0.0, 0.0),
                          (1.0, 1.0, 1.0),
                          (1.0, 1.0, 0.0),
                          (0.0, 0.0, 1.0),
                          (0.0, 1.0, 0.0),
                          (0.0, 1.0, 1.0),
                          (1.0, 0.0, 1.0),
                          (1.0, 0.0, 0.0)]
            for idx, point_corner in enumerate(map_corner):
                self.add_point(point_corner)
        # Edge
        map_edge = [(0.0, 0.0), (0.0, 1.0), (1.0, 0.0), (1.0, 1.0)]
        for i in range(3):  # 3 direction of edge node
            Edge_node = random.randint(0, 1)
            if Edge_node == 1:
                point_mod = 0.5 + random.uniform(0, 0.5) - 0.25
                for idx, point_edge in enumerate(map_edge):
                    point = list(point_edge)
                    point.insert(i, point_mod)
                    self.add_point(point)
        # Face
        map_face = [(0.25, 0.25), (0.75, 0.25), (0.5, 0.5), (0.25, 0.75), (0.75, 0.75)]
        for i in range(3):  # 3 direction of face node
            Face_node = [random.randint(0, 1) for _ in range(5)]
            for idx, point_face in enumerate(map_face):
                if Face_node[idx] == 1:
                    point = list(self.random_coordinate(map_face[idx], 0, 0.25))
                    point_sym = list(point)
                    point.insert(i, 0.0)
                    point_sym.insert(i, 1.0)
                    self.add_point(point)
                    self.add_point(point_sym)
        # Interior
        map_interior = [(0.25, 0.25, 0.25), (0.5, 0.25, 0.25), (0.75, 0.25, 0.25),
                        (0.25, 0.5, 0.25), (0.5, 0.5, 0.25), (0.75, 0.5, 0.25),
                        (0.25, 0.75, 0.25), (0.5, 0.75, 0.25), (0.75, 0.75, 0.25),

                        (0.25, 0.25, 0.5), (0.5, 0.25, 0.5), (0.75, 0.25, 0.5),
                        (0.25, 0.5, 0.5), (0.5, 0.5, 0.5), (0.75, 0.5, 0.5),
                        (0.25, 0.75, 0.5), (0.5, 0.75, 0.5), (0.75, 0.75, 0.5),

                        (0.25, 0.25, 0.75), (0.5, 0.25, 0.75), (0.75, 0.25, 0.75),
                        (0.25, 0.5, 0.75), (0.5, 0.5, 0.75), (0.75, 0.5, 0.75),
                        (0.25, 0.75, 0.75), (0.5, 0.75, 0.75), (0.75, 0.75, 0.75)]
        interior_node = [random.randint(0, 1) for _ in range(len(map_interior))]
        for idx, point_interior in enumerate(map_interior):
            if interior_node[idx] == 1:
                point = list(self.random_coordinate(map_interior[idx], 0, 0.1))
                self.add_point(point)

        self.remove_nodes_outside_unit_cube()
        print(len(self.nodes))

        self.remove_nodes_too_close_advanced()
        print(len(self.nodes))

        # Generate beam at least 2 beams per node
        for beam in range(1):
            beam = Beam(self.nodes[random.randint(0, len(self.nodes) - 1)],
                        self.nodes[random.randint(0, len(self.nodes) - 1)], beamRadius,
                        beamMaterial, 0)
            if self.beam_already_exist(beam):
                self.beams.append(beam)
        return self.nodes, self.beams

    def beam_already_exist(self, beam_test):
        # print(beam_test.point1, beam_test.point2)
        for beam in self.beams:
            if beam.point1 == beam_test.point1:
                if beam.point2 == beam_test.point2:
                    # print(beam.point1,beam.point2)
                    # print(False)
                    return False
            if beam.point1 == beam_test.point2:
                if beam.point2 == beam_test.point1:
                    # print(beam.point1,beam.point2)
                    # print(False)
                    return False
        # print(True)
        return True

    def remove_nodes_outside_unit_cube(self):
        """
        Remove nodes that are outside of the unit cube defined by the cell dimensions.
        """
        inside_nodes = []
        for point in self.nodes:
            if 0 <= point.x <= self.cellSizeX and 0 <= point.y <= self.cellSizeY and 0 <= point.z <= self.cellSizeZ:
                inside_nodes.append(point)
        self.nodes = inside_nodes

    def is_node_in_corner(self, node):
        """
        Check if a node is in a corner of the unit cube.

        :param node: The node to check.
        :return: True if the node is in a corner, False otherwise.
        """
        return ((node.x in [0, self.cellSizeX]) and
                (node.y in [0, self.cellSizeY]) and
                (node.z in [0, self.cellSizeZ]))

    def remove_nodes_too_close_advanced(self, min_distance=0.1):
        """
        Remove nodes that are too close together with a refined decision process.
        """
        to_remove = set()

        def is_surface_not_corner(node):
            return self.is_node_on_surface(node) and not self.is_node_in_corner(node)

        for i, node_i in enumerate(self.nodes):
            if i in to_remove:
                continue

            for j, node_j in enumerate(self.nodes[i + 1:], start=i + 1):
                if j in to_remove:
                    continue

                # Calculate Euclidean distance
                distance = math.sqrt((node_i.x - node_j.x) ** 2 +
                                     (node_i.y - node_j.y) ** 2 +
                                     (node_i.z - node_j.z) ** 2)

                if distance < min_distance:

                    if self.is_node_in_corner(node_i):
                        to_remove.add(j)
                        continue
                    elif self.is_node_in_corner(node_j):
                        to_remove.add(i)
                        continue

                    if self.is_node_on_surface(node_i) and not self.is_node_on_surface(node_j):
                        to_remove.add(i)
                    elif not self.is_node_on_surface(node_i) and self.is_node_on_surface(node_j):
                        to_remove.add(j)
                    else:

                        to_remove.add(j)  # Choix arbitraire pour supprimer
                        if is_surface_not_corner(node_j):
                            symmetric_node = self.find_symmetric_node(node_j)
                            if symmetric_node and not self.is_node_in_corner(symmetric_node):
                                symmetric_index = self.nodes.index(symmetric_node)
                                to_remove.add(symmetric_index)

        self.nodes = [node for index, node in enumerate(self.nodes) if index not in to_remove]

    def is_node_on_surface(self, node):
        """
        Check if a node is on the surface of the unit cube.

        :param node: The node to check.
        :return: True if the node is on the surface, False otherwise.
        """
        return (node.x in [0, self.cellSizeX] or
                node.y in [0, self.cellSizeY] or
                node.z in [0, self.cellSizeZ])

    def find_symmetric_node(self, node):
        """
        Find a symmetric node to a given node, based on the cube's faces. If a node is on the face X=0,
        its symmetric is on the face X=cellSizeX, and analogously for Y and Z axes.

        :param node: The node for which to find a symmetric counterpart.
        :return: The symmetric node, if found.
        """

        symmetric_x = self.cellSizeX - node.x if node.x in [0, self.cellSizeX] else node.x
        symmetric_y = self.cellSizeY - node.y if node.y in [0, self.cellSizeY] else node.y
        symmetric_z = self.cellSizeZ - node.z if node.z in [0, self.cellSizeZ] else node.z

        for n in self.nodes:
            if n.x == symmetric_x and n.y == symmetric_y and n.z == symmetric_z:
                return n
        return None

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
        color = ['blue', 'green', 'black', 'yellow', 'orange']
        for index, beam in enumerate(self.beams):
            point1 = beam.point1
            point2 = beam.point2
            ax.plot([point1.x, point2.x], [point1.y, point2.y], [point1.z, point2.z], color=color[beam.material])

    def getPointOnSurface(self, surfaceName):
        """
        Get the points on the surface specified in the global reference frame.

        Parameters:
        -----------
        surfaceName: str
            Name of the surface. Choose from 'Xmin', 'Xmax', 'Ymin', 'Ymax', 'Zmin', or 'Zmax'.

        Returns:
        --------
        list
           List of points on the specified surface.
        """
        surface_map = {
            "Xmin": self.coordinateCell[0],
            "Xmax": self.coordinateCell[0] + self.cellSize[0],
            "Ymin": self.coordinateCell[1],
            "Ymax": self.coordinateCell[1] + self.cellSize[1],
            "Zmin": self.coordinateCell[2],
            "Zmax": self.coordinateCell[2] + self.cellSize[2]
        }

        if surfaceName not in surface_map:
            raise ValueError(
                f"Surface '{surfaceName}' is not valid. Choose from 'Xmin', 'Xmax', 'Ymin', 'Ymax', 'Zmin', or 'Zmax'.")

        surface_value = surface_map[surfaceName]
        coord_index = {"Xmin": "x", "Xmax": "x", "Ymin": "y", "Ymax": "y", "Zmin": "z", "Zmax": "z"}

        return [point for beam in self.beams for point in [beam.point1, beam.point2] if
                getattr(point, coord_index[surfaceName]) == surface_value]

    def defineHybridRadius(self, hybridRadius):
        """
        Define hybrid radius for the cell

        Parameters:
        -----------
        hybridRadius: list of float dim 3
            Hybrid radius of the cell
        """
        self.hybridRadius = hybridRadius

    def getHybridRadius(self):
        """
        Return hybrid radius of the cell
        """
        return self.hybridRadius

    def getNodeOrderToSimulate(self):
        """
        Get the order of nodes to simulate in the cell
        """
        originalTags = [1000, 1001, 1002, 1003, 1004, 1005, 1006, 1007, 100, 101, 102, 103, 104, 105, 106, 107,
                        108, 109, 110, 111]
        tag_dict = {tag: None for tag in originalTags}
        for beam in self.beams:
            for point in [beam.point1, beam.point2]:
                if point.indexBoundary is not None:
                    tag = point.tagPoint(self.coordinateCell[0], self.coordinateCell[0] + self.cellSize[0],
                                         self.coordinateCell[1], self.coordinateCell[1] + self.cellSize[1],
                                         self.coordinateCell[2], self.coordinateCell[2] + self.cellSize[2])
                    if tag:  # Ensure tags is not an empty list
                        tag = tag[0]  # Take the first tag from the list
                        if tag in originalTags:
                            tag_dict[tag] = point
        return tag_dict

    def setDisplacementAtBoundaryNodes(self, displacementArray):
        """
        Set displacement at nodes.

        Parameters:
        ------------
        displacementArray: list or array-like
            Flattened array of displacement values.
        """
        for beam in self.beams:
            for point in [beam.point1, beam.point2]:
                if point.indexBoundary is not None:
                    index = point.indexBoundary * 6
                    displacement = displacementArray[index:index + 6]

                    # Filtrer les déplacements en fonction de fixedDOF
                    filtered_displacement = [
                        displacement[i] if point.fixedDOF[i] == 0 else 0
                        for i in range(6)
                    ]

                    point.setDisplacementValue(filtered_displacement)

    def getDisplacementAtBoundaryNodes(self, nodeList):
        """
        Get the displacement at nodes.

        Parameters:
        -----------
        nodeList: list of Point objects
            List of nodes to get the displacement.

        Returns:
        --------
        list
            A flattened list of displacement values.
        """
        displacementList = []
        for node in nodeList.values():
            if node:
                displacement = node.getDisplacementValue()
                displacementList.extend(displacement)  # Extend the list with the displacement values
            else:
                displacementList.extend([0, 0, 0, 0, 0, 0])  # Append zeros if the node is not found
        return displacementList

    def setReactionForceOnEachNodes(self, nodeList, reactionForce):
        """
        Set reaction force on each nodes.

        Parameters:
        -----------
        nodeList: list of Point objects
            List of nodes to set the reaction force.
        reactionForce: list
            List of reaction force values.
        """
        tagAlreadySet = []
        tagList = list(nodeList.keys())
        for beam in self.beams:
            for point in [beam.point1, beam.point2]:
                for tag, node in nodeList.items():
                    if node == point and tag not in tagAlreadySet:
                        point.setReactionForce(reactionForce[tagList.index(tag)])
                        tagAlreadySet.append(tag)
                        break
