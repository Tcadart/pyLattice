from Cell import *
import math
import random
import torch
from torch_geometric.data import Data
import matplotlib.pyplot as plt


class Lattice:
    """
    Generate lattice structures with a lot of different parameters
    """

    def __init__(self, cell_size_x: float, cell_size_y: float, cell_size_z: float,
                 num_cells_x: int, num_cells_y: int, num_cells_z: int,
                 Lattice_Type: int, Radius: float,
                 gradRadiusProperty, gradDimProperty, gradMatProperty,
                 simMethod: int = 0, uncertaintyNode: int = 0,
                 hybridLatticeData=None, periodicity: bool = 0):
        """
        Constructor general for the Lattice class.

        Parameter:
        -----------
        cellSizeX: integer
        cellSizeY: integer
        cellSizeZ: integer
            Dimension in each direction of the intial cell in the structure

        num_cells_x: integer
        num_cells_y: integer
        num_cells_z: integer
            Number of cells in each direction in the structure

        Lattice_Type: integer
            Geometry type the cell
                (-2 => Method random cell, -1 => Full random)
                (0 => BCC, 1 => Octet, 2 => OctetExt, 3 => OctetInt, 4 => BCCZ, 5 => Cubic, 6 => OctahedronZ,
                7 => OctahedronZcross, 8 => Kelvin, 9 => Cubic formulation 2 (centered), 10 => Cubic V3, 11 => Cubic V4,
                12 => New lattice (non connu) GPT generated, 13 => Diamond, 14 => Auxetic, 15 => Hichem, 16 => Hybrid1,
                17 => Hybrid2)
        Radius: float
            Initial radius geometry

        Gradient properties
        gradRadiusProperty: array of data as [GradDimRule,GradDimDirection,GradDimParameters]
            Radius gradient on the lattice structure
        gradDimProperty: array of data as [GradRadRule,GradRadDirection,GradRadParameters]
            Cell dimension gradient on the lattice structure
                GradRule => constant, linear, parabolic, sinusoide, exponential
                GradDirection => [bool,bool,bool] set integer to 1 to active gradient in direction [X,Y,Z] 0 inactive
                GradParameters => [float, float, float] variable in the gradient rule for each direction [X,Y,Z]
        gradMatProperty: array of data as [Multimat,GradMaterialDirection]
            Material gradient on the lattice structure
                Multimat => Type of multimaterial (0: inactive / 1: multimat by layer / -1: Full random)


        simMethod: integer (0: off / 1: on)
            Method of simulation with modification at node
        uncertaintyNode: integer (0: off / 1: on)
            Control if adding uncertainties on node position
        hybridLatticeData: array of 3 integer [RadiusOfGeometry1,RadiusOfGeometry2,RadiusOfGeometry3]
            Data of radius of each geometry on the hybrid lattice
        periodicity: boolean (0: off / 1: on)
            Applying periodicity on the outer box of the lattice structure to calculate penalization method
        """
        self.name = None
        self.yMin = None
        self.yMax = None
        self.xMax = None
        self.xMin = None
        self.zMax = None
        self.zMin = None

        self.cellSizeX = cell_size_x
        self.cellSizeY = cell_size_y
        self.cellSizeZ = cell_size_z
        self.numCellsX = num_cells_x
        self.numCellsY = num_cells_y
        self.numCellsZ = num_cells_z
        self.latticeType = Lattice_Type
        self.Radius = Radius
        self.gradRadius = self.getGradSettings(gradRadiusProperty)
        self.gradDim = self.getGradSettings(gradDimProperty)
        self.gradMat = self.gradMaterialSetting(gradMatProperty)
        self.simMethod = simMethod
        self.sizeX, self.sizeY, self.sizeZ = self.getSizeLattice()
        self.uncertaintyNode = uncertaintyNode
        self.hybridLatticeData = hybridLatticeData
        self.periodicity = periodicity  # Not finish to implemented
        self.penalizationCoefficient = 1.5  # Fixed with previous optimization

        self.cells = []
        self._nodes = []
        self._beams = []

        # Process
        self.generateLattice()
        self.getMinMaxValues()
        self.defineBeamNodeIndex()

        self.getAllAngles()

        # Case of penalization at beam near nodes
        if self.simMethod == 1:
            self.getBeamNodeMod()

        # Get some data about lattice structures
        self.getNodeData()
        self.getBeamData()

    @classmethod
    def simpleLattice(cls, cell_size_x, cell_size_y, cell_size_z, num_cells_x, num_cells_y, num_cells_z, Lattice_Type,
                      Radius):
        """
        Generate lattice structures with just simple parameters
        """
        # Define Default gradient properties
        GradDimRule = 'constant'
        GradDimDirection = [1, 0, 0]
        GradDimParameters = [0.0, 0.0, 0.0]
        GradRadRule = 'constant'
        GradRadDirection = [1, 0, 0]
        GradRadParameters = [0.0, 0.0, 0.0]
        Multimat = 0
        GradMaterialDirection = 1
        gradDimProperty = [GradDimRule, GradDimDirection, GradDimParameters]
        gradRadiusProperty = [GradRadRule, GradRadDirection, GradRadParameters]
        gradMatProperty = [Multimat, GradMaterialDirection]
        return cls(cell_size_x, cell_size_y, cell_size_z, num_cells_x, num_cells_y, num_cells_z, Lattice_Type,
                   Radius, gradRadiusProperty, gradDimProperty, gradMatProperty, 0, 0, None)

    @classmethod
    def hybridgeometry(cls, cell_size_x, cell_size_y, cell_size_z, simMethod, uncertaintyNode, hybridLatticeData):
        """
        Generate hybrid geometry structure with just some parameters
        """
        # Define Default gradient properties
        GradDimRule = 'constant'
        GradDimDirection = [1, 0, 0]
        GradDimParameters = [0.0, 0.0, 0.0]
        GradRadRule = 'constant'
        GradRadDirection = [1, 0, 0]
        GradRadParameters = [0.0, 0.0, 0.0]
        Multimat = 0
        GradMaterialDirection = 1
        gradDimProperty = [GradDimRule, GradDimDirection, GradDimParameters]
        gradRadiusProperty = [GradRadRule, GradRadDirection, GradRadParameters]
        gradMatProperty = [Multimat, GradMaterialDirection]
        return cls(cell_size_x, cell_size_y, cell_size_z, 1, 1, 1, 1000,
                   1, gradRadiusProperty, gradDimProperty, gradMatProperty, simMethod, uncertaintyNode,
                   hybridLatticeData, periodicity=1)

    @classmethod
    def latticeHybridForGraph(cls, hybridLatticeData):
        """
        Generate unit cell lattice structure with uniquely hybrid parameter for GNN dataset generation
        """
        cell_size_x = 1
        cell_size_y = 1
        cell_size_z = 1
        simMethod = 0
        uncertaintyNode = 0
        GradDimRule = 'constant'
        GradDimDirection = [1, 0, 0]
        GradDimParameters = [0.0, 0.0, 0.0]
        GradRadRule = 'constant'
        GradRadDirection = [1, 0, 0]
        GradRadParameters = [0.0, 0.0, 0.0]
        Multimat = 0
        GradMaterialDirection = 1
        gradDimProperty = [GradDimRule, GradDimDirection, GradDimParameters]
        gradRadiusProperty = [GradRadRule, GradRadDirection, GradRadParameters]
        gradMatProperty = [Multimat, GradMaterialDirection]
        return cls(cell_size_x, cell_size_y, cell_size_z, 1, 1, 1, 1000,
                   1, gradRadiusProperty, gradDimProperty, gradMatProperty, simMethod, uncertaintyNode,
                   hybridLatticeData)

    @property
    def nodes(self):
        return self._nodes

    @property
    def beams(self):
        return self._beams

    def getSizeLattice(self):
        """
        Computes the size of the lattice along each direction.

        Return:
        ---------
        sizeLattice: list of float in dim 3
            Length of the lattice in each direction
        """
        cellSize = [self.cellSizeX, self.cellSizeY, self.cellSizeZ]
        sizeLattice = [0, 0, 0]
        for direction in range(3):
            for dim in self.gradDim:
                sizeLattice[direction] += dim[direction] * cellSize[direction]
        return sizeLattice

    def getGradSettings(self, gradProperties):
        """
        Generate gradient settings based on the provided rule, direction, and parameters.

        Parameters:
        -----------
        Properties: list[Rule, Direction, Parameters]
            All types of properties for gradient definition

        Return:
        ---------
        gradientData: list[integer]
        Generated gradient settings (list of lists)
        """

        def applyRule(i, numberCell, dirValue, paramValue, rule):
            if i < numberCell:
                if i >= 1 and dirValue == 1:
                    if rule == 'constant':
                        return 1.0
                    elif rule == 'linear':
                        return i * paramValue
                    elif rule == 'parabolic':
                        return i * paramValue if i < numberCell / 2 else (numberCell - i - 1) * paramValue
                    elif rule == 'sinusoide':
                        if i < numberCell / 4:
                            return i * parameters
                        elif i < numberCell / 2:
                            return (numberCell / 2 - i) * parameters
                        elif i == numberCell / 2:
                            return 1.0
                        elif i < 3 / 4 * numberCell:
                            return (3 / 4 * numberCell - i) * (1 / parameters)
                    elif rule == 'exponential':
                        return math.exp(i * paramValue)
                return 1.0
            return 0.0

        # Extract gradient properties
        rule = gradProperties[0]
        direction = gradProperties[1]
        parameters = gradProperties[2]

        # Initialization matrix
        maxCells = max(self.numCellsX, self.numCellsY, self.numCellsZ)
        gradientData = [[0.0, 0.0, 0.0] for _ in range(maxCells)]

        # Processing multiple rules
        for i in range(maxCells):
            numberCells = [self.numCellsX, self.numCellsY, self.numCellsZ]
            for dimIndex in range(3):
                gradientData[i][dimIndex] = applyRule(i, numberCells[dimIndex], direction[dimIndex],
                                                      parameters[dimIndex],
                                                      rule)
        return gradientData

    def gradMaterialSetting(self, gradMatProperty):
        """
        Define gradient material settings

        Parameters:
        ------------
        gradMatProperty: list[Multimat, GradMaterialDirection]
            Set of properties for material gradient

        Returns:
        --------
        gradMat: list of integer
            list of material type in the structure
        """

        # Extract properties
        multimat = gradMatProperty[0]
        direction = gradMatProperty[1]

        if multimat == -1:  # Random
            gradMat = [[[random.randint(1, 3) for X in range(self.numCellsX)] for Y in range(self.numCellsY)] for Z in
                       range(self.numCellsZ)]
        if multimat == 0:  # Mono material
            gradMat = [[[1 for X in range(self.numCellsX)] for Y in range(self.numCellsY)] for Z in
                       range(self.numCellsZ)]
        elif multimat == 1:  # Graded material
            if direction == 1:
                gradMat = [[[X for X in range(self.numCellsX)] for Y in range(self.numCellsY)] for Z in
                           range(self.numCellsZ)]
            if direction == 2:
                gradMat = [[[Y for X in range(self.numCellsX)] for Y in range(self.numCellsY)] for Z in
                           range(self.numCellsZ)]
            if direction == 3:
                gradMat = [[[Z for X in range(self.numCellsX)] for Y in range(self.numCellsY)] for Z in
                           range(self.numCellsZ)]
        return gradMat

    def generateLattice(self):
        """
        Generates lattice structure based on specified parameters.
        """
        xCellStartInit = 0
        yCellStartInit = 0
        zCellStartInit = 0
        xCellStart = 0
        yCellStart = 0
        zCellStart = 0
        posCell = [0, 0, 0]
        for i in range(self.numCellsX):
            if i != 0:
                xCellStart += self.cellSizeX * self.gradDim[posCell[0]][0]
            else:
                xCellStart = xCellStartInit
            for j in range(self.numCellsY):
                if j != 0:
                    yCellStart += self.cellSizeY * self.gradDim[posCell[1]][1]
                else:
                    yCellStart = yCellStartInit
                for k in range(self.numCellsZ):
                    if k != 0:
                        zCellStart += self.cellSizeZ * self.gradDim[posCell[2]][2]
                    else:
                        zCellStart = zCellStartInit
                    posCell = [i, j, k]
                    initialCellSize = [self.cellSizeX, self.cellSizeY, self.cellSizeZ]
                    startCellPos = [xCellStart, yCellStart, zCellStart]
                    new_cell = []
                    # Case for normal lattice structures
                    if self.latticeType != -2 and self.latticeType < 1000:
                        new_cell = Cell(posCell, initialCellSize, startCellPos, self.latticeType, self.Radius,
                                        self.gradRadius, self.gradDim, self.gradMat)
                    # Case for hybrid lattice on 1 cell
                    elif self.latticeType == 1000:
                        latticeHybridType = [0, 16, 17]
                        for idx, radiusHybrid in enumerate(self.hybridLatticeData):
                            if radiusHybrid != 0.0:
                                if not new_cell:
                                    new_cell = Cell(posCell, initialCellSize, startCellPos,
                                                    latticeHybridType[idx], radiusHybrid, self.gradRadius,
                                                    self.gradDim, self.gradMat)
                                    for beam in new_cell.beams:
                                        beam.changeBeamType(idx + 100)
                                else:
                                    new_cell.generateBeamsInCell(latticeHybridType[idx], startCellPos, idx + 100)
                    # Case for randomized lattice
                    # else:
                    #     new_cell.generate_beams_random(self.Radius, self.gradRadius, self.gradDim, self.gradMat,
                    #                                    posCell)
                    self.cells.append(new_cell)
        if self.latticeType == 1000:
            self.checkHybridCollision()

    def defineBeamNodeIndex(self):
        """
        Define index at each beam and node
        """
        beamIndexed = {}
        nodeIndexed = {}
        nextBeamIndex = 0
        nextNodeIndex = 0
        # Define already indexed beam and node
        for cell in self.cells:
            for beam in cell.beams:
                if beam.index is not None:
                    beamIndexed[beam] = nextBeamIndex
                    nextBeamIndex += 1
                for node in [beam.point1, beam.point2]:
                    if node.index is not None:
                        nodeIndexed[node] = nextNodeIndex
                        nextNodeIndex += 1

        # Adding not indexed beam and node
        for cell in self.cells:
            for beam in cell.beams:
                if beam.index is None:
                    if beam not in beamIndexed:
                        beam.setIndex(nextBeamIndex)
                        beamIndexed[beam] = nextBeamIndex
                        nextBeamIndex += 1
                    else:
                        beam.setIndex(beamIndexed[beam])

                for node in [beam.point1, beam.point2]:
                    if node.index is None:
                        if node not in nodeIndexed:
                            node.setIndex(nextNodeIndex)
                            nodeIndexed[node] = nextNodeIndex
                            nextNodeIndex += 1
                        else:
                            node.setIndex(nodeIndexed[node])

    def getNodeData(self):
        """
        Retrieves node data for the lattice.
        data structure: each line represent a node with data [indexNode, X, Y, Z]
        """
        self._nodes = []
        for cell in self.cells:
            for beam in cell.beams:
                for node in [beam.point1, beam.point2]:
                    if node not in self._nodes:
                        self._nodes.append(node.getData())

    def getBeamData(self):
        """
        Retrieves beam data for the lattice.
        data structure: each line represent a beam with data [beamIndex, IndexPoint1, IndexPoint2, beamType]
        """
        self._beams = []
        for cell in self.cells:
            for beam in cell.beams:
                if beam not in self._beams:
                    self._beams.append(beam.getData())

    def visualizeLattice3D(self, beamColor: str = "Material"):
        """
        Visualizes the lattice in 3D using matplotlib.

        Parameter:
        -----------
        beamColor: string
            "Material" -> color by material
            "Type" -> color by type
        """
        fig = plt.figure()
        ax = fig.add_subplot(111, projection='3d')
        ax.set_title("Lattice généré")
        color = ['blue', 'green', 'red', 'yellow', 'orange']
        beamDraw = []
        nodeDraw = []
        for cell in self.cells:
            for beam in cell.beams:
                if beam not in beamDraw:
                    if beamColor == "Material":
                        colorBeam = color[beam.material]
                    elif beamColor == "Type":
                        colorBeam = color[int(str(beam.type)[0])]
                    point1 = beam.point1
                    point2 = beam.point2
                    ax.plot([point1.x, point2.x], [point1.y, point2.y], [point1.z, point2.z], color=colorBeam,
                            markersize=10)
                for node in [beam.point1, beam.point2]:
                    if node not in nodeDraw:
                        ax.scatter(node.x, node.y, node.z, c='black', s=5)

        ax.set_xlabel('X')
        ax.set_ylabel('Y')
        ax.set_zlabel('Z')
        ax.set_xlim3d(self.xMin, self.xMax)
        ax.set_ylim3d(self.yMin, self.yMax)
        ax.set_zlim3d(self.zMin, self.zMax)
        # ax.set_xlim3d(0, self.xMax)
        # ax.set_ylim3d(0, self.yMax)
        # ax.set_zlim3d(0, self.zMax)
        plt.show()

    def visualize_3d_random(self, ax):
        """
        Visualizes the lattice in 3D using matplotlib.

        Parameter:
        -----------
        ax: Axes3D object
        """

        def findColorPoint(x, y, z):
            """
            Define color node with position in lattice structure
            """
            # Coin
            if (x in [0, 1] and y in [0, 1] and z in [0, 1]):
                return 'red'
            elif x > 0 and x < 1 and y > 0 and y < 1 and z > 0 and z < 1:
                return 'black'
            else:
                return 'blue'

        for point in self.nodes_obj:
            x, y, z = point.x, point.y, point.z
            color_point = self.findColorPoint(x, y, z)
            ax.scatter(x, y, z, c=color_point, s=5)
        color = ['blue', 'green', 'red', 'yellow', 'orange']
        for beam in self.beams_obj:
            point1 = beam.point1
            point2 = beam.point2
            ax.plot([point1.x, point2.x], [point1.y, point2.y], [point1.z, point2.z], color=color[beam.material],
                    markersize=10)
        ax.set_xlabel('X')
        ax.set_ylabel('Y')
        ax.set_zlabel('Z')
        ax.set_xlim3d(0, self.xMax)
        ax.set_ylim3d(0, self.yMax)
        ax.set_zlim3d(0, self.zMax)


    def getListAngleBeam(self, beam, pointbeams):
        """
        Calculate angle between the considerate beam and beams contains in pointbeams

        Parameters:
        -----------
        beam: Beam object
            Beam where angle is computed on
        pointbeams: list of Beam object
            List of beam to calculate angle with considered beam

        Return:
        ---------
        non_zero_anglebeam: list of angle between considered beam and pointbeams beam list
        non_zero_radiusbeam: list of radius between considered beam and pointbeams beam list

        Special case when pointbeams is empty return max angle to minimize penalization zone
        """

        def getAngleBetweenBeams(beam1, beam2):
            """
            Calculates angle between 2 beams

            Return:
            --------
            Angle: float
                angle in degrees
            """
            u = beam1.point2 - beam1.point1
            if beam1.point1 == beam2.point1:
                v = beam2.point2 - beam2.point1
            else:
                v = beam2.point1 - beam2.point2
            dot_product = sum(a * b for a, b in zip(u, v))
            u_norm = math.sqrt(sum(a * a for a in u))
            v_norm = math.sqrt(sum(b * b for b in v))
            cos_theta = dot_product / (u_norm * v_norm)
            cos_theta = max(min(cos_theta, 1.0), -1.0)
            angle_rad = math.acos(cos_theta)
            angle_deg = math.degrees(angle_rad)
            return angle_deg

        anglebeam = []
        radiusBeam = []
        if len(pointbeams) > 1:
            for beampoint in pointbeams:
                radiusBeam.append(beampoint.radius)
                anglebeam.append(getAngleBetweenBeams(beam, beampoint))
        else:
            anglebeam.append(179.9)
            radiusBeam.append(beam.radius)
        non_zero_anglebeam = [angle for angle in anglebeam if angle >= 0.01]
        non_zero_radiusbeam = [radius for angle, radius in zip(anglebeam, radiusBeam) if angle >= 0.01]
        return non_zero_anglebeam, non_zero_radiusbeam

    def getConnectedBeams(self, beamList, beam):
        """
        get all beams connected to the interest beam

        Parameters:
        -----------
        beam: Beam object
            Beam of interest

        Returns:
        ---------
        point1beams: list of Beam Object
            list of beam connected to point1 of beam of interest
        point2beams: list of Beam Object
            list of beam connected to point2 of beam of interest
        """
        point1beams = []
        point2beams = []
        for beamidx in beamList:
            if beam.point1 == beamidx.point1 or beam.point1 == beamidx.point2:
                point1beams.append(beamidx)
            if beam.point2 == beamidx.point1 or beam.point2 == beamidx.point2:
                point2beams.append(beamidx)
            if self.periodicity:  # Periodicity is not finish
                tag1 = self.tagPoint(beam.point1)
                tag2 = self.tagPoint(beam.point2)
                tag1 = tag1[0] if len(tag1) == 1 else None
                tag2 = tag2[0] if len(tag2) == 1 else None
                point1_tag = self.tagPoint(beamidx.point1)
                point2_tag = self.tagPoint(beamidx.point2)
                if tag1 is not None and 999 < tag1 < 1008:  # Corner
                    if any(999 < tag < 1008 for tag in point1_tag):
                        point1beams.append(beamidx)
                    if any(999 < tag < 1008 for tag in point2_tag):
                        point1beams.append(beamidx)

                if tag2 is not None and 999 < tag2 < 1008:  # Corner
                    if any(999 < tag < 1008 for tag in point1_tag):
                        point2beams.append(beamidx)
                    if any(999 < tag < 1008 for tag in point2_tag):
                        point2beams.append(beamidx)
        return point1beams, point2beams

    def getAllAngles(self):
        """
        Calculates angles between beams in the lattice.

        Return:
        ---------
        angle:
            data structure => ((beam_index, Angle mininmum point 1, minRad1, Angle mininmum point 2, minRad2))
        """

        def findMinAngle(angles, radii):
            """
            Find Minimum angle between beams and radius connection to this particular beam
            """
            LValuesMax = 0
            LRadius = None
            LAngle = None
            for radius, angle in zip(radii, angles):
                L = self.functionPenalizationLzone(radius, angle)
                if L > LValuesMax:
                    LValuesMax = L
                    LRadius = radius
                    LAngle = angle
            return LAngle, LRadius

        # Create list of beam
        beamList = []
        for cell in self.cells:
            for beam in cell.beams:
                if beam not in beamList:
                    beamList.append(beam)

        angleList = {}
        for beam in beamList:
            # Determine beams on nodes
            point1beams, point2beams = self.getConnectedBeams(beamList, beam)

            # Determine angle for all beams connected at the node
            non_zero_anglebeam1, non_zero_radiusbeam1 = self.getListAngleBeam(beam, point1beams)
            non_zero_anglebeam2, non_zero_radiusbeam2 = self.getListAngleBeam(beam, point2beams)
            # Find the lowest angle
            LAngle1, LRadius1 = findMinAngle(non_zero_anglebeam1, non_zero_radiusbeam1)
            LAngle2, LRadius2 = findMinAngle(non_zero_anglebeam2, non_zero_radiusbeam2)
            angleList[beam.index] = (LRadius1, round(LAngle1, 2), LRadius2, round(LAngle2, 2))
            #BUG

        for cell in self.cells:
            for beam in cell.beams:
                beam.setAngle(angleList[beam.index])

    def getMinMaxValues(self):
        """
        Computes extremum values of coordinates in the lattice.

        Return:
        --------
        ExtrumumValues: tuple of floats (xMin, xMax, yMin, yMax, zMin, zMax)
        """
        if not self.cells:
            raise ValueError("No cells in the lattice.")

        # Flatten the list of nodes from all cells
        all_nodes = [point for cell in self.cells for beam in cell.beams for point in [beam.point1, beam.point2]]

        if not all_nodes:
            raise ValueError("No nodes in the cells of the lattice.")

        # Extract coordinates
        x_values = [node.x for node in all_nodes]
        y_values = [node.y for node in all_nodes]
        z_values = [node.z for node in all_nodes]

        self.xMin, self.xMax = min(x_values), max(x_values)
        self.yMin, self.yMax = min(y_values), max(y_values)
        self.zMin, self.zMax = min(z_values), max(z_values)

    def getBeamNodeMod(self):
        """
        Modifies beam and node data to model lattice structures for simulation with rigidity penalization at node
        """
        for cell in self.cells:
            beamsToRemove = []
            beamToAdd = []

            for beam in cell.beams:
                lengthMod = self.getLengthMod(beam)
                pointExt1 = beam.getPointOnBeamFromDistance(lengthMod[0], 1)
                pointExt1Obj = Point(pointExt1[0], pointExt1[1], pointExt1[2])
                pointExt2 = beam.getPointOnBeamFromDistance(lengthMod[1], 2)
                pointExt2Obj = Point(pointExt2[0], pointExt2[1], pointExt2[2])
                if beam.type == 0:
                    typeBeam = [1, 0, 1]
                else:
                    typeBeam = [beam.type + 100, beam.type, beam.type + 100]

                beamToAdd.append(
                    Beam(beam.point1, pointExt1Obj, beam.radius * self.penalizationCoefficient, beam.material,
                         typeBeam[0]))
                beamToAdd.append(Beam(pointExt1Obj, pointExt2Obj, beam.radius, beam.material, typeBeam[1]))
                beamToAdd.append(
                    Beam(pointExt2Obj, beam.point2, beam.radius * self.penalizationCoefficient, beam.material,
                         typeBeam[2]))

                beamsToRemove.append(beam)

            for addingBeam in beamToAdd:
                cell.addBeam(addingBeam)

            for beam in beamsToRemove:
                cell.removeBeam(beam)

        # Update index
        self.defineBeamNodeIndex()

    def functionPenalizationLzone(self, radius: float, angle: float):
        """
        Definition of the penalization function

        Parameters:
        ------------
        radius: float
            radius data
        angle: float
            angle data

        Returns:
        ---------
        L: float
            Length of the penalization zone
        """
        if angle > 170:
            L = 0.0000001
        else:
            L = radius / math.tan(math.radians(angle) / 2)
        return L

    def getLengthMod(self, beam):
        """
        Calculate and return length to modify in penalization method

        Parameter:
        ----------
        beam: Beam Object
            beam of interest

        Return:
        --------
        LenghtMod: list
            data structure: (Lmod point1,Lmod point2)
        """
        L1 = self.functionPenalizationLzone(beam.angle1[0], beam.angle1[1])
        L2 = self.functionPenalizationLzone(beam.angle2[0], beam.angle2[1])
        return L1, L2

    def removeCell(self, index):
        """
        Removes a cell from the lattice

        Parameters:
        ------------
        index: int
            index of the cell to remove
        """
        if 0 <= index < len(self.cells):
            del self.cells[index]
        else:
            raise IndexError("Invalid cell index.")

    def findMinimumBeamLength(self):
        """
        Find minimum beam length

        Returns:
        --------
        minLength: float
            Length of smallest beam in the lattice
        """
        minLength = 100000
        for cell in self.cells:
            for beam in cell.beams:
                if minLength > beam.getLength() > 0.0001 and (beam.type == 0 or beam.type == 2):
                    minLength = beam.getLength()
        return minLength

    def tagPoint(self, point):
        """
        Define standardized tags for a point

        Parameter:
        ----------
        point: Point Object

        Return:
        --------
        tags: array of integer
            List of tags of the point
        """
        tags = []

        # Faces
        if point.x == self.xMin and (self.yMin < point.y < self.yMax) and (
                self.zMin < point.z < self.zMax):
            tags.append(12)  # Face 1
        elif point.x == self.xMax and (self.yMin < point.y < self.yMax) and (
                self.zMin < point.z < self.zMax):
            tags.append(13)  # Face 2
        elif (self.xMin < point.x < self.xMax) and point.y == self.yMin and (
                self.zMin < point.z < self.zMax):
            tags.append(11)  # Face 3
        elif (self.xMin < point.x < self.xMax) and point.y == self.yMax and (
                self.zMin < point.z < self.zMax):
            tags.append(14)  # Face 4
        elif (self.xMin < point.x < self.xMax) and (
                self.yMin < point.y < self.yMax) and point.z == self.zMin:
            tags.append(10)  # Face 5
        elif (self.xMin < point.x < self.xMax) and (
                self.yMin < point.y < self.yMax) and point.z == self.zMax:
            tags.append(15)  # Face 6

        # Edge
        if point.x == self.xMin and point.y == self.yMin and (self.zMin < point.z < self.zMax):
            tags.append(102)  # Edge 0
        elif (self.xMin < point.x < self.xMax) and point.y == self.yMin and point.z == self.zMin:
            tags.append(100)  # Edge 1
        elif point.x == self.xMax and point.y == self.yMin and (self.zMin < point.z < self.zMax):
            tags.append(104)  # Edge 2
        elif (self.xMin < point.x < self.xMax) and point.y == self.yMin and point.z == self.zMax:
            tags.append(108)  # Edge 3
        elif point.x == self.xMin and (self.yMin < point.y < self.yMax) and point.z == self.zMin:
            tags.append(101)  # Edge 4
        elif point.x == self.xMax and (self.yMin < point.y < self.yMax) and point.z == self.zMin:
            tags.append(103)  # Edge 5
        elif point.x == self.xMin and point.y == self.yMax and (self.zMin < point.z < self.zMax):
            tags.append(106)  # Edge 6
        elif (self.xMin < point.x < self.xMax) and point.y == self.yMax and point.z == self.zMin:
            tags.append(105)  # Edge 7
        elif point.x == self.xMax and point.y == self.yMax and (self.zMin < point.z < self.zMax):
            tags.append(107)  # Edge 8
        elif (self.xMin < point.x < self.xMax) and point.y == self.yMax and point.z == self.zMax:
            tags.append(111)  # Edge 9
        elif point.x == self.xMin and (self.yMin < point.y < self.yMax) and point.z == self.zMax:
            tags.append(109)  # Edge 10
        elif point.x == self.xMax and (self.yMin < point.y < self.yMax) and point.z == self.zMax:
            tags.append(110)  # Edge 11

        # Corner
        if point.x == self.xMin and point.y == self.yMin and point.z == self.zMin:
            tags.append(1000)  # Corner 0
        elif point.x == self.xMax and point.y == self.yMin and point.z == self.zMin:
            tags.append(1001)  # Corner 1
        elif point.x == self.xMin and point.y == self.yMax and point.z == self.zMin:
            tags.append(1002)  # Corner 2
        elif point.x == self.xMax and point.y == self.yMax and point.z == self.zMin:
            tags.append(1003)  # Corner 3
        elif point.x == self.xMin and point.y == self.yMin and point.z == self.zMax:
            tags.append(1004)  # Corner 4
        elif point.x == self.xMax and point.y == self.yMin and point.z == self.zMax:
            tags.append(1005)  # Corner 5
        elif point.x == self.xMin and point.y == self.yMax and point.z == self.zMax:
            tags.append(1006)  # Corner 6
        elif point.x == self.xMax and point.y == self.yMax and point.z == self.zMax:
            tags.append(1007)  # Corner 7

        return tags

    def getTagList(self):
        """
        Get the tag for all points in lattice

        Returns:
        --------
        tagList: list of int
            List of all tags of each point in lattice
        """
        tagList = []
        nodeAlreadyAdded = []
        for cell in self.cells:
            for beam in cell.beams:
                for node in [beam.point1, beam.point2]:
                    if node not in nodeAlreadyAdded:
                        tagList.append(self.tagPoint(node))
                        nodeAlreadyAdded.append(node)
        return tagList

    def getConnectedNode(self, node):
        """
        Get all nodes connected to the input node with a beam

        Parameter:
        -----------
        node: point object

        Return:
        --------
        connectedNode: List of point object
        """
        connectedNode = []
        nodeIndexRef = node.index
        for cell in self.cells:
            for beam in cell.beams:
                if beam.point1.index == nodeIndexRef:
                    connectedNode.append(beam.point2)
                if beam.point2.index == nodeIndexRef:
                    connectedNode.append(beam.point1)
        return connectedNode

    def toucanLatticeModifier(self):
        """
        Fun lattice bio inspired toucan
        """
        toucanPourcent = 10
        for i in range(int(math.ceil(len(self.beams_obj) * toucanPourcent / 100))):
            nodeInit = random.choice(self.nodes_obj)
            node1 = random.choice(self.getConnectedNode(nodeInit))
            node2 = node1
            while (node2 == node1):
                node2 = random.choice(self.getConnectedNode(node1))
            self.toucanModifier.append(
                [[nodeInit.x, nodeInit.y, nodeInit.z], [node1.x, node1.y, node1.z], [node2.x, node2.y, node2.z]])
        return self.toucanModifier

    def findBoundaryBeams(self):
        """
        Find boundary beams and change type of beam

        Return:
        -------
        boundaryBeams: List of beam object
        """
        boundaryBeams = []
        for cell in self.cells:
            for beam in cell.beams:
                if self.isNodeOnBoundary(beam.point1) or self.isNodeOnBoundary(beam.point2):
                    beam.changeBeamType(2)
                    boundaryBeams.append(beam)
        return boundaryBeams

    def isNodeOnBoundary(self, node):
        return (node.x == self.xMin or node.x == self.xMax or node.y == self.yMin or node.y == self.yMax or
                            node.z == self.zMin or node.z == self.zMax)
    def findBoundaryNodes(self):
        """
        Find boundary nodes

        Returns:
        ---------
        boundaryNodes: List of point object
        """
        boundaryNodes = []
        for cell in self.cells:
            for beam in cell.beams:
                for node in [beam.point1, beam.point2]:
                    if (node.x == self.xMin or node.x == self.xMax or node.y == self.yMin or node.y == self.yMax or
                            node.z == self.zMin or node.z == self.zMax):
                        boundaryNodes.append(node)
        return boundaryNodes

    def getName(self):
        """
        Determine name of the lattice

        Returns:
        ---------
        name: string
        """
        nameList = [
            "BCC",  # 0
            "Octet",  # 1
            "OctetExt",  # 2
            "OctetInt",  # 3
            "BCCZ",  # 4
            "Cubic",  # 5
            "OctahedronZ",  # 6
            "OctahedronZcross",  # 7
            "Kelvin",  # 8
            "CubicV2",  # 9 (centered)
            "CubicV3",  # 10
            "CubicV4",  # 11
            "NewlatticeUnknown",  # 12 GPT generated
            "Diamond",  # 13
            "Auxetic"  # 14
        ]
        if self.latticeType == 1000:
            self.name = "Hybrid"
        else:
            self.name = nameList[self.latticeType]
        if self.simMethod == 1:
            self.name += "Mod"
        return self.name

    def checkHybridCollision(self):
        """
        Check if beam in hybrid configuration is cut by a point in the geometry
        Change the beam configuration of collisionned beams
        """
        for cell in self.cells:
            cellPoints = cell.getAllPoints()
            for node in cellPoints:
                for beam in cell.beams:
                    if beam.isPointOnBeam(node):
                        typeBeamToRemove = beam.type  # Get beam to remove type to apply in new separated beams
                        beam1 = Beam(beam.point1, node, beam.radius, beam.material, typeBeamToRemove)
                        beam2 = Beam(beam.point2, node, beam.radius, beam.material, typeBeamToRemove)
                        cell.removeBeam(beam)
                        cell.addBeam(beam1)
                        cell.addBeam(beam2)

    def getPosData(self):
        """
        Retrieves node position data for the lattice
        data structure : each line represent a node with data [X, Y, Z]
        """
        posData = []
        nodeAlreadyAdded = []
        for cell in self.cells:
            for beam in cell.beams:
                for node in [beam.point1, beam.point2]:
                    if node not in nodeAlreadyAdded:
                        posData.append([node.getPos()])
                        nodeAlreadyAdded.append(node)
        return posData

    def getEdgeIndex(self):
        """
        Retrieves edge index data for the lattice.
        data structure: each line represent a beam with data [IndexPoint1, IndexPoint2]
        """
        edgeIndex = []
        beamAlreadyAdded = []
        for cell in self.cells:
            for beam in cell.beams:
                if beam not in beamAlreadyAdded:
                    edgeIndex.append([beam.point1.index, beam.point2.index])
                    beamAlreadyAdded.append(beam)
        return edgeIndex

    def getBeamType(self): # A modifier
        """
        Retrieves beam type data for the lattice.
        data structure: each line is the type of the beam with index the line index
        """
        BeamType = []
        for index, beam in enumerate(self.beams_obj):
            BeamType.append([beam.type])
        return BeamType

    def changeHybridData(self, hybridRadiusData):
        """
        Change radius data of the hybrid lattice to efficiently generate graph from training neural network

        Parameters:
        ------------
        hybridRadiusData: array of dim 3
            Array of radius data of hybrid lattice
        """
        radiusDict = {
            100: hybridRadiusData[0],
            101: hybridRadiusData[1],
            102: hybridRadiusData[2],
        }
        for cell in self.cells:
            for beam in cell.beams:
                if beam.type in radiusDict:
                    beam.radius = radiusDict[beam.type]

    def attractorLattice(self, ):
        def distance(point1, point2):
            """
            Calculate distance between two points
            """
            return math.sqrt((point2.x - point1.x) ** 2 + (point2.y - point1.y) ** 2 + (point2.z - point1.z) ** 2)

        def movePointAttracted(point1, attractorPoint, alpha):
            """
            Move point1 relative from attractorPoint with coefficient alpha
            """
            Length = distance(point1, attractorPoint)
            DR = [(attractorPoint.x - point1.x) / Length, (attractorPoint.y - point1.y) / Length,
                  (attractorPoint.z - point1.z) / Length]
            factor = [(alpha) * dr for dr in DR]
            pointMod = [point1.x, point1.y, point1.z]
            pointMod = [p1 + p2 for p1, p2 in zip(pointMod, factor)]
            point1.movePoint(pointMod[0], pointMod[1], pointMod[2])

        pointAttractor = Point(5, 0.5, -2)
        alpha = 0.5
        for cell in self.cells:
            for beam in cell.beams:
                movePointAttracted(beam.point1, pointAttractor, alpha)
                movePointAttracted(beam.point2, pointAttractor, alpha)
        self.getMinMaxValues()
