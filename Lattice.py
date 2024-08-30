from __future__ import print_function, division

import numpy as np

from Cell import *
import math
import random
import sys
if sys.version_info[0] == 3:
    import matplotlib.pyplot as plt
    from mpl_toolkits.mplot3d.art3d import Line3DCollection


class Lattice(object):
    """
    Generate lattice structures with a lot of different parameters
    """

    def __init__(self, cell_size_x, cell_size_y, cell_size_z,
                 num_cells_x, num_cells_y, num_cells_z,
                 Lattice_Type, Radius,
                 gradRadiusProperty, gradDimProperty, gradMatProperty,
                 simMethod=0, uncertaintyNode=0,
                 hybridLatticeData=None, hybridGeomType=None, periodicity=0, erasedParts=None):
        """
        Constructor general for the Lattice class.

        Parameter:
        -----------
        cellSizeX: float
        cellSizeY: float
        cellSizeZ: float
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
        erasedParts: list of float in dim 6
            (xStart, yStart, zStart, xDim, yDim, zDim) of the erased region
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
        self.hybridGeomType = hybridGeomType
        self.periodicity = periodicity  # Not finish to implemented
        self.penalizationCoefficient = 1.5  # Fixed with previous optimization
        self.erasedParts = erasedParts

        self.cells = []
        self._nodes = []
        self._beams = []

        # Process
        self.generateLattice()
        self.getMinMaxValues()
        self.defineBeamNodeIndex()
        self.defineCellIndex()

        self.applyTagToAllPoint()
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
        return cls(cell_size_x, cell_size_y, cell_size_z, num_cells_x, num_cells_y, num_cells_z, Lattice_Type, Radius,
                   gradRadiusProperty, gradDimProperty, gradMatProperty)

    @classmethod
    def hybridgeometry(cls, cell_size_x, cell_size_y, cell_size_z, simMethod, uncertaintyNode, hybridLatticeData,
                       hybridGeomType=None):
        """
        Generate hybrid geometry structure with just some parameters
        """
        # Define Default gradient properties
        if hybridGeomType is None:
            hybridGeomType = [0, 16, 17]

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
                   hybridLatticeData, hybridGeomType=hybridGeomType, periodicity=True)

    @classmethod
    def latticeHybridForGraph(cls, hybridLatticeData, hybridGeomType):
        """
        Generate unit cell lattice structure with uniquely hybrid parameter for GNN dataset generation

        Parameters:
        -----------
        hybridLatticeData: list of float
            List of radius for each geometry
        hybridGeomType: list of integer
            List of geometry type for each radius
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
                   hybridLatticeData, hybridGeomType=hybridGeomType)

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
        gradProperties: list[Rule, Direction, Parameters]
            All types of properties for gradient definition

        Return:
        ---------
        gradientData: list[list[float]]
            Generated gradient settings (list of lists)
        """

        def apply_constant_rule(i, paramValue):
            return 1.0

        def apply_linear_rule(i, paramValue):
            return i * paramValue

        def apply_parabolic_rule(i, numberCell, paramValue):
            return i * paramValue if i < numberCell / 2 else (numberCell - i - 1) * paramValue

        def apply_sinusoide_rule(i, numberCell, paramValue):
            if i < numberCell / 4:
                return i * paramValue
            elif i < numberCell / 2:
                return (numberCell / 2 - i) * paramValue
            elif i == numberCell / 2:
                return 1.0
            elif i < 3 / 4 * numberCell:
                return (3 / 4 * numberCell - i) * (1 / paramValue)
            return 1.0

        def apply_exponential_rule(i, paramValue):
            return math.exp(i * paramValue)

        rule_functions = {
            'constant': apply_constant_rule,
            'linear': apply_linear_rule,
            'parabolic': apply_parabolic_rule,
            'sinusoide': apply_sinusoide_rule,
            'exponential': apply_exponential_rule
        }

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
                if i < numberCells[dimIndex] and direction[dimIndex] == 1:
                    rule_function = rule_functions.get(rule, apply_constant_rule)
                    if rule == 'parabolic' or rule == 'sinusoide':
                        gradientData[i][dimIndex] = rule_function(i, numberCells[dimIndex], parameters[dimIndex])
                    else:
                        gradientData[i][dimIndex] = rule_function(i, parameters[dimIndex])
                else:
                    gradientData[i][dimIndex] = 1.0
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

    def isNotInErasedRegion(self, startCellPos):
        """
        Check if the cell is not in the erased region

        Parameters:
        -----------
        startCellPos: list of float in dim 6
            (xStart, yStart, zStart, xDim, yDim, zDim) of the erased region
        """
        if self.erasedParts is None:
            return False
        counterIn = 0
        for delPart in self.erasedParts:
            for direction in range(3):
                if delPart[direction] <= startCellPos[direction] <= delPart[direction + 3] + delPart[direction]:
                    counterIn += 1
        return counterIn == 3

    def generateLattice(self):
        """
        Generates lattice structure based on specified parameters.

        Return:
        --------
        cells: list of Cell objects
            List of cells in the lattice
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
                    if not self.isNotInErasedRegion(startCellPos):
                        new_cell = []
                        # Case for normal lattice structures
                        if self.latticeType != -2 and self.latticeType < 1000:
                            new_cell = Cell(posCell, initialCellSize, startCellPos, self.latticeType, self.Radius,
                                            self.gradRadius, self.gradDim, self.gradMat)
                        # Case for hybrid lattice on 1 cell
                        elif self.latticeType == 1000:
                            for idx, radiusHybrid in enumerate(self.hybridLatticeData):
                                if radiusHybrid != 0.0:
                                    if not new_cell:
                                        new_cell = Cell(posCell, initialCellSize, startCellPos,
                                                        self.hybridGeomType[idx], radiusHybrid, self.gradRadius,
                                                        self.gradDim, self.gradMat)
                                        for beam in new_cell.beams:
                                            beam.changeBeamType(idx + 100)
                                    else:
                                        new_cell.getBeamRadius(self.gradRadius, radiusHybrid)
                                        new_cell.generateBeamsInCell(self.hybridGeomType[idx], startCellPos, idx + 100)
                            new_cell.defineHybridRadius(self.hybridLatticeData)
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

    def defineCellIndex(self):
        """
        Define index at each cell
        """
        cellIndexed = {}
        nextCellIndex = 0
        for cell in self.cells:
            if cell.index is not None:
                cellIndexed[cell] = nextCellIndex
                nextCellIndex += 1

        for cell in self.cells:
            if cell.index is None:
                if cell not in cellIndexed:
                    cell.setIndex(nextCellIndex)
                    cellIndexed[cell] = nextCellIndex
                    nextCellIndex += 1

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

    def getNodeObject(self):
        """
        Retrieves node object for the lattice.
        """
        nodeObjList = []
        for cell in self.cells:
            for beam in cell.beams:
                for node in [beam.point1, beam.point2]:
                    if node not in nodeObjList:
                        nodeObjList.append(node)
        return nodeObjList

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

    def getBeamObject(self):
        """
        Retrieves beam object for the lattice.
        """
        beamObjList = []
        self._beams = []
        for cell in self.cells:
            for beam in cell.beams:
                if beam not in beamObjList:
                    beamObjList.append(beam)
        return beamObjList

    def visualizeLattice3D(self, beamColor="Material", voxelViz=False, deformedForm=False):
        """
        Visualizes the lattice in 3D using matplotlib.

        Parameter:
        -----------
        beamColor: string (default: "Material")
            "Material" -> color by material
            "Type" -> color by type
        voxelViz: boolean (default: False)
            True -> voxel visualization
            False -> beam visualization
        deformedForm: boolean (default: False)
            True -> deformed form
        """
        fig = plt.figure()
        ax = fig.add_subplot(111, projection='3d')
        ax.set_title("Lattice generated")
        color = ['blue', 'green', 'red', 'yellow', 'orange']

        if not voxelViz:
            beamDraw = []
            lines = []
            colors = []
            nodeX = []
            nodeY = []
            nodeZ = []
            nodeDraw = set()
            for cell in self.cells:
                for beam in cell.beams:
                    if deformedForm:
                        node1 = beam.point1.getDeformedPos()
                        node2 = beam.point2.getDeformedPos()
                    else:
                        node1 = beam.point1.x, beam.point1.y, beam.point1.z
                        node2 = beam.point2.x, beam.point2.y, beam.point2.z
                    if beam not in beamDraw:
                        if beamColor == "Material":
                            colorBeam = color[beam.material]
                        elif beamColor == "Type":
                            colorBeam = color[int(str(beam.type)[0])]
                        lines.append([(node1[0], node1[1], node1[2]), (node2[0], node2[1], node2[2])])
                        colors.append(colorBeam)
                        beamDraw.append(beam)
                    for node in [node1, node2]:
                        if (node[0], node[1], node[2]) not in nodeDraw:
                            nodeDraw.add((node[0], node[1], node[2]))
                            nodeX.append(node[0])
                            nodeY.append(node[1])
                            nodeZ.append(node[2])

            line_collection = Line3DCollection(lines, colors=colors, linewidths=2)
            ax.add_collection3d(line_collection)

            ax.scatter(nodeX, nodeY, nodeZ, c='black', s=5)
        elif voxelViz:
            for cell in self.cells:
                x, y, z = cell.coordinateCell
                dx, dy, dz = cell.cellSize

                if beamColor == "Material":
                    colorCell = color[cell.beams[0].material]
                elif beamColor == "Type":
                    colorCell = color[int(str(cell.latticeType)[0])]
                ax.bar3d(x, y, z, dx, dy, dz, color=colorCell, alpha=1, shade=True, edgecolor='k')

        ax.set_xlabel('X')
        ax.set_ylabel('Y')
        ax.set_zlabel('Z')
        limMin = min(self.xMin, self.yMin, self.zMin)
        limMax = max(self.xMax, self.yMax, self.zMax)
        ax.set_xlim3d(limMin, limMax)
        ax.set_ylim3d(limMin, limMax)
        ax.set_zlim3d(limMin, limMax)
        # ax.set_xlim3d(self.xMin, self.zMax)
        # ax.set_ylim3d(self.yMin, self.yMax)
        # ax.set_zlim3d(self.zMin, self.zMax)
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
            if beam1.point1 == beam2.point1:
                u = beam1.point2 - beam1.point1
                v = beam2.point2 - beam2.point1
            elif beam1.point1 == beam2.point2:
                u = beam1.point2 - beam1.point1
                v = beam2.point1 - beam2.point2
            elif beam1.point2 == beam2.point1:
                u = beam1.point1 - beam1.point2
                v = beam2.point2 - beam2.point1
            elif beam1.point2 == beam2.point2:
                u = beam1.point1 - beam1.point2
                v = beam2.point1 - beam2.point2
            else:
                # This case should not occur if beams are connected at one point,
                # but you may handle it if needed
                raise ValueError("Beams are not connected at any point")

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
            radiusBeam.append(beam.radius)
            anglebeam.append(179.9)
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
                tag1 = beam.point1.tag
                tag2 = beam.point2.tag
                tag1 = tag1[0] if len(tag1) == 1 else None
                tag2 = tag2[0] if len(tag2) == 1 else None
                point1_tag = beamidx.point1.tag
                point2_tag = beamidx.point2.tag
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
                lengthMod = beam.getLengthMod()
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

    def functionPenalizationLzone(self, radius, angle):
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
                        tagList.append(node.tag)
                        nodeAlreadyAdded.append(node)
        return tagList

    def getTagListBoundary(self):
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
                    if node not in nodeAlreadyAdded and self.isNodeOnBoundary(node):
                        tagList.append(node.tag)
                        nodeAlreadyAdded.append(node)
        return tagList

    def applyTagToAllPoint(self):
        """
        Generate tag to all nodes in lattice
        """
        for cell in self.cells:
            for beam in cell.beams:
                for node in [beam.point1, beam.point2]:
                    tag = node.tagPoint(self.xMin, self.xMax, self.yMin, self.yMax, self.zMin, self.zMax)
                    node.setTag(tag)

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

    # def toucanLatticeModifier(self):
    #     """
    #     Fun lattice bio inspired toucan
    #     """
    #     toucanPourcent = 10
    #     for i in range(int(math.ceil(len(self.beams_obj) * toucanPourcent / 100))):
    #         nodeInit = random.choice(self.nodes_obj)
    #         node1 = random.choice(self.getConnectedNode(nodeInit))
    #         node2 = node1
    #         while (node2 == node1):
    #             node2 = random.choice(self.getConnectedNode(node1))
    #         self.toucanModifier.append(
    #             [[nodeInit.x, nodeInit.y, nodeInit.z], [node1.x, node1.y, node1.z], [node2.x, node2.y, node2.z]])
    #     return self.toucanModifier

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

    def getBeamType(self):
        """
        Retrieves beam type data for the lattice.
        data structure: each line is the type of the beam with index the line index
        """
        beamType = []
        beamAlreadyAdded = []
        for cell in self.cells:
            for beam in cell.beams:
                if beam not in beamAlreadyAdded:
                    beamType.append([beam.type])
        return beamType

    def getAllBeamLength(self):
        """
        Retrieves beam type data for the lattice.
        data structure: each line is the type of the beam with index the line index
        """
        beamLength = []
        beamAlreadyAdded = []
        for cell in self.cells:
            for beam in cell.beams:
                if beam not in beamAlreadyAdded:
                    beamLength.append([beam.length])
        return beamLength

    def changeHybridData(self, hybridRadiusData):
        """
        Change radius data of the hybrid lattice to efficiently generate graph from training neural network

        Parameters:
        ------------
        hybridRadiusData: array of dim 3
            Array of radius data of hybrid lattice
        """
        if len(hybridRadiusData) != 3:
            raise ValueError("Invalid hybrid radius data.")

        radiusDict = {
            100: hybridRadiusData[0],
            101: hybridRadiusData[1],
            102: hybridRadiusData[2],
        }
        for cell in self.cells:
            for beam in cell.beams:
                if beam.type in radiusDict:
                    beam.radius = radiusDict[beam.type]

    def attractorLattice(self, PointAttractorList=None, alpha=0.5, inverse=False):
        """
        Attract lattice to a specific point

        Parameters:
        -----------
        PointAttractor: list of float in dim 3
            Coordinates of the attractor point (default: None)
        alpha: float
            Coefficient of attraction (default: 0.5)
        inverse: bool
            If True, points farther away are attracted less (default: False)
        """

        def distance(point1, point2):
            """
            Calculate distance between two points
            """
            return math.sqrt((point2.x - point1.x) ** 2 + (point2.y - point1.y) ** 2 + (point2.z - point1.z) ** 2)

        def movePointAttracted(point, attractorPoint, alpha, inverse):
            """
            Move point1 relative from attractorPoint with coefficient alpha
            """
            Length = distance(point, attractorPoint)
            if inverse:
                factor = alpha / Length if Length != 0 else alpha
            else:
                factor = alpha * Length

            DR = [(attractorPoint.x - point.x) * factor, (attractorPoint.y - point.y) * factor,
                  (attractorPoint.z - point.z) * factor]

            pointMod = [point.x, point.y, point.z]
            pointMod = [p1 + p2 for p1, p2 in zip(pointMod, DR)]
            point.movePoint(pointMod[0], pointMod[1], pointMod[2])

        if PointAttractorList is None:
            pointAttractor = Point(5, 0.5, -2)
        else:
            pointAttractor = Point(PointAttractorList[0], PointAttractorList[1], PointAttractorList[2])

        for cell in self.cells:
            for beam in cell.beams:
                movePointAttracted(beam.point1, pointAttractor, alpha, inverse)
                movePointAttracted(beam.point2, pointAttractor, alpha, inverse)
        self.getMinMaxValues()

    def curveLattice(self, center_x, center_y, center_z, curvature_strength=0.1):
        """
        Curve the lattice structure around a given center.

        Parameters:
        -----------
        center_x: float
            The x-coordinate of the center of the curvature.
        center_y: float
            The y-coordinate of the center of the curvature.
        center_z: float
            The z-coordinate of the center of the curvature.
        curvature_strength: float (default: 0.1)
            The strength of the curvature applied to the lattice.
            Positive values curve upwards, negative values curve downwards.
        """
        for cell in self.cells:
            for beam in cell.beams:
                for node in [beam.point1, beam.point2]:
                    x, y, z = node.x, node.y, node.z
                    # Calculate the distance from the center of curvature
                    dx = x - center_x
                    dy = y - center_y
                    dz = z - center_z
                    new_z = z - curvature_strength * (dx ** 2 + dy ** 2 + dz ** 2)
                    node.movePoint(x, y, new_z)
        self.getMinMaxValues()

    def cylindrical_transform(self, radius):
        """
        Apply cylindrical transformation to the lattice structure.
        To create stent structures, 1 cell in the X direction is required and you can choose any number of cells in
        the Y and Z direction.

        Parameters:
        -----------
        radius: float
            Radius of the cylinder.
        """
        max_y = self.sizeY
        for cell in self.cells:
            for beam in cell.beams:
                for node in [beam.point1, beam.point2]:
                    x, y, z = node.x, node.y, node.z
                    # Convert Cartesian coordinates (x, y, z) to cylindrical coordinates (r, theta, z)
                    theta = (y / max_y) * 2 * math.pi  # theta = (y / total height) * 2 * pi
                    new_x = radius * math.cos(theta)
                    new_y = radius * math.sin(theta)
                    node.movePoint(new_x, new_y, z)
        self.getMinMaxValues()
        self.deleteDuplicatedBeams()

    def moveToCylinderForm(self, radius):
        """
        Move the lattice to a cylindrical form.

        Parameters:
        -----------
        radius: float
            Radius of the cylinder.
        """
        if radius <= self.xMax/2:
            raise ValueError("The radius of the cylinder is too small: minimum value = ", self.xMax/2)

        # Find moving distance
        def formula(x): # formula of arc of circle
            return radius - math.sqrt(radius ** 2 - (x - self.xMax/2) ** 2)

        for cell in self.cells:
            for beam in cell.beams:
                for node in [beam.point1, beam.point2]:
                    x, y, z = node.x, node.y, node.z
                    new_z = z - formula(x)
                    node.movePoint(x, y, new_z)
        self.getMinMaxValues()


    def deleteDuplicatedBeams(self):
        """
        Delete duplicated beams in the lattice
        """
        beamList = []
        for cell in self.cells:
            for beam in cell.beams:
                if beam not in beamList:
                    beamList.append(beam)
                else:
                    cell.removeBeam(beam)

    def getRelativeDensity(self):
        """
        Get relative density of the lattice
        """
        volumeLattice = self.getSizeLattice()[0] * self.getSizeLattice()[1] * self.getSizeLattice()[2]
        volumeBeams = 0
        for cell in self.cells:
            for beam in cell.beams:
                volumeBeams += beam.getVolume()
        return volumeBeams / volumeLattice

    def getNumberOfBeams(self):
        """
        Get number of beams in the lattice

        Returns:
        --------
        numBeams: int
            Number of beams in the lattice
        """
        numBeams = 0
        for cell in self.cells:
            numBeams += len(cell.beams)
        return numBeams

    def getNumberOfNodes(self):
        """
        Get number of nodes in the lattice

        Returns:
        --------
        numNodes: int
            Number of nodes in the lattice
        """
        numNodes = 0
        nodeIndexList = []
        for cell in self.cells:
            for beam in cell.beams:
                for node in [beam.point1, beam.point2]:
                    if node.index not in nodeIndexList:
                        nodeIndexList.append(node.index)
                        numNodes += 1
        return numNodes

    def latticeInfo(self):
        """
        Get information about the lattice
        """
        self.getName()
        print("Lattice name: ", self.name)
        latticeDim = self.getSizeLattice()
        print("Lattice size X: ", latticeDim[0])
        print("Lattice size Y: ", latticeDim[1])
        print("Lattice size Z: ", latticeDim[2])

        print("Number of beams: ", self.getNumberOfBeams())
        print("Number of nodes: ", self.getNumberOfNodes())

    def applyBoundaryConditionsOnSurface(self, cellList, surface, valueDisplacement):
        """
        Apply boundary conditions to the lattice

        Parameters:
        -----------
        cellList: list of int
            List of cell index to apply boundary conditions
        surface: str
            Surface to apply boundary conditions (Xmin, Xmax, Ymin, Ymax, Zmin, Zmax)
        valueDisplacement: np.array dim 6
            Displacement value to apply to the boundary conditions
        """
        if surface not in ["Xmin", "Xmax", "Ymin", "Ymax", "Zmin", "Zmax"]:
            raise ValueError("Invalid surface name.")
        if len(valueDisplacement) != 6:
            raise ValueError("Invalid displacement value, need dimension 6.")
        if self.cells[-1].index < max(cellList):
            raise ValueError("Invalid cell index, cell do not exist.")

        pointList = []
        for cell in self.cells:
            if cell.index in cellList:
                pointList.append(cell.getPointOnSurface(surface))
        pointList = [point for sublist in pointList for point in sublist]
        indexBoundaryList = []
        for point in pointList:
            indexBoundaryList.append(point.indexBoundary)
        for cell in self.cells:
            for beam in cell.beams:
                for node in [beam.point1, beam.point2]:
                    if node.indexBoundary in indexBoundaryList:
                        node.setDisplacementValue(valueDisplacement)

    def applyBoundaryConditionsOnNode(self, nodeList, valueDisplacement):
        """
        Apply boundary conditions to the lattice

        Parameters:
        -----------
        nodeList: list of int
            List of node index to apply boundary conditions
        valueDisplacement: np.array dim 6
            Displacement value to apply to the boundary conditions
        """
        if len(valueDisplacement) != 6:
            raise ValueError("Invalid displacement value, need dimension 6.")
        if self.getNumberOfNodes() < max(nodeList):
            raise ValueError("Invalid node index, node do not exist.")

        indexBoundaryList = []
        for node in nodeList:
            if node < 0 or node >= self.getNumberOfNodes():
                raise ValueError("Node index out of range.")
            indexBoundaryList.append(node)

        for cell in self.cells:
            for beam in cell.beams:
                for node in [beam.point1, beam.point2]:
                    if node.index in indexBoundaryList:
                        node.setDisplacementValue(valueDisplacement)

    def fixDOFOnSurface(self, cellList, surface, dofFixed):
        """
        Fix degree of freedom on the surface of the lattice

        Parameters:
        -----------
        cellList: list of int
            List of cell index to apply boundary conditions
        surface: str
            Surface to apply boundary conditions (Xmin, Xmax, Ymin, Ymax, Zmin, Zmax)
        dofFixed: list of int
            List of degree of freedom to fix (0: x, 1: y, 2: z, 3: Rx, 4: Ry, 5: Rz)
        """
        if surface not in ["Xmin", "Xmax", "Ymin", "Ymax", "Zmin", "Zmax"]:
            raise ValueError("Invalid surface name.")
        if self.cells[-1].index < max(cellList):
            raise ValueError("Invalid cell index, cell do not exist.")

        pointList = []
        for cell in self.cells:
            if cell.index in cellList:
                pointList.append(cell.getPointOnSurface(surface))
        pointList = [point for sublist in pointList for point in sublist]
        indexBoundaryList = []
        for point in pointList:
            indexBoundaryList.append(point.indexBoundary)
        for cell in self.cells:
            for beam in cell.beams:
                for node in [beam.point1, beam.point2]:
                    if node.indexBoundary in indexBoundaryList:
                        node.fixDOF(dofFixed)

    def getDisplacementGlobal(self):
        """
        Get global displacement of the lattice

        Returns:
        --------
        globalDisplacement: dict
            Dictionary of global displacement with indexBoundary as key and displacement vector as value
        """
        globalDisplacement = {i: [0, 0, 0, 0, 0, 0] for i in range(self.maxIndexBoundary + 1)}
        for cell in self.cells:
            for beam in cell.beams:
                for node in [beam.point1, beam.point2]:
                    if node.indexBoundary is not None:
                        globalDisplacement[node.indexBoundary] = node.getDisplacementValue()
        return np.concatenate(list(globalDisplacement.values()))

    def defineNodeIndexBoundary(self):
        """
        Define tag for all boundary nodes
        """
        IndexCounter = 0
        nodeAlreadyIndexed = {}
        for cell in self.cells:
            for beam in cell.beams:
                for node in [beam.point1, beam.point2]:
                    if node.tagPoint(cell.coordinateCell[0], cell.coordinateCell[0] + cell.cellSize[0],
                                     cell.coordinateCell[1], cell.coordinateCell[1] + cell.cellSize[1],
                                     cell.coordinateCell[2], cell.coordinateCell[2] + cell.cellSize[2]):
                        if node in nodeAlreadyIndexed:
                            node.setIndexBoundary(nodeAlreadyIndexed[node])
                        else:
                            nodeAlreadyIndexed[node] = IndexCounter
                            node.setIndexBoundary(IndexCounter)
                            IndexCounter += 1
        self.maxIndexBoundary = IndexCounter - 1

    def getGlobalReactionForce(self):
        """
        Get local reaction force of the lattice and sum if identical TagIndex
        """
        globalReactionForce = {i: [0, 0, 0, 0, 0, 0] for i in range(self.maxIndexBoundary + 1)}
        for cell in self.cells:
            for beam in cell.beams:
                for node in [beam.point1, beam.point2]:
                    if node.indexBoundary is not None:
                        globalReactionForce[node.indexBoundary] += node.getReactionForce()
        return globalReactionForce

    def getTotalDOF(self):
        """
        Get total number of degree of freedom in the lattice
        """
        totalDOF = 0
        for cell in self.cells:
            for beam in cell.beams:
                for node in [beam.point1, beam.point2]:
                    if node.indexBoundary is not None:
                        totalDOF += np.sum(node.fixDOF == 0)
        return totalDOF





