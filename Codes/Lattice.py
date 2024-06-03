from Cellule import *
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
                 hybridLatticeData = None, periodicity: bool = 0):
        """
        Constructor general for the Lattice class.

        Parameter:
        -----------
        cell_size_x: integer
        cell_size_y: integer
        cell_size_z: integer
            Dimension in each direction of the intial cell in the structure

        num_cells_x: integer
        num_cells_y: integer
        num_cells_z: integer
            Number of cell in each direction in the structure

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
        self.periodicity = periodicity # Not finish to implemented
        self.penalizationCoefficient = 1.5 # Fixed with previous optimization

        self.cells = []
        self.angles = []
        self._nodes = []
        self._beams = []
        self._radius = []
        self._material = []
        self._angle = []
        self.nodes_obj = []
        self.beams_obj = []
        self.posCell = []
        self.toucanModifier = []
        self.generateLattice()
        self.getNodesObj()
        self.getBeamsObj()

        self.getMinMaxValues()
        if self.latticeType == 1000: # case for hybrid lattice structures
            self.checkHybridCollision()
            # self.LatticeToPyGeometricData()
        self.getAllAngles()

        # Case of penalization at beam near nodes
        if self.simMethod == 1:
            self.getBeamNodeMod()

        # self.findBoundaryBeams()
        # self.attractorLattice()
        # Get some data about lattice structures
        self.getNodeData()
        self.getBeamData()
        self.getRadiusData()
        self.getMaterialData()


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
                 Radius,gradRadiusProperty,gradDimProperty,gradMatProperty,0,0, None)

    @classmethod
    def hybridgeometry(cls, cell_size_x, cell_size_y, cell_size_z,simMethod,uncertaintyNode,hybridLatticeData):
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
                 1,gradRadiusProperty,gradDimProperty,gradMatProperty,simMethod,uncertaintyNode,
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
                 1,gradRadiusProperty,gradDimProperty,gradMatProperty,simMethod,uncertaintyNode,
                 hybridLatticeData)

    @property
    def nodes(self):
        return self._nodes
    
    @property
    def beams(self):
        return self._beams

    @property
    def radius(self):
        return self._radius

    @property
    def material(self):
        return self._material
    
    @property
    def angle(self):
        return self._angle

    @property
    def centerCell(self):
        return self._centerCell

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
                gradientData[i][dimIndex] = applyRule(i, numberCells[dimIndex], direction[dimIndex], parameters[dimIndex],
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
            gradMat = [[[1 for X in range(self.numCellsX)] for Y in range(self.numCellsY)] for Z in range(self.numCellsZ)]
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
        xCellStart = 0
        yCellStart = 0
        zCellStart = 0
        self.cells = []
        for i in range(self.numCellsX):
            if i != 0:
                xCellStart += self.cellSizeX * self.gradDim[posCell[0]][0]
            else:
                xCellStart = 0
            for j in range(self.numCellsY):
                if j != 0:
                    yCellStart += self.cellSizeY * self.gradDim[posCell[1]][1]
                else:
                    yCellStart = 0
                for k in range(self.numCellsZ):
                    if k != 0:
                        zCellStart += self.cellSizeZ * self.gradDim[posCell[2]][2]
                    else:
                        zCellStart = 0
                    posCell = [i,j,k]
                    cellSizeX = self.cellSizeX * self.gradDim[posCell[0]][0]
                    cellSizeY = self.cellSizeY * self.gradDim[posCell[1]][1]
                    cellSizeZ = self.cellSizeZ * self.gradDim[posCell[2]][2]
                    new_cell = Cellule(cellSizeX, cellSizeY, cellSizeZ, xCellStart, yCellStart, zCellStart)
                    # Case for normal lattice structures
                    if self.latticeType != -2 and self.latticeType < 1000:
                        new_cell.generate_beams_from_given_point_list(self.latticeType, self.Radius, self.gradRadius, self.gradDim, self.gradMat, posCell)
                    # Case for hybrid lattice on 1 cell
                    elif self.latticeType == 1000:
                        latticeHybridType = [0,16,17]
                        for idx, radiusHybrid in enumerate(self.hybridLatticeData):
                            new_cell.generate_beams_from_given_point_list(latticeHybridType[idx], radiusHybrid,
                                                                          self.gradRadius, self.gradDim, self.gradMat, posCell)
                            for beam in new_cell.beams:
                                beam.changeBeamType(idx+100)
                            if idx < len(self.hybridLatticeData)-1:
                                self.cells.append(new_cell)
                                self.posCell.append([i, j, k])
                                new_cell = Cellule(cellSizeX, cellSizeY, cellSizeZ, xCellStart, yCellStart, zCellStart)
                    # Case for randomized lattice
                    else:
                        new_cell.generate_beams_random(self.Radius, self.gradRadius, self.gradDim, self.gradMat, posCell)
                    self.cells.append(new_cell)
                    self.posCell.append([i,j,k])

    def getNodesObj(self):
        """
        Retrieves the list of unique nodes objects from cells.
        """
        seen_coordinates = set()
        for cell in self.cells:
            for node in cell.nodes:
                node_coordinates = (node.x, node.y, node.z)
                if node_coordinates not in seen_coordinates:
                    seen_coordinates.add(node_coordinates)
                    self.nodes_obj.append(node)
    
    def getBeamsObj(self):
        """
        Retrieves the list of unique beams objects from cells.
        """
        seen_been = set()
        for cell in self.cells:
            for beam in cell.beams:
                beam_idx = (self.getPointIndex(beam.point1), self.getPointIndex(beam.point2))
                if beam_idx not in seen_been:
                    seen_been.add(beam_idx)
                    self.beams_obj.append(beam)

    def getNodeData(self):
        """
        Retrieves node data for the lattice.
        data structure: each line represent a node with data [indexNode, X, Y, Z]
        """
        self._nodes = []
        for index, point in enumerate(self.nodes_obj):
            self._nodes.append([index, point.x, point.y, point.z])

    def getPosData(self):
        """
        Retrieves node position data for the lattice
        data structure : each line represent a node with data [X, Y, Z]
        """
        posData = []
        for index, point in enumerate(self.nodes_obj):
            posData.append([point.x, point.y, point.z])
        return posData

    def getBeamData(self):
        """
        Retrieves beam data for the lattice.
        data structure: each line represent a beam with data [beamIndex, IndexPoint1, IndexPoint2, beamType]
        """
        self._beams = []
        for index, beam in enumerate(self.beams_obj):
            self._beams.append([index, self.getPointIndex(beam.point1), self.getPointIndex(beam.point2),beam.type])

    def getEdgeIndex(self):
        """
        Retrieves edge index data for the lattice.
        data structure: each line represent a beam with data [IndexPoint1, IndexPoint2]
        """
        edgeIndex = []
        for index, beam in enumerate(self.beams_obj):
            edgeIndex.append([self.getPointIndex(beam.point1), self.getPointIndex(beam.point2)])
        return edgeIndex

    def getBeamType(self):
        """
        Retrieves beam type data for the lattice.
        data structure: each line is the type of the beam with index the line index
        """
        BeamType = []
        for index, beam in enumerate(self.beams_obj):
            BeamType.append([beam.type])
        return BeamType

    def getPointIndex(self, point):
        """
        Retrieves the index of a point in the list of nodes.

        Parameter:
        ----------
        point: Point object

        Return:
        -------
        Node index: integer
        """
        for index, node in enumerate(self.nodes_obj):
            if point == node:
                return index

    def getBeamIndex(self, beam):
        """
        Retrieves the index of a beam in the list of nodes.

        Parameter:
        ----------
        beam: Beam object

        Return:
        -------
        Node index: integer
        """
        for index, beamlist in enumerate(self.beams_obj):
            if beam == beamlist:
                return index

    def getRadiusData(self):
        """
        Retrieves radius data for the beams.
        """
        for beam in self.beams_obj:
            self._radius.append(beam.radius)
    
    def getMaterialData(self):
        """
        Retrieves material data for the beams.
        """
        for beam in self.beams_obj:
            self._material.append(beam.material)

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
        for point in self.nodes_obj:
            x, y, z = point.x, point.y, point.z
            ax.scatter(x, y, z, c='black', s=5)
        color = ['blue','green','red','yellow','orange']
        for beam in self.beams_obj:
            point1 = beam.point1
            point2 = beam.point2
            if beamColor == "Material":
                ax.plot([point1.x, point2.x], [point1.y, point2.y], [point1.z, point2.z], color=color[beam.material], markersize=10)
            elif beamColor == "Type":
                ax.plot([point1.x, point2.x], [point1.y, point2.y], [point1.z, point2.z],
                        color=color[int(str(beam.type)[0])],
                        markersize=10)
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
        color = ['blue','green','red','yellow','orange']
        for beam in self.beams_obj:
            point1 = beam.point1
            point2 = beam.point2
            ax.plot([point1.x, point2.x], [point1.y, point2.y], [point1.z, point2.z], color=color[beam.material], markersize=10)
        ax.set_xlabel('X')
        ax.set_ylabel('Y')
        ax.set_zlabel('Z')
        ax.set_xlim3d(0, self.xMax)
        ax.set_ylim3d(0, self.yMax)
        ax.set_zlim3d(0, self.zMax)

    def getAngleBetweenBeams(self, beam1, beam2):
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
        anglebeam = []
        radiusBeam = []
        if len(pointbeams) > 1:
            for beampoint in pointbeams:
                radiusBeam.append(beampoint.radius)
                anglebeam.append(self.getAngleBetweenBeams(beam, beampoint))
        else:
            anglebeam.append(179.9)
            radiusBeam.append(beam.radius)
        non_zero_anglebeam = [angle for angle in anglebeam if angle >= 0.01]
        non_zero_radiusbeam = [radius for angle, radius in zip(anglebeam, radiusBeam) if angle >= 0.01]
        return non_zero_anglebeam, non_zero_radiusbeam

    def getConnectedBeams(self, beam):
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
        for beamidx in self.beams_obj:
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

        def findMinAngle(angles,radii):
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
            return LAngle,LRadius

        for index, beam in enumerate(self.beams_obj):
            # Determine beams on nodes
            point1beams, point2beams = self.getConnectedBeams(beam)

            # Determine angle for all beams connected at the node
            non_zero_anglebeam1,non_zero_radiusbeam1 = self.getListAngleBeam(beam, point1beams)
            non_zero_anglebeam2,non_zero_radiusbeam2 = self.getListAngleBeam(beam, point2beams)
            # Find the lowest angle
            LAngle1,LRadius1 = findMinAngle(non_zero_anglebeam1, non_zero_radiusbeam1)
            LAngle2,LRadius2 = findMinAngle(non_zero_anglebeam2,non_zero_radiusbeam2)
            self.angles.append((index, round(LAngle1,2), LRadius1, round(LAngle2,2), LRadius2))
        return self.angles


    def getMinMaxValues(self):
        """
        Computes extremum values of coordinates in the lattice.

        Return:
        --------
        ExtrumumValues: tuple of floats (xMin, xMax, yMin, yMax, zMin, zMax)
        """
        x_values = [node.x for cell in self.cells for node in cell.nodes]
        y_values = [node.y for cell in self.cells for node in cell.nodes]
        z_values = [node.z for cell in self.cells for node in cell.nodes]

        self.xMin = min(x_values)
        self.xMax = max(x_values)
        self.yMin = min(y_values)
        self.yMax = max(y_values)
        self.zMin = min(z_values)
        self.zMax = max(z_values)

        return self.xMin, self.xMax, self.yMin, self.yMax, self.zMin, self.zMax

    def getNbBeamCell(self):
        """
        Gets the number of beams in a cell based on lattice type.

        Return:
        ----------
        nbBeam: integer
            Number of beams in the lattice cell
        """
        self.nbBeam = len(Cellule.Lattice_geometry(self.cells[0], self.latticeType))
        return self.nbBeam
    
    def getBeamNodeMod(self):
        """
        Modifies beam and node data to model lattice structures for simulation with rigidity penalization at node
        """
        beamMod = []
        indexCell = -1
        self.getNbBeamCell()
        for index, beam in enumerate(self.beams_obj):
            lengthMod = self.getLengthMod(beam)
            if index%(self.nbBeam) == 0:
                indexCell = indexCell+1
            print("lengthMod", lengthMod)
            pointExt1 = beam.getPointOnBeamFromDistance(lengthMod[0], 1)
            pointExt1Obj = Point(pointExt1[0], pointExt1[1], pointExt1[2])
            pointExt2 = beam.getPointOnBeamFromDistance(lengthMod[1], 2)
            pointExt2Obj = Point(pointExt2[0], pointExt2[1], pointExt2[2])
            if beam.type == 0:
                typeBeam = [1,0,1]
            else:
                typeBeam = [beam.type+100,beam.type,beam.type+100]
            beamExt1 = Beam(beam.point1, pointExt1Obj, beam.radius * self.penalizationCoefficient, beam.material,
                            typeBeam[0])
            beamCenter = Beam(pointExt1Obj, pointExt2Obj, beam.radius, beam.material, typeBeam[1])
            beamExt2 = Beam(pointExt2Obj, beam.point2, beam.radius * self.penalizationCoefficient, beam.material,
                            typeBeam[2])
            print((beamExt1.getLength() / beam.getLength()) * 100,
                  (beamCenter.getLength() / beam.getLength()) * 100,
                  (beamExt2.getLength() / beam.getLength()) * 100)
            print((beamExt1.getLength() / beam.getLength()) * 100 +
                  (beamCenter.getLength() / beam.getLength()) * 100 +
                  (beamExt2.getLength() / beam.getLength()) * 100)

            self.nodes_obj.append(pointExt1Obj)
            self.nodes_obj.append(pointExt2Obj)
            beamMod.append(beamExt1)
            beamMod.append(beamExt2)
            beamMod.append(beamCenter)
        self.beams_obj = beamMod

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

    def getLengthMod2(self):
        """
        Calculate and return length to modify in penalization method

        Return:
        --------
        LenghtMod: list of tuples ((beam_index,Lmod point1,Lmod point2))
        """

        lengthMod = []
        for index, angle1, radius1, angle2, radius2 in self.angles:
            L1 = self.functionPenalizationLzone(radius1,angle1)
            L2 = self.functionPenalizationLzone(radius2,angle2)
            lengthMod.append((index,L1,L2))
        return lengthMod

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
        beamIndex = self.getBeamIndex(beam)
        L1 = self.functionPenalizationLzone(self.angles[beamIndex][2],self.angles[beamIndex][1])
        L2 = self.functionPenalizationLzone(self.angles[beamIndex][4],self.angles[beamIndex][3])
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

    def removeBeam(self, index):
        """
        Removes a beam from the lattice

        Parameters:
        ------------
        index: int
            index of the beam to remove
        """
        if 0 <= index < len(self.beams_obj):
            del self.beams_obj[index]
        else:
            raise IndexError("Invalid beam index.")

    def removeNode(self, index):
        """
        Removes a node from the lattice

        Parameters:
        ------------
        index: int
            index of the node to remove
        """
        if 0 <= index < len(self.nodes_obj):
            del self.nodes_obj[index]
        else:
            raise IndexError("Invalid node index.")

    def getCenterCells(self):
        """
        Determine center cell of the lattice
        """
        self._centerCell = []
        for index, cell in enumerate(self.cells):
            self._centerCell.append([cell.centerCell,self.posCell[index]])

    def findMinimumBeamLength(self):
        """
        Find minimum beam length

        Returns:
        --------
        minLength: float
            Length of smallest beam in the lattice
        """
        minLength = 1000
        for index, beam in enumerate(self.beams_obj):
            if beam.getLength() <minLength and beam.getLength() >0.0001 and (beam.type == 0 or beam.type == 2):
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
        if point.x == self.xMin and (point.y > self.yMin and point.y < self.yMax) and (
                point.z > self.zMin and point.z < self.zMax):
            tags.append(12)  # Face 1
        elif point.x == self.xMax and (point.y > self.yMin and point.y < self.yMax) and (
                point.z > self.zMin and point.z < self.zMax):
            tags.append(13)  # Face 2
        elif (point.x > self.xMin and point.x < self.xMax) and point.y == self.yMin and (
                point.z > self.zMin and point.z < self.zMax):
            tags.append(11)  # Face 3
        elif (point.x > self.xMin and point.x < self.xMax) and point.y == self.yMax and (
                point.z > self.zMin and point.z < self.zMax):
            tags.append(14)  # Face 4
        elif (point.x > self.xMin and point.x < self.xMax) and (
                point.y > self.yMin and point.y < self.yMax) and point.z == self.zMin:
            tags.append(10)  # Face 5
        elif (point.x > self.xMin and point.x < self.xMax) and (
                point.y > self.yMin and point.y < self.yMax) and point.z == self.zMax:
            tags.append(15)  # Face 6

        # Edge
        if point.x == self.xMin and point.y == self.yMin and (point.z > self.zMin and point.z < self.zMax):
            tags.append(102)  # Edge 0
        elif (point.x > self.xMin and point.x < self.xMax) and point.y == self.yMin and point.z == self.zMin:
            tags.append(100)  # Edge 1
        elif point.x == self.xMax and point.y == self.yMin and (point.z > self.zMin and point.z < self.zMax):
            tags.append(104)  # Edge 2
        elif (point.x > self.xMin and point.x < self.xMax) and point.y == self.yMin and point.z == self.zMax:
            tags.append(108)  # Edge 3
        elif point.x == self.xMin and (point.y > self.yMin and point.y < self.yMax) and point.z == self.zMin:
            tags.append(101)  # Edge 4
        elif point.x == self.xMax and (point.y > self.yMin and point.y < self.yMax) and point.z == self.zMin:
            tags.append(103)  # Edge 5
        elif point.x == self.xMin and point.y == self.yMax and (point.z > self.zMin and point.z < self.zMax):
            tags.append(106)  # Edge 6
        elif (point.x > self.xMin and point.x < self.xMax) and point.y == self.yMax and point.z == self.zMin:
            tags.append(105)  # Edge 7
        elif point.x == self.xMax and point.y == self.yMax and (point.z > self.zMin and point.z < self.zMax):
            tags.append(107)  # Edge 8
        elif (point.x > self.xMin and point.x < self.xMax) and point.y == self.yMax and point.z == self.zMax:
            tags.append(111)  # Edge 9
        elif point.x == self.xMin and (point.y > self.yMin and point.y < self.yMax) and point.z == self.zMax:
            tags.append(109)  # Edge 10
        elif point.x == self.xMax and (point.y > self.yMin and point.y < self.yMax) and point.z == self.zMax:
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
        for node in self.nodes_obj:
            tagList.append(self.tagPoint(node))
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
        for beam in self.beams_obj:
            if beam.point1 == node:
                connectedNode.append(beam.point2)
            if beam.point2 == node:
                connectedNode.append(beam.point1)
        return connectedNode

    def toucanLatticeModifier(self):
        toucanPourcent = 10
        for i in range(int(math.ceil(len(self.beams_obj)*toucanPourcent/100))):
            nodeInit = random.choice(self.nodes_obj)
            node1 = random.choice(self.getConnectedNode(nodeInit))
            node2 = node1
            while(node2 == node1):
                node2 = random.choice(self.getConnectedNode(node1))
            self.toucanModifier.append([[nodeInit.x,nodeInit.y,nodeInit.z],[node1.x,node1.y,node1.z],[node2.x,node2.y,node2.z]])
        return self.toucanModifier

    def findBoundaryBeams(self):
        """
        Find boundary beams and change type of beam

        Return:
        -------
        boundaryBeams: List of beam object
        """
        boundaryNodes = self.findBoundaryNodes()
        boundaryBeams = []
        for beam in self.beams_obj:
            if beam.point1 in boundaryNodes or beam.point2 in boundaryNodes:
                beam.changeBeamType(2)
                boundaryBeams.append(beam)
        return boundaryBeams

    def findBoundaryNodes(self):
        """
        Find boundary nodes

        Returns:
        ---------
        boundaryNodes: List of point object
        """
        boundaryNodes = []
        for node in self.nodes_obj:
            if node.x == self.xMin or node.x == self.xMax or node.y == self.yMin or node.y == self.yMax or node.z == self.zMin or node.z == self.zMax:
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
        for idxnode, node in enumerate(self.nodes_obj):
            for idxbeam, beam in enumerate(self.beams_obj):
                if self.isPointOnLine(beam.point1,beam.point2,node):
                    typeBeamToRemove = beam.type # Get beam to remove type to apply in new separated beams
                    self.removeBeam(idxbeam)
                    beam1= Beam(beam.point1, node, beam.radius, beam.material, typeBeamToRemove)
                    self.beams_obj.append(beam1)
                    beam2 = Beam(beam.point2, node, beam.radius, beam.material, typeBeamToRemove)
                    self.beams_obj.append(beam2)

    def isPointOnLine(self, point1, point2, node):
        """
        Find if input node if on the line formed by point1 and point2

        Return
        -------
        boolean: True => Point on line
        """
        vector1 = (point2.x - point1.x, point2.y - point1.y, point2.z - point1.z)
        vector2 = (node.x - point1.x, node.y - point1.y, node.z - point1.z)

        if (node.x == point1.x and node.y == point1.y and node.z == point1.z) or (node.x == point2.x and node.y == point2.y and node.z == point2.z):
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

    def LatticeToPyGeometricData(self):
        """
        Convert lattice data to PyTorch Geometric Data

        Return:
        ---------
        self.dataPyGeometric: Object Data Pytorch
            Data(x=x, edge_index=edge_index, edge_attr=edge_attr)
        """
        self.getNodeData()
        self.getBeamData()
        # Create tensor for pytorch
        # Node position
        x = torch.tensor([node[1:] for node in self._nodes], dtype=torch.float)
        # Beam connexion
        edge_index = torch.tensor([[beam[1], beam[2]] for beam in self._beams], dtype=torch.long).t().contiguous()
        # Beam features (Radius)
        edge_attr = torch.tensor([beam.radius for beam in self.beams_obj], dtype=torch.float)

        # Create Data object for PyTorch
        self.dataPyGeometric = Data(x=x, edge_index=edge_index, edge_attr=edge_attr)
        return self.dataPyGeometric

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
        for beam in self.beams_obj:
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
            factor = [(alpha)*dr for dr in DR]
            pointMod = [point1.x, point1.y, point1.z]
            pointMod = [p1 + p2 for p1, p2 in zip(pointMod, factor)]
            point1.movePoint(pointMod[0], pointMod[1], pointMod[2])

        pointAttractor = Point(5,0.5,-2)
        alpha = 0.5
        for beam in self.beams_obj:
            movePointAttracted(beam.point1, pointAttractor, alpha)
            movePointAttracted(beam.point2, pointAttractor, alpha)
        self.getMinMaxValues()


