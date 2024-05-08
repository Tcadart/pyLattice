from Cellule import *
import math
import random



class Lattice:
    """
    Generate lattice structures with parameters
    """
    def __init__(self, cell_size_x, cell_size_y, cell_size_z, num_cells_x, num_cells_y, num_cells_z, Lattice_Type,
                 Radius,gradRadiusProperty,gradDimProperty,gradMatProperty,simMethod = 0,uncertaintyNode = 0,
                 hybridLatticeData = None):
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
        """        
        self.cellSizeX = cell_size_x
        self.cellSizeY = cell_size_y
        self.cellSizeZ = cell_size_z
        self.numCellsX = num_cells_x
        self.numCellsY = num_cells_y
        self.numCellsZ = num_cells_z
        self.latticeType = Lattice_Type
        self.Radius = Radius
        self.gradRadius = self.gradSettings(gradRadiusProperty)
        self.gradDim = self.gradSettings(gradDimProperty)
        self.gradMat = self.gradMaterialSetting(gradMatProperty)
        self.simMethod = simMethod
        self.sizeX = self.getSize(0)
        self.sizeY = self.getSize(1)
        self.sizeZ = self.getSize(2)
        self.uncertaintyNode = uncertaintyNode
        self.hybridLatticeData = hybridLatticeData

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
        self.generate_custom_lattice()
        self.getNodesObj()
        self.getBeamsObj()
        if self.latticeType == 1000:
            self.checkHybridCollision()
        self.Getangle()
        if self.simMethod == 1:
            self.getBeamNodeMod()

        self.extremumFunction()
        self.findBoundaryBeams()

        # Get some data about lattice structures
        self.getNodeData()
        self.getBeamData()
        self.getRadiusData()
        self.getMaterialData()


    @classmethod
    def simpleLattice(cls, cell_size_x, cell_size_y, cell_size_z, num_cells_x, num_cells_y, num_cells_z, Lattice_Type,
                 Radius):
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

    def getSize(self, direction):
        """
        Computes the size of the lattice along a given direction.

        Parameter:
        ------------
        direction: integer
            0=X, 1=Y, 2=Z

        Return:
        ---------
        Length of the lattice in direction: float
        """        
        length = 0
        if direction == 0:
            for dim in self.gradDim:
                length += dim[0] * self.cellSizeX
        elif direction == 1:
            for dim in self.gradDim:
                length += dim[1] * self.cellSizeY
        elif direction == 2:
            for dim in self.gradDim:
                length += dim[2] * self.cellSizeZ
        return length

    def gradSettings(self, gradProperties):
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

        # def calculatePositionPointOnPlane(CoeffPlan, PositionPlan , x, y, z):
        #     value = CoeffPlan[0]*x + CoeffPlan[1]*y + CoeffPlan[2]*z + PositionPlan
        #     if value > 0:
        #         return 0
        #     else:
        #         return 1
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
        # elif Multimat == 2:
        #     gradMat = [[[0 for X in range(number_cell_X)] for Y in range(number_cell_Y)] for Z in range(number_cell_Z)]
        #     centerCell = []
        #     for x in range(number_cell_X):
        #         for y in range(number_cell_Y):
        #             for z in range(number_cell_Z):
        #                 centerCell.append([0.5+x,0.5+y,0.5+z])
        #     for PositionPlan in range(number_cell_Z):
        #         IdxMat = PositionPlan%2
        #         for IdxCenterCells in range(len(lattice.centerCell)):
        #             if calculatePositionPointOnPlane([1,0,0], PositionPlan , lattice.centerCell[IdxCenterCells,0], lattice.centerCell[IdxCenterCells,1], lattice.centerCell[IdxCenterCells,2]):
        #                 gradMat[lattice.centerCell[IdxCenterCells,3]][lattice.centerCell[IdxCenterCells,4]][lattice.centerCell[IdxCenterCells,5]] = IdxMat
        return gradMat


    
    def generate_custom_lattice(self):
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
                    if self.latticeType != -2 and self.latticeType < 1000:
                        new_cell.generate_beams_from_given_point_list(self.latticeType, self.Radius, self.gradRadius, self.gradDim, self.gradMat, posCell)
                    elif self.latticeType == 1000:
                        latticeHybridType = [0,16,17]
                        for idx, radiusHybrid in enumerate(self.hybridLatticeData):
                            new_cell.generate_beams_from_given_point_list(latticeHybridType[idx], radiusHybrid,
                                                                          self.gradRadius, self.gradDim, self.gradMat, posCell)
                            if idx < 2:
                                self.cells.append(new_cell)
                                self.posCell.append([i, j, k])
                                new_cell = Cellule(cellSizeX, cellSizeY, cellSizeZ, xCellStart, yCellStart, zCellStart)
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
        data structure : each line represent a node with data [indexNode, X, Y, Z]
        """
        for index, point in enumerate(self.nodes_obj):
            self._nodes.append([index, point.x, point.y, point.z])
    
    def getBeamData(self):
        """
        Retrieves beam data for the lattice.
        data structure: each line represent a beam with data [beamIndex, IndexPoint1, IndexPoint2, beamType]
        """
        for index, beam in enumerate(self.beams_obj):
            self._beams.append([index, self.getPointIndex(beam.point1), self.getPointIndex(beam.point2),beam.type])

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
            if point.x == node.x and point.y == node.y and point.z == node.z:
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

    def visualize_3d(self, ax):
        """
        Visualizes the lattice in 3D using matplotlib.

        Parameter:
        -----------
        ax: Axes3D object
        """
        for point in self.nodes_obj:
            x, y, z = point.x, point.y, point.z
            ax.scatter(x, y, z, c='black', s=5)
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

    def Getangle(self):
        """
        Calculates angles between beams in the lattice.

        Return:
        ---------
        angle:
            data structure => ((beam_index, Angle mininmum point 1, minRad1, Angle mininmum point 2, minRad2))
        """
        def calculate_angle(u, v):
            """
            Calculates angle between 2 beams

            Return:
            --------
            Angle: float
                angle in degrees
            """
            dot_product = sum(a * b for a, b in zip(u, v))
            u_norm = math.sqrt(sum(a * a for a in u))
            v_norm = math.sqrt(sum(b * b for b in v))
            cos_theta = dot_product / (u_norm * v_norm)
            cos_theta = max(min(cos_theta, 1.0), -1.0)
            angle_rad = math.acos(cos_theta)
            angle_deg = math.degrees(angle_rad)
            return angle_deg

        def getAngleBeam(pointIdx1,pointIdx2,pointbeams):
            anglebeam = []
            radiusBeam = []
            u = [pointIdx2.x - pointIdx1.x, pointIdx2.y - pointIdx1.y, pointIdx2.z - pointIdx1.z]
            if len(pointbeams) > 1:
                for beampoint in pointbeams:
                    radiusBeam.append(beampoint.radius)
                    if beampoint.point1.x == pointIdx1.x and beampoint.point1.y == pointIdx1.y and beampoint.point1.z == pointIdx1.z:
                        v = [beampoint.point2.x - beampoint.point1.x, beampoint.point2.y - beampoint.point1.y,
                                beampoint.point2.z - beampoint.point1.z]
                    else:
                        v = [beampoint.point1.x - beampoint.point2.x, beampoint.point1.y - beampoint.point2.y,
                                    beampoint.point1.z - beampoint.point2.z]
                    anglebeam.append(calculate_angle(u, v))
            else:
                anglebeam.append(179.9)
                radiusBeam.append(beam.radius)
            non_zero_anglebeam = [angle for angle in anglebeam if angle >= 0.01]
            non_zero_radiusbeam = [radius for angle, radius in zip(anglebeam, radiusBeam) if angle >= 0.01]
            return non_zero_anglebeam,non_zero_radiusbeam

        def findMinAngle(non_zero_anglebeam,non_zero_radiusbeam):
            minAngle = min(non_zero_anglebeam)
            minRad = min(non_zero_radiusbeam)
            return minAngle,minRad
        
        for index, beam in enumerate(self.beams_obj):
            point1beams = []
            point2beams = []
            pointIdx1 = beam.point1
            pointIdx2 = beam.point2

            # Determine beams on nodes
            for beamidx in self.beams_obj:
                if pointIdx1 in (beamidx.point1, beamidx.point2):
                    point1beams.append(beamidx)
                if pointIdx2 in (beamidx.point1, beamidx.point2):
                    point2beams.append(beamidx)

            # Determine angle for all beams connected at the node
            non_zero_anglebeam1,non_zero_radiusbeam1 = getAngleBeam(pointIdx1,pointIdx2,point1beams)
            non_zero_anglebeam2,non_zero_radiusbeam2 = getAngleBeam(pointIdx2,pointIdx1,point2beams)
            
            # Find the lowest angle
            minAngle1,minRad1 = findMinAngle(non_zero_anglebeam1,non_zero_radiusbeam1)
            minAngle2,minRad2 = findMinAngle(non_zero_anglebeam2,non_zero_radiusbeam2)
            self.angles.append((index, round(minAngle1,2), minRad1, round(minAngle2,2), minRad2))

        return self.angles


    def extremumFunction(self):
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
        nbBeam = len(Cellule.Lattice_geometry(self.cells[0],self.latticeType))
        return nbBeam
    
    def getBeamNodeMod(self):
        """
        Modifies beam and node data to model lattice structures for simulation with rigidity penalization at node
        """
        def distance(point1, point2):
            """
            Calculate distance between two points
            """
            return math.sqrt((point2.x - point1.x)**2 + (point2.y - point1.y)**2 + (point2.z - point1.z)**2)

        def findPointMod(point1, point2, lengthMod):
            DR = [(point2.x - point1.x)/distance(point1,point2), (point2.y - point1.y)/distance(point1,point2), (point2.z - point1.z)/distance(point1,point2)]
            factor = [dr * lengthMod for dr in DR]
            pointMod = [point1.x, point1.y, point1.z]
            pointMod = [p1 + p2 for p1, p2 in zip(pointMod, factor)]
            pointModObj = Point(pointMod[0], pointMod[1], pointMod[2])
            return pointModObj

        lengthMod = self.getLengthMod()
        beamMod = []
        indexCell = -1
        for index, beam in enumerate(self.beams_obj):
            if index%(self.getNbBeamCell()) == 0:
                indexCell = indexCell+1
            pointExt1Obj = findPointMod(beam.point1,beam.point2, lengthMod[index][1])
            pointExt2Obj = findPointMod(beam.point2,beam.point1, lengthMod[index][2])

            beamExt1 = Beam(beam.point1, pointExt1Obj, self.Radius * 1.5, self.cellSizeX, self.cellSizeY, self.cellSizeZ, self.gradRadius, self.gradMat, self.posCell[indexCell], 1)
            beamCenter = Beam(pointExt1Obj, pointExt2Obj, self.Radius, self.cellSizeX, self.cellSizeY, self.cellSizeZ, self.gradRadius, self.gradMat, self.posCell[indexCell], 0)
            beamExt2 = Beam(pointExt2Obj, beam.point2, self.Radius * 1.5, self.cellSizeX, self.cellSizeY, self.cellSizeZ, self.gradRadius, self.gradMat, self.posCell[indexCell], 1)

            self.nodes_obj.append(pointExt1Obj)
            self.nodes_obj.append(pointExt2Obj)
            beamMod.append(beamExt1)
            beamMod.append(beamExt2)
            beamMod.append(beamCenter)
        self.beams_obj = beamMod

    def getLengthMod(self):
        """
        Calculate and return length to modify in penalization method

        Return:
        --------
        LenghtMod: list of tuples ((beam_index,Lmod point1,Lmod point2))
        """
        def getlength(angle, radius):
            if angle > 170:
                L = 0.0001
            else:
                L = radius / math.tan(math.radians(angle) / 2)
            return L

        lengthMod = []
        for index, angle1, radius1, angle2, radius2 in self.angles:
            L1 = getlength(angle1,radius1)
            L2 = getlength(angle2,radius2)
            lengthMod.append((index,L1,L2))
        return lengthMod


    def removeCell(self, index):
        """
        Removes a cell from the lattice
        """
        if 0 <= index < len(self.cells):
            del self.cells[index]
        else:
            raise IndexError("Invalid cell index.")

    def removeBeam(self, index):
        """
        Removes a beam from the lattice
        """
        if 0 <= index < len(self.beams_obj):
            del self.beams_obj[index]
        else:
            raise IndexError("Invalid beam index.")

    def removeNode(self, index):
        """
        Removes a node from the lattice
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

    def getNumberBeamCell(self):
        """
        Calculate number of beams in a cell
        """
        cell = self.cells[0]
        self.lenghtLat = len(cell.beams)
        return self.lenghtLat

    def findMinimumBeamLength(self):
        """
        Find minimum beam length
        """
        minLength = 1000
        for index, beam in enumerate(self.beams_obj):
            if beam.get_length()<minLength and beam.get_length()>0.0001 and (beam.type == 0 or beam.type == 2):
                minLength = beam.get_length()
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
                    self.removeBeam(idxbeam)
                    beam1= Beam(beam.point1, node, self.Radius, self.cellSizeX, self.cellSizeY, self.cellSizeZ,
                                self.gradRadius, self.gradMat,[0,0,0],0)
                    self.beams_obj.append(beam1)
                    beam2 = Beam(beam.point2, node, self.Radius, self.cellSizeX, self.cellSizeY, self.cellSizeZ,
                                 self.gradRadius, self.gradMat, [0,0,0], 0)
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