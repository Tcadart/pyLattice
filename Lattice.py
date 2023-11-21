from Cellule import *
import math


class Lattice:
    """
    Represents a lattice of cells.
    """
    def __init__(self, cell_size_x, cell_size_y, cell_size_z, num_cells_x, num_cells_y, num_cells_z, Lattice_Type, Radius,gradRadius,gradDim,gradMat,simMethod):
        """
        Constructor for the Lattice class.

        :param cell_size_x: integer
        :param cell_size_y: integer
        :param cell_size_z: integer
        :param num_cells_x: integer
        :param num_cells_y: integer
        :param num_cells_z: integer
        :param Lattice_Type: string
        :param Radius: float
        :param gradRadius: list of integer
        :param gradDim: list of integer
        :param gradMat: list of integer
        :param simMethod: integer
        """        
        self.cell_size_x = cell_size_x
        self.cell_size_y = cell_size_y
        self.cell_size_z = cell_size_z
        self.num_cells_x = num_cells_x
        self.num_cells_y = num_cells_y
        self.num_cells_z = num_cells_z
        self.Lattice_Type = Lattice_Type
        self.Radius = Radius
        self.gradRadius = gradRadius
        self.gradDim = gradDim
        self.gradMat = gradMat
        self.simMethod = simMethod
        self.size_x = self.getSize(0)
        self.size_y = self.getSize(1)
        self.size_z = self.getSize(2)

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
        self.generate_custom_lattice()
        self.getNodesObj()
        self.getBeamsObj()
        self.Getangle()
        if simMethod == 1:
            self.getBeamNodeMod()
        self.getNodeData()
        self.getBeamData()
        self.getRadiusData()
        self.getMaterialData()
        self.extremumFunction()
        self.getCenterCells()
        self.getNumberBeamCell()

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

    # def generate_random_cells(self, Tmin, Tmax, padding):
    #     """
    #     Generates a random cell.

    #     :param Tmin: float
    #     :param Tmax: float
    #     :param padding: float
    #     :return: Cellule object
    #     """        
    #     cell = Cellule(self.cell_size_x, self.cell_size_y, self.cell_size_z, 0, 0, 0)
    #     cell.generate_nodes(Tmin, Tmax, padding)
    #     time.sleep(1)
    #     cell.merge_nodes()
    #     cell.generate_beams()
    #     cell.remove_unused_nodes()
    #     return cell

    def getSize(self, direction):
        """
        Computes the size of the lattice along a given direction.

        :param direction: integer
        :return: float
        """        
        length = 0
        if direction == 0:
            for dim in self.gradDim:
                length += dim[0] * self.cell_size_x
        elif direction == 1:
            for dim in self.gradDim:
                length += dim[1] * self.cell_size_y
        elif direction == 2:
            for dim in self.gradDim:
                length += dim[2] * self.cell_size_z
        return length

    # def generate_random_lattice(self, Tmin, Tmax, padding):
    #     """
    #     Generates a random lattice.

    #     :param Tmin: float
    #     :param Tmax: float
    #     :param padding: float
    #     :return: list of Cellule objects
    #     """        
    #     self.cells = []
    #     random_cell = self.generate_random_cells(Tmin, Tmax, padding)

    #     for i in range(self.num_cells_x):
    #         for j in range(self.num_cells_y):
    #             for k in range(self.num_cells_z):
    #                 new_cell = copy.deepcopy(random_cell)
    #                 new_nodes = new_cell.nodes
    #                 new_cell.nodes = new_nodes
    #                 new_cell.merge_nodes()
    #                 new_beams = new_cell.beams
    #                 new_cell.beams = new_beams
    #                 new_cell.remove_unused_nodes()
    #                 new_cell.translate((i, j, k))
    #                 self.cells.append(new_cell)
    #     return self.cells
    
    def generate_custom_lattice(self):
        """
        Generates a custom lattice based on specified parameters.

        :return: tuple containing list of Cellule objects and list of positions
        """
        xCellStart = 0
        yCellStart = 0
        zCellStart = 0
        self.cells = []
        for i in range(self.num_cells_x):
            if i != 0:
                xCellStart += self.cell_size_x*self.gradDim[posCell[0]][0]
            else:
                xCellStart = 0
            for j in range(self.num_cells_y):
                if j != 0:
                    yCellStart += self.cell_size_y*self.gradDim[posCell[1]][1]
                else:
                    yCellStart = 0
                for k in range(self.num_cells_z):
                    if k != 0:
                        zCellStart += self.cell_size_z*self.gradDim[posCell[2]][2]
                    else:
                        zCellStart = 0
                    posCell = [i,j,k]
                    cellSizeX = self.cell_size_x*self.gradDim[posCell[0]][0]
                    cellSizeY = self.cell_size_y*self.gradDim[posCell[1]][1]
                    cellSizeZ = self.cell_size_z*self.gradDim[posCell[2]][2]
                    new_cell = Cellule(cellSizeX, cellSizeY, cellSizeZ, xCellStart, yCellStart, zCellStart)
                    new_cell.generate_beams_from_given_point_list(self.Lattice_Type, self.Radius, self.gradRadius, self.gradDim, self.gradMat, posCell)
                    self.cells.append(new_cell)
                    self.posCell.append([i,j,k])
        return self.cells, self.posCell

    def getNodesObj(self):
        """
        Retrieves the list of unique nodes from cells.

        :return: list of Node objects
        """
        seen_coordinates = set()
        for cell in self.cells:
            for node in cell.nodes:
                node_coordinates = (node.x, node.y, node.z)
                if node_coordinates not in seen_coordinates:
                    seen_coordinates.add(node_coordinates)
                    self.nodes_obj.append(node)
        return self.nodes_obj
    
    def getBeamsObj(self):
        """
        Retrieves the list of unique beams from cells.

        :return: list of Beam objects
        """
        seen_been = set()
        for cell in self.cells:
            for beam in cell.beams:
                beam_idx = (self.getPointIndex(beam.point1), self.getPointIndex(beam.point2))
                if beam_idx not in seen_been:
                    seen_been.add(beam_idx)
                    self.beams_obj.append(beam)
        return self.beams_obj

    def getNodeData(self):
        """
        Retrieves node data for the lattice.

        :return: list of lists
        """
        for index, point in enumerate(self.nodes_obj):
            self._nodes.append([index, point.x, point.y, point.z])
        return self._nodes
    
    def getBeamData(self):
        """
        Retrieves beam data for the lattice.

        :return: list of lists
        """
        for index, beam in enumerate(self.beams_obj):
            self._beams.append([index, self.getPointIndex(beam.point1), self.getPointIndex(beam.point2),beam.type])
        return self._beams

    def getPointIndex(self, point):
        """
        Retrieves the index of a point in the list of nodes.

        :param point: Point object
        :return: integer or None
        """
        for index, node in enumerate(self.nodes_obj):
            if point.x == node.x and point.y == node.y and point.z == node.z:
                return index

    def getRadiusData(self):
        """
        Retrieves radius data for the beams.

        :return: list of floats
        """
        for beam in self.beams_obj:
            self._radius.append(beam.radius)
        return self._radius
    
    def getMaterialData(self):
        """
        Retrieves material data for the beams.

        :return: list of strings
        """
        for beam in self.beams_obj:
            self._material.append(beam.material)
        return self._material

    def visualize_3d(self, ax):
        """
        Visualizes the lattice in 3D using matplotlib.

        :param ax: Axes3D object
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

    def Getangle(self):
        """
        Calculates angles between beams in the lattice.

        :return: list of tuples ((beam_index, Angle mininmum point 1, minRad1, Angle mininmum point 2, minRad2))
        """
        def calculate_angle(u, v):
            """
            Calculates angle between 2 beams

            :return angle in degres
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

        :return: tuple of floats (xMin, xMax, yMin, yMax, zMin, zMax)
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

    def getNbBeamCell(self, latticeType):
        """
        Gets the number of beams in a cell based on lattice type.

        :param latticeType: string
        :return: integer
        """
        nbBeam = len(Cellule.Lattice_geometry(self.cells[0],latticeType))
        return nbBeam
    
    def getBeamNodeMod(self):
        """
        Gets modified beams and nodes for simulation method 1.

        :return: tuple containing list of Beam objects and list of Node objects
        """
        def distance(point1, point2):
            return math.sqrt((point2.x - point1.x)**2 + (point2.y - point1.y)**2 + (point2.z - point1.z)**2)

        lengthMod = self.getLengthMod()
        beamMod = []
        indexCell = -1
        for index, beam in enumerate(self.beams_obj):
            if index%(self.getNbBeamCell(self.Lattice_Type)) == 0:
                indexCell = indexCell+1
            DR = [x / math.sqrt(sum([(beam.point2.x-beam.point1.x)**2, (beam.point2.y-beam.point1.y)**2, (beam.point2.z-beam.point1.z)**2])) for x in [beam.point2.x-beam.point1.x, beam.point2.y-beam.point1.y, beam.point2.z-beam.point1.z]]
            factor1 = [dr * lengthMod[index][1] for dr in DR]
            pointExt1 = [beam.point1.x, beam.point1.y, beam.point1.z]
            pointExt1 = [p1 + p2 for p1, p2 in zip(pointExt1, factor1)]
            DR = [x / math.sqrt(sum([(beam.point2.x-beam.point1.x)**2, (beam.point2.y-beam.point1.y)**2, (beam.point2.z-beam.point1.z)**2])) for x in [beam.point1.x-beam.point2.x, beam.point1.y-beam.point2.y, beam.point1.z-beam.point2.z]]
            factor2 = [dr * lengthMod[index][2] for dr in DR]
            pointExt2 = [beam.point2.x, beam.point2.y, beam.point2.z]
            pointExt2 = [p1 + p2 for p1, p2 in zip(pointExt2, factor2)]
            pointExt1Obj = Point(pointExt1[0],pointExt1[1],pointExt1[2])
            pointExt2Obj = Point(pointExt2[0],pointExt2[1],pointExt2[2])
            beamExt1 = Beam(beam.point1,pointExt1Obj,self.Radius*1.5,self.cell_size_x, self.cell_size_y, self.cell_size_z, self.gradRadius, self.gradMat, self.posCell[indexCell],1)
            beamCenter = Beam(pointExt1Obj,pointExt2Obj,self.Radius,self.cell_size_x, self.cell_size_y, self.cell_size_z, self.gradRadius, self.gradMat, self.posCell[indexCell],0)
            beamExt2 = Beam(pointExt2Obj,beam.point2,self.Radius*1.5,self.cell_size_x, self.cell_size_y, self.cell_size_z, self.gradRadius, self.gradMat, self.posCell[indexCell],1)

            self.nodes_obj.append(pointExt1Obj)
            self.nodes_obj.append(pointExt2Obj)
            beamMod.append(beamExt1)
            beamMod.append(beamExt2)
            beamMod.append(beamCenter)
        self.beams_obj = beamMod
        return self.beams_obj, self.nodes_obj

    def getLengthMod(self):
        """
        Gets modified lengths of beams for simulation method 1.

        :return: list of tuples ((beam_index,Lmod point1,Lmod point2))
        """
        lengthMod = []
        for index, angle1, radius1, angle2, radius2 in self.angles:
            if angle1 > 170:
                L1 = 0.0001
            else:
                a = radius1/math.sin(math.radians(angle1))
                b = radius1/math.sin(math.radians(angle1))
                L1 = math.sqrt(math.pow(a,2)+math.pow(b,2)-2*a*b*math.cos(math.pi-math.radians(angle1)))
            if angle2 > 170:
                L2 = 0.0001
            else:
                a = radius2/math.sin(math.radians(angle2))
                b = radius2/math.sin(math.radians(angle2))
                L2 = math.sqrt(math.pow(a,2)+math.pow(b,2)-2*a*b*math.cos(math.pi-math.radians(angle2)))
            lengthMod.append((index,L1,L2))
        return lengthMod


    def remove_cell(self, index):
        """
        Removes a cell from the lattice.

        :param index: integer
        :raises IndexError: if index is out of bounds
        """
        if 0 <= index < len(self.cells):
            del self.cells[index]
        else:
            raise IndexError("Invalid cell index.")

    def getCenterCells(self):
        self._centerCell = []
        for index, cell in enumerate(self.cells):
            self._centerCell.append([cell.centerCell,self.posCell[index]])
        return self._centerCell

    def getNumberBeamCell(self):
        Lat = self.cells[0].Lattice_geometry(self.Lattice_Type)
        self.lenghtLat = len(Lat)
        return self.lenghtLat