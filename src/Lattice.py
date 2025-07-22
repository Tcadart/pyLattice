import json
import os
import pickle
from statistics import mean

import joblib
from colorama import Style, Fore
from matplotlib import pyplot as plt

from .Cell import *
import math
import random

import numpy as np
from scipy.sparse.linalg import splu
from scipy.linalg import pinvh
from scipy.sparse import coo_matrix
import trimesh
from trimesh.creation import cylinder
from .Timing import *

timing = Timing()


class Lattice(object):
    """
    Generate lattice structures with a lot of different parameters
    """

    def __init__(self, cell_size_x: float, cell_size_y: float, cell_size_z: float,
                 num_cells_x: int, num_cells_y: int, num_cells_z: int,
                 Lattice_Type: list[int], Radius: list[float], materialName: str,
                 gradRadiusProperty: list, gradDimProperty: list, gradMatProperty: list,
                 simMethod: int = 0, uncertaintyNode: float = 0.0,
                 periodicity: int = 0, erasedParts: list = None, randomHybrid: bool = False,
                 meshObject: "meshObject" = None, printing: bool = False, symmetryData: dict = None):
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

        Lattice_Type: list of integer
            Geometry type the cell
                (-2 => Method random cell, -1 => Full random)
                (0 => BCC, 1 => Octet, 2 => OctetExt, 3 => OctetInt, 4 => BCCZ, 5 => Cubic, 6 => OctahedronZ,
                7 => OctahedronZcross, 8 => Kelvin, 9 => Cubic formulation 2 (centered), 10 => Cubic V3, 11 => Cubic V4,
                12 => New lattice (non connu) GPT generated, 13 => Diamond, 14 => Auxetic, 15 => Hichem, 16 => Hybrid1,
                17 => Hybrid2)
        Radius: list of float
            Initial radius geometry
        materialName: string
            Name of the default material in the lattice structure ('Ti-6Al-4V', 'VeroClear'...)
            Possible to add more material in the Material.py file

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
        randomHybrid: boolean
            Randomize the hybrid lattice structure
        meshObject: meshObject
            Mesh object to check if the lattice structure is inside the mesh
        printing: boolean
            Print information about the lattice structure
        symmetryData: dictionary {"symPlane": string, "symPoint": tuple}
            Data to apply symmetry on the lattice structure
        """
        self._validate_inputs(cell_size_x, cell_size_y, cell_size_z, num_cells_x, num_cells_y, num_cells_z,
                              Lattice_Type, Radius, materialName, gradRadiusProperty, gradDimProperty, gradMatProperty,
                              simMethod, uncertaintyNode, periodicity, erasedParts,
                              randomHybrid)

        self.name = None
        self.yMin = None
        self.yMax = None
        self.xMax = None
        self.xMin = None
        self.zMax = None
        self.zMin = None
        self.latticeDimensionsDict = None
        self.dictSchurComplement = None
        self.objectifData = None

        self.cellSizeX = cell_size_x
        self.cellSizeY = cell_size_y
        self.cellSizeZ = cell_size_z
        self.numCellsX = num_cells_x
        self.numCellsY = num_cells_y
        self.numCellsZ = num_cells_z
        self.latticeType = Lattice_Type
        self.Radius = Radius
        self.materialName = materialName
        self.gradRadius = self.getGradSettings(gradRadiusProperty)
        self.gradDim = self.getGradSettings(gradDimProperty)
        self.gradMat = self.gradMaterialSetting(gradMatProperty)
        self.simMethod = simMethod
        self.sizeX, self.sizeY, self.sizeZ = self.getSizeLattice()
        self.uncertaintyNode = uncertaintyNode
        self.periodicity = periodicity  # Not finish to implemented
        self.erasedParts = erasedParts
        self.randomHybrid = randomHybrid
        self.meshObject = meshObject
        self.printing = printing

        self.cells = []
        self._nodes = []
        self._beams = []

        # Simulation necessary
        self.freeDOF = None  # Free DOF gradient conjugate gradient method
        self.maxIndexBoundary = None
        self.globalDisplacementIndex = None
        self.initialValueObjective = None
        self.initialRelativeDensityConstraint = None
        self.initialContinuityConstraint = None
        self.relativeDensityPoly = []
        self.relativeDensityPolyDeriv = []
        self.nDOFperNode = 6
        self.parameterOptimization = []
        self.krigingModelRelativeDensity = None
        self.penalizationCoefficient = 1.5  # Fixed with previous optimization
        self.meshLattice = None

        # Process
        self.generateLattice()
        if symmetryData is not None:
            self.applySymmetry(symmetryData["symPlane"], symmetryData["symPoint"])
        self.getMinMaxValues()
        self.defineBeamNodeIndex()
        self.defineCellIndex()
        self.defineCellNeighbours()
        self.setPointLocalTag()

        self.applyTagToAllPoint()

        # Case of penalization at beam near nodes
        if self.simMethod == 1:
            self.getAllAngles()
            self.setBeamNodeMod()

        # Get some data about lattice structures
        # self.getNodeData()
        # self.getBeamData()
        # print("Data retrieved")

        # Simulation FenicsX necessaries
        # Define global indexation
        self.defineNodeIndexBoundary()
        if self.printing:
            self.printStatistics()

        # Optimization necessary
        self.loadRelativeDensityModel()

    @classmethod
    def simpleLattice(cls, cell_size_x: float, cell_size_y: float, cell_size_z: float,
                      num_cells_x: int, num_cells_y: int, num_cells_z: int,
                      Lattice_Type: int, Radius: float) -> "Lattice":
        """
        Generate lattice structures with just simple parameters
        """
        Lattice_Type = [Lattice_Type]
        Radius = [Radius]
        nameMaterial = "VeroClear"
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
                   nameMaterial, gradRadiusProperty, gradDimProperty, gradMatProperty)

    @classmethod
    def hybridgeometry(cls, cell_size_x: float, cell_size_y: float, cell_size_z: float,
                       simMethod: int, Radius: list[float], latticeType: list[int], uncertaintyNode: int,
                       periodicity: bool = True) -> "Lattice":
        """
        Generate hybrid geometry structure with just some parameters
        """
        nameMaterial = "VeroClear"
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
        return cls(cell_size_x, cell_size_y, cell_size_z, 1, 1, 1, latticeType,
                   Radius, nameMaterial, gradRadiusProperty, gradDimProperty, gradMatProperty, simMethod,
                   uncertaintyNode, periodicity=periodicity)

    @classmethod
    def loadLatticeObject(cls, file_name: str = "LatticeObject", folder: str = "Saved_Lattice") -> "Lattice":
        """
        Load a lattice object from a file.

        Parameters:
        -----------
        file_name: str
            Name of the file to load (with or without the '.pkl' extension).
        folder: str
            Folder where the file is located.

        Returns:
        --------
        Lattice
            The loaded lattice object.
        """
        if not file_name.endswith(".pkl"):
            file_name += ".pkl"

        file_path = os.path.join(folder, file_name)

        if not os.path.exists(file_path):
            raise FileNotFoundError(f"The file {file_path} does not exist.")

        with open(file_path, "rb") as file:
            lattice = pickle.load(file)

        print(f"Lattice loaded successfully from {file_path}")
        return lattice

    def saveLatticeObject(self, file_name: str = "LatticeObject") -> None:
        """
        Save the lattice object to a file.

        Parameters:
        -----------
        file_name: str
            Name of the file to save (with or without the '.pkl' extension).
        """
        folder = "Saved_Lattice"
        os.makedirs(folder, exist_ok=True)

        if not file_name.endswith(".pkl"):
            file_name += ".pkl"

        file_path = os.path.join(folder, file_name)

        with open(file_path, "wb") as file:
            pickle.dump(self, file)

        print(f"Lattice saved successfully to {file_path}")

    @property
    def nodes(self):
        return self._nodes

    @property
    def beams(self):
        return self._beams

    def __repr__(self) -> str:
        string = f"Lattice name: {self.name}\n"
        string += f"Dimensions: {self.sizeX} x {self.sizeY} x {self.sizeZ}\n"
        string += f"Number of cells: {self.numCellsX} x {self.numCellsY} x {self.numCellsZ}\n"
        string += f"Cell size: {self.cellSizeX} x {self.cellSizeY} x {self.cellSizeZ}\n"
        string += f"Material: {self.materialName}\n"
        if self.randomHybrid:
            string += f"Random radius hybrid lattice\n"
        else:
            string += f"Radius: {self.Radius}\n"
        return string

    def _validate_inputs(self, cell_size_x, cell_size_y, cell_size_z,
                         num_cells_x, num_cells_y, num_cells_z,
                         Lattice_Type, Radius, materialName, gradRadiusProperty, gradDimProperty, gradMatProperty,
                         simMethod, uncertaintyNode, periodicity, erasedParts,
                         randomHybrid):

        # Check cell sizes
        assert isinstance(cell_size_x, (int, float)) and cell_size_x > 0, "cell_size_x must be a positive number"
        assert isinstance(cell_size_y, (int, float)) and cell_size_y > 0, "cell_size_y must be a positive number"
        assert isinstance(cell_size_z, (int, float)) and cell_size_z > 0, "cell_size_z must be a positive number"

        # Check number of cells
        assert isinstance(num_cells_x, int) and num_cells_x > 0, "num_cells_x must be a positive integer"
        assert isinstance(num_cells_y, int) and num_cells_y > 0, "num_cells_y must be a positive integer"
        assert isinstance(num_cells_z, int) and num_cells_z > 0, "num_cells_z must be a positive integer"

        # Check lattice type
        assert isinstance(Lattice_Type, list), "Lattice_Type must be a list"
        assert all(isinstance(lt, int) for lt in Lattice_Type), "All elements of Lattice_Type must be integers"

        # Check radius
        assert isinstance(Radius, list), "Radius must be a list"
        assert all(isinstance(r, float) for r in Radius), "All radius values must be floats"
        assert len(Radius) == len(Lattice_Type), "The number of radius must be equal to the number of lattice types"

        # Check material name
        assert isinstance(materialName, str), "materialName must be a string"

        # Check gradient properties
        assert isinstance(gradRadiusProperty, list), "gradRadiusProperty must be a list"
        assert isinstance(gradDimProperty, list), "gradDimProperty must be a list"
        assert isinstance(gradMatProperty, list), "gradMatProperty must be a list"

        # Check optional parameters
        assert isinstance(simMethod, int), "simMethod must be an integer"
        assert isinstance(uncertaintyNode, float), "uncertaintyNode must be a float"

        assert isinstance(periodicity, int), "periodicity must be an integer"

        if erasedParts is not None:
            for erasedPart in erasedParts:
                assert len(erasedPart) == 6 and all(
                    isinstance(x, float) for x in erasedPart), "erasedParts must be a list of 6 floats"
        assert isinstance(randomHybrid, bool), "randomHybrid must be a boolean"

    def getSizeLattice(self) -> list[float]:
        """
        Computes the size of the lattice along each direction.

        Return:
        ---------
        sizeLattice: list of float in dim 3
            Length of the lattice in each direction
        """
        sizeLattice = [0.0, 0.0, 0.0]
        for direction in range(3):
            total_length = 0.0
            # Get the number of cells in the current direction
            num_cells = [self.numCellsX, self.numCellsY, self.numCellsZ][direction]
            cell_size = [self.cellSizeX, self.cellSizeY, self.cellSizeZ][direction]
            gradient_factors = [grad[direction] for grad in self.gradDim[:num_cells]]

            # Calculate the total length by summing cell sizes multiplied by their respective gradient factors
            for factor in gradient_factors:
                total_length += factor * cell_size
            sizeLattice[direction] = total_length

        return sizeLattice

    @timing.timeit
    def getGradSettings(self, gradProperties: list) -> list[list[float]]:
        """
        Generate gradient settings based on the provided rule, direction, and parameters.

        Parameters:
        -----------
        gradProperties: list[Rule, Direction, Parameters]
            All types of properties for gradient definition.

        Return:
        ---------
        gradientData: list[list[float]]
            Generated gradient settings (list of lists).
        """

        def apply_rule(i: int, total_cells: int, param_value: float, rule: str) -> float:
            """
            Apply a specific gradient rule to calculate the gradient factor.

            Parameters:
            -----------
            i : int
                Current cell index in active direction.
            total_cells : int
                Total number of cells in the active direction.
            param_value : float
                Gradient parameter value.
            rule : str
                The gradient rule to apply ('constant', 'linear', etc.).

            Returns:
            --------
            float
                Calculated gradient factor.
            """
            mid = total_cells / 2
            match rule:
                case 'constant':
                    return 1.0
                case 'linear':
                    return 1.0 + i * param_value
                case 'parabolic':
                    if i < mid:
                        return 1.0 + (i / mid) * param_value
                    else:
                        return 1.0 + ((total_cells - i - 1) / mid) * param_value
                case 'sinusoide':
                    return 1.0 + param_value * math.sin((i / total_cells) * math.pi)
                case 'exponential':
                    return 1.0 + math.exp(i * param_value)
                case _:
                    raise ValueError(f"Unknown gradient rule: {rule}")

        # Extract gradient properties
        rule, direction, parameters = gradProperties

        # Determine the number of cells in each direction
        number_cells = [self.numCellsX, self.numCellsY, self.numCellsZ]

        indices = [0, 0, 0]

        gradientData = []

        for _ in range(max(number_cells)):
            gradientData.append([
                apply_rule(indices[dim], number_cells[dim], parameters[dim], rule) if direction[dim] == 1 else 1.0
                for dim in range(3)
            ])

            for dim in range(3):
                if direction[dim] == 1 and indices[dim] < number_cells[dim] - 1:
                    indices[dim] += 1
        return gradientData

    @timing.timeit
    def gradMaterialSetting(self, gradMatProperty: list) -> list:
        """
        Define gradient material settings.

        Parameters:
        ------------
        gradMatProperty: list[Multimat, GradMaterialDirection]
            Set of properties for material gradient.

        Returns:
        --------
        gradMat: list
            3D list representing the material type in the structure.
        """
        multimat, direction = gradMatProperty

        # Initialize gradMat based on `multimat` value
        if multimat == -1:  # Random materials
            return [[[random.randint(1, 3) for _ in range(self.numCellsX)] for _ in range(self.numCellsY)] for _ in
                    range(self.numCellsZ)]

        if multimat == 0:  # Single material
            return [[[1 for _ in range(self.numCellsX)] for _ in range(self.numCellsY)] for _ in range(self.numCellsZ)]

        if multimat == 1:  # Graded materials
            # Generate gradient based on the direction
            return [
                [
                    [
                        X if direction == 1 else Y if direction == 2 else Z
                        for X in range(self.numCellsX)
                    ]
                    for Y in range(self.numCellsY)
                ]
                for Z in range(self.numCellsZ)
            ]

        # Default case: return an empty gradMat if no valid `multimat` is provided
        return []

    @timing.timeit
    def isNotInErasedRegion(self, startCellPos: list[float]) -> bool:
        """
        Check if the cell is not in the erased region or inside the mesh.

        Parameters:
        -----------
        startCellPos: list of float
            (xStart, yStart, zStart) position of the cell to check.

        Returns:
        --------
        bool:
            True if the cell should be removed.
        """
        # Vérifier si le point est dans `erasedParts`
        if self.erasedParts is not None:
            for delPart in self.erasedParts:
                inside_erased = all(
                    delPart[direction] <= startCellPos[direction] <= delPart[direction] + delPart[direction + 3]
                    for direction in range(3)
                )

                if inside_erased:
                    return True  # La cellule est supprimée si elle est dans `erasedParts`

        return False  # cell removed

    def isPointInMesh(self, point):
        """
        Check if the point is inside the mesh.
        """
        if self.meshObject is not None:
            return self.meshObject.is_inside_mesh(point)

    def isCellInMesh(self, cell) -> bool:
        """
        Check if the cell is in the erased region.
        """
        cellBoundaryPoint = cell.getCellCornerCoordinates()
        for point in cellBoundaryPoint:
            self.isPointInMesh(point)
            if self.isPointInMesh(point):
                return True
        return False

    @timing.timeit
    def generateLattice(self) -> None:
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
                        if self.randomHybrid:
                            radius = [random.uniform(0.01, 0.1) for _ in self.latticeType]
                            while sum(radius) == 0:
                                radius = [random.uniform(0.01, 0.1) for _ in self.latticeType]
                        else:
                            radius = self.Radius
                        new_cell = Cell(posCell, initialCellSize, startCellPos, self.latticeType,
                                        radius, self.gradRadius, self.gradDim, self.gradMat,
                                        self.uncertaintyNode)
                        if self.meshObject is not None and self.isCellInMesh(new_cell):
                            self.cells.append(new_cell)
                        elif self.meshObject is None:
                            self.cells.append(new_cell)
                        else:
                            del new_cell
        if len(self.latticeType) > 1:
            self.checkHybridCollision()

    @timing.timeit
    def defineBeamNodeIndex(self) -> None:
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

    @timing.timeit
    def defineCellIndex(self) -> None:
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

    @timing.timeit
    def defineCellNeighbours(self) -> None:
        """
        Define neighbours for each cell in the lattice, with periodic boundaries if enabled.
        """
        cell_dict = {tuple(cell.posCell): cell for cell in self.cells}

        neighbor_offsets = [
            (-self.cellSizeX, 0, 0), (self.cellSizeX, 0, 0),
            (0, -self.cellSizeY, 0), (0, self.cellSizeY, 0),
            (0, 0, -self.cellSizeZ), (0, 0, self.cellSizeZ)
        ]

        for cell in self.cells:
            cell.neighbourCells = []
            for offset in neighbor_offsets:
                raw_pos = (
                    cell.posCell[0] + offset[0],
                    cell.posCell[1] + offset[1],
                    cell.posCell[2] + offset[2]
                )

                if self.periodicity:
                    # périodicité X
                    if raw_pos[0] < self.xMin:
                        neighbor_x = self.xMax + offset[0]
                    elif raw_pos[0] >= self.xMax:
                        neighbor_x = self.xMin
                    else:
                        neighbor_x = raw_pos[0]
                    # périodicité Y
                    if raw_pos[1] < self.yMin:
                        neighbor_y = self.yMax + offset[1]
                    elif raw_pos[1] >= self.yMax:
                        neighbor_y = self.yMin
                    else:
                        neighbor_y = raw_pos[1]
                    # périodicité Z
                    if raw_pos[2] < self.zMin:
                        neighbor_z = self.zMax + offset[2]
                    elif raw_pos[2] >= self.zMax:
                        neighbor_z = self.zMin
                    else:
                        neighbor_z = raw_pos[2]

                    neighbor_pos = (neighbor_x, neighbor_y, neighbor_z)
                else:
                    if not (self.xMin <= raw_pos[0] <= self.xMax and
                            self.yMin <= raw_pos[1] <= self.yMax and
                            self.zMin <= raw_pos[2] <= self.zMax):
                        continue
                    neighbor_pos = raw_pos
                if neighbor_pos in cell_dict:
                    cell.addCellNeighbour(cell_dict[neighbor_pos])

    def getNodeData(self) -> None:
        """
        Retrieves node data for the lattice.
        data structure: each line represents a node with data [indexNode, X, Y, Z]
        """
        self._nodes = []
        for cell in self.cells:
            for beam in cell.beams:
                for node in [beam.point1, beam.point2]:
                    if node not in self._nodes:
                        self._nodes.append(node.getData())

    def getNodeObject(self) -> list:
        """
        Retrieves node object for the lattice.

        Returns:
        --------
        nodeObjList: list
            List of node objects in the lattice.
        """
        nodeObjList = []
        for cell in self.cells:
            for beam in cell.beams:
                for node in [beam.point1, beam.point2]:
                    if node not in nodeObjList:
                        nodeObjList.append(node)
        return nodeObjList

    def getBeamData(self) -> None:
        """
        Retrieves beam data for the lattice.
        data structure: each line represents a beam with data [indexBeam, indexNode1, indexNode2, type]
        """
        self._beams = []
        for cell in self.cells:
            for beam in cell.beams:
                if beam not in self._beams:
                    self._beams.append(beam.getData())

    def getBeamObject(self) -> list:
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

    @timing.timeit
    def getListAngleBeam(self, beam: "Beam", pointbeams: list["Beam"]) -> tuple[list[float], list[float]]:
        """
        Calculate an angle between the considerate beam and beams contains in pointbeams

        Parameters:
        -----------
        beam: Beam object
            Beam where an angle is computed on
        pointbeams: list of Beam object
            List of beam to calculate an angle with considered beam

        Return:
        ---------
        non_zero_anglebeam: list of an angle between considered beam and pointbeams beam list
        non_zero_radiusbeam: list of radius between a considered beam and pointbeams beam list

        Special case when pointbeams is an empty return max angle to minimize penalization zone
        """

        def getAngleBetweenBeams(beam1, beam2):
            """
            Calculates angle between 2 beams

            Return:
            --------
            Angle: float
                angle in degrees
            """
            p1, p2 = None, None
            if self.periodicity:
                # Vérification pour les coins (tags entre 1000 et 1007)
                for idx1, p1_candidate in enumerate([beam1.point1, beam1.point2]):
                    if p1_candidate.tag and 1000 <= p1_candidate.tag[0] <= 1007:
                        p1 = idx1  # Enregistrer l'indice pour beam1
                        for idx2, p2_candidate in enumerate([beam2.point1, beam2.point2]):
                            if p2_candidate.tag and 1000 <= p2_candidate.tag[0] <= 1007:
                                p2 = idx2  # Enregistrer l'indice pour beam2
                                break  # Sortir de la boucle si une correspondance est trouvée

                # Vérification pour les arêtes (tags spécifiques)
                if p1 is None and p2 is None:  # Ne vérifie les arêtes que si aucun coin n'est trouvé
                    list_tag_edge = [[102, 104, 106, 107], [100, 108, 105, 111], [101, 109, 103, 110]]
                    for tag_list in list_tag_edge:
                        for idx1, p1_candidate in enumerate([beam1.point1, beam1.point2]):
                            if p1_candidate.tag and p1_candidate.tag[0] in tag_list:
                                p1 = idx1  # Enregistrer l'indice pour beam1
                                for idx2, p2_candidate in enumerate([beam2.point1, beam2.point2]):
                                    if p2_candidate.tag and p2_candidate.tag[0] in tag_list:
                                        p2 = idx2  # Enregistrer l'indice pour beam2
                                        break  # Sortir de la boucle si une correspondance est trouvée

                if p1 is None and p2 is None:
                    list_face_tag = [[10, 15], [11, 14], [12, 13]]
                    for face_tag in list_face_tag:
                        for idx1, p1_candidate in enumerate([beam1.point1, beam1.point2]):
                            if p1_candidate.tag and p1_candidate.tag[0] in face_tag:
                                p1 = idx1
                                for idx2, p2_candidate in enumerate([beam2.point1, beam2.point2]):
                                    if p2_candidate.tag and p2_candidate.tag[0] in face_tag:
                                        p2 = idx2
                                        break

            if beam1.point1 == beam2.point1 or (p1 == 0 and p2 == 0):
                u = beam1.point2 - beam1.point1
                v = beam2.point2 - beam2.point1
            elif beam1.point1 == beam2.point2 or (p1 == 0 and p2 == 1):
                u = beam1.point2 - beam1.point1
                v = beam2.point1 - beam2.point2
            elif beam1.point2 == beam2.point1 or (p1 == 1 and p2 == 0):
                u = beam1.point1 - beam1.point2
                v = beam2.point2 - beam2.point1
            elif beam1.point2 == beam2.point2 or (p1 == 1 and p2 == 1):
                u = beam1.point1 - beam1.point2
                v = beam2.point1 - beam2.point2
            else:
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

    @timing.timeit
    def getConnectedBeams(self, beamList: list["Beam"], beam: "Beam") -> tuple[list["Beam"], list["Beam"]]:
        """
        Get all beams connected to the interest beam.

        Parameters:
        -----------
        beam: Beam
            Beam of interest

        Returns:
        --------
        point1beams: list[Beam]
            Beams connected to point1
        point2beams: list[Beam]
            Beams connected to point2
        """
        point1beams = []
        point2beams = []

        for beamidx in beamList:
            if beam.point1 == beamidx.point1 or beam.point1 == beamidx.point2:
                if beamidx not in point1beams:
                    point1beams.append(beamidx)
            if beam.point2 == beamidx.point1 or beam.point2 == beamidx.point2:
                if beamidx not in point2beams:
                    point2beams.append(beamidx)

            # Gestion de la périodicité
            if self.periodicity and (beam.point1.tag is not [] or beam.point2.tag is not []):
                def check_periodic_connection(point, beam_idx, connected_beams, tag_range):
                    periodic = False
                    if point.tag and point.tag[0] in tag_range:
                        if any(tag in tag_range for tag in beam_idx.point1.tag) or \
                                any(tag in tag_range for tag in beam_idx.point2.tag):
                            if any(tag in range(1000, 1008) for tag in point.localTag) and (
                                    any(tag in range(1000, 1008) for tag in beam_idx.point1.localTag) or
                                    any(tag in range(1000, 1008) for tag in beam_idx.point2.localTag)
                            ):
                                periodic = True
                            edge_tags_list = [[102, 104, 106, 107], [100, 108, 105, 111], [101, 109, 103, 110]]
                            for tags in edge_tags_list:
                                if any(tag in tags for tag in point.localTag) and (
                                        any(tag in tags for tag in beam_idx.point1.localTag) or
                                        any(tag in tags for tag in beam_idx.point2.localTag)
                                ):
                                    periodic = True
                                    break
                            face_tags = [[10, 15], [11, 14], [12, 13]]
                            for tags in face_tags:
                                if any(tag in tags for tag in point.localTag) and (
                                        any(tag in tags for tag in beam_idx.point1.localTag) or
                                        any(tag in tags for tag in beam_idx.point2.localTag)
                                ):
                                    periodic = True
                                    break
                    if beam_idx not in connected_beams and periodic:
                        connected_beams.append(beam_idx)

                # Coins
                check_periodic_connection(beam.point1, beamidx, point1beams, range(1000, 1008))
                check_periodic_connection(beam.point2, beamidx, point2beams, range(1000, 1008))

                # Arêtes
                edge_tags_list = [[102, 104, 106, 107], [100, 108, 105, 111], [101, 109, 103, 110]]
                for tags in edge_tags_list:
                    check_periodic_connection(beam.point1, beamidx, point1beams, tags)
                    check_periodic_connection(beam.point2, beamidx, point2beams, tags)

                # Faces
                face_tags = [[10, 15], [11, 14], [12, 13]]
                for tags in face_tags:
                    check_periodic_connection(beam.point1, beamidx, point1beams, tags)
                    check_periodic_connection(beam.point2, beamidx, point2beams, tags)

        return point1beams, point2beams

    @timing.timeit
    def getAllAngles(self) -> None:
        """
        Calculates angles between beams in the lattice.

        Return:
        ---------
        angle:
            data structure => ((beam_index, Angle mininmum point 1, minRad1, Angle mininmum point 2, minRad2))
        """

        @timing.timeit
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
        for cell in self.cells:
            beamList = []
            cellListNeighbours = cell.getNeighbourCells()
            cellListNeighbours.append(cell)  # Include the cell itself
            for neighbour in cellListNeighbours:
                for beam in neighbour.beams:
                    if beam not in beamList:
                        beamList.append(beam)
            angleList = {}
            for beam in cell.beams:
                # Determine beams on nodes
                point1beams, point2beams = self.getConnectedBeams(beamList, beam)
                # Determine angle for all beams connected at the node
                non_zero_anglebeam1, non_zero_radiusbeam1 = self.getListAngleBeam(beam, point1beams)
                non_zero_anglebeam2, non_zero_radiusbeam2 = self.getListAngleBeam(beam, point2beams)
                # Find the lowest angle
                LAngle1, LRadius1 = findMinAngle(non_zero_anglebeam1, non_zero_radiusbeam1)
                LAngle2, LRadius2 = findMinAngle(non_zero_anglebeam2, non_zero_radiusbeam2)
                angleList[beam.index] = (LRadius1, round(LAngle1, 2), LRadius2, round(LAngle2, 2))
                beam.setAngle(angleList[beam.index])

    @timing.timeit
    def getMinMaxValues(self) -> None:
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

        self.setLatticeDimensionsDict()

    def setLatticeDimensionsDict(self):
        """
        Set lattice dimensions in a dictionary format
        """
        self.latticeDimensionsDict = {
            "xMin": self.xMin,
            "xMax": self.xMax,
            "yMin": self.yMin,
            "yMax": self.yMax,
            "zMin": self.zMin,
            "zMax": self.zMax
        }

    @timing.timeit
    def setBeamNodeMod(self) -> None:
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
                pointExt1Obj.setNodeMod(True)
                pointExt2 = beam.getPointOnBeamFromDistance(lengthMod[1], 2)
                pointExt2Obj = Point(pointExt2[0], pointExt2[1], pointExt2[2])
                pointExt2Obj.setNodeMod(True)

                b1 = Beam(beam.point1, pointExt1Obj, beam.radius, beam.material,
                          beam.type)
                b1.setBeamMod()
                b2 = Beam(pointExt1Obj, pointExt2Obj, beam.radius, beam.material, beam.type)
                b3 = Beam(pointExt2Obj, beam.point2, beam.radius, beam.material,
                          beam.type)
                b3.setBeamMod()

                beamToAdd.append((b1, b2, b3))

                beamsToRemove.append(beam)

            for addingBeam in beamToAdd:
                cell.addBeam(addingBeam)

            for beam in beamsToRemove:
                cell.removeBeam(beam)

        # Update index
        self.defineBeamNodeIndex()

    def functionPenalizationLzone(self, radius: float, angle: float) -> float:
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

    def removeCell(self, index: int) -> None:
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

    def findMinimumBeamLength(self) -> float:
        """
        Find minimum beam length

        Returns:
        --------
        minLength: float
            Length of the smallest beam in the lattice
        """
        minLength = 100000
        for cell in self.cells:
            for beam in cell.beams:
                if minLength > beam.getLength() > 0.0001 and (beam.type == 0 or beam.type == 2):
                    minLength = beam.getLength()
        return minLength

    def getTagList(self) -> list[int]:
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

    def getTagListBoundary(self) -> list[int]:
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

    @timing.timeit
    def applyTagToAllPoint(self) -> None:
        """
        Generate tag to all nodes in lattice
        """
        for cell in self.cells:
            for beam in cell.beams:
                for node in [beam.point1, beam.point2]:
                    tag = node.tagPoint(self.getLatticeBoundaryBox())
                    node.setTag(tag)

    def getLatticeBoundaryBox(self) -> list[float]:
        """
        Get the boundary box of the lattice
        """
        return [self.xMin, self.xMax, self.yMin, self.yMax, self.zMin, self.zMax]

    def getConnectedNode(self, node: "Point") -> list["Point"]:
        """
        Get all nodes connected to the input node with a beam

        Parameter:
        -----------
        node: point object

        Return:
        --------
        connectedNode: List of a point object
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

    def findBoundaryBeams(self) -> list["Beam"]:
        """
        Find boundary beams and change the type of beam

        Return:
        -------
        boundaryBeams: List of a beam object
        """
        boundaryBeams = []
        for cell in self.cells:
            for beam in cell.beams:
                if self.isNodeOnBoundary(beam.point1) or self.isNodeOnBoundary(beam.point2):
                    beam.changeBeamType(2)
                    boundaryBeams.append(beam)
        return boundaryBeams

    def isNodeOnBoundary(self, node: "Point") -> bool:
        """
        Get boolean that give information of boundary node

        Parameters:
        -----------
        node : Point object

        Returns:
        ----------
        boolean: (True if node on boundary)
        """
        return (node.x == self.xMin or node.x == self.xMax or node.y == self.yMin or node.y == self.yMax or
                node.z == self.zMin or node.z == self.zMax)

    def findBoundaryNodes(self) -> list["Point"]:
        """
        Find boundary nodes

        Returns:
        ---------
        boundaryNodes: List of a point object
        """
        boundaryNodes = []
        for cell in self.cells:
            for beam in cell.beams:
                for node in [beam.point1, beam.point2]:
                    if (node.x == self.xMin or node.x == self.xMax or node.y == self.yMin or node.y == self.yMax or
                            node.z == self.zMin or node.z == self.zMax):
                        boundaryNodes.append(node)
        return boundaryNodes

    def getName(self) -> str:
        """
        Determine the name of the lattice

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

    def checkHybridCollision(self) -> None:
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

    def getPosData(self) -> list[list[float]]:
        """
        Retrieves position data for the lattice.

        Returns:
        --------
        posData: list of list of float
            List of node positions
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

    def getEdgeIndex(self) -> list[list[int]]:
        """
        Retrieves edge index data for the lattice.

        Returns:
        --------
        edgeIndex: list of list of int
            List of edge index
        """
        edgeIndex = []
        beamAlreadyAdded = []
        for cell in self.cells:
            for beam in cell.beams:
                if beam not in beamAlreadyAdded:
                    edgeIndex.append([beam.point1.index, beam.point2.index])
                    beamAlreadyAdded.append(beam)
        return edgeIndex

    def getBeamType(self) -> list[int]:
        """
        Retrieves beam type data for the lattice.

        Returns:
        --------
        beamType: list of int
            List of beam types
        """
        beamType = []
        beamAlreadyAdded = []
        for cell in self.cells:
            for beam in cell.beams:
                if beam not in beamAlreadyAdded:
                    beamType.append([beam.type])
        return beamType

    def getAllBeamLength(self) -> list[float]:
        """
        Retrieves beam length data for the lattice.

        Returns:
        --------
        beamLength: list of float
            List of beam lengths
        """
        beamLength = []
        beamAlreadyAdded = []
        for cell in self.cells:
            for beam in cell.beams:
                if beam not in beamAlreadyAdded:
                    beamLength.append([beam.length])
        return beamLength

    def changeHybridData(self, hybridRadiusData: list[float]) -> None:
        """
        Change radius data for hybrid lattice

        Parameters:
        ------------
        hybridRadiusData: list of float
            List of radius data for hybrid lattice
        """
        if len(hybridRadiusData) != len(self.Radius):
            raise ValueError("Invalid hybrid radius data.")
        for cell in self.cells:
            for beam in cell.beams:
                if beam.modBeam:
                    beam.radius = hybridRadiusData[beam.type] * beam.penalizationCoefficient
                else:
                    beam.radius = hybridRadiusData[beam.type]

    def attractorLattice(self, PointAttractorList: list[float] = None, alpha: float = 0.5,
                         inverse: bool = False) -> None:
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

        def movePointAttracted(point, attractorPoint, alpha_coeff, inverse_bool):
            """
            Move point1 relative from attractorPoint with coefficient alpha

            Parameters:
            -----------
            point: Point object
                Point to move
            attractorPoint: Point object
                Attractor point
            alpha_coeff: float
                Coefficient of attraction
            inverse: bool
                If True, points farther away are attracted less
            """
            Length = distance(point, attractorPoint)
            if inverse_bool:
                factor = alpha_coeff / Length if Length != 0 else alpha_coeff
            else:
                factor = alpha_coeff * Length

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

    def curveLattice(self, center_x: float, center_y: float, center_z: float,
                     curvature_strength: float = 0.1) -> None:
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

    def cylindrical_transform(self, radius: float) -> None:
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

    def moveToCylinderForm(self, radius: float) -> None:
        """
        Move the lattice to a cylindrical form.

        Parameters:
        -----------
        radius: float
            Radius of the cylinder.
        """
        if radius <= self.xMax / 2:
            raise ValueError("The radius of the cylinder is too small: minimum value = ", self.xMax / 2)

        # Find moving distance
        def formula(x_coords):
            """
            Formula to calculate the new z-coordinate of the node.

            Parameters:
            -----------
            x: float
                x-coordinate of the node.
            """
            return radius - math.sqrt(radius ** 2 - (x_coords - self.xMax / 2) ** 2)

        for cell in self.cells:
            for beam in cell.beams:
                for node in [beam.point1, beam.point2]:
                    x, y, z = node.x, node.y, node.z
                    new_z = z - formula(x)
                    node.movePoint(x, y, new_z)
        self.getMinMaxValues()

    def fitToSurface(self, equation: callable, mode: str = "z", params: dict = None):
        """
        Adjust the lattice nodes to follow a surface defined by an equation.

        Parameters:
        -----------
        equation : callable
            Function representing the surface. For example, a lambda function or a normal function.
            Example: lambda x, y: x**2 + y**2 (for a paraboloid).
        mode : str
            Adjustment mode:
            - "z": Adjust nodes on a surface (z = f(x, y)).
            - "z_plan": Adjust nodes on a plan (z = f(x, y)) without changing the z-coordinate.
        params : dict
            Additional parameters for the equation or mode (e.g., radius, angle, etc.).
        """
        if params is None:
            params = {}
        nodeAlreadyChanged = []
        for cell in self.cells:
            for beam in cell.beams:
                for node in [beam.point1, beam.point2]:
                    x, y, z = node.x, node.y, node.z
                    if node not in nodeAlreadyChanged:
                        nodeAlreadyChanged.append(node)
                        # Adjust for a surface \( z = f(x, y) \)
                        if mode == "z":
                            new_z = equation(x, y, **params)
                            new_z = z + new_z
                            node.movePoint(x, y, new_z)
                        elif mode == "z_plan":
                            new_z = equation(x, y, **params)
                            node.movePoint(x, y, new_z)

                        # Other modes can be added here (e.g. cylindrical, spherical)
                        else:
                            raise ValueError(f"Mode '{mode}' non supporté.")

        # Update lattice limits after adjustment
        self.getMinMaxValues()

    def deleteDuplicatedBeams(self) -> None:
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

    def getRelativeDensityConstraint(self, relativeDensityMax, geomScheme) -> float:
        """
        Get relative density of the lattice
        """
        relativeDensity = self.getRelativeDensity(geomScheme)
        print("Relative density: ", relativeDensity)
        error = relativeDensity - relativeDensityMax
        print("Relative density maximum: ", relativeDensityMax)
        print("Relative density error: ", error)
        return error

    def getRelativeDensity(self, geomScheme=None) -> float:
        """
        Get mean relative density of all cells in lattice

        Returns:
        --------
        meanRelDens: float
            Mean relative density of the lattice
        """
        cellRelDens = []
        for cell in self.cells:
            if self.krigingModelRelativeDensity is not None:
                cellRelDens.append(cell.getRelativeDensityKriging(self.krigingModelRelativeDensity, geomScheme))
            else:
                cellRelDens.append(cell.getRelativeDensityCell())
        meanRelDens = mean(cellRelDens)
        return meanRelDens

    def defineRelativeDensityFunction(self, degree: int = 3) -> None:
        """
        Define relative density function
        Possible to define a more complex function with dependency on hybrid cells

        Parameters:
        -----------
        degree: int
            Degree of the polynomial function
        """
        if len(self.relativeDensityPoly) == 0:
            fictiveCell = Cell([0, 0, 0], [self.cellSizeX, self.cellSizeY, self.cellSizeZ],
                               [0, 0, 0], self.latticeType, self.Radius, self.gradRadius, self.gradDim, self.gradMat,
                               self.uncertaintyNode)
            domainRadius = np.linspace(0.01, 0.1, 10)
            for idxRad in range(len(self.Radius)):
                radius = np.zeros(len(self.Radius))
                relativeDensity = []
                for domainIdx in domainRadius:
                    radius[idxRad] = domainIdx
                    fictiveCell.changeBeamRadius([radius], self.gradRadius)
                    relativeDensity.append(fictiveCell.getRelativeDensityCell())
                poly_coeffs = np.polyfit(domainRadius, relativeDensity, degree).flatten()
                poly = np.poly1d(poly_coeffs)
                self.relativeDensityPoly.append(poly)
                self.relativeDensityPolyDeriv.append(poly.deriv())

    def getRelativeDensityGradient(self) -> list[float]:
        """
        Get relative density gradient of the lattice

        Returns:
        --------
        grad: list of float
            Gradient of relative density
        """
        if len(self.relativeDensityPoly) == 0:
            self.defineRelativeDensityFunction()
        if len(self.cells[0].radius) != len(self.relativeDensityPoly):
            raise ValueError("Invalid radius data.")

        grad = []
        for cell in self.cells:
            grad.append(cell.getRelativeDensityGradient(self.relativeDensityPolyDeriv))
        return grad

    def getRelativeDensityGradientKriging(self, geomScheme=None) -> list[float]:
        """
        Get relative density gradient of the lattice using kriging model

        Returns:
        --------
        grad: list of float
            Gradient of relative density
        """
        grad = []
        numberOfCells = len(self.cells)
        if geomScheme is None or len(geomScheme) != 3:
            geomScheme = [i < len(self.Radius) for i in range(3)]

        for cell in self.cells:
            gradient3Geom = cell.getRelativeDensityGradientKrigingCell(self.krigingModelRelativeDensity,
                                                                       geomScheme) / numberOfCells
            grad.extend(gradient3Geom[geomScheme])
        return grad

    def getRadiusContinuityDifference(self, delta: float = 0.01) -> list[float]:
        """
        Get the difference in radius between connected beams in the lattice

        Parameters:
        -----------
        delta: float
            Minimum difference in radius between connected cells
        """
        radiusContinuityDifference = []
        for cell in self.cells:
            radiusCell = cell.getRadius()
            for neighbours in cell.neighbourCells:
                for rad in range(len(radiusCell)):
                    radiusContinuityDifference.append((radiusCell[rad] - neighbours.getRadius()[rad]) ** 2 - delta ** 2)
        return radiusContinuityDifference

    def getRadiusContinuityJacobian(self) -> np.ndarray:
        """
        Compute the Jacobian of the radius continuity constraint.

        Returns:
        --------
        np.ndarray
            Jacobian matrix of shape (num_constraints, num_radii)
        """
        rows = []
        cols = []
        values = []
        constraint_index = 0

        for cell in self.cells:
            radiusCell = cell.getRadius()
            for neighbour in cell.neighbourCells:
                radiusNeighbour = neighbour.getRadius()
                for rad in range(len(radiusCell)):
                    i = cell.index * len(radiusCell) + rad
                    j = neighbour.index * len(radiusCell) + rad
                    diff = radiusCell[rad] - radiusNeighbour[rad]

                    rows.append(constraint_index)
                    cols.append(i)
                    values.append(2 * diff)

                    rows.append(constraint_index)
                    cols.append(j)
                    values.append(-2 * diff)

                    constraint_index += 1

        jacobian = np.zeros((constraint_index, self.getNumberParametersOptimization()))
        for r, c, v in zip(rows, cols, values):
            jacobian[r, c] = v

        return jacobian

    def changeBeamRadiusForType(self, typeToChange: int, newRadius: float) -> None:
        """
        Change radius of beam for specific type

        Parameters:
        -----------
        typeToChange: int
            Type of beam to change
        newRadius: float
            New radius of beam
        """
        for cell in self.cells:
            for beam in cell.beams:
                if beam.type == typeToChange:
                    beam.radius = newRadius

    def getNumberOfBeams(self) -> int:
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

    def getNumberOfNodes(self) -> int:
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

    def latticeInfo(self) -> None:
        """
        Print information about the lattice
        """
        self.getName()
        print("Lattice name: ", self.name)
        latticeDim = self.getSizeLattice()
        print("Lattice size X: ", latticeDim[0])
        print("Lattice size Y: ", latticeDim[1])
        print("Lattice size Z: ", latticeDim[2])

        print("Number of beams: ", self.getNumberOfBeams())
        print("Number of nodes: ", self.getNumberOfNodes())

    def applyBoundaryConditionsOnSurface(self, surfaceNames: list[str], valueDisplacement: list[float],
                                         DOF: list[int]) -> None:
        """
        Apply boundary conditions to the lattice

        Parameters:
        -----------
        surfaceNames: list[str]
            List of surfaces to apply boundary conditions (e.g., ["Xmin", "Xmax", "Ymin"])
        valueDisplacement: list of float
            Displacement value to apply to the boundary conditions
        DOF: list of int
            Degree of freedom to fix (0: x, 1: y, 2: z, 3: Rx, 4: Ry, 5: Rz)
        """
        self.applyAllConstraintsOnNodes(surfaceNames, valueDisplacement, DOF, "Displacement")

    def applyAllConstraintsOnNodes(self, surfaceNames: list[str], value: list[float], DOF: list[int],
                                   type: str = "Displacement", surfaceNamePoint: list[str] = None) -> None:
        """
        Apply boundary conditions to the lattice

        Parameters:
        -----------
        surfaceNames: list[str]
            List of surfaces to apply constraint (e.g., ["Xmin", "Xmax", "Ymin"])
        value: list of float
            Values to apply to the constraint
        DOF: list of int
            Degree of freedom to apply constraint (0: x, 1: y, 2: z, 3: Rx, 4: Ry, 5: Rz)
        type: str
            Type of constraint (Displacement, Force)

        """
        if surfaceNamePoint is None:
            pointSet = self.findPointOnLatticeSurface(surfaceNames)
        else:
            pointSet = self.findPointOnLatticeSurfaceComplex(surfaceNames, surfaceNamePoint)

        indexBoundaryList = {point.indexBoundary for point in pointSet}

        for cell in self.cells:
            for beam in cell.beams:
                for node in [beam.point1, beam.point2]:
                    if node.indexBoundary in indexBoundaryList:
                        for val, DOFi in zip(value, DOF):
                            if type == "Displacement":
                                node.setDisplacementValue(val, DOFi)
                                node.fixDOF([DOFi])
                            elif type == "Force":
                                node.setForceValue(val, DOFi)

    def findPointOnLatticeSurface(self, surfaceNames: list[str]) -> set["Point"]:
        """
        Find points on the surface of the lattice

        Parameters:
        -----------
        surfaceNames: list[str]
            List of surfaces to find points on (e.g., ["Xmin", "Xmax", "Ymin"])

        Returns:
        --------
        pointSet: set of Point objects
            Set of points found on the specified surfaces
        """
        valid_surfaces = {"Xmin", "Xmax", "Ymin", "Ymax", "Zmin", "Zmax"}

        if not all(surface in valid_surfaces for surface in surfaceNames):
            raise ValueError("Invalid surface name(s).")

        cellLists = [set(self.getCellSurface(surface)) for surface in surfaceNames]
        cellList = set.intersection(*cellLists)  # Union of all cell indices from given surfaces

        if self.cells[-1].index < max(cellList, default=-1):
            raise ValueError("Invalid cell index, some cells do not exist.")

        pointSet = None
        for cell in self.cells:
            if cell.index in cellList:
                cellPointSets = [set(cell.getPointOnSurface(surface)) for surface in surfaceNames]
                if cellPointSets:
                    if pointSet is None:
                        pointSet = set.intersection(*cellPointSets)
                    else:
                        pointSet.update(set.intersection(*cellPointSets))
        pointSet = pointSet if pointSet is not None else set()

        if pointSet == set():
            raise ValueError("No points found on the specified surfaces.")

        return pointSet

    def findPointOnLatticeSurfaceComplex(self, surfaceNamesCell: list[str], surfaceNamePoint: list[str]) \
            -> set["Point"]:
        """
        Find points on the surface of the lattice

        Parameters:
        -----------
        surfaceNames: list[str]
            List of surfaces to find points on (e.g., ["Xmin", "Xmax", "Ymin"])

        Returns:
        --------
        pointSet: set of Point objects
            Set of points found on the specified surfaces
        """
        valid_surfaces = {"Xmin", "Xmax", "Ymin", "Ymax", "Zmin", "Zmax", "Xmid", "Ymid", "Zmid"}

        if not all(surface in valid_surfaces for surface in surfaceNamesCell):
            raise ValueError("Invalid surface name(s).")

        cellLists = [set(self.getCellSurface(surface)) for surface in surfaceNamesCell]
        cellList = set.intersection(*cellLists)  # Union of all cell indices from given surfaces

        if self.cells[-1].index < max(cellList, default=-1):
            raise ValueError("Invalid cell index, some cells do not exist.")

        pointSet = None
        for cell in self.cells:
            if cell.index in cellList:
                cellPointSets = [set(cell.getPointOnSurface(surface)) for surface in surfaceNamePoint]
                if cellPointSets:
                    if pointSet is None:
                        pointSet = set.intersection(*cellPointSets)
                    else:
                        pointSet.update(set.intersection(*cellPointSets))
        pointSet = pointSet if pointSet is not None else set()

        if pointSet == set():
            raise ValueError("No points found on the specified surfaces.")

        return pointSet

    def applyBoundaryConditionsOnNode(self, nodeList: list[int], valueDisplacement: list[float],
                                      DOF: list[int]) -> None:
        """
        Apply boundary conditions to the lattice

        Parameters:
        -----------
        nodeList: list of int
            List of node index to apply boundary conditions
        valueDisplacement: float
            Displacement value to apply to the boundary conditions
        DOF: int
            Degree of freedom to fix (0: x, 1: y, 2: z, 3: Rx, 4: Ry, 5: Rz)
        """
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
                        for val, DOFi in zip(valueDisplacement, DOF):
                            node.setDisplacementValue(val, DOFi)
                            node.fixDOF([DOFi])

    def applyForceOnSurface(self, surfaceName: list[str], valueForce: list[float], DOF: list[int]) -> None:
        """
        Apply force to the lattice

        Parameters:
        -----------
        surface: str
            Surface to apply force (Xmin, Xmax, Ymin, Ymax, Zmin, Zmax)
        valueForce: list of float
            Force value to apply to the boundary conditions
        DOF: list of int
            List of degree of freedom to fix (0: x, 1: y, 2: z, 3: Rx, 4: Ry, 5: Rz)
        """
        self.applyAllConstraintsOnNodes(surfaceName, valueForce, DOF, "Force")

    def fixDOFOnSurface(self, surfaceName: list[str], dofFixed: list[int]) -> None:
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
        self.applyAllConstraintsOnNodes(surfaceName, [0.0 for _ in dofFixed], dofFixed, "Displacement")

    def fixDOFOnNode(self, nodeList: list[int], dofFixed: list[int]) -> None:
        """
        Fix degree of freedom on the surface of the lattice

        Parameters:
        -----------
        nodeList: list of int
            List of node index to apply boundary conditions
        dofFixed: list of int
            List of degree of freedom to fix (0: x, 1: y, 2: z, 3: Rx, 4: Ry, 5: Rz)
        """
        if self.getNumberOfNodes() < max(nodeList):
            raise ValueError("Invalid node index, node do not exist.")

        for node in nodeList:
            if node < 0 or node >= self.getNumberOfNodes():
                raise ValueError("Node index out of range.")

        for cell in self.cells:
            for beam in cell.beams:
                for node in [beam.point1, beam.point2]:
                    if node.index in nodeList:
                        node.fixDOF(dofFixed)

    def setRandomDisplacementOnCell(self) -> None:
        """
        Set random displacement on the lattice cells
        """
        for cell in self.cells:
            for beam in cell.beams:
                for node in [beam.point1, beam.point2]:
                    if node.indexBoundary is not None:
                        for dof in range(6):
                            if dof < 3:
                                node.setDisplacementValue(random.uniform(-0.1, 0.1), dof)
                            else:
                                node.setDisplacementValue(random.uniform(-0.01, 0.01), dof)
                            node.fixDOF([dof])

    def setDisplacementWithVector(self, displacementMatrix: list[float]) -> None:
        """
        Set displacement on the lattice with vector

        Parameters:
        -----------
        displacementMatrix: list of float of dim n_nodes*n_dofperNode
            Displacement matrix to apply to the lattice
        """
        for cell in self.cells:
            nodeInOrder = cell.getNodeOrderToSimulate()
            idxNode = 0
            for idx, node in enumerate(nodeInOrder.values()):
                if node is not None:
                    node.setDisplacementVector(displacementMatrix[idxNode])
                    node.fixDOF([i for i in range(6)])
                    idxNode += 1

    def getDisplacementGlobal(self, withFixed: bool = False, OnlyImposed: bool = False, printLevel=0) \
            -> tuple[list[float], list[int]]:
        """
        Get global displacement of the lattice

        Parameters:
        -----------
        withFixed: bool
            If True, return displacement of all nodes, else return only free degree of freedom

        Returns:
        --------
        globalDisplacement: dict
            Dictionary of global displacement with indexBoundary as key and displacement vector as value
        globalDisplacementIndex: list of int
            List of indexBoundary of the lattice
        """
        globalDisplacement = []
        globalDisplacementIndex = []
        processed_nodes = set()
        for cell in self.cells:
            for beam in cell.beams:
                for node in [beam.point1, beam.point2]:
                    if node.indexBoundary is not None and node.indexBoundary not in processed_nodes:
                        for i in range(6):
                            if node.fixedDOF[i] == 0 and not OnlyImposed:
                                globalDisplacement.append(node.displacementValue[i])
                                globalDisplacementIndex.append(node.indexBoundary)
                            elif node.fixedDOF[i] == 0 and node.appliedForce[i] == 0:
                                globalDisplacement.append(0)
                            elif withFixed or OnlyImposed:
                                globalDisplacement.append(node.displacementValue[i])
                                globalDisplacementIndex.append(node.indexBoundary)

                        processed_nodes.add(node.indexBoundary)
        if not OnlyImposed:
            self.globalDisplacementIndex = globalDisplacementIndex
        if printLevel > 2:
            print("globalDisplacement: ", globalDisplacement)
            print("globalDisplacementIndex: ", globalDisplacementIndex)
        return globalDisplacement, globalDisplacementIndex

    @timing.timeit
    def defineNodeIndexBoundary(self) -> None:
        """
        Define boundary tag for all boundary nodes and calculate the total number of boundary nodes
        """
        IndexCounter = 0
        nodeAlreadyIndexed = {}
        for cell in self.cells:
            for beam in cell.beams:
                for node in [beam.point1, beam.point2]:
                    localTag = node.tagPoint(cell.getCellBoundaryBox())
                    node.setLocalTag(localTag)
                    if localTag:
                        if node in nodeAlreadyIndexed:
                            node.setIndexBoundary(nodeAlreadyIndexed[node])
                        else:
                            nodeAlreadyIndexed[node] = IndexCounter
                            node.setIndexBoundary(IndexCounter)
                            IndexCounter += 1
        self.maxIndexBoundary = IndexCounter - 1

    def setPointLocalTag(self) -> None:
        """
        Set local tag for all points in the lattice based on their position within the cell boundary box.
        """
        for cell in self.cells:
            for beam in cell.beams:
                for node in [beam.point1, beam.point2]:
                    localTag = node.tagPoint(cell.getCellBoundaryBox())
                    node.setLocalTag(localTag)

    def getGlobalReactionForce(self, appliedForceAdded: bool = False) -> dict:
        """
        Get local reaction force of the lattice and sum if identical TagIndex

        Returns:
        --------
        globalReactionForce: dict
            Dictionary of global reaction force with indexBoundary as key and reaction force vector as value
        """
        globalReactionForce = {i: [0, 0, 0, 0, 0, 0] for i in range(self.maxIndexBoundary + 1)}
        for cell in self.cells:
            nodeIndexProcessed = set()
            for beam in cell.beams:
                for node in [beam.point1, beam.point2]:
                    if node.indexBoundary is not None and node.index not in nodeIndexProcessed:
                        globalReactionForce[node.indexBoundary] = [
                            x + y for x, y in zip(globalReactionForce[node.indexBoundary], node.getReactionForce())
                        ]
                        if appliedForceAdded:
                            for i in range(6):
                                if node.appliedForce[i] != 0:
                                    globalReactionForce[node.indexBoundary][i] = node.appliedForce[i]
                        nodeIndexProcessed.add(node.index)
        return globalReactionForce

    def getGlobalReactionForceWithoutFixedDOF(self, globalReactionForce: dict, rightHandSide: bool = False) \
            -> np.ndarray:
        """
        Get global reaction force of free degree of freedom

        Parameters:
        -----------
        globalReactionForce: dict
            Dictionary of global reaction force with indexBoundary as key and reaction force vector as value

        Returns:
        --------
        globalReactionForceWithoutFixedDOF: np.ndarray
            Array of global reaction force without fixed degree of freedom
        """
        globalReactionForceWithoutFixedDOF = []
        processed_nodes = set()
        for cell in self.cells:
            for beam in cell.beams:
                for node in [beam.point1, beam.point2]:
                    if node.indexBoundary is not None and node.indexBoundary not in processed_nodes:
                        # Append reaction force components where fixedDOF is 0
                        RFToAdd = []
                        for i in range(6):
                            if node.appliedForce[i] != 0 and rightHandSide:
                                RFToAdd.append(-node.appliedForce[i])
                                # Add a sign minus because right-hand side already with a sign minus see (b = -b)
                            elif node.fixedDOF[i] == 0:
                                RFToAdd.append(globalReactionForce[node.indexBoundary][i])
                        globalReactionForceWithoutFixedDOF.append(RFToAdd)
                        # globalReactionForceWithoutFixedDOF.append([
                        #     v1 for v1, v2 in zip(globalReactionForce[node.indexBoundary], node.fixedDOF)
                        #     if v2 == 0])
                        # print(globalReactionForceWithoutFixedDOF[-1])
                        # Mark this node as processed
                        processed_nodes.add(node.indexBoundary)
        return np.concatenate(globalReactionForceWithoutFixedDOF)

    def getFreeDOF(self):
        """
        Get total number of degrees of freedom in the lattice
        """
        self.freeDOF = 0
        processed_nodes = set()
        for cell in self.cells:
            for beam in cell.beams:
                for node in [beam.point1, beam.point2]:
                    if node.indexBoundary is not None and node.indexBoundary not in processed_nodes:
                        self.freeDOF += node.fixedDOF.count(0)
                        processed_nodes.add(node.indexBoundary)

    def setGlobalFreeDOFIndex(self) -> None:
        """
        Set global free degree of freedom index for all nodes in boundary
        """
        counter = 0
        processed_nodes = {}
        for cell in self.cells:
            for beam in cell.beams:
                for node in [beam.point1, beam.point2]:
                    if node.indexBoundary is not None:
                        if node.indexBoundary not in processed_nodes.keys():
                            for i in np.where(np.array(node.fixedDOF) == 0)[0]:
                                node.globalFreeDOFIndex[i] = counter
                                counter += 1
                            processed_nodes[node.indexBoundary] = node.globalFreeDOFIndex
                        else:
                            node.globalFreeDOFIndex[:] = processed_nodes[node.indexBoundary]

    def initializeReactionForce(self) -> None:
        """
        Initialize reaction force of all nodes to 0 on each DOF
        """
        for cell in self.cells:
            for beam in cell.beams:
                for node in [beam.point1, beam.point2]:
                    node.initializeReactionForce()

    def initializeDisplacementToZero(self) -> None:
        """
        Initialize displacement of all nodes to zero on each DOF
        """
        for cell in self.cells:
            for beam in cell.beams:
                for node in [beam.point1, beam.point2]:
                    node.initializeDisplacementToZero()

    def buildCouplingOperatorForEachCells(self) -> None:
        """
        Build coupling operator for each cell in the lattice
        """
        for cell in self.cells:
            cell.buildCouplingOperator(self.freeDOF)

    def buildLUSchurComplement(self, dictSchurComplement: dict = None) -> splu:
        """
        Build LU decomposition of the Schur complement matrix for the lattice

        Parameters:
        -----------
        schurComplementMatrix: coo_matrix
            Schur complement matrix of the lattice
        """
        from ConjugateGradientMethod.Utils_Schur import loadSchurComplement, getSref_nearest

        if dictSchurComplement is None:
            nameFileSchur = "ConjugateGradientMethod/schurComplement/Hybrid_" + str(
                len(self.latticeType)) + ".npz"
            dictSchurComplement = loadSchurComplement(nameFileSchur)

        self.buildCouplingOperatorForEachCells()
        globalSchurComplement = coo_matrix((self.freeDOF, self.freeDOF))
        for cell in self.cells:
            schurComplementMatrix = coo_matrix(getSref_nearest(cell.radius, SchurDict=dictSchurComplement,
                                                               printing=False))
            globalSchurComplement += cell.buildPreconditioner(schurComplementMatrix)

        if np.any(globalSchurComplement.sum(axis=1) == 0):
            print("Attention : There are some rows with all zeros in the Schur complement matrix.")
        cond_number = np.linalg.cond(globalSchurComplement.toarray())
        print("Condition number of the Schur complement matrix: ", cond_number)

        # Factorize preconditioner
        LUSchurComplement = None
        inverseSchurComplement = None
        if cond_number > 1e15:
            inverseSchurComplement = np.linalg.pinv(globalSchurComplement.toarray())
            # inverseSchurComplement = None
            print("Using pseudo-inverse of the Schur complement matrix.")
        else:
            globalSchurComplement = globalSchurComplement.tocsc()
            LUSchurComplement = splu(globalSchurComplement)
            print("Using LU decomposition of the Schur complement matrix.")
        return LUSchurComplement, inverseSchurComplement

    def getCellSurface(self, surface: str) -> list[int]:
        """
        Get a cell list on the surface of the lattice.

        Parameters:
        -----------
        surface: str
            Surface to get points ("Xmin", "Xmax", "Ymin", "Ymax", "Zmin", "Zmax", "Xmid", "Ymid", "Zmid")

        Returns:
        --------
        cellTagList: list of cell index
        """
        valid_surfaces = ["Xmin", "Xmax", "Ymin", "Ymax", "Zmin", "Zmax", "Xmid", "Ymid", "Zmid"]
        if surface not in valid_surfaces:
            raise ValueError("Invalid surface name.")

        mid_dict = {
            "Xmid": 0.5 * (self.xMin + self.xMax),
            "Ymid": 0.5 * (self.yMin + self.yMax),
            "Zmid": 0.5 * (self.zMin + self.zMax)
        }

        surface_dict = {
            "Xmin": ("x", self.xMin),
            "Xmax": ("x", self.xMax),
            "Ymin": ("y", self.yMin),
            "Ymax": ("y", self.yMax),
            "Zmin": ("z", self.zMin),
            "Zmax": ("z", self.zMax),
            "Xmid": ("x", mid_dict["Xmid"]),
            "Ymid": ("y", mid_dict["Ymid"]),
            "Zmid": ("z", mid_dict["Zmid"])
        }

        axis, valueSurface = surface_dict[surface]

        cellTagList = []
        for cell in self.cells:
            listPointOnSurface = cell.getPointOnSurface(surface)
            for point in listPointOnSurface:
                if getattr(point, axis) == valueSurface:
                    cellTagList.append(cell.index)
                    break

        return cellTagList

    def getRadius(self) -> float:
        """
        ################### TEMPORARY FUNCTION ###################
        Get the radius of the lattice

        Returns:
        --------
        radius: float
            Radius of the lattice
        """
        for cell in self.cells:
            for beam in cell.beams:
                return beam.radius

    def unnormalize_r(self, r_norm):
        """
        Denormalize optimization parameters

        Parameters:
        -----------
        r_norm: float
            Normalized optimization parameter
        """
        r_min, r_max = 0.01, 0.1
        return r_min + (r_max - r_min) * r_norm  # Dé-normalisation

    def setOptimizationParameters(self, optimizationParameters: list[float], geomScheme: list[bool]) -> None:
        """
        Set optimization parameters for the lattice

        Parameters:
        -----------
        optimizationParameters: list of float
            List of optimization parameters
        geomScheme: list of bool
            List of N boolean values indicating the scheme of geometry to optimize
        """
        if len(optimizationParameters) != self.getNumberParametersOptimization(geomScheme):
            raise ValueError("Invalid number of optimization parameters.")

        if geomScheme is None:
            numberOfParametersPerCell = len(self.latticeType)
        else:
            numberOfParametersPerCell = sum(geomScheme)

        for cell in self.cells:
            startIdx = cell.index * numberOfParametersPerCell
            endIdx = (cell.index + 1) * numberOfParametersPerCell
            radius = optimizationParameters[startIdx:endIdx]

            if len(radius) != len(cell.radius):
                # Reconstruct the full radius vector based on geomScheme
                full_radius = []
                i = 0  # index for radius (optimization vector)
                for keep, old in zip(geomScheme, cell.radius):
                    if keep:
                        full_radius.append(radius[i])
                        i += 1
                    else:
                        full_radius.append(old)
                radius = full_radius

            cell.changeBeamRadius(radius, self.gradRadius)

    def calculateObjective(self, typeObjective: str) -> float:
        """
        Calculate objective function for the lattice optimization

        Parameters
        ----------
        typeObjective: str
            Type of objective function to calculate (Compliance...)

        Returns
        -------
        objectiveValue: float
            Objective function value
        """
        if typeObjective == "Compliance":
            reactionForce = self.getGlobalReactionForce(appliedForceAdded=True)
            reaction_force_array = np.array(list(reactionForce.values())).flatten()
            displacement = np.array(self.getDisplacementGlobal(OnlyImposed=True)[0])
            objective = 0.5 * np.dot(reaction_force_array, displacement)
            if self.printing > 2:
                np.set_printoptions(threshold=np.inf)
                print("Reaction force: ", reaction_force_array[displacement != 0])
                print("Displacement: ", displacement[displacement != 0])
                print("Compliance: ", objective)
        elif typeObjective == "Displacement":
            setNode = self.findPointOnLatticeSurface(surfaceNames=self.objectifData["surface"])
            displacements = []
            for node in setNode:
                for dof in self.objectifData["DOF"]:
                    if dof < 0 or dof > 5:
                        raise ValueError("Invalid degree of freedom index.")
                    displacements.append(node.displacementValue[dof])
            displacements = np.array(displacements)
            objective = sum(abs(displacements)) / len(displacements)

        elif typeObjective == "Stiffness":
            pass
        return objective

    def getNumberParametersOptimization(self, geomScheme) -> int:
        """
        Get number of parameters for optimization

        Returns:
        --------
        numParameters: int
            Number of parameters for optimization
        """
        numParameters = 0
        for cell in self.cells:
            if geomScheme is None:
                numParameters += len(cell.radius)
            else:
                numParameters += sum(geomScheme)
        return numParameters

    def applyReactionForceOnNodeList(self, reactionForce: list, nodeCoordinatesList: list):
        """
        Apply reaction force on node list

        Parameters:
        -----------
        reactionForce: list of float
            Reaction force to apply
        nodeCoordinatesList: list of float
            Coordinates of the node
        """
        nodeCoordinatesArray = np.array(nodeCoordinatesList)

        for cell in self.cells:
            for beam in cell.beams:
                for node in [beam.point1, beam.point2]:
                    nodeCoord = np.array([node.x, node.y, node.z])
                    match = np.all(nodeCoordinatesArray == nodeCoord, axis=1)
                    if np.any(match):
                        index = np.where(match)[0][0]
                        node.setReactionForce(reactionForce[index])

    def applyDisplacementOnNodeList(self, displacement: list, nodeCoordinatesList: list):
        """
        Apply displacement on node list

        Parameters:
        -----------
        displacement: list of float
            Displacement to apply
        nodeCoordinatesList: list of float
            Coordinates of the node
        """
        nodeCoordinatesArray = np.array(nodeCoordinatesList)

        for cell in self.cells:
            for beam in cell.beams:
                for node in [beam.point1, beam.point2]:
                    nodeCoord = np.array([node.x, node.y, node.z])
                    match = np.all(nodeCoordinatesArray == nodeCoord, axis=1)
                    if np.any(match):
                        index = np.where(match)[0][0]
                        node.setDisplacementVector(displacement[index])
                        node.fixDOF([i for i in range(6)])

    def expandSchurToFullBasis(self, SchurReduced, nodeInOrder):
        """
        Expand the reduced Schur complement matrix into the full 156-DOF boundary space.

        Parameters:
        -----------
        SchurReduced : np.ndarray
            The Schur complement matrix computed only for the active boundary nodes.
        nodeInOrder : dict
            Dictionary with node tags as keys and corresponding node objects as values.

        Returns:
        --------
        SchurFull : np.ndarray
            The expanded Schur complement matrix in the full 156-DOF boundary space.
        """

        num_total_boundary_nodes = len(nodeInOrder)  # Nombre total de nœuds frontière
        total_dofs = num_total_boundary_nodes * 6  # Chaque nœud a 6 DOFs
        SchurFull = np.zeros((total_dofs, total_dofs))  # Matrice complète initialisée à zéro

        # Construire un mapping entre les indices de SchurReduced et les indices globaux
        boundary_dof_map = {}  # Associe index local -> index global dans SchurFull
        dof_counter = 0  # Compteur pour attribuer des indices locaux à la matrice réduite

        for node_idx, node in enumerate(nodeInOrder.values()):
            if node is not None:  # Le nœud est utilisé
                for dof in range(6):  # Chaque nœud a 6 DOFs
                    boundary_dof_map[dof_counter] = node_idx * 6 + dof
                    dof_counter += 1

        # Remplir SchurFull en utilisant les indices mappés
        for i in range(SchurReduced.shape[0]):
            global_i = boundary_dof_map[i]  # Trouver l'index global correspondant
            for j in range(SchurReduced.shape[1]):
                global_j = boundary_dof_map[j]
                SchurFull[global_i, global_j] = SchurReduced[i, j]

        return SchurFull

    def printStatistics(self):
        """
        Print statistics about the lattice
        """
        print("Number of cells: ", len(self.cells))
        # print("Number of beams: ", self.getNumberOfBeams())
        # print("Number of nodes: ", self.getNumberOfNodes())
        print("Relative density: ", self.getRelativeDensity())
        radMax, radMin = self.getRadiusMinMax()
        print("Radius max: ", radMax)
        print("Radius min: ", radMin)

    def getRadiusMinMax(self):
        """
        Get the maximum and minimum radius of the lattice

        Returns:
        --------
        radMax: float
            Maximum radius of the lattice
        radMin: float
            Minimum radius of the lattice
        """
        radMin = 1000000
        radMax = 0
        for cell in self.cells:
            for beam in cell.beams:
                if not beam.modBeam:
                    radMin = min(radMin, beam.radius)
                    radMax = max(radMax, beam.radius)
        return radMax, radMin

    def addMeshObject(self, meshObject):
        """
        Add a mesh object to the lattice

        Parameters:
        -----------
        meshObject: MeshObject
            Mesh object to add to the lattice
        """
        self.meshObject = meshObject

    def cutBeamsAtMeshIntersection(self):
        """
        Cut beams at the intersection with the mesh
        """
        if self.meshObject is None:
            raise ValueError("A mesh object must be assigned to the lattice before cutting beams.")

        new_beams = []
        beams_to_remove = []

        for cell in self.cells:
            for beam in cell.beams:
                if not beam.modBeam:
                    p1_inside = self.meshObject.is_inside_mesh([beam.point1.x, beam.point1.y, beam.point1.z])
                    p2_inside = self.meshObject.is_inside_mesh([beam.point2.x, beam.point2.y, beam.point2.z])

                    if not p1_inside and not p2_inside:
                        # The Beam is outside the mesh, remove it
                        beams_to_remove.append(beam)
                    elif not p1_inside or not p2_inside:
                        # The Beam intersects the mesh, cut it
                        intersection_point = beam.findIntersectionWithMesh(self.meshObject)
                        if intersection_point is not None:
                            new_point = Point(intersection_point[0], intersection_point[1], intersection_point[2])

                            if not p1_inside:
                                new_beam = Beam(new_point, beam.point2, beam.radius, beam.material, beam.type)
                            else:
                                new_beam = Beam(beam.point1, new_point, beam.radius, beam.material, beam.type)

                            new_beams.append(new_beam)
                            beams_to_remove.append(beam)
                else:
                    raise ValueError("Cutting is only available for non modified lattice.")
            # Apply changes
            for beam in beams_to_remove:
                cell.removeBeam(beam)
            for beam in new_beams:
                cell.addBeam(beam)

            new_beams = []
            beams_to_remove = []

    @timing.timeit
    def applySymmetry(self, symmetry_plane: str, reference_point: tuple = (0, 0, 0)) -> None:
        """
        Apply symmetry to the lattice structure based on a reference point.

        Parameters:
        -----------
        symmetry_plane : str
            The plane of symmetry, can be "XY", "XZ", "YZ" (default symmetries),
            or "X", "Y", "Z" for symmetry about a specific coordinate.
        reference_point : tuple (x_ref, y_ref, z_ref), optional
            The reference point for the symmetry. Defaults to (0,0,0).
        """

        if symmetry_plane not in ["XY", "XZ", "YZ", "X", "Y", "Z"]:
            raise ValueError("Invalid symmetry plane. Choose from 'XY', 'XZ', 'YZ', 'X', 'Y', or 'Z'.")

        x_ref, y_ref, z_ref = reference_point
        new_cells = []
        node_map = {}

        for cell in self.cells:
            new_pos = list(cell.posCell)
            new_start_pos = list(cell.coordinateCell)
            mirrored_beams = []

            for beam in cell.beams:
                new_point1 = Point(beam.point1.x, beam.point1.y, beam.point1.z)
                new_point2 = Point(beam.point2.x, beam.point2.y, beam.point2.z)

                # Apply symmetry transformation based on the selected plane
                if symmetry_plane == "XY":
                    new_point1.z = 2 * z_ref - new_point1.z
                    new_point2.z = 2 * z_ref - new_point2.z
                    new_start_pos[2] = 2 * z_ref - new_start_pos[2]
                elif symmetry_plane == "XZ":
                    new_point1.y = 2 * y_ref - new_point1.y
                    new_point2.y = 2 * y_ref - new_point2.y
                    new_start_pos[1] = 2 * y_ref - new_start_pos[1]
                elif symmetry_plane == "YZ":
                    new_point1.x = 2 * x_ref - new_point1.x
                    new_point2.x = 2 * x_ref - new_point2.x
                    new_start_pos[0] = 2 * x_ref - new_start_pos[0]
                elif symmetry_plane == "X":
                    new_point1.x = 2 * x_ref - new_point1.x
                    new_point2.x = 2 * x_ref - new_point2.x
                    new_start_pos[0] = 2 * x_ref - new_start_pos[0]
                elif symmetry_plane == "Y":
                    new_point1.y = 2 * y_ref - new_point1.y
                    new_point2.y = 2 * y_ref - new_point2.y
                    new_start_pos[1] = 2 * y_ref - new_start_pos[1]
                elif symmetry_plane == "Z":
                    new_point1.z = 2 * z_ref - new_point1.z
                    new_point2.z = 2 * z_ref - new_point2.z
                    new_start_pos[2] = 2 * z_ref - new_start_pos[2]

                # Ensure uniqueness of nodes
                if new_point1 not in node_map:
                    node_map[new_point1] = new_point1
                if new_point2 not in node_map:
                    node_map[new_point2] = new_point2

                mirrored_beams.append(
                    Beam(node_map[new_point1], node_map[new_point2], beam.radius, beam.material, beam.type))

            # Create a new mirrored cell
            new_cell = Cell(new_pos, cell.cellSize, new_start_pos, cell.latticeType,
                            cell.radius, cell.gradRadius, cell.gradDim, cell.gradMat, cell.uncertaintyNode)

            new_cell.beams = mirrored_beams
            new_cells.append(new_cell)

        self.cells.extend(new_cells)
        self.getMinMaxValues()  # Recalculate the lattice boundaries

    @timing.timeit
    def loadRelativeDensityModel(self, model_path="Lattice/Saved_Lattice/RelativeDensityKrigingModel.pkl"):
        """
        Load the relative density model from a file

        Parameters:
        -----------
        model_path: str
            Path to the model file

        Returns:
        --------
        model: Kriging
            The loaded model
        """
        if not os.path.exists(model_path):
            print(f"Model file not found: {model_path}")
        else:
            gpr = joblib.load(model_path)
            self.krigingModelRelativeDensity = gpr

    def deleteBeamsUnderThreshold(self, threshold: float = 0.01) -> None:
        """
        Delete beams with radius under a certain threshold

        Parameters:
        -----------
        threshold: float
            Threshold value for beam radius
        """
        for cell in self.cells:
            beamsToRemove = []
            for beam in cell.beams:
                if beam.radius <= threshold:
                    beamsToRemove.append(beam)
            for beam in beamsToRemove:
                cell.removeBeam(beam)

    def deleteBeamsGeomScheme(self, geomScheme: list[bool]) -> None:
        """
        Delete beams based on the geometry scheme.
        If geomScheme[i] is False, the beam of type i will be removed.
        Usefull for hybrid lattices geometry.

        Parameters:
        -----------
        geomScheme: list of bool
            List of N boolean values indicating the scheme of geometry to optimize
            If geomScheme[i] is False, the beam of type i will be removed.
        """
        for cell in self.cells:
            beamsToRemove = []
            for beam in cell.beams:
                if not geomScheme[beam.type]:
                    beamsToRemove.append(beam)
            for beam in beamsToRemove:
                cell.removeBeam(beam)

    def setObjectiveData(self, objectifData: dict) -> None:
        """
        Add objective data to the lattice

        Parameters:
        -----------
        objectifData: dict
            Dictionary containing objective data
        """
        self.objectifData = objectifData

    @timing.timeit
    def generateMeshLattice(self, sectionPrecision: int = 8, cutMeshAtBoundary: bool = False):
        """
        Generate a mesh representation of the lattice structure.

        Parameters:
        -----------
        sectionPrecision: int
            Number of sections for the cylindrical mesh representation of beams.
        cutMeshAtBoundary: bool
            If True, cut the mesh at the bounding box of the lattice.
        """
        mesh_list = []

        bounds_min = np.array([self.xMin, self.yMin, self.zMin])
        bounds_max = np.array([self.xMax, self.yMax, self.zMax])
        center = (bounds_min + bounds_max) / 2
        size = bounds_max - bounds_min

        bounding_box = trimesh.creation.box(extents=size, transform=trimesh.transformations.translation_matrix(center))

        for cell in self.cells:
            for beam in cell.beams:
                p1 = np.array(beam.point1.getPos())
                p2 = np.array(beam.point2.getPos())
                height = beam.length

                direction = p2 - p1
                direction /= height

                cyl = cylinder(radius=beam.radius, height=height, sections=sectionPrecision, cap_ends=True)

                R = trimesh.geometry.align_vectors([0, 0, 1], direction)
                T = trimesh.transformations.translation_matrix((p1 + p2) / 2)

                cyl.apply_transform(R)
                cyl.apply_transform(T)

                if cutMeshAtBoundary:
                    cyl = trimesh.boolean.intersection([cyl, bounding_box])
                mesh_list.append(cyl)

        self.meshLattice = trimesh.boolean.union(mesh_list)
        print(f"Mesh lattice has been constructed")

    def are_cells_identical(self) -> bool:
        """
        Check if all cells in the list are identical based on their attributes and beams.
        Print the result.
        Possible upgrade could be to use a more sophisticated comparison method (Only beam length is checked for now).
        """
        if len(self.cells) < 2:
            print(Fore.GREEN + "Only one or no cell: considered identical." + Style.RESET_ALL)
            return True

        reference = self.cells[0]
        attrs_to_check = [
            "latticeType",
            "radius",
            "cellSize",
            "gradRadius",
            "gradDim",
        ]

        for i, cell in enumerate(self.cells[1:], start=1):
            for attr in attrs_to_check:
                if not np.array_equal(getattr(reference, attr), getattr(cell, attr)):
                    print(
                        Fore.RED + f"Difference found in attribute '{attr}' between cell 0 and cell {i}" + Style.RESET_ALL)
                    return False

            if len(reference.beams) != len(cell.beams):
                print(Fore.RED + f"Different number of beams between cell 0 and cell {i}" + Style.RESET_ALL)
                return False

            for j, (b1, b2) in enumerate(zip(reference.beams, cell.beams)):
                if not b1.is_identical_to(b2):
                    print(Fore.RED + f"Difference found in beam {j} between cell 0 and cell {i}" + Style.RESET_ALL)
                    return False

        print(Fore.GREEN + "All cells are identical." + Style.RESET_ALL)
        return True
