import numpy as np

from .Point import *
from .Beam import *
from .Geometry_Lattice import Lattice_geometry

from scipy.sparse import coo_matrix


class Cell(object):
    """
    Define Cell data for lattice structure
    """

    def __init__(self, posCell: list, initialCellSize: list, startCellPos: list, latticeType: list[int],
                 Radius: list[float], gradRadius: list, gradDim: list, gradMat: list, uncertaintyNode: float = 0.0):
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
        uncertaintyNode: float
            Standard deviation for adding uncertainty to node coordinates. Defaults to 0.0.
        """
        if len(Radius) == 1:
            if latticeType[0] == 0:
                self.originalTags = [1000, 1001, 1002, 1003, 1004, 1005, 1006, 1007]
                self.originalCellGeom = [0, 0, 0, 0, 0, 0, 0, 0]
            elif latticeType[0] == 16:
                self.originalTags = [100, 101, 102, 103, 104, 105, 106, 107, 108, 109, 110, 111]
                self.originalCellGeom = [1, 1, 1, 1, 1, 1, 1, 1, 1, 1, 1, 1]
            elif latticeType[0] == 19:
                self.originalTags = [10, 11, 12, 13, 14, 15]
                self.originalCellGeom = [2, 2, 2, 2, 2, 2,]
            else:
                raise ValueError("Lattice type not recognized")
        elif len(Radius) == 2:
            if latticeType[0] == 0 and latticeType[1] == 16:
                self.originalTags = [1000, 1001, 1002, 1003, 1004, 1005, 1006, 1007,
                                     100, 101, 102, 103, 104, 105, 106, 107, 108, 109, 110, 111]
                self.originalCellGeom = [0, 0, 0, 0, 0, 0, 0, 0,
                                         1, 1, 1, 1, 1, 1, 1, 1, 1, 1, 1, 1]
            else:
                raise ValueError("Lattice type not recognized")
        else:
            self.originalTags = [1000, 1001, 1002, 1003, 1004, 1005, 1006, 1007,
                                 10, 11, 12, 13, 14, 15,
                                 100, 101, 102, 103, 104, 105, 106, 107, 108, 109, 110, 111]
            self.originalCellGeom = [0, 0, 0, 0, 0, 0, 0, 0,
                                     2, 2, 2, 2, 2, 2,
                                     1, 1, 1, 1, 1, 1, 1, 1, 1, 1, 1, 1]

        self.centerPoint = None
        self._beamMaterial = None
        self.cellSize = None
        self.posCell = posCell
        self.coordinateCell = startCellPos
        self.beams = []
        self.index = None
        self.latticeType = latticeType
        self.radius = Radius
        self.matB = None  # B matrix (Coupling matrix)
        self.uncertaintyNode = uncertaintyNode
        self.gradRadius = gradRadius
        self.gradMat = gradMat
        self.gradDim = gradDim
        self.neighbourCells = []

        idxCell = 0
        for idx, rad in enumerate(self.radius):
            if rad > 0.0:
                if idxCell == 0:
                    self.getBeamMaterial(gradMat)
                    beamRadius = self.getBeamRadius(gradRadius, rad)
                    self.getCellSize(initialCellSize, gradDim)
                    self.generateBeamsInCell(self.latticeType[idx], startCellPos, beamRadius, idx)
                    self.getCellCenter(startCellPos)
                else:
                    hybridRadius = self.getBeamRadius(gradRadius, rad)
                    self.generateBeamsInCell(self.latticeType[idx], startCellPos, hybridRadius, idx)
                idxCell += 1

    def __repr__(self) -> str:
        return f"Cell(Coordinates:{self.coordinateCell}, Size: {self.cellSize}, Index:{self.index})"


    def generateBeamsInCell(self, latticeType: int, startCellPos: list, beamRadius: float, beamType: int = 0) -> None:
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
        pointDict = {}
        for line in Lattice_geometry(latticeType):
            x1, y1, z1, x2, y2, z2 = map(float, line)
            if (x1, y1, z1) in pointDict:
                point1 = pointDict[(x1, y1, z1)]
            else:
                point1 = Point(x1 * self.cellSize[0] + startCellPos[0], y1 * self.cellSize[1] + startCellPos[1],
                               z1 * self.cellSize[2] + startCellPos[2], self.uncertaintyNode)
                pointDict[(x1, y1, z1)] = point1
            if (x2, y2, z2) in pointDict:
                point2 = pointDict[(x2, y2, z2)]
            else:
                point2 = Point(x2 * self.cellSize[0] + startCellPos[0], y2 * self.cellSize[1] + startCellPos[1],
                               z2 * self.cellSize[2] + startCellPos[2], self.uncertaintyNode)
                pointDict[(x2, y2, z2)] = point2
            beam = Beam(point1, point2, beamRadius, self._beamMaterial, beamType)
            self.beams.append(beam)

    def getBeamMaterial(self, gradMat: list) -> None:
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

    def getBeamRadius(self, gradRadius: list, BaseRadius: float) -> float:
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
        beamRadius = (BaseRadius * gradRadius[self.posCell[0]][0] * gradRadius[self.posCell[1]][1] *
                      gradRadius[self.posCell[2]][2])
        return beamRadius

    def getCellSize(self, initialCellSize: list, gradDim: list) -> None:
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

    def getCellCenter(self, startCellPos: list) -> None:
        """
        Calculate the center point of the cell
        """
        self.centerPoint = [startCellPos[i] + self.cellSize[i] / 2 for i in range(3)]

    def getAllPoints(self) -> list:
        """
        Determine list of points in cell
        """
        pointList = []
        for beam in self.beams:
            for point in [beam.point1, beam.point2]:
                if point not in pointList:
                    pointList.append(point)
        return pointList

    def removeBeam(self, beamToDelete: "Beam") -> None:
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

    def addBeam(self, beamToAdd: "Beam") -> None:
        """
        Adding beam to cell
        """
        if isinstance(beamToAdd, Beam):
            self.beams.append(beamToAdd)
        elif isinstance(beamToAdd, tuple):
            for beam in beamToAdd:
                self.beams.append(beam)
        else:
            raise ValueError("Invalid beam type")

    def setIndex(self, index: int) -> None:
        """
        Set cell index
        """
        self.index = index

    def getPointOnSurface(self, surfaceName: str) -> list:
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
                "Surface " + str(surfaceName) + " is not valid. Choose from 'Xmin', 'Xmax', 'Ymin', 'Ymax', 'Zmin', "
                                                "or 'Zmax'.")

        surface_value = surface_map[surfaceName]
        coord_index = {"Xmin": "x", "Xmax": "x", "Ymin": "y", "Ymax": "y", "Zmin": "z", "Zmax": "z"}

        return [point for beam in self.beams for point in [beam.point1, beam.point2] if
                getattr(point, coord_index[surfaceName]) == surface_value]

    def getRadius(self) -> list[float]:
        """
        Get the radius of the beam
        """
        return self.radius

    def getNodeOrderToSimulate(self) -> dict:
        """
        Get the order of nodes to simulate in the cell
        """
        tag_dict = {tag: None for tag in self.originalTags}

        for beam in self.beams:
            if beam.radius > 0:
                for point in [beam.point1, beam.point2]:
                    if point.indexBoundary is not None:
                        tag = point.localTag
                        if tag:  # Ensure tags is not an empty list
                            if tag[0] in self.originalTags:
                                tag_dict[tag[0]] = point
        return tag_dict

    def getNodesOrderNN(self, nodeInOrder: dict, numberRadiusNN) -> dict:
        """
        Get the nodes order for the neural network

        Parameters:
        -----------
        nodeInOrder: dict
            Dictionary of nodes in order
        originalCellGeom: list
            Original cell geometry
        """
        tag_dictNN = nodeInOrder.copy()
        idx = 0
        for i, key in nodeInOrder.items():
            if key:
                tag_dictNN[i] = 1
            # elif self.originalCellGeom[idx] < numberRadiusNN:
            #     tag_dictNN[i] = 1
            else:
                tag_dictNN[i] = 0
            idx += 1
        return tag_dictNN


    def setDisplacementAtBoundaryNodes(self, displacementArray: list, displacementIndex: list, printLevel = 0) -> None:
        """
        Set displacement at nodes.

        Parameters:
        ------------
        displacementArray: list or array-like
            Flattened array of displacement values.
        displacementIndex: array of int
            Boundary node index of each displacement value.
        """
        if printLevel > 1:
            print("Displacement array", displacementArray)
            print("Displacement index", displacementIndex)
            print("Non-zero displacements:", displacementArray[displacementArray != 0])
        for beam in self.beams:
            for point in [beam.point1, beam.point2]:
                if point.indexBoundary is not None and point.indexBoundary in displacementIndex:
                    index = displacementIndex.index(point.indexBoundary)
                    indexActual = 0
                    for i in range(6):
                        if point.fixedDOF[i] == 0:  # Filter out the fixed DOF
                            point.setDisplacementValue(displacementArray[index + indexActual], i)
                            indexActual += 1

    def getDisplacementAtNodes(self, nodeList: dict, numberRadiusNN: int) -> list:
        """
        Get the displacement at nodes.

        Parameters:
        -----------
        nodeList: list of Point objects
            List of nodes to get the displacement.
        numberRadiusNN: int
            Number of radius for the neural network

        Returns:
        --------
        list
            A flattened list of displacement values.
        """
        nodeListNN = self.getNodesOrderNN(nodeList, numberRadiusNN)

        displacementList = []
        nullDisplacement = [0.0,0.0,0.0,0.0,0.0,0.0]
        for key, node in nodeList.items():
            if node:
                displacement = node.getDisplacementValue()
                displacementList.append(displacement)
            elif nodeListNN[key] == 1:
                displacementList.append(nullDisplacement)
        return displacementList

    def setReactionForceOnEachNodes(self, nodeList: dict, reactionForce: list) -> None:
        """
        Set reaction force on each node.

        Parameters:
        -----------
        nodeList: dict
            Dictionary mapping node tags to Point objects.
        reactionForce: list
            List of reaction force values corresponding to the nodes.
        """
        # Vérification que la longueur des listes correspond
        if len([v for v in nodeList.values() if v is not None]) != len(reactionForce):
            print("Lenght nodeList ",len([v for v in nodeList.values() if v is not None])
                  , "Lenght reactionForce ", len(reactionForce))
            print("nodeList", nodeList)
            print("reactionForce", reactionForce)
            raise ValueError("Mismatch: nodeList and reactionForce must have the same length.")

        # for beam in self.beams:
        #     for point in [beam.point1, beam.point2]:
        #         print(point, point.indexBoundary)
        #         if point.indexBoundary is not None:
        #             index = list(nodeList.keys()).index(point.localTag[0])
        #             point.setReactionForce(reactionForce[index])
        idx = 0
        for node in nodeList:
            if nodeList[node]:
                nodeList[node].setReactionForce(reactionForce[idx])
                idx += 1

    def getNumberOfBoundaryNodes(self) -> int:
        """
        Get the number of boundary nodes in the cell
        """
        return len(
            [beam for beam in self.beams for point in [beam.point1, beam.point2] if point.indexBoundary is not None])

    def buildCouplingOperator(self, nbFreeDOF: int) -> None:
        """
        Build the coupling operator for the cell

        Parameters:
        -----------
        nbFreeDOF: int
            Number of free degrees of freedom
        """
        data = []
        row, col = [], []
        listBndNodes = []
        for beam in self.beams:
            for point in [beam.point1, beam.point2]:
                if point.indexBoundary is not None and point.indexBoundary not in listBndNodes:
                    localNodeIndex = self.originalTags.index(point.localTag[0])
                    listBndNodes.append(point.indexBoundary)
                    for i in range(6):
                        if point.globalFreeDOFIndex[i] is not None:
                            data.append(1)
                            col.append(localNodeIndex * 6 + i)
                            row.append(point.globalFreeDOFIndex[i])
        nbBndDOFloc = len(listBndNodes) * 6
        shapeB = (nbFreeDOF, nbBndDOFloc)
        self.matB = coo_matrix((data, (row, col)), shape=shapeB)

    def buildPreconditioner(self, SchurMatrix: "coo_matrix") -> "coo_matrix":
        """
        Build the preconditioner part for the cell

        Parameters:
        -----------
        SchurMatrix: coo_matrix
            Schur matrix
        """
        if self.matB is None:
            raise ValueError("Coupling matrix has not been built yet. Please build it first.")
        if self.matB.shape[1] != SchurMatrix.shape[0]:
            print("Shape of B matrix", self.matB.shape)
            print("Shape of Schur matrix", SchurMatrix.shape)
            raise ValueError("Incompatible dimensions between the coupling matrix and the Schur matrix.")

        return self.matB @ SchurMatrix @ self.matB.transpose()

    def getInternalEnergy(self) -> float:
        """
        Get the internal energy of the cell
        """
        internalEnergy = 0
        for beam in self.beams:
            for point in [beam.point1, beam.point2]:
                if point.indexBoundary is not None:
                    pointEnergy = point.calculatePointEnergy()
                    # if pointEnergy < 0:
                    #     print("Negative energy", pointEnergy)
                    internalEnergy += pointEnergy
        return internalEnergy

    def getDisplacementData(self) -> list:
        """
        Build and return displacement data on cell for dataset generation
        """
        allBoundaryDisplacementData = []
        for beam in self.beams:
            for point in [beam.point1, beam.point2]:
                if point.indexBoundary is not None:
                    allBoundaryDisplacementData.append(point.getDisplacementValue())
        return allBoundaryDisplacementData

    def changeBeamRadius(self, newRadius: list, gradRadius: list = None) -> None:
        """
        ATTENTION /!\ : BEAM MOD IS NOT WORKING
        Change beam radius in the cell

        Parameters:
        -----------
        newRadius: list
            beam radius wanted to assign
        hybridData: list
            Hybrid data type
        gradRadius: list
            Gradient of the radius
        penalizationCoeff: float
        """
        assert len(newRadius) == len(self.radius), ("Length of new radius vector and already cell radius vector needs "
                                                    "to be equal ")
        beamRadius = []
        for rad in newRadius:
            beamRadius.append(self.getBeamRadius(gradRadius, rad))

        for beam in self.beams:
            if beam.modBeam:
                beam.setRadius(beamRadius[beam.type] * beam.penalizationCoefficient)
            else:
                beam.setRadius(beamRadius[beam.type])

        self.radius = newRadius

    def getVolumeCell(self) -> float:
        """
        Get the volume of the cell
        """
        return self.cellSize[0] * self.cellSize[1] * self.cellSize[2]

    def getRelativeDensityCell(self) -> float:
        """
        Get the relative density of the cell
        """
        volumeBeams = 0
        for beam in self.beams:
            volumeBeams += beam.getVolume()
        return volumeBeams / self.getVolumeCell()

    def getVolumeGeomSeparated(self) -> list:
        volumes = np.zeros(len(self.radius))
        for beam in self.beams:
            if not beam.modBeam:
                volumeBeam = beam.getVolume()
                volumes[beam.type] += volumeBeam
        return volumes

    def getRelativeDensityKriging(self, krigingModel, geomScheme=None) -> float:
        """
        Get the relative density of the cell using kriging model

        Parameters:
        -----------
        krigingModel: Kriging
            Kriging model to use for prediction
        """
        radii = np.zeros(3)
        if geomScheme is not None:
            true_indices = [i for i, val in enumerate(geomScheme) if val]
            for idx, val in zip(true_indices, self.radius):
                radii[idx] = val
        else:
            for idx, rad in enumerate(self.radius):
                radii[idx] = rad
        radii = np.array(radii).reshape(-1, 3)
        relativeDensity = krigingModel.predict(radii)[0]
        return relativeDensity

    def getRelativeDensityGradient(self, relativeDensityPolyDeriv) -> float:
        """
        Get the gradient of the relative density

        Parameters:
        -----------
        relativeDensityPolyDeriv: list
            List of polynomial derivative functions

        Returns:
        --------
        deriv: float
            Derivative of the relative density
        """
        deriv = 0
        for idx, polyDeriv in enumerate(relativeDensityPolyDeriv):
            deriv += polyDeriv(self.radius[idx])
        return deriv

    def getRelativeDensityGradientKrigingCell(self, gpr,geomScheme=None) -> np.ndarray:
        """
        Retourne le gradient de la fonction volume par rapport aux rayons (dérivée partielle).

        Paramètres :
        ------------
        gpr : GaussianProcessRegressor
            Modèle de Kriging entraîné.
        radii : np.ndarray
            Tableau de taille (n_samples, 3) contenant les rayons.

        Retourne :
        ----------
        gradients : np.ndarray
            Gradient du volume par rapport aux rayons (shape: (n_samples, 3)).
        """
        epsilon = 1e-3
        radii = np.zeros(3)
        if geomScheme is not None:
            true_indices = [i for i, val in enumerate(geomScheme) if val]
            for idx, val in zip(true_indices, self.radius):
                radii[idx] = val
        else:
            for idx, rad in enumerate(self.radius):
                radii[idx] = rad
        radii = np.array(radii).reshape(-1, 3)
        grad = np.zeros(3)

        for idx, rad in enumerate(grad):
            perturbed_radii = radii.copy()
            perturbed_radii[0][idx] += epsilon
            grad[idx] = (gpr.predict(perturbed_radii) - gpr.predict(radii)) / epsilon
        return grad

    def getNumberNodesAtBoundary(self):
        """
        Get the number of nodes at the boundary

        Returns:
        --------
        int
            Number of nodes at the boundary
        """
        counterNodes = 0
        nodeAlreadyCounted = []
        for beam in self.beams:
            for point in [beam.point1, beam.point2]:
                if point.indexBoundary is not None and point.indexBoundary not in nodeAlreadyCounted:
                    counterNodes += 1
                    nodeAlreadyCounted.append(point.indexBoundary)
        return counterNodes

    def getCellBoundaryBox(self):
        """
        Get the boundary box of the cell

        Returns:
        --------
        list
            List of the boundary box coordinates
        """
        xMin = self.coordinateCell[0]
        xMax = self.coordinateCell[0] + self.cellSize[0]
        yMin = self.coordinateCell[1]
        yMax = self.coordinateCell[1] + self.cellSize[1]
        zMin = self.coordinateCell[2]
        zMax = self.coordinateCell[2] + self.cellSize[2]
        return [xMin, xMax, yMin, yMax, zMin, zMax]

    def getRGBcolorDependingOfRadius(self):
        return tuple(r / 0.1 for r in self.radius)

    def getCellCornerCoordinates(self) -> list:
        """
        Get the corner coordinates of the cell.

        Returns:
        --------
        list of tuples
            List of (x, y, z) coordinates of the corner points.
        """
        x0, y0, z0 = self.coordinateCell  # Position de départ de la cellule
        dx, dy, dz = self.cellSize  # Dimensions de la cellule

        # Liste des 8 coins de la cellule
        corners = [
            (x0, y0, z0),
            (x0 + dx, y0, z0),
            (x0, y0 + dy, z0),
            (x0 + dx, y0 + dy, z0),
            (x0, y0, z0 + dz),
            (x0 + dx, y0, z0 + dz),
            (x0, y0 + dy, z0 + dz),
            (x0 + dx, y0 + dy, z0 + dz),
        ]

        return corners

    def addCellNeighbour(self, neighbourCell: "Cell") -> None:
        """
        Add a neighbour cell to the current cell

        Parameters:
        -----------
        neighbourCell: Cell
            Neighbour cell to add
        """
        self.neighbourCells.append(neighbourCell)

    def getNeighbourCells(self) -> list:
        """
        Get the list of neighbour cells

        Returns:
        --------
        list
            List of neighbour cells
        """
        return self.neighbourCells

    def printCellData(self):
        print("Cell position: ", self.posCell)
        print("Cell coordinates: ", self.coordinateCell)
        print("Cell size: ", self.cellSize)
        print("Lattice type: ", self.latticeType)
        print("Beam radius: ", self.radius)
        print("Beam material: ", self._beamMaterial)
        print("Beams in cell: ", self.beams)
        print("Cell center point: ", self.centerPoint)
        print("Cell index: ", self.index)
        print("Beam material: ", self._beamMaterial)
        print("Coupling matrix: ", self.matB)
        print("Number of beams: ", len(self.beams))
        print("Volume of the cell: ", self.getVolumeCell())
        print("Relative density: ", self.getRelativeDensityCell())
        print("Number of nodes at boundary: ", self.getNumberNodesAtBoundary())

    def getTranslationRigidBody(self):
        """
        Get the translation of the rigid body
        """
        translation = np.zeros(3)
        for beam in self.beams:
            for point in [beam.point1, beam.point2]:
                if point.indexBoundary is not None:
                    translation += point.getDisplacementValue()[:3]
        return translation / self.getNumberNodesAtBoundary()

    def getRotationRigidBody(self):
        """
        Get the rotation matrix of the rigid body using SVD.
        """
        all_points = self.getAllPoints()
        initial_positions = np.array([point.getPos() for point in all_points])  # P_i
        final_positions = np.array([point.getDeformedPos() for point in all_points])  # P_i'

        # Soustraction du centre de gravité
        center_initial = np.mean(initial_positions, axis=0)
        center_final = np.mean(final_positions, axis=0)
        P = initial_positions - center_initial
        P_prime = final_positions - center_final

        # Matrice de covariance
        H = P.T @ P_prime

        # Décomposition SVD
        U, _, Vt = np.linalg.svd(H)
        R = Vt.T @ U.T

        # Correction si nécessaire (assurer que R est une rotation propre)
        if np.linalg.det(R) < 0:
            Vt[-1, :] *= -1
            R = Vt.T @ U.T

        return R