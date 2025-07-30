"""
cell.py
"""
import json

import numpy as np
from colorama import Fore, Style

from beam import *
from geometries.geometries_utils import *

# from scipy.sparse import coo_matrix


class Cell(object):
    """
    Define Cell data for lattice structure
    """

    def __init__(self, pos_cell: list, initial_cell_size: list, coordinate_cell: list, geom_types: list[str],
                 radii: list[float], grad_radius: list, grad_dim: list, grad_mat: list, uncertainty_node: float = 0.0,
                 verbose: int = 0):
        """
        Initialize a Cell with its dimensions and position

        Parameters:
        -----------
        pos_cell: list
            Position of the cell in the lattice
        initial_cell_size: list
            Initial size of the cell
        coordinate_cell: list
            Coordinates of the cell minimum corner in the lattice
        geom_types: list[str]
            Type of lattice geometry
        radii: float
            Base radii of the beam
        grad_radius: list
            Gradient of the radii
        grad_dim: list
            Gradient of the dimensions
        grad_mat: list
            Gradient of the material
        uncertainty_node: float
            Standard deviation for adding uncertainty to node coordinates. Defaults to 0.0.
        _verbose: bool
            If True, prints additional information during initialization. Defaults to False.
        """
        self.original_cell_geom = None
        self.original_tags = None
        self.center_point = None
        self._beamMaterial = None
        self.cell_size = None
        self.pos_cell: list[int] = pos_cell
        self.coordinate_cell: list[float] = coordinate_cell
        self.beams: Optional[list] = []
        self.index: Optional[int] = None
        self.geom_types: list[str] = geom_types
        self.radii: list[float] = radii
        self.coupling_matrix_B: Optional = None  # B matrix (Coupling matrix)
        self.uncertainty_node: float = uncertainty_node
        self.grad_radius: list = grad_radius
        self.grad_mat: list = grad_mat
        self.grad_dim: list = grad_dim
        self._verbose: int = verbose
        self.neighbour_cells: Optional = []

        self.define_original_tags()
        self.generate_cell_properties(initial_cell_size)

    def __repr__(self) -> str:
        return f"Cell(Coordinates:{self.coordinate_cell}, Size: {self.cell_size}, Index:{self.index})"

    @property
    def volume(self):
        """ Calculate the volume of the cell."""
        return self.cell_size[0] * self.cell_size[1] * self.cell_size[2]

    @property
    def relative_density(self) -> float:
        """
        Calculate the relative density of the cell based on the volume of beams and the cell volume.
        """
        volumeBeams = 0
        for beam in self.beams:
            volumeBeams += beam.volume
        return volumeBeams / self.volume

    @property
    def volume_each_geom(self) -> np.ndarray:
        """
        Get the volume of the cell separated by geometry type_beam.
        """
        volumes = np.zeros(len(self.radii))
        for beam in self.beams:
            if not beam.beam_mod:
                volumeBeam = beam.volume
                volumes[beam.type_beam] += volumeBeam
        return volumes

    @property
    def boundary_box(self) -> list:
        """
        Get the boundary box of the cell

        Returns:
        --------
        list
            List of the boundary box coordinates
        """
        xMin = self.coordinate_cell[0]
        xMax = self.coordinate_cell[0] + self.cell_size[0]
        yMin = self.coordinate_cell[1]
        yMax = self.coordinate_cell[1] + self.cell_size[1]
        zMin = self.coordinate_cell[2]
        zMax = self.coordinate_cell[2] + self.cell_size[2]
        return [xMin, xMax, yMin, yMax, zMin, zMax]

    @property
    def corner_coordinates(self) -> list:
        """
        Get the corner coordinates of the cell.

        Returns:
        --------
        list of tuples
            List of (x, y, z) coordinates of the corner points.
        """
        x0, y0, z0 = self.coordinate_cell
        dx, dy, dz = self.cell_size

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

    def generate_cell_properties(self, initialCellSize):
        """
        Generate a cell object with beams and nodes based on the lattice type_beam and radii.

        Parameters:
        -----------
        initialCellSize: list
            Initial size of the cell without modification
        """
        idxCell = 0
        for idx, radius in enumerate(self.radii):
            if radius > 0.0:
                if idxCell == 0:
                    self.get_beam_material()
                    beamRadius = self.get_radius(radius)
                    self.get_cell_size(initialCellSize)
                    self.generate_beams(self.geom_types[idx], beamRadius, idx)
                    self.get_cell_center()
                else:
                    hybridRadius = self.get_radius(radius)
                    self.generate_beams(self.geom_types[idx], hybridRadius, idx)
                idxCell += 1

    def define_original_tags(self) -> None:
        """
        Define original tags and cell geometry based on the lattice type_beam and radii.

        Parameters:
        -----------
        geom_types: list[int]
            List of lattice types
        radii: list[float]
            List of beam radii
        """
        if len(self.radii) == 1:
            if self.geom_types[0] == 0:
                self.original_tags = [1000, 1001, 1002, 1003, 1004, 1005, 1006, 1007]
                self.original_cell_geom = [0, 0, 0, 0, 0, 0, 0, 0]
            elif self.geom_types[0] == 16:
                self.original_tags = [100, 101, 102, 103, 104, 105, 106, 107, 108, 109, 110, 111]
                self.original_cell_geom = [1, 1, 1, 1, 1, 1, 1, 1, 1, 1, 1, 1]
            elif self.geom_types[0] == 19:
                self.original_tags = [10, 11, 12, 13, 14, 15]
                self.original_cell_geom = [2, 2, 2, 2, 2, 2]
            else:
                pass
        elif len(self.radii) == 2:
            if self.geom_types[0] == 0 and self.geom_types[1] == 16:
                self.original_tags = [1000, 1001, 1002, 1003, 1004, 1005, 1006, 1007,
                                      100, 101, 102, 103, 104, 105, 106, 107, 108, 109, 110, 111]
                self.original_cell_geom = [0, 0, 0, 0, 0, 0, 0, 0,
                                           1, 1, 1, 1, 1, 1, 1, 1, 1, 1, 1, 1]
            else:
                raise ValueError("Lattice type_beam not recognized")
        else:
            self.original_tags = [1000, 1001, 1002, 1003, 1004, 1005, 1006, 1007,
                                  10, 11, 12, 13, 14, 15,
                                  100, 101, 102, 103, 104, 105, 106, 107, 108, 109, 110, 111]
            self.original_cell_geom = [0, 0, 0, 0, 0, 0, 0, 0,
                                       2, 2, 2, 2, 2, 2,
                                       1, 1, 1, 1, 1, 1, 1, 1, 1, 1, 1, 1]

    def generate_beams(self, latticeType: int, beamRadius: float, beamType: int = 0) -> None:
        """
        Generate beams and nodes using a given lattice type_beam and parameters.

        Parameters:
        -----------
        geom_types: int
            Type of lattice geometry
        startCellPos: list
            Position of the start of the cell
        beamType: int
            Type of beam
        """
        pointDict = {}
        for line in get_beam_structure(latticeType):
            x1, y1, z1, x2, y2, z2 = map(float, line)
            if (x1, y1, z1) in pointDict:
                point1 = pointDict[(x1, y1, z1)]
            else:
                point1 = Point(x1 * self.cell_size[0] + self.coordinate_cell[0], y1 * self.cell_size[1] + self.coordinate_cell[1],
                               z1 * self.cell_size[2] + self.coordinate_cell[2], self.uncertainty_node)
                pointDict[(x1, y1, z1)] = point1
            if (x2, y2, z2) in pointDict:
                point2 = pointDict[(x2, y2, z2)]
            else:
                point2 = Point(x2 * self.cell_size[0] + self.coordinate_cell[0], y2 * self.cell_size[1] + self.coordinate_cell[1],
                               z2 * self.cell_size[2] + self.coordinate_cell[2], self.uncertainty_node)
                pointDict[(x2, y2, z2)] = point2
            beam = Beam(point1, point2, beamRadius, self._beamMaterial, beamType)
            self.beams.append(beam)


    def get_beam_material(self) -> None:
        """
        Get the material of the beam based on the gradient and position.

        Parameters:
        -----------
        grad_mat: list
            Gradient of the material

        Returns:
        ---------
        materialType: int
            Material index of the beam
        """
        self._beamMaterial = self.grad_mat[self.pos_cell[2]][self.pos_cell[1]][self.pos_cell[0]]

    def get_radius(self, base_radius: float) -> float:
        """
        Calculate and return the beam radii

        Parameters:
        -----------
        grad_radius: list
            Gradient of the radii
        BaseRadius: float
            Base radii of the beam

        Returns:
        ---------
        actualBeamRadius: float
            Calculated beam radii
        """
        beamRadius = (base_radius * self.grad_radius[self.pos_cell[0]][0] *
                      self.grad_radius[self.pos_cell[1]][1] *
                      self.grad_radius[self.pos_cell[2]][2])
        return beamRadius

    def get_cell_size(self, initial_cell_size: list) -> None:
        """
        Calculate and return the cell size

        Parameters:
        -----------
        initialCellSize: 3-array
            Dimension of the initial cell without modification
        grad_dim:

        Returns:
        ---------
        cell_size : float
            Calculated beam radii
        """
        self.cell_size = [initial_size * self.grad_dim[pos][i] for i, (initial_size, pos) in
                          enumerate(zip(initial_cell_size, self.pos_cell))]

    def get_cell_center(self) -> None:
        """
        Calculate the center point of the cell
        """
        self.center_point = [self.coordinate_cell[i] + self.cell_size[i] / 2 for i in range(3)]

    def get_list_points(self) -> list:
        """
        Determine a list of points in cell
        """
        pointList = []
        for beam in self.beams:
            for point in [beam.point1, beam.point2]:
                if point not in pointList:
                    pointList.append(point)
        return pointList

    def remove_beam(self, beam_to_delete: "Beam") -> None:
        """
        Removes a beam from the lattice

        Parameters:
        ------------
        beamToDelete: beam Object
            Beam to remove
        """
        try:
            self.beams.remove(beam_to_delete)
        except ValueError:
            if self._verbose > 0:
                print("Beam not found in the list for delete:", beam_to_delete)

    def add_beam(self, beam_to_add: "Beam") -> None:
        """
        Adding beam to cell
        """
        if isinstance(beam_to_add, Beam):
            self.beams.append(beam_to_add)
        elif isinstance(beam_to_add, tuple):
            for beam in beam_to_add:
                self.beams.append(beam)
        else:
            raise ValueError("Invalid beam type_beam")

    def get_point_on_surface(self, surfaceName: str) -> list:
        """
        Get the points on the surface specified in the global reference frame.

        Parameters:
        -----------
        surfaceName: str
            Name of the surface. Choose from 'Xmin', 'Xmax', 'Ymin', 'Ymax', 'Zmin', 'Zmax', 'Xmid', 'Ymid', 'Zmid'.
            If 'Xmid', 'Ymid', or 'Zmid' is specified, it returns the points at the bottom of the cell

        Returns:
        --------
        list
           List of points on the specified surface.
        """
        boundaryBox = self.boundary_box
        surface_map = {
            "Xmin": boundaryBox[0],
            "Xmax": boundaryBox[1],
            "Ymin": boundaryBox[2],
            "Ymax": boundaryBox[3],
            "Zmin": boundaryBox[4],
            "Zmax": boundaryBox[5],
            "Xmid": self.coordinate_cell[0],
            "Ymid": self.coordinate_cell[1],
            "Zmid": self.coordinate_cell[2]
        }

        coord_index = {
            "Xmin": "x", "Xmax": "x", "Xmid": "x",
            "Ymin": "y", "Ymax": "y", "Ymid": "y",
            "Zmin": "z", "Zmax": "z", "Zmid": "z"
        }

        if surfaceName not in surface_map:
            raise ValueError(
                f"Surface '{surfaceName}' is not valid. Choose from 'Xmin', 'Xmax', 'Ymin', 'Ymax', 'Zmin', 'Zmax', "
                f"'Xmid', 'Ymid', 'Zmid'.")

        surface_value = surface_map[surfaceName]
        axis = coord_index[surfaceName]

        return list({
            point for beam in self.beams
            for point in [beam.point1, beam.point2]
            if getattr(point, axis) == surface_value
        })

    def get_node_order_to_simulate(self) -> dict:
        """
        Get the order of nodes to simulate in the cell
        """
        tag_dict = {tag: None for tag in self.original_tags}

        for beam in self.beams:
            if beam.radii > 0:
                for point in [beam.point1, beam.point2]:
                    if point.index_boundary is not None:
                        tag = point.local_tag
                        if tag:  # Ensure tags is not an empty list
                            if tag[0] in self.original_tags:
                                tag_dict[tag[0]] = point
        return tag_dict

    def getNodesOrderNN(self, nodeInOrder: dict, numberRadiusNN) -> dict:
        """
        Get the nodes order for the neural network

        Parameters:
        -----------
        nodeInOrder: dict
            Dictionary of nodes in order
        original_cell_geom: list
            Original cell geometry
        """
        # TODO : Check utility
        tag_dictNN = nodeInOrder.copy()
        idx = 0
        for i, key in nodeInOrder.items():
            if key:
                tag_dictNN[i] = 1
            # elif self.original_cell_geom[idx] < numberRadiusNN:
            #     tag_dictNN[i] = 1
            else:
                tag_dictNN[i] = 0
            idx += 1
        return tag_dictNN

    def set_displacement_at_boundary_nodes(self, displacementArray: list, displacementIndex: list, printLevel=0) -> \
            None:
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
                if point.index_boundary is not None and point.index_boundary in displacementIndex:
                    index = displacementIndex.index(point.index_boundary)
                    indexActual = 0
                    for i in range(6):
                        if point.fixed_DOF[i] == 0:  # Filter out the fixed DOF
                            point.displacement_vector[i] = displacementArray[index + indexActual]
                            indexActual += 1

    def get_displacement_at_nodes(self, nodeList: dict, numberRadiusNN: int) -> list:
        """
        Get the displacement at nodes.

        Parameters:
        -----------
        nodeList: list of Point objects
            List of nodes to get the displacement.
        numberRadiusNN: int
            Number of radii for the neural network

        Returns:
        --------
        list
            A flattened list of displacement values.
        """
        nodeListNN = self.getNodesOrderNN(nodeList, numberRadiusNN)

        displacementList = []
        nullDisplacement = [0.0, 0.0, 0.0, 0.0, 0.0, 0.0]
        for key, node in nodeList.items():
            if node:
                displacement = node.displacement_vector
                displacementList.append(displacement)
            elif nodeListNN[key] == 1:
                displacementList.append(nullDisplacement)
        return displacementList

    def set_reaction_force_on_nodes(self, nodeList: dict, reactionForce: list) -> None:
        """
        Set reaction force on each node.

        Parameters:
        -----------
        nodeList: dict
            Dictionary mapping node tags to Point objects.
        reactionForce: list
            List of reaction force values corresponding to the nodes.
        """
        #TODO : Verify if used

        # Vérification que la longueur des listes correspond
        if len([v for v in nodeList.values() if v is not None]) != len(reactionForce):
            print("Lenght nodeList ", len([v for v in nodeList.values() if v is not None])
                  , "Lenght reactionForce ", len(reactionForce))
            print("nodeList", nodeList)
            print("reactionForce", reactionForce)
            raise ValueError("Mismatch: nodeList and reactionForce must have the same length.")

        # for beam in self.beams:
        #     for point in [beam.point1, beam.point2]:
        #         if point.index_boundary is not None:
        #             index = list(nodeList.keys()).index(point.local_tag[0])
        #             point.setReactionForce(reactionForce[index])
        idx = 0
        for node in nodeList:
            if nodeList[node]:
                nodeList[node].set_reaction_force(reactionForce[idx])
                idx += 1

    def get_number_boundary_nodes(self) -> int:
        """
        Get the number of boundary nodes in the cell
        """
        return len(
            [beam for beam in self.beams for point in [beam.point1, beam.point2] if point.index_boundary is not None])

    def build_coupling_operator(self, nb_free_DOF: int) -> None:
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
                if point.index_boundary is not None and point.index_boundary not in listBndNodes:
                    localNodeIndex = self.original_tags.index(point.local_tag[0])
                    listBndNodes.append(point.index_boundary)
                    for i in range(6):
                        if point.global_free_DOF_index[i] is not None:
                            data.append(1)
                            col.append(localNodeIndex * 6 + i)
                            row.append(point.global_free_DOF_index[i])
        nbBndDOFloc = len(listBndNodes) * 6
        shapeB = (nb_free_DOF, nbBndDOFloc)
        self.coupling_matrix_B = coo_matrix((data, (row, col)), shape=shapeB)

    def build_preconditioner(self, SchurMatrix: "coo_matrix") -> "coo_matrix":
        """
        Build the preconditioner part for the cell

        Parameters:
        -----------
        SchurMatrix: coo_matrix
            Schur matrix
        """
        if self.coupling_matrix_B is None:
            raise ValueError("Coupling matrix has not been built yet. Please build it first.")
        if self.coupling_matrix_B.shape[1] != SchurMatrix.shape[0]:
            print("Shape of B matrix", self.coupling_matrix_B.shape)
            print("Shape of Schur matrix", SchurMatrix.shape)
            raise ValueError("Incompatible dimensions between the coupling matrix and the Schur matrix.")

        return self.coupling_matrix_B @ SchurMatrix @ self.coupling_matrix_B.transpose()

    def get_internal_energy(self) -> float:
        """
        Get the internal energy of the cell
        """
        internalEnergy = 0
        for beam in self.beams:
            for point in [beam.point1, beam.point2]:
                if point.index_boundary is not None:
                    pointEnergy = point.calculate_point_energy()
                    if pointEnergy < 0:
                        raise ValueError("Negative energy detected at point with index " + str(point.index_boundary))
                    internalEnergy += pointEnergy
        return internalEnergy

    def get_displacement_data(self) -> list:
        """
        Build and return displacement data on cell for dataset generation
        """
        allBoundaryDisplacementData = []
        for beam in self.beams:
            for point in [beam.point1, beam.point2]:
                if point.index_boundary is not None:
                    allBoundaryDisplacementData.append(point.displacement_vector)
        return allBoundaryDisplacementData

    def change_beam_radius(self, newRadius: list) -> None:
        """
        ATTENTION: BEAM MOD IS NOT WORKING
        Change beam radii in the cell

        Parameters:
        -----------
        newRadius: list
            beam radii wanted to assign
        hybridData: list
            Hybrid data type_beam
        """
        print(Fore.RED + "WARNING: Beam modification is not implemented yet. " + Style.RESET_ALL)
        assert len(newRadius) == len(self.radii), ("Length of new radii vector and already cell radii vector needs "
                                                    "to be equal ")
        beamRadius = []
        for rad in newRadius:
            beamRadius.append(self.get_radius(rad))

        for beam in self.beams:
            if beam.beam_mod:
                beam.radii = beamRadius[beam.type_beam] * beam.penalization_coefficient
            else:
                beam.radii = beamRadius[beam.type_beam]

        self.radii = newRadius

    def get_relative_density_kriging(self, krigingModel, geomScheme=None) -> float:
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
            for idx, val in zip(true_indices, self.radii):
                radii[idx] = val
        else:
            for idx, rad in enumerate(self.radii):
                radii[idx] = rad
        radii = np.array(radii).reshape(-1, 3)
        relativeDensity = krigingModel.predict(radii)[0]
        return relativeDensity

    def get_relative_density_gradient(self, relativeDensityPolyDeriv) -> float:
        """
        Get the gradient of the relative density

        Parameters:
        -----------
        relative_density_poly_deriv: list
            List of polynomial derivative functions

        Returns:
        --------
        deriv: float
            Derivative of the relative density
        """
        deriv = 0
        for idx, polyDeriv in enumerate(relativeDensityPolyDeriv):
            deriv += polyDeriv(self.radii[idx])
        return deriv

    def get_relative_density_gradient_kriging(self, gpr, geomScheme=None) -> np.ndarray:
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
            for idx, val in zip(true_indices, self.radii):
                radii[idx] = val
        else:
            for idx, rad in enumerate(self.radii):
                radii[idx] = rad
        radii = np.array(radii).reshape(-1, 3)
        grad = np.zeros(3)

        for idx, rad in enumerate(grad):
            perturbed_radii = radii.copy()
            perturbed_radii[0][idx] += epsilon
            grad[idx] = (gpr.predict(perturbed_radii) - gpr.predict(radii)) / epsilon
        return grad

    def get_number_nodes_at_boundary(self):
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
                if point.index_boundary is not None and point.index_boundary not in nodeAlreadyCounted:
                    counterNodes += 1
                    nodeAlreadyCounted.append(point.index_boundary)
        return counterNodes

    def get_RGBcolor_depending_of_radius(self):
        """
        Get the RGB color of the cell depending on the radii.
        """
        return tuple(r / 0.1 for r in self.radii)

    def add_cell_neighbour(self, neighbour_cell: "Cell") -> None:
        """
        Add a neighbour cell to the current cell if it's not already present.

        Parameters:
        -----------
        neighbourCell: Cell
            Neighbour cell to add
        """
        if neighbour_cell not in self.neighbour_cells:
            self.neighbour_cells.append(neighbour_cell)

    def print_data(self):
        """
        Print the data of the cell for debugging purposes.
        """
        print("Cell position: ", self.pos_cell)
        print("Cell coordinates: ", self.coordinate_cell)
        print("Cell size: ", self.cell_size)
        print("Lattice type_beam: ", self.geom_types)
        print("Beam radii: ", self.radii)
        print("Beam material: ", self._beamMaterial)
        print("Beams in cell: ", self.beams)
        print("Cell center point: ", self.center_point)
        print("Cell index: ", self.index)
        print("Beam material: ", self._beamMaterial)
        print("Coupling matrix: ", self.coupling_matrix_B)
        print("Number of beams: ", len(self.beams))
        print("Volume of the cell: ", self.volume)
        print("Relative density: ", self.relative_density)
        print("Number of nodes at boundary: ", self.get_number_nodes_at_boundary())

    def get_translation_rigid_body(self):
        """
        Get the translation of the rigid body
        """
        translation = np.zeros(3)
        for beam in self.beams:
            for point in [beam.point1, beam.point2]:
                if point.index_boundary is not None:
                    translation += point.displacement_vector[:3]
        return translation / self.get_number_nodes_at_boundary()

    def get_rotation_rigid_body(self):
        """
        Get the rotation matrix of the rigid body using SVD.
        """
        all_points = self.get_list_points()
        initial_positions = np.array([point.coordinates for point in all_points])  # P_i
        final_positions = np.array([point.deformed_coordinates for point in all_points])  # P_i'

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
