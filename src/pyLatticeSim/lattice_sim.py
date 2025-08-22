from typing import TYPE_CHECKING
import numpy as np
from scipy.sparse import coo_matrix

from pyLattice.lattice import Lattice
from pyLattice.utils import open_lattice_parameters
from pyLatticeSim.utils_simulation import get_schur_complement


if TYPE_CHECKING:
    from mesh_file.mesh_trimmer import MeshTrimmer

from pyLattice.timing import *
timing = Timing()


class LatticeSim(Lattice):
    def __init__(self, name_file: str, mesh_trimmer: "MeshTrimmer" = None, verbose: int = 0):
        super().__init__(name_file, mesh_trimmer, verbose)
        self._simulation_flag = True

        self.enable_periodicity = None  # Warning not working for graded structures
        self.boundary_conditions = None
        self.enable_simulation_properties = None
        self.material_name = None

        self.define_simulation_parameters(name_file)
        assert self.material_name is not None, "Material name_lattice must be defined for simulation properties."

        self.free_DOF = None  # Free DOF gradient conjugate gradient method
        self.max_index_boundary = None
        self.global_displacement_index = None
        self.n_DOF_per_node: int = 6  # Number of DOF per node (3 translation + 3 rotation)
        self.penalization_coefficient: float = 1.5  # Fixed with previous optimization

        self.get_all_angles()
        self.set_beam_node_mod()
        # Define global indexation
        self.define_node_index_boundary()
        self.set_boundary_conditions()

    def define_simulation_parameters(self, name_file: str):
        """
        Define simulation parameters from the input file.

        Parameters
        ----------
        name_file : str
            Name of the input file
        """
        lattice_parameters = open_lattice_parameters(name_file)

        sim_params = lattice_parameters.get("simulation_parameters", {})
        self.enable_simulation_properties = bool(sim_params.get("enable", False))
        self.material_name = sim_params.get("material", "VeroClear")
        self.enable_periodicity = sim_params.get("periodicity", False)

        self.boundary_conditions = lattice_parameters.get("boundary_conditions", {})


    def apply_displacement_surface(self, surfaceNames: list[str], valueDisplacement: list[float],
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
        self.apply_constraints_nodes(surfaceNames, valueDisplacement, DOF, "Displacement")

    def apply_constraints_nodes(self, surfaces: list[str], value: list[float], DOF: list[int],
                                type_constraint: str = "Displacement", surface_cells: list[str] = None) -> None:
        """
        Apply boundary conditions to the lattice

        Parameters:
        -----------
        surfaces: list[str]
            List of surfaces to apply constraint (e.g., ["Xmin", "Xmax", "Ymin"])
        value: list of float
            Values to apply to the constraint
        DOF: list of int
            Degree of freedom to apply constraint (0: x, 1: y, 2: z, 3: Rx, 4: Ry, 5: Rz)
        type_beam: str
            Type of constraint (Displacement, Force)
        surface_cells: list[str], optional
            List of surfaces to find points on cells (e.g., ["Xmin", "Xmax", "Ymin"]). If None, uses surfaceNames.
        """
        pointSet = self.find_point_on_lattice_surface(surfaces, surface_cells)

        indexBoundaryList = {point.index_boundary for point in pointSet}

        for cell in self.cells:
            for beam in cell.beams:
                for node in [beam.point1, beam.point2]:
                    if node.index_boundary in indexBoundaryList:
                        for val, DOFi in zip(value, DOF):
                            if type_constraint == "Displacement":
                                node.displacement_vector[DOFi] = val
                                node.fix_DOF([DOFi])
                            elif type_constraint == "Force":
                                node.applied_force[DOFi] = val
                            else:
                                raise ValueError("Invalid type of constraint. Use 'Displacement' or 'Force'.")

    def apply_force_surface(self, surfaceName: list[str], valueForce: list[float], DOF: list[int]) -> None:
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
        self.apply_constraints_nodes(surfaceName, valueForce, DOF, "Force")

    def fix_DOF_on_surface(self, surfaceName: list[str], dofFixed: list[int]) -> None:
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
        self.apply_constraints_nodes(surfaceName, [0.0 for _ in dofFixed], dofFixed, "Displacement")

    def fix_DOF_on_node(self, nodeList: list[int], dofFixed: list[int]) -> None:
        """
        Fix degree of freedom on the surface of the lattice

        Parameters:
        -----------
        nodeList: list of int
            List of node index to apply boundary conditions
        dofFixed: list of int
            List of degree of freedom to fix (0: x, 1: y, 2: z, 3: Rx, 4: Ry, 5: Rz)
        """
        if self.get_number_nodes() < max(nodeList):
            raise ValueError("Invalid node index, node do not exist.")

        for node in nodeList:
            if node < 0 or node >= self.get_number_nodes():
                raise ValueError("Node index out of range.")

        for cell in self.cells:
            for beam in cell.beams:
                for node in [beam.point1, beam.point2]:
                    if node.index in nodeList:
                        node.fix_DOF(dofFixed)

    def get_global_displacement(self, withFixed: bool = False, OnlyImposed: bool = False) \
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
            Dictionary of global displacement with index_boundary as key and displacement vector as value
        global_displacement_index: list of int
            List of index_boundary of the lattice
        """
        globalDisplacement = []
        globalDisplacementIndex = []
        processed_nodes = set()
        for cell in self.cells:
            for beam in cell.beams:
                for node in [beam.point1, beam.point2]:
                    if node.index_boundary is not None and node.index_boundary not in processed_nodes:
                        for i in range(6):
                            if node.fixed_DOF[i] == 0 and not OnlyImposed:
                                globalDisplacement.append(node.displacement_vector[i])
                                globalDisplacementIndex.append(node.index_boundary)
                            elif node.fixed_DOF[i] == 0 and node.applied_force[i] == 0:
                                globalDisplacement.append(0)
                            elif withFixed or OnlyImposed:
                                globalDisplacement.append(node.displacement_vector[i])
                                globalDisplacementIndex.append(node.index_boundary)

                        processed_nodes.add(node.index_boundary)
        if not OnlyImposed:
            self.global_displacement_index = globalDisplacementIndex
        if self._verbose > 2:
            print("globalDisplacement: ", globalDisplacement)
            print("global_displacement_index: ", globalDisplacementIndex)
        return globalDisplacement, globalDisplacementIndex

    @timing.timeit
    def define_node_index_boundary(self) -> None:
        """
        Define boundary tag for all boundary nodes and calculate the total number of boundary nodes
        """
        IndexCounter = 0
        nodeAlreadyIndexed = {}
        for cell in self.cells:
            for beam in cell.beams:
                for node in [beam.point1, beam.point2]:
                    localTag = node.tag_point(cell.boundary_box)
                    node.set_local_tag(localTag)
                    if localTag:
                        if node in nodeAlreadyIndexed:
                            node.index_boundary = nodeAlreadyIndexed[node]
                        else:
                            nodeAlreadyIndexed[node] = IndexCounter
                            node.index_boundary = IndexCounter
                            IndexCounter += 1
        self.max_index_boundary = IndexCounter - 1


    def get_global_reaction_force(self, appliedForceAdded: bool = False) -> dict:
        """
        Get local reaction force of the lattice and sum if identical TagIndex

        Returns:
        --------
        globalReactionForce: dict
            Dictionary of global reaction force with index_boundary as key and reaction force vector as value
        """
        globalReactionForce = {i: [0, 0, 0, 0, 0, 0] for i in range(self.max_index_boundary + 1)}
        for cell in self.cells:
            nodeIndexProcessed = set()
            for beam in cell.beams:
                for node in [beam.point1, beam.point2]:
                    if node.index_boundary is not None and node.index not in nodeIndexProcessed:
                        globalReactionForce[node.index_boundary] = [
                            x + y for x, y in zip(globalReactionForce[node.index_boundary], node.reaction_force_vector)
                        ]
                        if appliedForceAdded:
                            for i in range(6):
                                if node.applied_force[i] != 0:
                                    globalReactionForce[node.index_boundary][i] = node.applied_force[i]
                        nodeIndexProcessed.add(node.index)
        return globalReactionForce

    def get_global_reaction_force_without_fixed_DOF(self, globalReactionForce: dict, rightHandSide: bool = False) \
            -> np.ndarray:
        """
        Get global reaction force of free degree of freedom

        Parameters:
        -----------
        globalReactionForce: dict
            Dictionary of global reaction force with index_boundary as key and reaction force vector as value

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
                    if node.index_boundary is not None and node.index_boundary not in processed_nodes:
                        # Append reaction force components where fixed_DOF is 0
                        RFToAdd = []
                        for i in range(6):
                            if node.applied_force[i] != 0 and rightHandSide:
                                RFToAdd.append(-node.applied_force[i])
                                # Add a sign minus because right-hand side already with a sign minus see (b = -b)
                            elif node.fixed_DOF[i] == 0:
                                RFToAdd.append(globalReactionForce[node.index_boundary][i])
                        globalReactionForceWithoutFixedDOF.append(RFToAdd)
                        # globalReactionForceWithoutFixedDOF.append([
                        #     v1 for v1, v2 in zip(globalReactionForce[node.index_boundary], node.fixed_DOF)
                        #     if v2 == 0])
                        # print(globalReactionForceWithoutFixedDOF[-1])
                        # Mark this node as processed
                        processed_nodes.add(node.index_boundary)
        return np.concatenate(globalReactionForceWithoutFixedDOF)

    def define_free_DOF(self):
        """
        Get total number of degrees of freedom in the lattice
        """
        self.free_DOF = 0
        processed_nodes = set()
        for cell in self.cells:
            for beam in cell.beams:
                for node in [beam.point1, beam.point2]:
                    if node.index_boundary is not None and node.index_boundary not in processed_nodes:
                        self.free_DOF += node.fixed_DOF.count(0)
                        processed_nodes.add(node.index_boundary)

    def set_global_free_DOF_index(self) -> None:
        """
        Set global free degree of freedom index for all nodes in boundary
        """
        counter = 0
        processed_nodes = {}
        for cell in self.cells:
            for beam in cell.beams:
                for node in [beam.point1, beam.point2]:
                    if node.index_boundary is not None:
                        if node.index_boundary not in processed_nodes.keys():
                            for i in np.where(np.array(node.fixed_DOF) == 0)[0]:
                                node.global_free_DOF_index[i] = counter
                                counter += 1
                            processed_nodes[node.index_boundary] = node.global_free_DOF_index
                        else:
                            node.global_free_DOF_index[:] = processed_nodes[node.index_boundary]

    def initialize_reaction_force(self) -> None:
        """
        Initialize reaction force of all nodes to 0 on each DOF
        """
        for cell in self.cells:
            for beam in cell.beams:
                for node in [beam.point1, beam.point2]:
                    node.initialize_reaction_force()

    def initialize_displacement(self) -> None:
        """
        Initialize displacement of all nodes to zero on each DOF
        """
        for cell in self.cells:
            for beam in cell.beams:
                for node in [beam.point1, beam.point2]:
                    node.initialize_displacement()

    def build_coupling_operator_cells(self) -> None:
        """
        Build coupling operator for each cell in the lattice
        """
        for cell in self.cells:
            cell.build_coupling_operator(self.free_DOF)


    def get_radius_temp(self) -> float:
        """
        ################### TEMPORARY FUNCTION ###################
        Get the radii of the lattice

        Returns:
        --------
        radii: float
            radii of the lattice
        """
        for cell in self.cells:
            for beam in cell.beams:
                return beam.radius

    def apply_reaction_force_on_node_list(self, reactionForce: list, nodeCoordinatesList: list):
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
                        node.set_reaction_force(reactionForce[index])

    def apply_displacement_on_node_list(self, displacement: list, nodeCoordinatesList: list):
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
                        node.displacement_vector = displacement[index]
                        node.fix_DOF([i for i in range(6)])

    def set_boundary_conditions(self) -> None:
        """
        Set boundary conditions on the lattice.
        """
        def check_data_boundary_condition_validity(data_dict_valid: dict) -> None:
            """
            Check if the data of the boundary condition is valid
            """
            if "Surface" not in data_dict_valid or "Value" not in data_dict_valid or "DOF" not in data_dict_valid:
                raise ValueError("Invalid boundary condition data. 'Surface', 'Value' and 'DOF' are required.")
            if not isinstance(data_dict_valid["Surface"], list):
                raise ValueError("Surface must be a list of strings.")
            if not isinstance(data_dict_valid["Value"], list):
                raise ValueError("Value must be a list of floats.")
            if not isinstance(data_dict_valid["DOF"], list):
                raise ValueError("DOF must be a list of strings.")
            if len(data_dict_valid["Value"]) != len(data_dict_valid["DOF"]):
                raise ValueError("Value and DOF must have the same length.")
            if not all(dof in ["X", "Y", "Z", "RX", "RY", "RZ"] for dof in data_dict_valid["DOF"]):
                raise ValueError("DOF must be one of 'X', 'Y', 'Z', 'RX', 'RY', 'RZ'.")
            if not all(surface in ["Xmin", "Xmax", "Ymin", "Ymax", "Zmin", "Zmax", "Xmid", "Ymid", "Zmid"]
                       for surface in data_dict_valid["Surface"]):
                raise ValueError("Surface must be one of 'Xmin', 'Xmax', 'Ymin', 'Ymax', 'Zmin', 'Zmax', "
                                 "'Xmid', 'Ymid', 'Zmid'.")


        DOF_map = {"X": 0, "Y": 1, "Z": 2, "RX": 3, "RY": 4, "RZ": 5}
        for key, dict_data in self.boundary_conditions.items():
            if key not in ["Force", "Displacement"]:
                raise ValueError(f"Invalid boundary condition type: {key}. Must be 'Force' or 'Displacement'.")
            for name_condition, data in dict_data.items():
                check_data_boundary_condition_validity(data)
                numeric_DOFs = [DOF_map[dof] for dof in data["DOF"]]
                surface_cells = data.get("SurfaceCells", None)
                self.apply_constraints_nodes(data["Surface"], data["Value"], numeric_DOFs, key, surface_cells)

    def calculate_schur_complement_cells(self):
        """
        Calculate the Schur complement for each cell in the lattice.
        Save the result in the cell.schur_complement attribute.
        """
        schur_complement_calculated: dict = {}

        for cell in self.cells:
            geom_key = tuple(cell.geom_types) if isinstance(cell.geom_types, list) else cell.geom_types
            radius_key = tuple(round(float(r), 8) for r in cell.radii)

            if geom_key not in schur_complement_calculated:
                schur_complement_calculated[geom_key] = {}

            if radius_key not in schur_complement_calculated[geom_key]:
                schur_complement_cell = get_schur_complement(self, cell.index)
                schur_complement_calculated[geom_key][radius_key] = schur_complement_cell
                if self._verbose > 1:
                    print(f"Schur complement calculated for geometry {geom_key} with radii {radius_key}.")
            else:
                schur_complement_cell = schur_complement_calculated[geom_key][radius_key]

            cell.schur_complement = schur_complement_cell


