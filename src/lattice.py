"""
lattice.py

Generate lattice structures with various parameters and properties.

This module provides the Lattice class, which allows for the creation and manipulation of lattice structures,
including the definition of cell dimensions, material properties, and gradient settings. It also supports simulation methods and uncertainty handling.

Created in 2023 by Cadart Thomas, University of technology Belfort Montbéliard.
"""
import os
import pickle
from statistics import mean

import joblib
from scipy.sparse.linalg import splu
import gmsh

from cell import *
from timing import *
from utils import _validate_inputs
from gradient_properties import get_grad_settings, grad_material_setting, grad_settings_constant
from mesh_file.mesh_trimmer import MeshTrimmer

timing = Timing()


class Lattice(object):
    """
    Generate lattice structures with a lot of different parameters
    """

    def __init__(self, cell_size_x: float, cell_size_y: float, cell_size_z: float,
                 num_cells_x: int, num_cells_y: int, num_cells_z: int,
                 geom_types: list[str], radii: list[float], material_name: str = "",
                 grad_radius_property: list = None, grad_dim_property: list = None, grad_mat_property: list = None,
                 uncertainty_node: float = 0.0, enable_periodicity: bool = False, eraser_blocks: list = None,
                 mesh_trimmer: "MeshTrimmer" = None, symmetry_lattice: dict = None,
                 enable_simulation_properties: bool = False, verbose: int = 0) -> None:
        """
        Constructor general for the Lattice class.

        Parameter:
        -----------
        cell_size_x: float
        cell_size_y: float
        cell_size_z: float
            Dimension in each direction of the intial cell in the structure

        num_cells_x: integer
        num_cells_y: integer
        num_cells_z: integer
            Number of cells in each direction in the structure

        geom_types: list of string
            Name of geometry types of the cell structure.
        radii: list of float
            Initial radii geometry
        material_name: string
            Name of the default material in the lattice structure ('Ti-6Al-4V', 'VeroClear'...)
            Possible to add more material in the Material.py file

        Gradient properties
        grad_radius_property: array of data as [GradDimRule,GradDimDirection,GradDimParameters]
            radii gradient on the lattice structure
        grad_dim_property: array of data as [GradRadRule,GradRadDirection,GradRadParameters]
            Cell dimension gradient on the lattice structure
                GradRule => constant, linear, parabolic, sinusoide, exponential
                GradDirection => [bool,bool,bool] set integer to True to active gradient in direction [X,Y,Z] False inactive
                GradParameters => [float, float, float] variable in the gradient rule for each direction [X,Y,Z]
        grad_mat_property: array of data as [Multimat,GradMaterialDirection]
            Material gradient on the lattice structure
                Multimat => Type of multimaterial (0: inactive / 1: multimat by layer / -1: Full random)

        uncertainty_node: float
            Control if adding uncertainties on node position
        enable_periodicity: boolean
            Applying enable_periodicity on the outer box of the lattice structure to calculate penalization method
        eraser_blocks: list of float in dim 6
            (xStart, yStart, zStart, xDim, yDim, zDim) of the erased region
        mesh_trimmer: class object MeshTrimmer
            Class to trim the lattice structure with a mesh.
        symmetry_lattice: dictionary {"sym_plane": string, "sym_point": tuple}
            Data to apply symmetry on the lattice structure
        enable_simulation_properties: boolean
            If True, the lattice will generate properties necessary for simulation and optimization
            And joint beam penalization method will be applied
        _verbose: boolean
            If True, print statistics and information during the lattice generation process
        """
        _validate_inputs(cell_size_x, cell_size_y, cell_size_z, num_cells_x, num_cells_y, num_cells_z, geom_types,
                         radii, material_name, grad_radius_property, grad_dim_property, grad_mat_property,
                         uncertainty_node, enable_periodicity, eraser_blocks)

        self.name_lattice: str = "Lattice"
        self.x_min, self.y_min, self.z_min = None, None, None
        self.x_max, self.y_max, self.z_max = None, None, None
        self.edge_tags = None
        self.face_tags = None
        self.corner_tags = None
        self.lattice_dimension_dict = None
        self.dict_schur_complement = None
        self.objectif_data = None
        self.occupancy_matrix = None

        self.cell_size_x = cell_size_x
        self.cell_size_y = cell_size_y
        self.cell_size_z = cell_size_z
        self.num_cells_x = num_cells_x
        self.num_cells_y = num_cells_y
        self.num_cells_z = num_cells_z
        self.geom_types = geom_types
        self.radii = radii
        self.material_name = material_name
        if grad_radius_property is not None:
            self.grad_radius = get_grad_settings(self.num_cells_x, self.num_cells_y, self.num_cells_z,
                                                 grad_radius_property)
        else:
            self.grad_radius = grad_settings_constant(self.num_cells_x, self.num_cells_y, self.num_cells_z)
        if grad_dim_property is not None:
            self.grad_dim = get_grad_settings(self.num_cells_x, self.num_cells_y, self.num_cells_z, grad_dim_property)
        else:
            self.grad_dim = grad_settings_constant(self.num_cells_x, self.num_cells_y, self.num_cells_z)
        if grad_mat_property is not None:
            self.grad_mat = grad_material_setting(self.num_cells_x, self.num_cells_y, self.num_cells_z,
                                                  grad_mat_property)
        else:
            self.grad_mat = grad_settings_constant(self.num_cells_x, self.num_cells_y, self.num_cells_z,
                                                   material_gradient=True)
        self.enable_simulation_properties = enable_simulation_properties
        self.size_x, self.size_y, self.size_z = self.get_size_lattice()
        self.uncertainty_node = uncertainty_node
        self.enable_periodicity = enable_periodicity  # Warning not working for graded structures
        self.eraser_blocks = eraser_blocks
        self.mesh_trimmer = mesh_trimmer
        self._verbose: int = verbose

        self.cells = []

        # Simulation necessary
        self.free_DOF = None  # Free DOF gradient conjugate gradient method
        self.max_index_boundary = None
        self.global_displacement_index = None
        self.initial_value_objective = None
        self.initial_relative_density_constraint = None
        self.initial_continuity_constraint = None
        self.relative_density_poly = []
        self.relative_density_poly_deriv = []
        self.n_DOF_per_node: int = 6  # Number of DOF per node (3 translation + 3 rotation)
        self.parameter_optimization = []
        self.kriging_model_relative_density = None
        self.penalization_coefficient: float = 1.5  # Fixed with previous optimization
        self.mesh_lattice = None
        self.timing = timing

        # Generate global structure
        self.generate_lattice()

        # Generate important data for the lattice structure
        self.set_tag_classification()
        self.get_lattice_dimensions()
        self.define_beam_node_index()
        self.define_cell_index()
        self.define_cell_neighbours()
        self.set_point_local_tag()
        self.apply_tag_all_point()

        # Simulation necessaries
        if self.enable_simulation_properties == 1:
            self.get_all_angles()
            self.set_beam_node_mod()

            assert self.material_name is not None, "Material name_lattice must be defined for simulation properties."
            # Define global indexation
            self.define_node_index_boundary()
            # Optimization necessary
            self.load_relative_density_model()

        if symmetry_lattice is not None:
            self.apply_symmetry(symmetry_lattice["sym_plane"], symmetry_lattice["sym_point"])

        if self._verbose > 0:
            self.are_cells_identical()
            self.print_statistics_lattice()
            self.timing.summary()

    @classmethod
    def from_json(cls, name_file: str, mesh_trimmer: "MeshTrimmer" = None, verbose: int = 0) -> "Lattice":
        """
        Load a lattice object from a JSON file.

        Parameters:
        -----------
        file_path: str
            Path to the JSON file containing lattice data.

        Returns:
        --------
        Lattice
            The loaded lattice object.
        """
        project_root = Path(__file__).resolve().parent.parent
        json_path = project_root / "preset_lattice" / name_file
        if json_path.suffix != ".json":
            json_path = json_path.with_suffix('.json')

        with open(json_path, 'r') as file:
            data = json.load(file)

        # Geometry
        geometry = data.get("geometry", {})
        cell_size = geometry.get("cell_size", {})
        number_of_cells = geometry.get("number_of_cells", {})

        cell_size_x = cell_size.get("x")
        cell_size_y = cell_size.get("y")
        cell_size_z = cell_size.get("z")
        num_cells_x = number_of_cells.get("x")
        num_cells_y = number_of_cells.get("y")
        num_cells_z = number_of_cells.get("z")
        radii = geometry.get("radii")
        geom_types = geometry.get("geom_types")

        if None in [cell_size_x, cell_size_y, cell_size_z, num_cells_x, num_cells_y, num_cells_z, radii, geom_types]:
            raise ValueError("Missing geometry parameters in JSON file.")

        # Gradient
        gradient = data.get("gradient", {})
        radius_grad = gradient.get("radii", {})
        dim_grad = gradient.get("cell_dimension", {})
        mat_grad = gradient.get("material", {})

        grad_radius_property = [
            radius_grad.get("rule", "constant"),
            [radius_grad.get("direction_x", False),
             radius_grad.get("direction_y", False),
             radius_grad.get("direction_z", False)],
            [radius_grad.get("parameter_x", 0.0),
             radius_grad.get("parameter_y", 0.0),
             radius_grad.get("parameter_z", 0.0)]
        ]

        grad_dim_property = [
            dim_grad.get("rule", "constant"),
            [dim_grad.get("direction_x", False),
             dim_grad.get("direction_y", False),
             dim_grad.get("direction_z", False)],
            [dim_grad.get("parameter_x", 0.0),
             dim_grad.get("parameter_y", 0.0),
             dim_grad.get("parameter_z", 0.0)]
        ]

        grad_mat_property = [
            mat_grad.get("type_beam", 0),
            mat_grad.get("direction", 0)
        ]

        # Supplementary
        supplementary = data.get("suplementary", {})
        uncertainty_node = supplementary.get("node_uncertainty", 0.0)

        # Erased blocks
        erased_blocks_json = supplementary.get("erased_blocks", {})
        erased_blocks = []
        for block in erased_blocks_json.values():
            start = block.get("start_point", {})
            dim = block.get("dimensions_block", {})
            erased_blocks.append([
                start.get("x", 0.0), start.get("y", 0.0), start.get("z", 0.0),
                dim.get("x", 0.0), dim.get("y", 0.0), dim.get("z", 0.0)
            ])

        if len(erased_blocks) == 0:
            erased_blocks = None

        # Symmetry
        symmetries = supplementary.get("symmetries", {})
        symmetry_lattice = None
        if symmetries:
            sym_plane = symmetries.get("plane", None)
            sym_point = symmetries.get("reference_point", {})
            symmetry_lattice = {
                "sym_plane": sym_plane,
                "sym_point": (sym_point.get("x", 0.0),
                              sym_point.get("y", 0.0),
                              sym_point.get("z", 0.0))
            }

        # Simulation activation
        sim_params = data.get("simulation_parameters", {})
        enable_simulation_properties = bool(sim_params.get("enable", False))
        material_name = sim_params.get("material", "VeroClear")
        periodicity = sim_params.get("periodicity", False)

        return cls(cell_size_x, cell_size_y, cell_size_z,
                   num_cells_x, num_cells_y, num_cells_z,
                   geom_types, radii,
                   material_name=material_name,
                   grad_radius_property=grad_radius_property,
                   grad_dim_property=grad_dim_property,
                   grad_mat_property=grad_mat_property,
                   uncertainty_node=uncertainty_node,
                   enable_periodicity=periodicity,
                   eraser_blocks=erased_blocks,
                   symmetry_lattice=symmetry_lattice,
                   enable_simulation_properties=enable_simulation_properties,
                   mesh_trimmer=mesh_trimmer,
                   verbose=verbose)

    @classmethod
    def open_pickle_lattice(cls, file_name: str = "LatticeObject") -> "Lattice":
        """
        Load a lattice pickle from a file.

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
        project_root = Path(__file__).resolve().parent.parent
        path = project_root / "saved_lattice_file" / file_name
        if path.suffix != ".pkl":
            path = path.with_suffix('.pkl')

        if not os.path.exists(path):
            raise FileNotFoundError(f"The file {path} does not exist.")

        with open(path, "rb") as file:
            lattice = pickle.load(file)

        print(f"Lattice loaded successfully from {path}")
        return lattice

    def __repr__(self) -> str:
        string = f"Lattice name_lattice: {self.name_lattice}\n"
        string += f"Dimensions: {self.size_x} x {self.size_y} x {self.size_z}\n"
        string += f"Number of cells: {self.num_cells_x} x {self.num_cells_y} x {self.num_cells_z}\n"
        string += f"Cell size: {self.cell_size_x} x {self.cell_size_y} x {self.cell_size_z}\n"
        string += f"Material: {self.material_name}\n"
        string += f"radii: {self.radii}\n"
        return string

    @timing.timeit
    def generate_lattice(self):
        """
        Generate cells in the lattice structure based on cell size, number of cells, geometry types, and radii.
        Gradient informations and erased regions are also considered during cell generation.
        """
        x_cell_start_init = 0
        y_cell_start_init = 0
        z_cell_start_init = 0
        x_cell_start = 0
        y_cell_start = 0
        z_cell_start = 0
        pos_cell = [0, 0, 0]
        for i in range(self.num_cells_x):
            if i != 0:
                x_cell_start += self.cell_size_x * self.grad_dim[pos_cell[0]][0]
            else:
                x_cell_start = x_cell_start_init
            for j in range(self.num_cells_y):
                if j != 0:
                    y_cell_start += self.cell_size_y * self.grad_dim[pos_cell[1]][1]
                else:
                    y_cell_start = y_cell_start_init
                for k in range(self.num_cells_z):
                    if k != 0:
                        z_cell_start += self.cell_size_z * self.grad_dim[pos_cell[2]][2]
                    else:
                        z_cell_start = z_cell_start_init
                    pos_cell = [i, j, k]
                    initial_cell_size = [self.cell_size_x, self.cell_size_y, self.cell_size_z]
                    start_cell_pos = [x_cell_start, y_cell_start, z_cell_start]
                    if not self.is_not_in_erased_region(start_cell_pos):
                        radius = self.radii
                        new_cell = Cell(pos_cell, initial_cell_size, start_cell_pos, self.geom_types, radius,
                                        self.grad_radius, self.grad_dim, self.grad_mat, self.uncertainty_node,
                                        self._verbose)
                        if self.mesh_trimmer is not None and self.mesh_trimmer.is_cell_in_mesh(new_cell):
                            self.cells.append(new_cell)
                        elif self.mesh_trimmer is None:
                            self.cells.append(new_cell)
                        else:
                            del new_cell
        if len(self.geom_types) > 1:
            self.check_hybrid_collision()

    def set_tag_classification(self) -> None:
        """
        Define tag list classification.
        """
        self.edge_tags = [[102, 104, 106, 107], [100, 108, 105, 111], [101, 109, 103, 110]]
        self.face_tags = [[10, 15], [11, 14], [12, 13]]
        self.corner_tags = [[1000, 1001, 1002, 1003, 1004, 1005, 1006, 1007]]

    def get_size_lattice(self) -> list[float]:
        """
        Computes the size of the lattice along each direction.

        Return:
        ---------
        size_lattice: list of float in dim 3
            Length of the lattice in each direction
        """
        size_lattice = [0.0, 0.0, 0.0]
        for direction in range(3):
            total_length = 0.0
            num_cells = [self.num_cells_x, self.num_cells_y, self.num_cells_z][direction]
            cell_size = [self.cell_size_x, self.cell_size_y, self.cell_size_z][direction]
            gradient_factors = [grad[direction] for grad in self.grad_dim[:num_cells]]

            for factor in gradient_factors:
                total_length += factor * cell_size
            size_lattice[direction] = total_length

        return size_lattice

    @timing.timeit
    def is_not_in_erased_region(self, start_cell_pos: list[float]) -> bool:
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
        if self.eraser_blocks is not None:
            for delPart in self.eraser_blocks:
                inside_erased = all(
                    delPart[direction] <= start_cell_pos[direction] <= delPart[direction] + delPart[direction + 3]
                    for direction in range(3)
                )

                if inside_erased:
                    return True

        return False  # cell removed

    @timing.timeit
    def define_beam_node_index(self) -> None:
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
                        beam.index = nextBeamIndex
                        beamIndexed[beam] = nextBeamIndex
                        nextBeamIndex += 1
                    else:
                        beam.index = beamIndexed[beam]

                for node in [beam.point1, beam.point2]:
                    if node.index is None:
                        if node not in nodeIndexed:
                            node.index = nextNodeIndex
                            nodeIndexed[node] = nextNodeIndex
                            nextNodeIndex += 1
                        else:
                            node.index = nodeIndexed[node]

    @timing.timeit
    def define_cell_index(self) -> None:
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
                    cell.index = nextCellIndex
                    cellIndexed[cell] = nextCellIndex
                    nextCellIndex += 1

    @timing.timeit
    def define_cell_neighbours(self) -> None:
        """
        Define neighbours for each cell in the lattice, with periodic boundaries if enabled.
        """
        cell_dict = {tuple(cell.pos_cell): cell for cell in self.cells}

        neighbor_offsets = [
            (-self.cell_size_x, 0, 0), (self.cell_size_x, 0, 0),
            (0, -self.cell_size_y, 0), (0, self.cell_size_y, 0),
            (0, 0, -self.cell_size_z), (0, 0, self.cell_size_z)
        ]
        localBoundaryBox = False
        if self.occupancy_matrix is None and self.enable_periodicity and self.eraser_blocks is not None:
            self.get_cell_occupancy_matrix()
            localBoundaryBox = True

        for cell in self.cells:
            cell.neighbour_cells = []
            boundaryBox = self.get_relative_boundary_box(cell) if localBoundaryBox else self.get_lattice_boundary_box()

            for offset in neighbor_offsets:
                raw_pos = (
                    cell.pos_cell[0] + offset[0],
                    cell.pos_cell[1] + offset[1],
                    cell.pos_cell[2] + offset[2]
                )

                if self.enable_periodicity:
                    # enable_periodicity X
                    if raw_pos[0] < boundaryBox[0]:
                        neighbor_x = boundaryBox[1] + offset[0]
                    elif raw_pos[0] >= boundaryBox[1]:
                        neighbor_x = boundaryBox[0]
                    else:
                        neighbor_x = raw_pos[0]
                    # enable_periodicity Y
                    if raw_pos[1] < boundaryBox[2]:
                        neighbor_y = boundaryBox[3] + offset[1]
                    elif raw_pos[1] >= boundaryBox[3]:
                        neighbor_y = boundaryBox[2]
                    else:
                        neighbor_y = raw_pos[1]
                    # enable_periodicity Z
                    if raw_pos[2] < boundaryBox[4]:
                        neighbor_z = boundaryBox[5] + offset[2]
                    elif raw_pos[2] >= boundaryBox[5]:
                        neighbor_z = boundaryBox[4]
                    else:
                        neighbor_z = raw_pos[2]

                    neighbor_pos = (neighbor_x, neighbor_y, neighbor_z)
                else:
                    if not (boundaryBox[0] <= raw_pos[0] <= boundaryBox[1] and
                            boundaryBox[2] <= raw_pos[1] <= boundaryBox[3] and
                            boundaryBox[4] <= raw_pos[2] <= boundaryBox[5]):
                        continue
                    neighbor_pos = raw_pos
                if neighbor_pos in cell_dict:
                    cell.add_cell_neighbour(cell_dict[neighbor_pos])

    @timing.timeit
    def get_list_angle_beam(self, beam: "Beam", pointbeams: list["Beam"]) -> tuple[list[float], list[float]]:
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
        non_zero_radiusbeam: list of radii between a considered beam and pointbeams beam list

        Special case when pointbeams is an empty return max angle to minimize penalization zone
        """
        anglebeam = []
        radiusBeam = []
        if len(pointbeams) > 1:
            for beampoint in pointbeams:
                radiusBeam.append(beampoint.radius)
                anglebeam.append(beam.get_angle_between_beams(beampoint, self.enable_periodicity))
        else:  # Not connected beam
            radiusBeam.append(beam.radius)
            anglebeam.append(179.9)
        non_zero_anglebeam = [angle for angle in anglebeam if angle >= 0.01]
        non_zero_radiusbeam = [radius for angle, radius in zip(anglebeam, radiusBeam) if angle >= 0.01]
        return non_zero_anglebeam, non_zero_radiusbeam

    @timing.timeit
    def get_connected_beams(self, beamList: list["Beam"], beam: "Beam") -> tuple[list["Beam"], list["Beam"]]:
        """
        Get all beams connected to the interest beam.

        Parameters:
        -----------
        beamList: list of Beam objects
            List of all beams in the lattice.
        beam: Beam
            Beam of interest.
        """
        point1beams = set()
        point2beams = set()

        tag_checks = [*[(tags, 'corner') for tags in self.corner_tags], *[(tags, 'edge') for tags in self.edge_tags],
                      *[(tags, 'face') for tags in self.face_tags]]

        def is_periodic_connected(p, b_idx, tags_range):
            """
            Check if the point p is periodic connected to the beam index b_idx
            """
            if not p.tag:
                return False
            p_tag_set = set(p.tag)
            if not p_tag_set.intersection(tags_range):
                return False

            bp1_tag = set(b_idx.point1.tag or [])
            bp2_tag = set(b_idx.point2.tag or [])
            if not (bp1_tag.union(bp2_tag)).intersection(tags_range):
                return False

            p_local = set(p.local_tag)
            bp1_local = set(b_idx.point1.local_tag or [])
            bp2_local = set(b_idx.point2.local_tag or [])

            if tags_range == self.corner_tags:
                return bool(p_local.intersection(tags_range) and (bp1_local.union(bp2_local)).intersection(tags_range))
            else:
                return bool(p_local.intersection(tags_range) and (bp1_local.union(bp2_local)).intersection(tags_range))

        for beamidx in beamList:
            if beam.point1 in [beamidx.point1, beamidx.point2]:
                point1beams.add(beamidx)
            if beam.point2 in [beamidx.point1, beamidx.point2]:
                point2beams.add(beamidx)

            if not self.enable_periodicity:
                continue

            for tags_range, _ in tag_checks:
                if is_periodic_connected(beam.point1, beamidx, tags_range):
                    point1beams.add(beamidx)
                if is_periodic_connected(beam.point2, beamidx, tags_range):
                    point2beams.add(beamidx)

        return list(point1beams), list(point2beams)

    @timing.timeit
    def get_all_angles(self) -> None:
        """
        Calculates angles between beams in the lattice.

        Return:
        ---------
        angle:
            data structure => ((beam_index, Angle mininmum point 1, minRad1, Angle mininmum point 2, minRad2))
        """

        @timing.timeit
        def find_min_angle(angles, radii):
            """
            Find the Minimum angle between beams and radii connection to this particular beam
            """
            LValuesMax = 0
            LRadius = None
            LAngle = None
            for radius, angle in zip(radii, angles):
                L = function_penalization_Lzone((radius, angle))
                if L > LValuesMax:
                    LValuesMax = L
                    LRadius = radius
                    LAngle = angle
            return LAngle, LRadius

        # Create the list of beam objects for each cell with neighbors cells
        for cell in self.cells:
            beamList = []
            cellListNeighbours = cell.neighbour_cells
            cellListNeighbours.append(cell)  # Include the cell itself
            for neighbour in cellListNeighbours:
                for beam in neighbour.beams:
                    if beam not in beamList:
                        beamList.append(beam)
            angleList = {}
            for beam in cell.beams:
                # Determine beams on nodes
                point1beams, point2beams = self.get_connected_beams(beamList, beam)
                # Determine angles for all beams connected at the node
                non_zero_anglebeam1, non_zero_radiusbeam1 = self.get_list_angle_beam(beam, point1beams)
                non_zero_anglebeam2, non_zero_radiusbeam2 = self.get_list_angle_beam(beam, point2beams)
                # Find the lowest angle
                LAngle1, LRadius1 = find_min_angle(non_zero_anglebeam1, non_zero_radiusbeam1)
                LAngle2, LRadius2 = find_min_angle(non_zero_anglebeam2, non_zero_radiusbeam2)
                angleList[beam.index] = (LRadius1, round(LAngle1, 2), LRadius2, round(LAngle2, 2))
                beam.set_angle(angleList[beam.index])

    @timing.timeit
    def get_lattice_dimensions(self) -> None:
        """
        Computes extremum values of coordinates in the lattice.

        Return:
        --------
        ExtrumumValues: tuple of floats (x_min, x_max, y_min, y_max, z_min, z_max)
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

        self.x_min, self.x_max = min(x_values), max(x_values)
        self.y_min, self.y_max = min(y_values), max(y_values)
        self.z_min, self.z_max = min(z_values), max(z_values)

        self.set_lattice_dimensions_dict()

    def set_lattice_dimensions_dict(self):
        """
        Set lattice dimensions in a dictionary format
        """
        self.lattice_dimension_dict = {
            "x_min": self.x_min,
            "x_max": self.x_max,
            "y_min": self.y_min,
            "y_max": self.y_max,
            "z_min": self.z_min,
            "z_max": self.z_max
        }

    @timing.timeit
    def set_beam_node_mod(self) -> None:
        """
        Modifies beam and node data to model lattice structures for simulation with rigidity penalization at node
        """
        for cell in self.cells:
            beamsToRemove = []
            beamToAdd = []

            for beam in cell.beams:
                lengthMod = beam.get_length_mod()
                pointExt1 = beam.get_point_on_beam_at_distance(lengthMod[0], 1)
                pointExt1.node_mod = True
                pointExt2 = beam.get_point_on_beam_at_distance(lengthMod[1], 2)
                pointExt2.node_mod = True

                b1 = Beam(beam.point1, pointExt1, beam.radius, beam.material, beam.type_beam)
                b1.set_beam_mod()
                b2 = Beam(pointExt1, pointExt2, beam.radius, beam.material, beam.type_beam)
                b3 = Beam(pointExt2, beam.point2, beam.radius, beam.material, beam.type_beam)
                b3.set_beam_mod()

                beamToAdd.append((b1, b2, b3))

                beamsToRemove.append(beam)

            for addingBeam in beamToAdd:
                cell.add_beam(addingBeam)

            for beam in beamsToRemove:
                cell.remove_beam(beam)

        # Update index
        self.define_beam_node_index()

    def remove_cell(self, index: int) -> None:
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

    def find_minimum_beam_length(self) -> float:
        """
        Find minimum beam length

        Returns:
        --------
        minLength: float
            Length of the smallest beam in the lattice
        """
        minLength = float('inf')
        for cell in self.cells:
            for beam in cell.beams:
                if minLength > beam.length > 0.0001 and (beam.type_beam == 0 or beam.type_beam == 2):
                    minLength = beam.length
        return minLength

    def get_tag_list(self) -> list[int]:
        """
        Get the tag for all unique points in the lattice.

        Returns:
        --------
        tagList: list of int
            List of all tags of each unique point in the lattice.
        """
        tagList = []
        seen_nodes = set()
        for cell in self.cells:
            for beam in cell.beams:
                for node in (beam.point1, beam.point2):
                    if node not in seen_nodes:
                        tagList.append(node.tag)
                        seen_nodes.add(node)
        return tagList

    def get_tag_list_boundary(self) -> list[int]:
        """
        Get the tag for boundary points in the lattice.

        Returns:
        --------
        tagList: list of int
            List of all tags of each point in lattice
        """
        boundary_box_lattice = self.get_lattice_boundary_box()
        tagList = []
        seen_nodes = set()
        for cell in self.cells:
            for beam in cell.beams:
                for node in [beam.point1, beam.point2]:
                    if node not in seen_nodes and node.is_on_boundary(boundary_box_lattice):
                        tagList.append(node.tag)
                        seen_nodes.add(node)
        return tagList

    @timing.timeit
    def apply_tag_all_point(self) -> None:
        """
        Assign a tag to all nodes in the lattice structure.
        Tags are assigned relative to either the global bounding box
        or a local (cell-relative) bounding box if erased parts are used.
        """
        use_local_box = self.eraser_blocks is not None

        if self.occupancy_matrix is None and use_local_box:
            self.get_cell_occupancy_matrix()

        global_box = None if use_local_box else self.get_lattice_boundary_box()

        for cell in self.cells:
            local_box = self.get_relative_boundary_box(cell) if use_local_box else None
            for beam in cell.beams:
                for node in [beam.point1, beam.point2]:
                    tag = node.tag_point(local_box if use_local_box else global_box)
                    node.tag = tag

    def get_cell_occupancy_matrix(self):
        """
        Generate a 3D boolean matrix indicating presence of a cell at each (i, j, k) position.

        Returns:
        --------
        occupancy_matrix: np.ndarray of shape (num_cells_x, num_cells_y, num_cells_z)
            True if a cell is present at the corresponding position, False otherwise.
        """
        self.occupancy_matrix = np.empty((self.num_cells_x, self.num_cells_y, self.num_cells_z), dtype=object)
        for cell in self.cells:
            i, j, k = cell.pos_cell
            self.occupancy_matrix[i, j, k] = cell

    def get_cells_at_index(self, axis: str, index: int) -> list:
        """
        Get all cells at a specific index along a specified axis.

        Parameters:
        -----------
        axis: str
            Axis to query ('x', 'y', or 'z').
        index: int
            Index along the specified axis.
        """
        if axis == 'x':
            return [c for c in self.occupancy_matrix[index, :, :].flatten() if c is not None]
        elif axis == 'y':
            return [c for c in self.occupancy_matrix[:, index, :].flatten() if c is not None]
        elif axis == 'z':
            return [c for c in self.occupancy_matrix[:, :, index].flatten() if c is not None]
        else:
            raise ValueError("Axis must be 'x', 'y', or 'z'")

    def get_relative_boundary_box(self, cell) -> list[float]:
        """
        Get the relative boundary box of a cell in the lattice.
        It corresponds to the minimum and maximum dimension of the lattice for each axis with cell continuity.
        Useful for structures with erased parts or periodic boundaries.

        Parameters:
        -----------
        cell: Cell object
            The cell for which the boundary box is computed.
        """
        x, y, z = cell.pos_cell
        bbox = []

        for axis, index in zip(
                ['x', 'y', 'z'],
                [x, y, z],
        ):
            arrayCell = self.get_cells_at_index(axis, index)
            min_val, max_val = np.inf, -np.inf

            for cellIn in arrayCell:
                cellInBoundaryBox = cellIn.boundary_box
                if axis == 'x':
                    bounds = cellInBoundaryBox[0:2]
                elif axis == 'y':
                    bounds = cellInBoundaryBox[2:4]
                elif axis == 'z':
                    bounds = cellInBoundaryBox[4:6]

                if bounds[0] < min_val:
                    min_val = bounds[0]
                if bounds[1] > max_val:
                    max_val = bounds[1]

            bbox.append(min_val)
            bbox.append(max_val)

        return bbox

    def get_lattice_boundary_box(self) -> list[float]:
        """
        Get the boundary box of the lattice
        """
        return [self.x_min, self.x_max, self.y_min, self.y_max, self.z_min, self.z_max]

    def get_connected_node(self, node: "Point") -> list["Point"]:
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

    def find_boundary_beams(self) -> list["Beam"]:
        """
        Find boundary beams and change the type_beam of beam

        Return:
        -------
        boundaryBeams: List of a beam object
        """
        boundary_box_lattice = self.get_lattice_boundary_box()
        boundaryBeams = []
        for cell in self.cells:
            for beam in cell.beams:
                if beam.point1.is_on_boundary(boundary_box_lattice) or beam.point2.is_on_boundary(boundary_box_lattice):
                    beam.type_beam = 2
                    boundaryBeams.append(beam)
        return boundaryBeams

    def find_boundary_nodes(self) -> list["Point"]:
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
                    if (node.x == self.x_min or node.x == self.x_max or node.y == self.y_min or node.y == self.y_max or
                            node.z == self.z_min or node.z == self.z_max):
                        boundaryNodes.append(node)
        return boundaryNodes

    def get_name_lattice(self) -> str:
        """
        Determine the name_lattice of the lattice

        Returns:
        ---------
        name_lattice: string
        """
        if len(self.geom_types) == 1:
            self.name_lattice = self.geom_types[0]
        else:
            self.name_lattice = "Hybrid_"
            for i, geom_type in enumerate(self.geom_types):
                self.name_lattice += geom_type + "_"
        if self.enable_simulation_properties == 1:
            self.name_lattice += "_Mod"
        return self.name_lattice

    def check_hybrid_collision(self) -> None:
        """
        Check if beam in hybrid configuration is cut by a point in the geometry
        Change the beam configuration of collisionned beams
        """
        for cell in self.cells:
            cellPoints = cell.get_list_points()
            for node in cellPoints:
                for beam in cell.beams:
                    if beam.is_point_on_beam(node):
                        typeBeamToRemove = beam.type_beam  # Get beam to remove type_beam to apply in new separated beams
                        beam1 = Beam(beam.point1, node, beam.radius, beam.material, typeBeamToRemove)
                        beam2 = Beam(beam.point2, node, beam.radius, beam.material, typeBeamToRemove)
                        cell.remove_beam(beam)
                        cell.add_beam(beam1)
                        cell.add_beam(beam2)

    def get_node_coordinates_data(self) -> list[list[float]]:
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
                        posData.append([node.coordinates])
                        nodeAlreadyAdded.append(node)
        return posData

    def get_edge_index(self) -> list[list[int]]:
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

    def get_beam_type(self) -> list:
        """
        Retrieves beam type_beam data for the lattice.

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
                    beamType.append([beam.type_beam])
        return beamType

    def get_all_beam_length(self) -> list[list[float]]:
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

    def change_beam_radius_data_hybrid_lattice(self, hybridRadiusData: list[float]) -> None:
        """
        Change radii data for hybrid lattice

        Parameters:
        ------------
        hybridRadiusData: list of float
            List of radii data for hybrid lattice
        """
        if len(hybridRadiusData) != len(self.radii):
            raise ValueError("Invalid hybrid radii data.")
        for cell in self.cells:
            for beam in cell.beams:
                if beam.beam_mod:
                    beam.radius = hybridRadiusData[beam.type_beam] * beam.penalization_coefficient
                else:
                    beam.radius = hybridRadiusData[beam.type_beam]

    @timing.timeit
    def delete_duplicated_beams(self) -> None:
        """
        Delete duplicated beams in the lattice
        """
        beamList = []
        for cell in self.cells:
            for beam in cell.beams:
                if beam not in beamList:
                    beamList.append(beam)
                else:
                    cell.remove_beam(beam)

    def get_relative_density_constraint(self, relativeDensityMax, geomScheme) -> float:
        """
        Get relative density of the lattice
        """
        relativeDensity = self.get_relative_density(geomScheme)
        print("Relative density: ", relativeDensity)
        error = relativeDensity - relativeDensityMax
        print("Relative density maximum: ", relativeDensityMax)
        print("Relative density error: ", error)
        return error

    def get_relative_density(self, geom_scheme=None) -> float:
        """
        Get mean relative density of all cells in lattice

        Returns:
        --------
        meanRelDens: float
            Mean relative density of the lattice
        """
        cellRelDens = []
        for cell in self.cells:
            if self.kriging_model_relative_density is not None:
                cellRelDens.append(cell.get_relative_density_kriging(self.kriging_model_relative_density, geom_scheme))
            else:
                cellRelDens.append(cell.relative_density)
        meanRelDens = mean(cellRelDens)
        return meanRelDens

    def define_relative_density_function(self, degree: int = 3) -> None:
        """
        Define relative density function
        Possible to define a more complex function with dependency on hybrid cells

        Parameters:
        -----------
        degree: int
            Degree of the polynomial function
        """
        if len(self.relative_density_poly) == 0:
            fictiveCell = Cell([0, 0, 0], [self.cell_size_x, self.cell_size_y, self.cell_size_z], [0, 0, 0],
                               self.geom_types, self.radii, self.grad_radius, self.grad_dim, self.grad_mat,
                               self.uncertainty_node, self._verbose)
            domainRadius = np.linspace(0.01, 0.1, 10)
            for idxRad in range(len(self.radii)):
                radius = np.zeros(len(self.radii))
                relativeDensity = []
                for domainIdx in domainRadius:
                    radius[idxRad] = domainIdx
                    fictiveCell.change_beam_radius([radius], self.grad_radius)
                    relativeDensity.append(fictiveCell.relative_density())
                poly_coeffs = np.polyfit(domainRadius, relativeDensity, degree).flatten()
                poly = np.poly1d(poly_coeffs)
                self.relative_density_poly.append(poly)
                self.relative_density_poly_deriv.append(poly.deriv())

    def get_relative_density_gradient(self) -> list[float]:
        """
        Get relative density gradient of the lattice

        Returns:
        --------
        grad: list of float
            Gradient of relative density
        """
        if len(self.relative_density_poly) == 0:
            self.define_relative_density_function()
        if len(self.cells[0].radii) != len(self.relative_density_poly):
            raise ValueError("Invalid radii data.")

        grad = []
        for cell in self.cells:
            grad.append(cell.get_relative_density_gradient())
        return grad

    def get_relative_density_gradient_kriging(self, geom_scheme=None) -> list[float]:
        """
        Get relative density gradient of the lattice using kriging model

        Returns:
        --------
        grad: list of float
            Gradient of relative density
        """
        grad = []
        numberOfCells = len(self.cells)
        if geom_scheme is None or len(geom_scheme) != 3:
            geom_scheme = [i < len(self.radii) for i in range(3)]

        for cell in self.cells:
            gradient3Geom = cell.get_relative_density_gradient_kriging(self.kriging_model_relative_density,
                                                                       geom_scheme) / numberOfCells
            grad.extend(gradient3Geom[geom_scheme])
        return grad

    def get_radius_continuity_difference(self, delta: float = 0.01) -> list[float]:
        """
        Get the difference in radii between connected beams in the lattice

        Parameters:
        -----------
        delta: float
            Minimum difference in radii between connected cells
        """
        radiusContinuityDifference = []
        for cell in self.cells:
            radiusCell = cell.radii
            for neighbours in cell.neighbour_cells:
                for rad in range(len(radiusCell)):
                    radiusContinuityDifference.append((radiusCell[rad] - neighbours.radii[rad]) ** 2 - delta ** 2)
        return radiusContinuityDifference

    def get_radius_continuity_jacobian(self) -> np.ndarray:
        """
        Compute the Jacobian of the radii continuity constraint.

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
            radiusCell = cell.radii
            for neighbour in cell.neighbour_cells:
                radiusNeighbour = neighbour.radii
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

        jacobian = np.zeros((constraint_index, self.get_number_parameters_optimization()))
        for r, c, v in zip(rows, cols, values):
            jacobian[r, c] = v

        return jacobian

    def change_beam_radius_depending_type(self, typeToChange: int, newRadius: float) -> None:
        """
        Change radii of beam for specific type_beam

        Parameters:
        -----------
        typeToChange: int
            Type of beam to change
        newRadius: float
            New radii of beam
        """
        for cell in self.cells:
            for beam in cell.beams:
                if beam.type_beam == typeToChange:
                    beam.radius = newRadius

    def get_number_beams(self) -> int:
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

    def get_number_nodes(self) -> int:
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

    def apply_boundary_conditions_surface(self, surfaceNames: list[str], valueDisplacement: list[float],
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

    def apply_constraints_nodes(self, surfaceNames: list[str], value: list[float], DOF: list[int],
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
        type_beam: str
            Type of constraint (Displacement, Force)

        """
        if surfaceNamePoint is None:
            pointSet = self.find_point_on_lattice_surface(surfaceNames)
        else:
            pointSet = self.find_point_on_lattice_surface_complex(surfaceNames, surfaceNamePoint)

        indexBoundaryList = {point.index_boundary for point in pointSet}

        for cell in self.cells:
            for beam in cell.beams:
                for node in [beam.point1, beam.point2]:
                    if node.index_boundary in indexBoundaryList:
                        for val, DOFi in zip(value, DOF):
                            if type == "Displacement":
                                node.displacement_vector[DOFi] = val
                                node.fix_DOF([DOFi])
                            elif type == "Force":
                                node.applied_force[DOFi] = val

    def find_point_on_lattice_surface(self, surfaceNames: list[str]) -> set["Point"]:
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
            raise ValueError("Invalid surface name_lattice(s).")

        cellLists = [set(self.get_cell_on_surface(surface)) for surface in surfaceNames]
        cellList = set.intersection(*cellLists)  # Union of all cell indices from given surfaces

        if self.cells[-1].index < max(cellList, default=-1):
            raise ValueError("Invalid cell index, some cells do not exist.")

        pointSet = None
        for cell in self.cells:
            if cell.index in cellList:
                cellPointSets = [set(cell.get_point_on_surface(surface)) for surface in surfaceNames]
                if cellPointSets:
                    if pointSet is None:
                        pointSet = set.intersection(*cellPointSets)
                    else:
                        pointSet.update(set.intersection(*cellPointSets))
        pointSet = pointSet if pointSet is not None else set()

        if pointSet == set():
            raise ValueError("No points found on the specified surfaces.")

        return pointSet

    def find_point_on_lattice_surface_complex(self, surfaceNamesCell: list[str], surfaceNamePoint: list[str]) \
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
            raise ValueError("Invalid surface name_lattice(s).")

        cellLists = [set(self.get_cell_on_surface(surface)) for surface in surfaceNamesCell]
        cellList = set.intersection(*cellLists)  # Union of all cell indices from given surfaces

        if self.cells[-1].index < max(cellList, default=-1):
            raise ValueError("Invalid cell index, some cells do not exist.")

        pointSet = None
        for cell in self.cells:
            if cell.index in cellList:
                cellPointSets = [set(cell.get_point_on_surface(surface)) for surface in surfaceNamePoint]
                if cellPointSets:
                    if pointSet is None:
                        pointSet = set.intersection(*cellPointSets)
                    else:
                        pointSet.update(set.intersection(*cellPointSets))
        pointSet = pointSet if pointSet is not None else set()

        if pointSet == set():
            raise ValueError("No points found on the specified surfaces.")

        return pointSet

    def apply_boundary_conditions_node(self, nodeList: list[int], valueDisplacement: list[float],
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
        #TODO: Check if used
        if self.get_number_nodes() < max(nodeList):
            raise ValueError("Invalid node index, node do not exist.")

        indexBoundaryList = []
        for node in nodeList:
            if node < 0 or node >= self.get_number_nodes():
                raise ValueError("Node index out of range.")
            indexBoundaryList.append(node)

        for cell in self.cells:
            for beam in cell.beams:
                for node in [beam.point1, beam.point2]:
                    if node.index in indexBoundaryList:
                        for val, DOFi in zip(valueDisplacement, DOF):
                            node.displacement_vector[DOFi] = val
                            node.fix_DOF([DOFi])

    def apply_force_on_surface(self, surfaceName: list[str], valueForce: list[float], DOF: list[int]) -> None:
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

    def set_displacement_with_vector(self, displacementMatrix: list[float]) -> None:
        """
        Set displacement on the lattice with vector

        Parameters:
        -----------
        displacementMatrix: list of float of dim n_nodes*n_dofperNode
            Displacement matrix to apply to the lattice
        """
        for cell in self.cells:
            nodeInOrder = cell.get_node_order_to_simulate()
            idxNode = 0
            for idx, node in enumerate(nodeInOrder.values()):
                if node is not None:
                    node.displacement_vector = displacementMatrix[idxNode]
                    node.fix_DOF([i for i in range(6)])
                    idxNode += 1

    def get_global_displacement(self, withFixed: bool = False, OnlyImposed: bool = False, printLevel=0) \
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
        if printLevel > 2:
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

    def set_point_local_tag(self) -> None:
        """
        Set local tag for all points in the lattice based on their position within the cell boundary box.
        """
        for cell in self.cells:
            for beam in cell.beams:
                for node in [beam.point1, beam.point2]:
                    localTag = node.tag_point(cell.boundary_box)
                    node.set_local_tag(localTag)

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

    def getGlobalReactionForceWithoutFixedDOF(self, globalReactionForce: dict, rightHandSide: bool = False) \
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

    def set_free_DOF(self):
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

    def build_LU_schur_complement(self, dictSchurComplement: dict = None) -> splu:
        """
        Build LU decomposition of the Schur complement matrix for the lattice

        Parameters:
        -----------
        schurComplementMatrix: coo_matrix
            Schur complement matrix of the lattice
        """
        #TODO : Voir quoi faire avec ça
        from ConjugateGradientMethod.Utils_Schur import loadSchurComplement, getSref_nearest

        if dictSchurComplement is None:
            nameFileSchur = "ConjugateGradientMethod/schurComplement/Hybrid_" + str(
                len(self.geom_types)) + ".npz"
            dictSchurComplement = loadSchurComplement(nameFileSchur)

        self.build_coupling_operator_cells()
        globalSchurComplement = coo_matrix((self.free_DOF, self.free_DOF))
        for cell in self.cells:
            schurComplementMatrix = coo_matrix(getSref_nearest(cell.radii, SchurDict=dictSchurComplement,
                                                               printing=False))
            globalSchurComplement += cell.build_preconditioner(schurComplementMatrix)

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

    def get_cell_on_surface(self, surface: str) -> list[int]:
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
            raise ValueError("Invalid surface name_lattice.")

        mid_dict = {
            "Xmid": 0.5 * (self.x_min + self.x_max),
            "Ymid": 0.5 * (self.y_min + self.y_max),
            "Zmid": 0.5 * (self.z_min + self.z_max)
        }

        surface_dict = {
            "Xmin": ("x", self.x_min),
            "Xmax": ("x", self.x_max),
            "Ymin": ("y", self.y_min),
            "Ymax": ("y", self.y_max),
            "Zmin": ("z", self.z_min),
            "Zmax": ("z", self.z_max),
            "Xmid": ("x", mid_dict["Xmid"]),
            "Ymid": ("y", mid_dict["Ymid"]),
            "Zmid": ("z", mid_dict["Zmid"])
        }

        axis, valueSurface = surface_dict[surface]

        cellTagList = []
        for cell in self.cells:
            listPointOnSurface = cell.get_point_on_surface(surface)
            for point in listPointOnSurface:
                if getattr(point, axis) == valueSurface:
                    cellTagList.append(cell.index)
                    break

        return cellTagList

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

    def unnormalize_r(self, r_norm):
        """
        Denormalize optimization parameters

        Parameters:
        -----------
        r_norm: float
            Normalized optimization parameter
        """
        #TODO : To delete
        r_min, r_max = 0.01, 0.1
        return r_min + (r_max - r_min) * r_norm  # Dé-normalisation

    def set_optimization_parameters(self, optimizationParameters: list[float], geomScheme: list[bool]) -> None:
        """
        Set optimization parameters for the lattice

        Parameters:
        -----------
        optimizationParameters: list of float
            List of optimization parameters
        geomScheme: list of bool
            List of N boolean values indicating the scheme of geometry to optimize
        """
        if len(optimizationParameters) != self.get_number_parameters_optimization(geomScheme):
            raise ValueError("Invalid number of optimization parameters.")

        if geomScheme is None:
            numberOfParametersPerCell = len(self.geom_types)
        else:
            numberOfParametersPerCell = sum(geomScheme)

        for cell in self.cells:
            startIdx = cell.index * numberOfParametersPerCell
            endIdx = (cell.index + 1) * numberOfParametersPerCell
            radius = optimizationParameters[startIdx:endIdx]

            if len(radius) != len(cell.radii):
                # Reconstruct the full radii vector based on geomScheme
                full_radius = []
                i = 0  # index for radii (optimization vector)
                for keep, old in zip(geomScheme, cell.radii):
                    if keep:
                        full_radius.append(radius[i])
                        i += 1
                    else:
                        full_radius.append(old)
                radius = full_radius

            cell.change_beam_radius(radius, self.grad_radius)

    def calculate_objective(self, typeObjective: str) -> float:
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
            reactionForce = self.get_global_reaction_force(appliedForceAdded=True)
            reaction_force_array = np.array(list(reactionForce.values())).flatten()
            displacement = np.array(self.get_global_displacement(OnlyImposed=True)[0])
            objective = 0.5 * np.dot(reaction_force_array, displacement)
            if self._verbose > 2:
                np.set_printoptions(threshold=np.inf)
                print("Reaction force: ", reaction_force_array[displacement != 0])
                print("Displacement: ", displacement[displacement != 0])
                print("Compliance: ", objective)
        elif typeObjective == "Displacement":
            setNode = self.find_point_on_lattice_surface(surfaceNames=self.objectif_data["surface"])
            displacements = []
            for node in setNode:
                for dof in self.objectif_data["DOF"]:
                    if dof < 0 or dof > 5:
                        raise ValueError("Invalid degree of freedom index.")
                    displacements.append(node.displacement_vector[dof])
            displacements = np.array(displacements)
            objective = sum(abs(displacements)) / len(displacements)

        elif typeObjective == "Stiffness":
            pass
        return objective

    def get_number_parameters_optimization(self, geomScheme) -> int:
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
                numParameters += len(cell.radii)
            else:
                numParameters += sum(geomScheme)
        return numParameters

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

    def expand_schur_to_full_basis(self, SchurReduced, nodeInOrder):
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
        #TODO : Check if used
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

    def print_statistics_lattice(self):
        """
        Print statistics about the lattice
        """
        print("Lattice name: ", self.get_name_lattice())
        print("Number of cells: ", len(self.cells))
        print("Number of beams: ", self.get_number_beams())
        print("Number of nodes: ", self.get_number_nodes())
        print("Relative density approximate: ", self.get_relative_density())
        radMax, radMin = self.get_beam_radius_min_max()
        print("radii max: ", radMax)
        print("radii min: ", radMin)
        print("Lattice dimensions: ", self.get_size_lattice())

    def get_beam_radius_min_max(self):
        """
        Get the maximum and minimum radii of the lattice

        Returns:
        --------
        radMax: float
            Maximum radii of the lattice
        radMin: float
            Minimum radii of the lattice
        """
        radMin = 1000000
        radMax = 0
        for cell in self.cells:
            for beam in cell.beams:
                if not beam.beam_mod:
                    radMin = min(radMin, beam.radius)
                    radMax = max(radMax, beam.radius)
        return radMax, radMin

    @timing.timeit
    def cut_beam_with_mesh_trimmer(self):
        """
        Cut beams in the lattice using the mesh trimmer.
        """
        if self.mesh_trimmer is None:
            raise ValueError("A mesh object must be assigned to the lattice before cutting beams.")
        self.mesh_trimmer.cut_beams_at_mesh_intersection(self.cells)

    @timing.timeit
    def apply_symmetry(self, symmetry_plane: str, reference_point: tuple = (0, 0, 0)) -> None:
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
        symmetry_plane = symmetry_plane.upper()
        if symmetry_plane not in ["XY", "XZ", "YZ", "X", "Y", "Z"]:
            raise ValueError("Invalid symmetry plane. Choose from 'XY', 'XZ', 'YZ', 'X', 'Y', or 'Z'.")

        x_ref, y_ref, z_ref = reference_point
        new_cells = []
        node_map = {}

        for cell in self.cells:
            new_pos = list(cell.pos_cell)
            new_start_pos = list(cell.coordinate_cell)
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
                    Beam(node_map[new_point1], node_map[new_point2], beam.radius, beam.material, beam.type_beam))

            # Create a new mirrored cell
            new_cell = Cell(new_pos, cell.cell_size, new_start_pos, cell.geom_types, cell.radii, cell.grad_radius,
                            cell.grad_dim, cell.grad_mat, cell.uncertainty_node, self._verbose)

            new_cell.beams = mirrored_beams
            new_cells.append(new_cell)

        self.cells.extend(new_cells)
        self.get_lattice_dimensions()  # Recalculate the lattice boundaries

    @timing.timeit
    def load_relative_density_model(self, model_path="Lattice/saved_lattice_file/RelativeDensityKrigingModel.pkl"):
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
            self.kriging_model_relative_density = gpr

    def delete_beams_under_radius_threshold(self, threshold: float = 0.01) -> None:
        """
        Delete beams with radii under a certain threshold

        Parameters:
        -----------
        threshold: float
            Threshold value for beam radii
        """
        for cell in self.cells:
            beamsToRemove = []
            for beam in cell.beams:
                if beam.radius <= threshold:
                    beamsToRemove.append(beam)
            for beam in beamsToRemove:
                cell.remove_beam(beam)

    def delete_beams_with_geom_scheme(self, geomScheme: list[bool]) -> None:
        """
        Delete beams based on the geometry scheme.
        If geomScheme[i] is False, the beam of type_beam i will be removed.
        Usefull for hybrid lattices geometry.

        Parameters:
        -----------
        geomScheme: list of bool
            List of N boolean values indicating the scheme of geometry to optimize
            If geomScheme[i] is False, the beam of type_beam i will be removed.
        """
        for cell in self.cells:
            beamsToRemove = []
            for beam in cell.beams:
                if not geomScheme[beam.type_beam]:
                    beamsToRemove.append(beam)
            for beam in beamsToRemove:
                cell.remove_beam(beam)

    def set_objective_data(self, objectifData: dict) -> None:
        """
        Add objective data to the lattice

        Parameters:
        -----------
        objectif_data: dict
            Dictionary containing objective data
        """
        #TODO : Voir la gestion des objectifs
        self.objectif_data = objectifData

    @timing.timeit
    def generate_mesh_lattice_Gmsh(self, cutMeshAtBoundary: bool = False, meshSize: float = 0.05,
                                   nameMesh: str = "Lattice", runGmshApp: bool = False,
                                   saveMesh: bool = False, saveSTL: bool = True, volume_computation: bool = False):
        """
        Generate a mesh representation of the lattice structure using GMSH.

        Parameters:
        -----------
        cutMeshAtBoundary: bool
            If True, cut the mesh at the bounding box of the lattice.
        meshSize: float
            Size of the mesh elements.
        runGmshApp: bool
            If True, run the GMSH application to visualize the mesh.
        saveMesh: bool
            If True, save the mesh to a file.
        nameMesh: str
            Name of the mesh file to save.
        """
        if self.enable_simulation_properties == 1:
            raise ValueError("mesh_file generation is not available for the current simulation method.")

        gmsh.initialize()
        gmsh.option.setNumber("General.Verbosity", 1)
        gmsh.model.add(nameMesh)
        dim = 3  # Dimension of the mesh

        all_tags = []
        for cell in self.cells:
            for beam in cell.beams:
                p1 = np.array(beam.point1.coordinates)
                p2 = np.array(beam.point2.coordinates)
                direction = p2 - p1
                if beam.index not in all_tags:
                    beamMesh = gmsh.model.occ.addCylinder(*p1, *direction, beam.radius, tag=beam.index)
                    all_tags.append(beamMesh)

        # Merge all beams into a single entity
        beam_entities = [(dim, tag) for tag in all_tags]
        lattice = gmsh.model.occ.fragment(beam_entities, [])

        if cutMeshAtBoundary:
            # Bounding box definition
            x0, y0, z0 = self.x_min, self.y_min, self.z_min
            dx, dy, dz = self.x_max - x0, self.y_max - y0, self.z_max - z0
            box = gmsh.model.occ.addBox(x0, y0, z0, dx, dy, dz)
            gmsh.model.occ.intersect(lattice[0], [(dim, box)], removeObject=True)

        gmsh.model.occ.synchronize()

        # Define mesh size
        points = gmsh.model.getEntities(dim=0)
        gmsh.model.mesh.setSize(points, meshSize)
        gmsh.model.occ.synchronize()

        if volume_computation:
            self.get_volume_mesh()


        gmsh.model.mesh.generate(dim)

        if runGmshApp:
            gmsh.fltk.run()

        project_root = Path(__file__).resolve().parent.parent
        if saveMesh:
            path = project_root / "mesh_file" / f"{nameMesh}.msh"
            gmsh.write(str(path))
            print("mesh_file saved at ", path)

        if saveSTL:
            path = project_root / "mesh_file" / f"{nameMesh}.stl"
            gmsh.write(str(path))
            print("mesh_file saved at ", path)

        gmsh.finalize()

    @timing.timeit
    def get_volume_mesh(self):
        """
        Compute the total volume of the mesh and the relative density of the lattice.
        """
        vol = 0.0
        for dim, tag in gmsh.model.getEntities(dim=3):
            vol += gmsh.model.occ.getMass(dim, tag)
        print(f"Total volume: {vol}")
        bb = self.get_lattice_boundary_box()
        volume_lattice = bb[1] - bb[0] * bb[3] - bb[2] * bb[5] - bb[4]
        relative_density = vol / volume_lattice
        print(f"Relative density: {relative_density:.4f}")

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
            "geom_types",
            "radii",
            "cell_size",
            "grad_radius",
            "grad_dim",
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
                if not b1.is_identical_to(b2, reference.cell_size):
                    print(b1, b2)
                    print(Fore.RED + f"Difference found in beam {j} between cell 0 and cell {i}" + Style.RESET_ALL)
                    return False

        print(Fore.GREEN + "All cells are identical." + Style.RESET_ALL)
        return True
