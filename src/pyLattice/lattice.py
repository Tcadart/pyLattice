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
from typing import TYPE_CHECKING

import gmsh

from .cell import *
from .timing import *
from .utils import _validate_inputs_lattice, open_lattice_parameters
from .gradient_properties import get_grad_settings, grad_material_setting, grad_settings_constant

if TYPE_CHECKING:
    from data.inputs.mesh_file.mesh_trimmer import MeshTrimmer

timing = Timing()


class Lattice(object):
    """
    Generate lattice structures with a lot of different parameters
    """

    def __init__(self, name_file: str, mesh_trimmer: "MeshTrimmer" = None, verbose: int = 0):
        """
        Constructor from a JSON file to create a lattice object.

        Parameters:
        -----------
        name_file: str
            Name of the JSON file containing lattice parameters.
        mesh_trimmer: MeshTrimmer, optional
            MeshTrimmer object to trim the lattice.
        _verbose: int, optional
            Verbosity level for logging.
        """
        self.mesh_trimmer = mesh_trimmer
        self._verbose: int = verbose
        self.timing = timing

        self.name_lattice: str = "Lattice"
        self.x_min, self.y_min, self.z_min = None, None, None
        self.x_max, self.y_max, self.z_max = None, None, None
        self.cell_size_x, self.cell_size_y, self.cell_size_z = 0, 0, 0
        self.num_cells_x, self.num_cells_y, self.num_cells_z = 0, 0, 0
        self.size_x, self.size_y, self.size_z = 0, 0, 0
        self.radii = None
        self.geom_types = None
        self.grad_dim, self.grad_radius, self.grad_mat = None, None, None
        self.symmetry_lattice = None
        self.uncertainty_node = None
        self.eraser_blocks = None
        self.enable_periodicity = False # Warning not working for graded structures
        self.enable_simulation_properties = False
        self._simulation_flag = False
        self._optimization_flag = False

        self.cells = []
        self.beams = set()
        self.nodes = set()
        self.edge_tags = None
        self.face_tags = None
        self.corner_tags = None
        self.lattice_dimension_dict = None
        self.occupancy_matrix = None
        self.mesh_lattice = None

        # Generate global structure
        self.extract_parameters_from_json(name_file)
        self.generate_lattice()

        if self.symmetry_lattice is not None:
            self.apply_symmetry()

        self.define_size_lattice()

        # Generate important data for the lattice structure
        self.set_tag_classification()
        self.define_lattice_dimensions()
        self.define_beam_node_index()
        self.define_cell_index()
        self.define_cell_neighbours()
        self.define_connected_beams_for_all_nodes()
        self.apply_tag_all_point()

        if self._verbose > 0:
            self.are_cells_identical()
            self.print_statistics_lattice()
            self.timing.summary()

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
        project_root = Path(__file__).resolve().parents[2]
        path = project_root / "data" / "outputs" / "saved_lattice_file" / file_name
        if path.suffix != ".pkl":
            path = path.with_suffix('.pkl')

        if not os.path.exists(path):
            raise FileNotFoundError(f"The file {path} does not exist.")

        try:
            with open(path, "rb") as file:
                lattice = pickle.load(file)
        except Exception as e:
            raise RuntimeError(f"Failed to load lattice from {path}: {e}")

        print(f"Lattice loaded successfully from {path}")
        return lattice

    def __repr__(self) -> str:
        string = f"Lattice name_lattice: {self.name_lattice}\n"
        string += f"Dimensions: {self.size_x} x {self.size_y} x {self.size_z}\n"
        string += f"Number of cells: {self.num_cells_x} x {self.num_cells_y} x {self.num_cells_z}\n"
        string += f"Cell size: {self.cell_size_x} x {self.cell_size_y} x {self.cell_size_z}\n"
        string += f"radii: {self.radii}\n"
        return string

    def extract_parameters_from_json(self, name_file: str):
        """
        Extract lattice parameters from a JSON file and set the corresponding attributes.

        Parameters:
        -----------
        name_file: str
            Name of the JSON file containing lattice parameters.
        """
        lattice_parameters = open_lattice_parameters(name_file)

        # Geometry
        geometry = lattice_parameters.get("geometry", {})
        cell_size = geometry.get("cell_size", {})
        number_of_cells = geometry.get("number_of_cells", {})

        self.cell_size_x = cell_size.get("x")
        self.cell_size_y = cell_size.get("y")
        self.cell_size_z = cell_size.get("z")
        self.num_cells_x = number_of_cells.get("x")
        self.num_cells_y = number_of_cells.get("y")
        self.num_cells_z = number_of_cells.get("z")
        self.radii = geometry.get("radii")
        self.geom_types = geometry.get("geom_types")

        if None in [self.cell_size_x, self.cell_size_y, self.cell_size_z,
                    self.num_cells_x, self.num_cells_y, self.num_cells_z,
                    self.radii, self.geom_types]:
            raise ValueError("Missing geometry parameters in JSON file.")

        # Gradient
        gradient = lattice_parameters.get("gradient", {})
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
        supplementary = lattice_parameters.get("suplementary", {})
        self.uncertainty_node = supplementary.get("node_uncertainty", 0.0)

        # Erased blocks
        erased_blocks_json = supplementary.get("erased_blocks", {})
        self.eraser_blocks = []
        for block in erased_blocks_json.values():
            start = block.get("start_point", {})
            dim = block.get("dimensions_block", {})
            self.eraser_blocks.append([
                start.get("x", 0.0), start.get("y", 0.0), start.get("z", 0.0),
                dim.get("x", 0.0), dim.get("y", 0.0), dim.get("z", 0.0)
            ])

        if len(self.eraser_blocks) == 0:
            self.eraser_blocks = None

        # Symmetry
        symmetries = supplementary.get("symmetries", {})
        if symmetries:
            sym_plane = symmetries.get("plane", None)
            sym_point = symmetries.get("reference_point", {})
            self.symmetry_lattice = {
                "sym_plane": sym_plane,
                "sym_point": (sym_point.get("x", 0.0),
                              sym_point.get("y", 0.0),
                              sym_point.get("z", 0.0))
            }

        _validate_inputs_lattice(self.cell_size_x, self.cell_size_y, self.cell_size_z, self.num_cells_x,
                                 self.num_cells_y, self.num_cells_z, self.geom_types, self.radii,
                                 grad_radius_property, grad_dim_property, grad_mat_property,
                                 self.uncertainty_node, self.eraser_blocks)

        self.define_gradient(grad_radius_property, grad_dim_property, grad_mat_property)


    def define_gradient(self, grad_radius_property, grad_dim_property, grad_mat_property):
        """
        Define gradient settings for radii, dimensions, and materials based on provided properties.
        """
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
            self.grad_mat = grad_settings_constant(self.num_cells_x, self.num_cells_y, self.num_cells_z, True)

    @timing.timeit
    def generate_lattice(self):
        """
        Generate cells in the lattice structure based on cell size, number of cells, geometry types, and radii.
        Gradient information and erased regions are also considered during cell generation.
        """
        nx, ny, nz = self.num_cells_x, self.num_cells_y, self.num_cells_z
        csx, csy, csz = self.cell_size_x, self.cell_size_y, self.cell_size_z
        grad = self.grad_dim
        initial_cell_size = [csx, csy, csz]
        added_nodes = {}
        added_beams = {}

        # Precompute cumulative start positions along each axis
        x_starts = [0.0] * nx
        y_starts = [0.0] * ny
        z_starts = [0.0] * nz
        for i in range(1, nx):
            x_starts[i] = x_starts[i - 1] + csx * grad[i - 1][0]
        for j in range(1, ny):
            y_starts[j] = y_starts[j - 1] + csy * grad[j - 1][1]
        for k in range(1, nz):
            z_starts[k] = z_starts[k - 1] + csz * grad[k - 1][2]

        append_cell = self.cells.append
        mesh_trimmer = self.mesh_trimmer

        for i in range(nx):
            xi = x_starts[i]
            for j in range(ny):
                yj = y_starts[j]
                for k in range(nz):
                    zk = z_starts[k]
                    start_cell_pos = [xi, yj, zk]
                    if self.is_not_in_erased_region(start_cell_pos):
                        continue  # skip erased regions

                    pos_cell = [i, j, k]
                    new_cell = Cell(
                        pos_cell, initial_cell_size, start_cell_pos,
                        self.geom_types, self.radii,
                        self.grad_radius, self.grad_dim, self.grad_mat,
                        self.uncertainty_node, self._verbose,
                        beams_already_defined = added_beams,
                        nodes_already_defined = added_nodes
                    )

                    if mesh_trimmer and not mesh_trimmer.is_cell_in_mesh(new_cell):
                        continue
                    append_cell(new_cell)
                    self.beams.update(new_cell.beams_cell)
                    self.nodes.update(new_cell.points_cell)

        if len(self.geom_types) > 1:
            self.check_hybrid_collision()

    def set_tag_classification(self) -> None:
        """
        Define tag list classification.
        """
        self.edge_tags = [[102, 104, 106, 107], [100, 108, 105, 111], [101, 109, 103, 110]]
        self.face_tags = [[10, 15], [11, 14], [12, 13]]
        self.corner_tags = [[1000, 1001, 1002, 1003, 1004, 1005, 1006, 1007]]

    def define_size_lattice(self):
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

        self.size_x, self.size_y, self.size_z = size_lattice

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
        # Beams
        existing_beam_idx = [b.index for b in self.beams if getattr(b, "index", None) is not None]
        next_beam_idx = (max(existing_beam_idx) + 1) if existing_beam_idx else 0

        def beam_key(b):
            """ Key for sorting beams by their endpoints and radius. """
            p1, p2 = b.point1, b.point2
            return (min((p1.x, p1.y, p1.z), (p2.x, p2.y, p2.z)),
                    max((p1.x, p1.y, p1.z), (p2.x, p2.y, p2.z)),
                    getattr(b, "radius", 0.0))

        for b in sorted(self.beams, key=beam_key):
            if getattr(b, "index", None) is None:
                b.index = next_beam_idx
                next_beam_idx += 1

        # Nodes
        existing_node_idx = [n.index for n in self.nodes if getattr(n, "index", None) is not None]
        next_node_idx = (max(existing_node_idx) + 1) if existing_node_idx else 0

        def node_key(n):
            """ Key for sorting nodes by their coordinates. """
            return n.x, n.y, n.z

        for n in sorted(self.nodes, key=node_key):
            if getattr(n, "index", None) is None:
                n.index = next_node_idx
                next_node_idx += 1

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

    def define_cell_neighbours(self) -> None:
        """
        Neighbour assignment using integer grid indices + optional periodic wrap.
        """
        # Build integer grid indices from positions
        xmin, xmax, ymin, ymax, zmin, zmax = self.get_lattice_boundary_box()
        inv_dx = 1.0 / self.cell_size_x
        inv_dy = 1.0 / self.cell_size_y
        inv_dz = 1.0 / self.cell_size_z

        def to_idx(p):
            """
            Convert position to integer grid index.
            """
            # round to be robust to tiny float noise
            return (
                int(round((p[0] - xmin) * inv_dx)),
                int(round((p[1] - ymin) * inv_dy)),
                int(round((p[2] - zmin) * inv_dz)),
            )

        idx_to_cell = {}
        for c in self.cells:
            idx_to_cell[to_idx(c.pos)] = c

        # Grid extents for wrapping / bounds
        xs = sorted({i for (i, _, _) in idx_to_cell})
        ys = sorted({j for (_, j, _) in idx_to_cell})
        zs = sorted({k for (_, _, k) in idx_to_cell})
        ix0, nx = xs[0], len(xs)
        iy0, ny = ys[0], len(ys)
        iz0, nz = zs[0], len(zs)

        # Offsets and their (axis, sign) labels
        neighbor_steps = [
            ((-1, 0, 0), ("x", "negatif")),
            ((+1, 0, 0), ("x", "positif")),
            ((0, -1, 0), ("y", "negatif")),
            ((0, +1, 0), ("y", "positif")),
            ((0, 0, -1), ("z", "negatif")),
            ((0, 0, +1), ("z", "positif")),
        ]

        # Assign neighbours
        for c in self.cells:
            # reset dict container
            c.neighbour_cells = {}

            ix, iy, iz = to_idx(c.pos)
            for (dx, dy, dz), (axis, sign) in neighbor_steps:
                nix, niy, niz = ix + dx, iy + dy, iz + dz

                if self.enable_periodicity:
                    # periodic wrap (relative to min index of each axis)
                    if nx > 0:
                        nix = ix0 + ((nix - ix0) % nx)
                    if ny > 0:
                        niy = iy0 + ((niy - iy0) % ny)
                    if nz > 0:
                        niz = iz0 + ((niz - iz0) % nz)
                else:
                    # hard bounds (skip if outside)
                    if not (ix0 <= nix < ix0 + nx and iy0 <= niy < iy0 + ny and iz0 <= niz < iz0 + nz):
                        continue

                neigh = idx_to_cell.get((nix, niy, niz))
                if neigh is not None:
                    c.add_cell_neighbour(axis, sign, neigh)

    @timing.timeit
    def define_connected_beams_for_all_nodes(self, merge_tol: float = 1e-9) -> None:
        """
        Populate `Point.connected_beams` for all nodes.
        Optionally merges connectivity for coincident/periodic nodes via a spatial hash.

        Parameters
        ----------
        merge_tol : float
            Quantization used to bucket nearly-identical coordinates (and periodic wraps).
        """
        # Reset connected_beams
        for p in self.nodes:
            p.connected_beams.clear()

        # Single pass: attach each beam to its two endpoints
        for b in self.beams:
            a, c = b.point1, b.point2
            if b not in a.connected_beams: a.connected_beams.append(b)
            if b not in c.connected_beams: c.connected_beams.append(b)

        #  Optional stitching of coincident & periodic nodes via spatial hashing
        # Build spatial buckets
        xmin, xmax, ymin, ymax, zmin, zmax = self.get_lattice_boundary_box()

        def qkey(x, y, z):
            """Quantized key for spatial hashing."""
            return round(x / merge_tol), round(y / merge_tol), round(z / merge_tol)

        buckets = {}  # geometric buckets (exact/coincident)
        periodic_buckets = {}  # buckets after wrapping (only if periodicity enabled)

        def add_to_bucket(d, p):
            """ Add point to geometric bucket."""
            k = qkey(p.x, p.y, p.z)
            d.setdefault(k, []).append(p)

        def add_to_periodic_bucket(d, p):
            """ Add point to periodic bucket (wrap max-face to min-face)."""
            wx = xmin if abs(p.x - xmax) <= merge_tol else p.x
            wy = ymin if abs(p.y - ymax) <= merge_tol else p.y
            wz = zmin if abs(p.z - zmax) <= merge_tol else p.z
            k = qkey(wx, wy, wz)
            d.setdefault(k, []).append(p)

        # Fill buckets
        for p in self.nodes:
            add_to_bucket(buckets, p)
            if getattr(self, "enable_periodicity", False):
                add_to_periodic_bucket(periodic_buckets, p)


        def merge_bucket_connectivity(b):
            """ Merge connectivity for all points in each bucket."""
            for pts in b.values():
                if len(pts) <= 1:
                    continue
                # union of beams connected to all points in the bucket
                merged = []
                seen = set()
                for pt in pts:
                    for bm in pt.connected_beams:
                        if id(bm) not in seen:
                            seen.add(id(bm))
                            merged.append(bm)
                # assign back
                for pt in pts:
                    pt.connected_beams = merged.copy()

        merge_bucket_connectivity(buckets)
        if getattr(self, "enable_periodicity", False):
            merge_bucket_connectivity(periodic_buckets)

    @timing.timeit
    def define_angles_between_beams(self) -> None:
        """
        Compute, for each beam, the penalization-optimal angle at its two endpoints
        using node-level connectivity (Point.connected_beams). Assumes
        `define_connected_beams_for_all_nodes()` has been called.
        """
        @timing.timeit
        def _best_by_penalization(angles: list[float], radii: list[float]) -> tuple[float, float]:
            """Return (angle, radius) that maximizes L = function_penalization_Lzone((radius, angle))."""
            if not angles:
                return 0.0, 0.0
            Lmax, best = -1.0, (0.0, 0.0)
            for r, a in zip(radii, angles):
                L = function_penalization_Lzone(r, a)
                if L > Lmax:
                    Lmax, best = L, (a, r)
            return best

        # Compute angles using local node neighborhoods only
        for b in self.beams:
            for p in (b.point1, b.point2):
                a1_list, r1_list = [], []
                for nb in p.connected_beams:
                    if nb is b: continue
                    try:
                        ang = b.get_angle_between_beams(nb, self.enable_periodicity)
                    except ValueError:
                        continue
                    if ang > 1e-12:
                        a1_list.append(ang);
                        r1_list.append(nb.radius)
                ang1, rad1 = _best_by_penalization(a1_list, r1_list)
                b.set_angle(rad1, ang1, p)

    @timing.timeit
    def define_lattice_dimensions(self) -> dict:
        """
        Computes extremum values of coordinates in the lattice.

        Return:
        --------
        lattice_dimension_dict: dict
            Dictionary containing min and max values for x, y, z coordinates.
        """
        if not self.nodes:
            raise ValueError("No nodes in the lattice.")
        xs = [n.x for n in self.nodes];
        ys = [n.y for n in self.nodes];
        zs = [n.z for n in self.nodes]
        self.x_min, self.x_max = min(xs), max(xs)
        self.y_min, self.y_max = min(ys), max(ys)
        self.z_min, self.z_max = min(zs), max(zs)
        self.lattice_dimension_dict = {
            "x_min": self.x_min, "x_max": self.x_max,
            "y_min": self.y_min, "y_max": self.y_max,
            "z_min": self.z_min, "z_max": self.z_max,
        }
        return self.lattice_dimension_dict

    def get_lattice_boundary_box(self) -> list[float]:
        """
        Get the boundary box of the lattice
        """
        return [self.x_min, self.x_max, self.y_min, self.y_max, self.z_min, self.z_max]

    @timing.timeit
    def set_penalized_beams(self) -> None:
        """
        Modifies beam and node data to model lattice structures for simulation with rigidity penalization at node.
        If L_zone at one end is 0 (or <= 0), we do not create a penalized segment at that end.
        If both are 0, the beam is left unchanged.
        """
        for cell in self.cells:
            beams_to_remove = []
            beams_to_add = []
            points_to_add = []

            for beam in list(cell.beams_cell):
                L1 = beam.angle_point_1.get("L_zone", 0)
                L2 = beam.angle_point_2.get("L_zone", 0)

                # No modification if both are zero or negative
                if L1 <= 0 and L2 <= 0:
                    continue

                # Start from original endpoints
                start = beam.point1
                end = beam.point2

                # Left/end-1 modification
                if L1 > 0:
                    pointExt1 = beam.get_point_on_beam_at_distance(L1, 1)
                    pointExt1.node_mod = True
                    points_to_add.append(pointExt1)
                    b1 = Beam(start, pointExt1, beam.radius, beam.material, beam.type_beam)
                    b1.set_beam_mod()
                    beams_to_add.append(b1)
                    start = pointExt1  # middle starts here

                # Right/end-2 modification
                if L2 > 0:
                    pointExt2 = beam.get_point_on_beam_at_distance(L2, 2)
                    pointExt2.node_mod = True
                    points_to_add.append(pointExt2)
                    mid_end = pointExt2
                else:
                    mid_end = end

                # Middle (unmodified) segment
                b_mid = Beam(start, mid_end, beam.radius, beam.material, beam.type_beam)
                beams_to_add.append(b_mid)

                # Final penalized segment at end-2 if needed
                if L2 > 0:
                    b3 = Beam(mid_end, end, beam.radius, beam.material, beam.type_beam)
                    b3.set_beam_mod()
                    beams_to_add.append(b3)

                beams_to_remove.append(beam)

            cell.add_beam(beams_to_add)
            cell.add_point(points_to_add)
            cell.remove_beam(beams_to_remove)

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
            for beam in cell.beams_cell:
                if minLength > beam.length > 0.0001 and (beam.type_beam == 0 or beam.type_beam == 2):
                    minLength = beam.length
        return minLength

    def get_tag_list(self) -> list[int]:
        """Get the tag for all unique points in the lattice."""
        return [n.tag for n in self.nodes]

    def get_tag_list_boundary(self) -> list[int]:
        """Get the tag for boundary points in the lattice."""
        bbox = self.get_lattice_boundary_box()
        return [n.tag for n in self.nodes if n.is_on_boundary(bbox)]

    @timing.timeit
    def apply_tag_all_point(self) -> None:
        """
        Assign a tag to all nodes in the lattice structure.
        Tags are assigned relative to either the global bounding box
        or a local (cell-relative) bounding box if erased parts are used.
        """
        use_local_box = self.eraser_blocks is not None
        if not use_local_box:
            bbox = self.get_lattice_boundary_box()
            for n in self.nodes:
                n.tag = n.tag_point(bbox)
            return

        if self.occupancy_matrix is None:
            self.get_cell_occupancy_matrix()
        for cell in self.cells:
            local_box = self.get_relative_boundary_box(cell)
            for n in cell.nodes_cell:
                n.tag = n.tag_point(local_box)

    def get_cell_occupancy_matrix(self):
        """
        Build a 3D matrix storing Cell objects (or None) at each (i, j, k).

        Returns
        -------
        occupancy_matrix : np.ndarray, shape (num_cells_x, num_cells_y, num_cells_z), dtype=object
            occupancy_matrix[i, j, k] is the Cell at that grid position, or None if empty.
        """
        nx, ny, nz = self.num_cells_x, self.num_cells_y, self.num_cells_z
        self.occupancy_matrix = np.empty((nx, ny, nz), dtype=object)
        self.occupancy_matrix[:] = None
        for cell in self.cells:
            i, j, k = cell.pos  # integer grid indices
            self.occupancy_matrix[i, j, k] = cell
        return self.occupancy_matrix

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
        x, y, z = cell.pos
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
                else:
                    raise ValueError("Axis must be 'x', 'y', or 'z'")

                if bounds[0] < min_val:
                    min_val = bounds[0]
                if bounds[1] > max_val:
                    max_val = bounds[1]

            bbox.append(min_val)
            bbox.append(max_val)

        return bbox

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
            for beam in cell.beams_cell:
                if beam.point1.index == nodeIndexRef:
                    connectedNode.append(beam.point2)
                if beam.point2.index == nodeIndexRef:
                    connectedNode.append(beam.point1)
        return connectedNode

    def find_boundary_beams(self) -> list["Beam"]:
        """Find boundary beams and change the type_beam of beam"""
        bbox = self.get_lattice_boundary_box()
        boundary = [b for b in self.beams if b.point1.is_on_boundary(bbox) or b.point2.is_on_boundary(bbox)]
        for b in boundary:
            b.type_beam = 2
        return boundary

    def find_boundary_nodes(self) -> list["Point"]:
        """Find boundary nodes"""
        return [n for n in self.nodes if n.x in (self.x_min, self.x_max)
                or n.y in (self.y_min, self.y_max)
                or n.z in (self.z_min, self.z_max)]

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
                for beam in cell.beams_cell:
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
        return [[n.coordinates] for n in self.nodes]

    def get_edge_index(self) -> list[list[int]]:
        """
        Retrieves edge index data for the lattice.

        Returns:
        --------
        edgeIndex: list of list of int
            List of edge index
        """
        return [[b.point1.index, b.point2.index] for b in self.beams]

    def get_all_beam_type(self) -> list:
        """
        Retrieves beam type_beam data for the lattice.

        Returns:
        --------
        beamType: list of int
            List of beam types
        """
        return [[b.type_beam] for b in self.beams]

    def get_all_beam_length(self) -> list[list[float]]:
        """
        Retrieves beam length data for the lattice.

        Returns:
        --------
        beamLength: list of float
            List of beam lengths
        """
        return [[b.length] for b in self.beams]

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
            for beam in cell.beams_cell:
                if beam.beam_mod:
                    beam.radius = hybridRadiusData[beam.type_beam] * beam.penalization_coefficient
                else:
                    beam.radius = hybridRadiusData[beam.type_beam]

    def get_relative_density(self) -> float:
        """
        Get mean relative density of all cells in lattice

        Returns:
        --------
        meanRelDens: float
            Mean relative density of the lattice
        """
        cellRelDens = []
        for cell in self.cells:
            cellRelDens.append(cell.relative_density)
        meanRelDens = mean(cellRelDens)
        return meanRelDens

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
            for beam in cell.beams_cell:
                if beam.type_beam == typeToChange:
                    beam.radius = newRadius

    def get_number_cells(self) -> int:
        """Get number of cells in the lattice"""
        return len(self.cells)

    def get_number_beams(self) -> int:
        """Get number of beams in the lattice"""
        return len(self.beams)

    def get_number_nodes(self) -> int:
        """Get number of nodes in the lattice"""
        return len(self.nodes)

    def find_point_on_lattice_surface(self, surfaceNames: list[str], surface_cells: list[str] = None) -> set["Point"]:
        """
        Find points on the surface of the lattice

        Parameters:
        -----------
        surfaceNames: list[str]
            List of surfaces to find points on (e.g., ["Xmin", "Xmax", "Ymin"])
        surface_cells: list[str], optional
            List of surfaces to find points on cells (e.g., ["Xmin", "Xmax", "Ymin"]). If None, use surfaceNames.

        Returns:
        --------
        pointSet: set of point objects
            Set of points found on the specified surfaces
        """
        valid_surfaces = {"Xmin", "Xmax", "Ymin", "Ymax", "Zmin", "Zmax", "Xmid", "Ymid", "Zmid"}

        if not all(surface in valid_surfaces for surface in surfaceNames):
            raise ValueError("Invalid surface name_lattice(s).")

        cell_list = self.get_cells_on_surfaces(surfaceNames)
        cell_indices = {c.index for c in cell_list}

        if self.cells[-1].index < max(cell_indices, default=-1):
            raise ValueError("Invalid cell index, some cells do not exist.")

        surface_on_cells = surface_cells if surface_cells is not None else surfaceNames
        pointSet = set()
        for cell in self.cells:
            if cell.index in cell_indices:
                cellPointSets = [set(cell.get_point_on_surface(surface)) for surface in surface_on_cells]
                if cellPointSets:
                    pointSet.update(set.intersection(*cellPointSets))

        if pointSet == set():
            raise ValueError("No points found on the specified surfaces.")

        return pointSet

    def get_cells_on_surfaces(self, surfaces: list[str]) -> list:
        """
        Return the list of Cell objects matching ordered extrema constraints like ["Xmin"], ["Xmin","Zmax"], etc.
        Filtering is iterative: e.g., first keep all cells at X minimum, then among those keep Z maximum.

        Parameters
        ----------
        surfaces : list[str]
            Each item is one of {"Xmin","Xmax","Ymin","Ymax","Zmin","Zmax"} (case-insensitive).

        Returns
        -------
        list
            List of Cell objects (possibly multiple if several share the same extreme index).
        """
        # Ensure occupancy matrix is available
        if getattr(self, "occupancy_matrix", None) is None or self.occupancy_matrix.size == 0:
            self.get_cell_occupancy_matrix()

        # Start from all occupied positions (use existing cells to avoid None handling)
        candidates = [tuple(cell.pos) for cell in self.cells]
        if not candidates:
            return []

        axis_map = {"X": 0, "Y": 1, "Z": 2}

        for token in surfaces:
            t = token.strip().lower()
            if not t:
                continue
            ax_char = t[0].upper()
            if ax_char not in axis_map:
                raise ValueError(f"Invalid axis in constraint '{token}', expected X/Y/Z with min/max.")
            ax = axis_map[ax_char]

            if "min" in t:
                extreme = min(idx[ax] for idx in candidates)
            elif "max" in t:
                extreme = max(idx[ax] for idx in candidates)
            else:
                raise ValueError(f"Invalid extrema in constraint '{token}', expected 'min' or 'max'.")

            candidates = [idx for idx in candidates if idx[ax] == extreme]
            if not candidates:
                return []

        # Map positions back to Cell objects via the occupancy matrix (ignore any None just in case)
        occ = self.occupancy_matrix
        return [occ[i, j, k] for (i, j, k) in candidates if occ[i, j, k] is not None]


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
        print("Lattice dimensions: ", self.size_x, self.size_y, self.size_z)

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
            for beam in cell.beams_cell:
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
    def apply_symmetry(self) -> None:
        """
        Apply symmetry to the lattice structure based on a reference point.
        """
        symmetry_plane, reference_point = self.symmetry_lattice["sym_plane"], self.symmetry_lattice["sym_point"]
        symmetry_plane = symmetry_plane.upper()
        if symmetry_plane not in ["XY", "XZ", "YZ", "X", "Y", "Z"]:
            raise ValueError("Invalid symmetry plane. Choose from 'XY', 'XZ', 'YZ', 'X', 'Y', or 'Z'.")

        x_ref, y_ref, z_ref = reference_point
        new_cells = []
        node_map = {}

        for cell in self.cells:
            new_pos = list(cell.pos)
            new_start_pos = np.array(cell.coordinate)
            mirrored_beams = []

            for beam in cell.beams_cell:
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
            new_cell = Cell(new_pos, cell.size, list(new_start_pos), cell.geom_types, cell.radii, cell.grad_radius,
                            cell.grad_dim, cell.grad_mat, cell.uncertainty_node, self._verbose)

            new_cell.beams_cell = mirrored_beams
            new_cells.append(new_cell)

        self.cells.extend(new_cells)
        self.define_lattice_dimensions()  # Recalculate the lattice boundaries

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
            for beam in cell.beams_cell:
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
            for beam in cell.beams_cell:
                if not geomScheme[beam.type_beam]:
                    beamsToRemove.append(beam)
            for beam in beamsToRemove:
                cell.remove_beam(beam)

    @timing.timeit
    def generate_mesh_lattice_Gmsh(self, cut_mesh_at_boundary: bool = False, mesh_size: float = 0.05,
                                   name_mesh: str = "Lattice",save_mesh: bool = False, save_STL: bool = True,
                                   volume_computation: bool = False):
        """
        Generate a 2D mesh representation of the lattice structure using GMSH.
        Generating 3D mesh for simulation is not currently supported, but will be in the future.

        Parameters:
        -----------
        cut_mesh_at_boundary: bool
            If True, the mesh will be cut at the boundary of the lattice.
        mesh_size: float
            Size of the mesh elements.
        name_mesh: str
            Name of the mesh to be generated.
        save_mesh: bool
            If True, the mesh will be saved to a file.
        save_STL: bool
            If True, the mesh will be saved in STL format.
        volume_computation: bool
            If True, the volume of the mesh will be computed and printed.
        """
        if self.enable_simulation_properties == 1:
            raise ValueError("mesh_file generation is not available for the current simulation method.")

        gmsh.initialize()
        gmsh.option.setNumber("General.Verbosity", 1)
        gmsh.model.add(name_mesh)
        dim = 3  # Dimension of the mesh

        all_tags = []
        for beam in self.beams:
            p1 = np.array(beam.point1.coordinates)
            p2 = np.array(beam.point2.coordinates)
            direction = p2 - p1
            if beam.index not in all_tags:
                tag = gmsh.model.occ.addCylinder(*p1, *direction, beam.radius, tag=beam.index)
                all_tags.append(tag)

        # Merge all beams into a single entity
        beam_entities = [(dim, tag) for tag in all_tags]
        lattice = gmsh.model.occ.fragment(beam_entities, [])

        if cut_mesh_at_boundary:
            # Bounding box definition
            x0, y0, z0 = self.x_min, self.y_min, self.z_min
            dx, dy, dz = self.x_max - x0, self.y_max - y0, self.z_max - z0
            box = gmsh.model.occ.addBox(x0, y0, z0, dx, dy, dz)
            gmsh.model.occ.intersect(lattice[0], [(dim, box)], removeObject=True)

        gmsh.model.occ.synchronize()

        # Define mesh size
        points = gmsh.model.getEntities(dim=0)
        gmsh.model.mesh.setSize(points, mesh_size)
        gmsh.model.occ.synchronize()

        if volume_computation:
            self.get_volume_mesh()

        gmsh.model.mesh.generate(2)

        project_root = Path(__file__).resolve().parents[2]
        if save_mesh:
            path = project_root / "data" / "outputs" / "mesh_file" / f"{name_mesh}.msh"
            gmsh.write(str(path))
            print("mesh_file saved at ", path)

        if save_STL:
            path = project_root / "mesh_file" / f"{name_mesh}.stl"
            gmsh.write(str(path))
            print("mesh_file saved at ", path)

        gmsh.finalize()

    @timing.timeit
    def get_volume_mesh(self) -> float:
        """
        Compute the total volume of the mesh and the relative density of the lattice.

        Returns:
        --------
        vol: float
            Total volume of the mesh.
        """
        vol = 0.0
        for dim, tag in gmsh.model.getEntities(dim=3):
            vol += gmsh.model.occ.getMass(dim, tag)
        if self._verbose > 0:
            print(f"Total volume: {vol}")
        return vol

    def get_relative_density_mesh(self) -> float:
        """
        Compute the relative density of the lattice based on the mesh volume.

        Returns:
        --------
        relative_density: float
            Relative density of the lattice.
        """
        bb = self.get_lattice_boundary_box()
        volume_lattice = bb[1] - bb[0] * bb[3] - bb[2] * bb[5] - bb[4]
        relative_density = self.get_volume_mesh() / volume_lattice
        if self._verbose > 0:
            print(f"Relative density: {relative_density:.4f}")
        return relative_density

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
            "size",
            "grad_radius",
            "grad_dim",
        ]

        for i, cell in enumerate(self.cells[1:], start=1):
            for attr in attrs_to_check:
                if not np.array_equal(getattr(reference, attr), getattr(cell, attr)):
                    print(
                        Fore.RED + f"Difference found in attribute '{attr}' between cell 0 and cell {i}" + Style.RESET_ALL)
                    return False

            if len(reference.beams_cell) != len(cell.beams_cell):
                print(Fore.RED + f"Different number of beams between cell 0 and cell {i}" + Style.RESET_ALL)
                return False

            for j, (b1, b2) in enumerate(zip(sorted(reference.beams_cell, key=id),
                                             sorted(cell.beams_cell, key=id))):
                if not b1.is_identical_to(b2, reference.size):
                    print(b1, b2)
                    print(Fore.RED + f"Difference found in beam {j} between cell 0 and cell {i}" + Style.RESET_ALL)
                    return False

        print(Fore.GREEN + "All cells are identical." + Style.RESET_ALL)
        return True

