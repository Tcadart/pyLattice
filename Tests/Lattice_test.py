from src.lattice import *


def test_simple_lattice():
    latticeTest = Lattice.simpleLattice(2, 2, 2,
                                        2, 2, 2, 0, 0.05)
    assert len(latticeTest.cells) == 8, "Lattice should have 8 cells"
    assert latticeTest.x_min == 0, "x_min should be 0"
    assert latticeTest.x_max == 4, "x_max should be 4"
    assert latticeTest.y_min == 0, "y_min should be 0"
    assert latticeTest.y_max == 4, "y_max should be 4"
    assert latticeTest.z_min == 0, "z_min should be 0"
    assert latticeTest.z_max == 4, "z_max should be 4"


def test_create_simple_lattice():
    lattice = Lattice.simpleLattice(1, 1, 1, 2, 2, 2, 0, 0.05)
    assert len(lattice.cells) == 8
    assert lattice.getNumberOfBeams() > 0
    assert lattice.get_number_nodes() > 0


def test_lattice_dimensions_and_bounds():
    lattice = Lattice.simpleLattice(1, 1, 1, 3, 1, 1, 0, 0.05)
    assert lattice.size_x == 3.0
    assert lattice.size_y == 1.0
    assert lattice.size_z == 1.0
    assert lattice.x_min == 0.0
    assert lattice.x_max == 3.0


def test_beam_and_node_counts():
    lattice = Lattice.simpleLattice(1, 1, 1, 2, 2, 2, 0, 0.05)
    num_beams = lattice.getNumberOfBeams()
    num_nodes = lattice.get_number_nodes()
    assert num_beams > 0
    assert num_nodes > 0


def test_relative_density():
    lattice = Lattice.simpleLattice(1, 1, 1, 2, 2, 2, 0, 0.05)
    rd = lattice.get_relative_density()
    assert 0 < rd < 1


def test_generate_mesh():
    lattice = Lattice.simpleLattice(1, 1, 1, 1, 1, 1, 0, 0.05)
    lattice.generateMeshLattice()
    assert lattice.mesh_lattice is not None
    assert lattice.mesh_lattice.faces.shape[0] > 0


def test_delete_small_beams():
    lattice = Lattice.simpleLattice(1, 1, 1, 2, 2, 2, 0, 0.0001)
    initial_beam_count = lattice.getNumberOfBeams()
    lattice.delete_beams_under_radius_threshold(threshold=0.001)
    final_beam_count = lattice.getNumberOfBeams()
    assert final_beam_count < initial_beam_count


def test_apply_symmetry():
    lattice = Lattice.simpleLattice(1, 1, 1, 1, 1, 1, 0, 0.05)
    initial_num_cells = len(lattice.cells)
    lattice.apply_symmetry("XY", (0, 0, 0))
    assert len(lattice.cells) == initial_num_cells * 2


def test_erased_cells():
    erased_zone = [[0.0, 0.0, 0.0, 0.0, 2.0, 2.0]]
    lattice = Lattice(1, 1, 1, 2, 2, 2, [0], [0.05], "VeroClear", ["constant", [0, 0, 0], [0, 0, 0]],
                      ["constant", [0, 0, 0], [0, 0, 0]], [0, 1], eraser_blocks=erased_zone)
    assert len(lattice.cells) == 4
