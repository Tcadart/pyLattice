from src.Lattice import *


def test_simple_lattice():
    latticeTest = Lattice.simpleLattice(2, 2, 2,
                                        2, 2, 2, 0, 0.05)
    assert len(latticeTest.cells) == 8, "Lattice should have 8 cells"
    assert latticeTest.xMin == 0, "xMin should be 0"
    assert latticeTest.xMax == 4, "xMax should be 4"
    assert latticeTest.yMin == 0, "yMin should be 0"
    assert latticeTest.yMax == 4, "yMax should be 4"
    assert latticeTest.zMin == 0, "zMin should be 0"
    assert latticeTest.zMax == 4, "zMax should be 4"


def test_create_simple_lattice():
    lattice = Lattice.simpleLattice(1, 1, 1, 2, 2, 2, 0, 0.05)
    assert len(lattice.cells) == 8
    assert lattice.getNumberOfBeams() > 0
    assert lattice.getNumberOfNodes() > 0


def test_lattice_dimensions_and_bounds():
    lattice = Lattice.simpleLattice(1, 1, 1, 3, 1, 1, 0, 0.05)
    assert lattice.size_x == 3.0
    assert lattice.size_y == 1.0
    assert lattice.size_z == 1.0
    assert lattice.xMin == 0.0
    assert lattice.xMax == 3.0


def test_beam_and_node_counts():
    lattice = Lattice.simpleLattice(1, 1, 1, 2, 2, 2, 0, 0.05)
    num_beams = lattice.getNumberOfBeams()
    num_nodes = lattice.getNumberOfNodes()
    assert num_beams > 0
    assert num_nodes > 0


def test_relative_density():
    lattice = Lattice.simpleLattice(1, 1, 1, 2, 2, 2, 0, 0.05)
    rd = lattice.getRelativeDensity()
    assert 0 < rd < 1


def test_generate_mesh():
    lattice = Lattice.simpleLattice(1, 1, 1, 1, 1, 1, 0, 0.05)
    lattice.generateMeshLattice()
    assert lattice.mesh_lattice is not None
    assert lattice.mesh_lattice.faces.shape[0] > 0


def test_delete_small_beams():
    lattice = Lattice.simpleLattice(1, 1, 1, 2, 2, 2, 0, 0.0001)
    initial_beam_count = lattice.getNumberOfBeams()
    lattice.deleteBeamsUnderThreshold(threshold=0.001)
    final_beam_count = lattice.getNumberOfBeams()
    assert final_beam_count < initial_beam_count


def test_apply_symmetry():
    lattice = Lattice.simpleLattice(1, 1, 1, 1, 1, 1, 0, 0.05)
    initial_num_cells = len(lattice.cells)
    lattice.applySymmetry("XY", (0, 0, 0))
    assert len(lattice.cells) == initial_num_cells * 2


def test_erased_cells():
    erased_zone = [[0.0, 0.0, 0.0, 0.0, 2.0, 2.0]]
    lattice = Lattice(1, 1, 1, 2, 2, 2, [0], [0.05], "VeroClear", ["constant", [0, 0, 0], [0, 0, 0]],
                      ["constant", [0, 0, 0], [0, 0, 0]], [0, 1], eraser_blocks=erased_zone)
    assert len(lattice.cells) == 4
