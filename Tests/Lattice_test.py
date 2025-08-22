# from src.pyLattice.lattice import Lattice
#
#
# def test_simple_lattice():
#     latticeTest = Lattice(2, 2, 2,
#                                   2, 2, 2, ["BCC"], [0.05])
#     assert len(latticeTest.cells) == 8, "Lattice should have 8 cells"
#     assert latticeTest.x_min == 0, "x_min should be 0"
#     assert latticeTest.x_max == 4, "x_max should be 4"
#     assert latticeTest.y_min == 0, "y_min should be 0"
#     assert latticeTest.y_max == 4, "y_max should be 4"
#     assert latticeTest.z_min == 0, "z_min should be 0"
#     assert latticeTest.z_max == 4, "z_max should be 4"
#
#
# def test_create_simple_lattice():
#     latticeTest = Lattice(1, 1, 1,
#                                   2, 2, 2, ["BCC"], [0.05])
#     assert len(latticeTest.cells) == 8
#     assert latticeTest.get_number_beams() > 0
#     assert latticeTest.get_number_nodes() > 0
#
#
# def test_lattice_dimensions_and_bounds():
#     latticeTest = Lattice(1, 1, 1,
#                                   3, 1, 1, ["BCC"], [0.05])
#     assert latticeTest.size_x == 3.0
#     assert latticeTest.size_y == 1.0
#     assert latticeTest.size_z == 1.0
#     assert latticeTest.x_min == 0.0
#     assert latticeTest.x_max == 3.0
#
#
# def test_beam_and_node_counts():
#     latticeTest = Lattice(1, 1, 1,
#                                   2, 2, 2, ["BCC"], [0.05])
#     num_beams = latticeTest.get_number_beams()
#     num_nodes = latticeTest.get_number_nodes()
#     assert num_beams > 0
#     assert num_nodes > 0
#
#
# def test_relative_density():
#     latticeTest = Lattice(1, 1, 1,
#                                   2, 2, 2, ["BCC"], [0.05])
#     rd = latticeTest.get_relative_density()
#     assert 0 < rd < 1
#
#
# def test_generate_mesh():
#     latticeTest = Lattice(1, 1, 1,
#                                   1, 1, 1, ["BCC"], [0.05])
#     latticeTest.generate_mesh_lattice_Gmsh()
#     assert latticeTest.mesh_lattice is not None
#     assert latticeTest.mesh_lattice.faces.shape[0] > 0
#
#
# def test_delete_small_beams():
#     latticeTest = Lattice(1, 1, 1,
#                                   2, 2, 2, ["BCC"], [0.00001])
#     initial_beam_count = latticeTest.get_number_beams()
#     latticeTest.delete_beams_under_radius_threshold(threshold=0.001)
#     final_beam_count = latticeTest.get_number_beams()
#     assert final_beam_count < initial_beam_count
#
#
# def test_apply_symmetry():
#     latticeTest = Lattice(1, 1, 1,
#                                   1, 1, 1, ["BCC"], [0.05])
#     initial_num_cells = len(latticeTest.cells)
#     latticeTest.apply_symmetry("XY", (0, 0, 0))
#     assert len(latticeTest.cells) == initial_num_cells * 2
#
#
# def test_erased_cells():
#     erased_zone = [[0.0, 0.0, 0.0, 0.0, 2.0, 2.0]]
#     lattice = Lattice(1, 1, 1, 2, 2, 2, ["BCC"], [0.05], "VeroClear", ["constant", [0, 0, 0], [0, 0, 0]],
#                       ["constant", [0, 0, 0], [0, 0, 0]], [0, 1], eraser_blocks=erased_zone)
#     assert len(lattice.cells) == 4
