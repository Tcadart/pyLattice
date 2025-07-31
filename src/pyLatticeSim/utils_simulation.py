from mpi4py import MPI


from .beam_model import *
from .full_scale_lattice_simulation import *
from pyLattice.lattice import Lattice
from .homogenization_cell import HomogenizedCell

def solve_FEM_FenicsX(lattice : "Lattice"):
    """
    Solve the finite element method problem using FenicsX for a given lattice.
    """
    LatticeModel = BeamModel(MPI.COMM_SELF, lattice=lattice)

    simulationModel = FullScaleLatticeSimulation(LatticeModel)
    simulationModel.prepare_simulation()
    simulationModel.apply_displacement_all_nodes_with_lattice_data()
    simulationModel.apply_force_on_all_nodes_with_lattice_data()

    simulationModel.define_L_form_null()
    simulationModel.solve_problem()
    assign_displacements_to_lattice_nodes(simulationModel)
    xsol, globalDisplacementIndex = lattice.get_global_displacement()
    return xsol, globalDisplacementIndex


def assign_displacements_to_lattice_nodes(simulation):
    """
    Assigns the displacement and rotation values from the simulation to the lattice nodes.
    """
    # Displacement
    displacement_fem = simulation.u.sub(0).collapse()
    coords_disp = np.round(displacement_fem.function_space.tabulate_dof_coordinates(), 5)
    values_disp = displacement_fem.x.array.reshape((-1, 3))

    # Rotations
    rotation_fem = simulation.u.sub(1).collapse()
    coords_rot = np.round(rotation_fem.function_space.tabulate_dof_coordinates(), 5)
    values_rot = rotation_fem.x.array.reshape((-1, 3))

    # Mapping dictionaries
    pos_to_disp = {tuple(coord): disp for coord, disp in zip(coords_disp, values_disp)}
    pos_to_rot = {tuple(coord): rot for coord, rot in zip(coords_rot, values_rot)}

    # Node assignment
    for cell in simulation.BeamModel.lattice.cells:
        for beam in cell.beams:
            for node in [beam.point1, beam.point2]:
                pos = tuple(np.round([node.x, node.y, node.z], 5))
                if pos in pos_to_disp:
                    node.displacement_vector[:3] = pos_to_disp[pos]
                else:
                    print(f"⚠️ Missing displacement for {pos}")
                if pos in pos_to_rot:
                    node.displacement_vector[3:] = pos_to_rot[pos]
                else:
                    print(f"⚠️ Missing rotation for {pos}")


def get_homogenized_properties(lattice: "Lattice") -> np.ndarray:
    """
    Perform homogenization analysis on a lattice structure.

    Parameters:
    -----------
    lattice: Lattice object
        The lattice structure to be homogenized.
    """
    if lattice.get_number_cells() > 1:
        raise ValueError("The lattice must contain only one cell for homogenization.")

    cell_model = BeamModel(MPI.COMM_SELF, lattice=lattice)

    # Initialization simulation
    homogenization_analysis = HomogenizedCell(cell_model)
    homogenization_analysis.prepare_simulation()

    homogenization_analysis.apply_dirichlet_for_homogenization()
    homogenization_analysis.periodic_boundary_condition()

    homogenization_analysis.solve_full_homogenization()

    homogenization_analysis.print_homogenized_matrix()
    homogenization_analysis.print_errors()

    # exportData = exportBeamData("Result/LatticeHomogenization", homogenization_analysis)
    # exportData.exportDataHomogenization()

    mat_Sorthotropic = homogenization_analysis.get_S_orthotropic()

    return mat_Sorthotropic

