from mpi4py import MPI
from typing import TYPE_CHECKING


from .beam_model import *
from .full_scale_lattice_simulation import *
from .homogenization_cell import HomogenizedCell
from .schur_complement import SchurComplement

if TYPE_CHECKING:
    from pyLatticeSim.lattice_sim import LatticeSim


def solve_FEM_FenicsX(lattice : "LatticeSim"):
    """
    Solve the finite element method problem using FenicsX for a given lattice.

    Parameters:
    -----------
    lattice: Lattice object
        The lattice structure to be simulated.

    Returns:
    --------
    xsol: numpy.ndarray
        The solution vector containing displacements.
    simulationModel: FullScaleLatticeSimulation
        The simulation model containing the results of the simulation.
    """
    # Generate the lattice model and mesh
    LatticeModel = BeamModel(MPI.COMM_SELF, lattice=lattice)

    # Define the FE model and apply boundary conditions
    simulationModel = FullScaleLatticeSimulation(LatticeModel)
    simulationModel.apply_displacement_all_nodes_with_lattice_data()
    simulationModel.apply_force_on_all_nodes_with_lattice_data()
    simulationModel.define_L_form_null()

    # Solve the problem
    simulationModel.solve_problem()

    # Assign results to lattice object
    simulationModel.set_result_diplacement_on_lattice_object()

    # Get results to return
    xsol, globalDisplacementIndex = lattice.get_global_displacement()
    return xsol, simulationModel

def get_homogenized_properties(lattice: "LatticeSim"):
    """
    Perform homogenization analysis on a lattice structure.

    Parameters:
    -----------
    lattice: Lattice object
        The lattice structure to be homogenized.

    Returns:
    --------
    mat_Sorthotropic: numpy.ndarray
        The homogenized orthotropic stiffness matrix.
    homogenization_analysis: HomogenizedCell
        The homogenization analysis object containing results and methods.
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

    mat_Sorthotropic = homogenization_analysis.get_S_orthotropic()

    return mat_Sorthotropic, homogenization_analysis

def get_schur_complement(lattice: "LatticeSim", cell_index: int = None):
    """
    Calculate the Schur complement of the stiffness matrix for a given lattice.

    Parameters:
    -----------
    lattice: Lattice object
        The lattice structure to be analyzed.
    cell_index: int, optional
        The index of the cell to be used for the Schur complement calculation.
        If None, the first cell is used.
    """
    if cell_index is None and lattice.get_number_cells() > 1:
        raise ValueError("The lattice must contain only one cell for Schur complement calculation or specify a cell_index.")

    cell_model = BeamModel(MPI.COMM_SELF, lattice=lattice, cell_index=cell_index)

    # Initialization simulation
    schur_complement_analysis = SchurComplement(cell_model)

    tags_nodes_boundary = lattice.cells[0].get_node_order_to_simulate() if cell_index is None else (
                            lattice.cells[cell_index].get_node_order_to_simulate())

    schur_complement, _ = schur_complement_analysis.calculate_schur_complement(tags_nodes_boundary)

    return schur_complement


