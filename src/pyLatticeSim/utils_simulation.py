from mpi4py import MPI


from .beam_model import *
from .full_scale_lattice_simulation import *
from pyLattice.lattice import Lattice
from .homogenization_cell import HomogenizedCell

def solve_FEM_FenicsX(lattice : "Lattice"):
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

def get_homogenized_properties(lattice: "Lattice"):
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

