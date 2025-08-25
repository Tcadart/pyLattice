"""
Domain Decomposition Solver for Full-Scale Lattice Simulation using FenicsX
"""
from colorama import Fore, Style
from typing import TYPE_CHECKING
from scipy.sparse.linalg import LinearOperator, splu
from scipy.sparse import coo_matrix, lil_matrix
import numpy as np

from pyLatticeSim.lattice_sim import LatticeSim
from pyLatticeSim.conjugate_gradient_solver import conjugate_gradient_solver

if TYPE_CHECKING:
    from mesh_file.mesh_trimmer import MeshTrimmer

class DomainDecompositionSolver(LatticeSim):
    """
    A class to handle domain decomposition simulations for lattice structures using FenicsX.
    """

    def __init__(self, name_file: str, mesh_trimmer: "MeshTrimmer" = None, verbose: int = 0):
        super().__init__(name_file, mesh_trimmer, verbose)
        self.domain_decomposition_solver = True

        self.calculate_schur_complement_cells()

        self.enable_preconditioner = None
        self.numberIterationMax = 0
        self._parameters_define = False
        self.preconditioner = None
        self.iteration = 0
        self.residuals = []


    def define_parameters(self, enable_precondioner: bool = True, numberIterationMax: int = 1000):
        """
        Define parameters for the domain decomposition solver.

        Parameters:
        enable_preconditioner: bool
            Enable the preconditioner for the conjugate gradient solver.
        numberIterationMax: int
            Maximum number of iterations for the conjugate gradient solver.
        """
        self.enable_preconditioner = enable_precondioner
        self.numberIterationMax = numberIterationMax
        self._parameters_define = True

    def _check_parameters_defined(self):
        if not self._parameters_define:
            raise ValueError(Fore.RED + "Parameters not defined. Please call define_parameters() before solving." + Style.RESET_ALL)

    def solve_DDM(self):
        """
        Solve the problem with the domain decomposition method.

        Parameters:
        enable_preconditioner: bool
            Enable the preconditioner for the conjugate gradient solver.
        """
        self._check_parameters_defined()

        # Free DOF
        self.define_free_DOF()
        if self.free_DOF == 0:
            raise ValueError(Fore.RED + "No free DOF in the lattice. Process aborted."+ Style.RESET_ALL)
        if self._verbose > 0:
            print(Fore.GREEN + "Free DOF", self.free_DOF, Style.RESET_ALL)

        # Calculate b
        if self._verbose > -1:
            print(Fore.GREEN + "Calculate right-hand side" + Style.RESET_ALL)
        globalDisplacement, _ = self.get_global_displacement()

        b = self.calculate_reaction_force_global(globalDisplacement, rightHandSide=True)
        b = -b # Change sign to have the right-hand side

        self.set_global_free_DOF_index()

        # Initialize local displacement to zero
        self.initialize_displacement()

        # Define the preconditioner
        self.define_preconditioner()

        A_operator = LinearOperator(shape=(self.free_DOF, self.free_DOF), matvec=self.calculate_reaction_force_global)

        print(Fore.GREEN + "Conjugate Gradient started."+ Style.RESET_ALL)

        tol = 1e-6
        mintol = 1e-5
        restart_every = 50
        alpha_max = 100
        xsol, info = conjugate_gradient_solver(A_operator, b, M = self.preconditioner, maxiter=self.numberIterationMax,
                                               tol = tol, mintol = mintol, restart_every = restart_every, alpha_max= alpha_max,
                        callback=lambda xk: self.cg_progress(xk, b, A_operator))

        if self._verbose > -1:
            if info == 0:
                print(Fore.GREEN + "Conjugate Gradient converged."+ Style.RESET_ALL)
            else:
                print(Fore.RED + "Conjugate Gradient did not converge."+ Style.RESET_ALL)

        # Reset boundary conditions
        self.set_boundary_conditions()
        return xsol, info, self.global_displacement_index, b

    def calculate_reaction_force_global(self, globalDisplacement, rightHandSide:bool = False):
        """
        Calculate RF global

        Parameters:
        """
        # Calculate RF local
        self.update_reaction_force_each_cell(globalDisplacement)

        # Calculate the RF global
        globalReactionForce = self.get_global_reaction_force()

        globalReactionForceWithoutFixedDOF = self.get_global_reaction_force_without_fixed_DOF(globalReactionForce,
                                                                                              rightHandSide)

        return globalReactionForceWithoutFixedDOF

    def update_reaction_force_each_cell(self, global_displacement):
        """
        Update RF local

        Parameters:
        -----------
        global_displacement: np.ndarray
            The global displacement vector.
        """
        # print("Update RF local with FenicsX")
        self.initialize_reaction_force()
        datasetDataCell = []
        for cell in self.cells:
            # Set Displacement on nodes (Global to Local)
            cell.set_displacement_at_boundary_nodes(global_displacement, self.global_displacement_index)

            reaction_force_cell = self.solve_sub_problem(cell)

            # Update the RF local in the cell
            cell.set_reaction_force_on_nodes(reaction_force_cell)
        return datasetDataCell

    def solve_sub_problem(self, cell):
        """
        Solve the subproblem on a cell to get the reaction force
        """
        if cell.node_in_order_simulation is None:
            node_in_order = cell.define_node_order_to_simulate()
        if self._verbose > 1:
            print("self.node_in_order", cell.node_in_order_simulation)

        # Check displacement null
        displacement_cell = cell.get_displacement_at_nodes(cell.node_in_order_simulation)
        if np.sum(displacement_cell) != 0:
            # Solve the local problem
            displacement = np.array(displacement_cell).flatten()
            reaction_force_cell = np.dot(cell.schur_complement, displacement)
            reaction_force_cell = reaction_force_cell.reshape(-1, 6)
        else:
            # If displacement is null, set reaction force to zero
            reaction_force_cell = displacement_cell
            if self._verbose > 1:
                print("Displacement is null")
        return reaction_force_cell

    def cg_progress(self, xk, b, A_operator):
        """
        Callback function to track and plot progress of the CG method.

        Parameters:
        -----------
        xk: np.array
            The current solution vector at the k-th iteration.
        b: np.array
            The right-hand side vector of the system.
        A_operator: callable or matrix
            The operator or matrix for the system Ax = b.
        iteration: list
            A list to keep track of the iteration number.
        residuals: list
            A list to store the residual norms for plotting.
        """
        if np.isnan(xk).any():
            raise ValueError(Fore.RED + "NaN detected in the Conjugate Gradient solution. Process aborted."+ Style.RESET_ALL)

        plotting = False
        # Increment iteration count
        self.iteration += 1

        if plotting:
            # Calculate the residual norm
            residual = b - A_operator @ xk
            residual_norm = np.linalg.norm(residual)

            # Append the residual norm for tracking
            self.residuals.append(residual_norm)
            # Save the data at the last iteration
            if residual_norm < 1e-8 or self.iteration == len(b):  # Last iteration condition
                save_file = "ConjugateGradientMethod/ProgressFile/progress_cg.txt"
                with open(save_file, "w") as f:
                    f.write("Conjugate Gradient Progress\n")
                    f.write("Iteration, Residual Norm\n")
                    for i, res in enumerate(self.residuals, 1):
                        f.write(f"{i}, {res:.6e}\n")
                print(f"Progress saved to {save_file}")
        else:
            if self._verbose > 0:
                print(Fore.BLUE + f"Iteration {self.iteration}" + Style.RESET_ALL)
            if self._verbose > 1:
                # Calculate the residual norm
                residual = b - A_operator @ xk
                residual_norm = np.linalg.norm(residual) / np.linalg.norm(b)

                # Append the residual norm for tracking
                self.residuals.append(residual_norm)
                print(f"Residual norm: {residual_norm:.6e}")

    def define_preconditioner(self):
        """
        Define the preconditioner for the conjugate gradient solver.
        Returns:

        """
        if self.enable_preconditioner:
            if self._verbose > 0:
                print(Fore.GREEN + "Define the preconditioner" + Style.RESET_ALL)
            LUSchurComplement, inverseSchurComplement = self.build_preconditioner()
            if LUSchurComplement is not None:
                self.preconditioner = LinearOperator(shape = (self.free_DOF, self.free_DOF),
                                                matvec=lambda x: LUSchurComplement.solve(x))
            elif inverseSchurComplement is not None:
                self.preconditioner = LinearOperator(shape = (self.free_DOF, self.free_DOF),
                                                matvec=lambda x: inverseSchurComplement @ x)

    def build_preconditioner(self):
        """
        Build LU decomposition of the Schur complement matrix for the lattice
        """
        self.build_coupling_operator_cells()
        global_schur_complement = lil_matrix((self.free_DOF, self.free_DOF))

        for cell in self.cells:
            local = cell.build_local_preconditioner()
            global_schur_complement += local

        global_schur_complement = coo_matrix(global_schur_complement)

        if np.any(global_schur_complement.sum(axis=1) == 0):
            print("Attention : There are some rows with all zeros in the Schur complement matrix.")

        # condition number
        cond_number = np.linalg.cond(global_schur_complement.toarray())
        if self._verbose > 0:
            print("Condition number of the Schur complement matrix: ", cond_number)

        # Factorize preconditioner
        LUSchurComplement = None
        inverseSchurComplement = None
        if cond_number > 1e15:
            inverseSchurComplement = np.linalg.pinv(global_schur_complement.toarray())
            if self._verbose > 0:
                print("Using pseudo-inverse of the Schur complement matrix.")
        else:
            global_schur_complement = global_schur_complement.tocsc()
            LUSchurComplement = splu(global_schur_complement)
            if self._verbose > 0:
                print("Using LU decomposition of the Schur complement matrix.")

        return LUSchurComplement, inverseSchurComplement