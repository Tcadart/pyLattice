import numpy as np
from dolfinx import fem
from petsc4py import PETSc
from dolfinx.fem import petsc


from .simulation_base import SimulationBase


class SchurComplement(SimulationBase):
    """
    Class to perform Schur Complement calculations on a given BeamModel.
    """

    def __init__(self, BeamModel):
        super().__init__(BeamModel)
        self.prepare_simulation()
        self.construct_K()


    def construct_K(self):
        """
        Construct K matrix from variational form
        """
        k = fem.form(self._k_form)
        self._K = petsc.assemble_matrix(k)
        self._K.assemble()

    def calculate_schur_complement(self, tags_nodes_boundary):
        """
        Calculate the Schur complement of the stiffness matrix.

        Parameters:
        -----------
        tags_nodes_boundary: array of int
            Tags identifying the nodes on the boundary
        """
        # Collect all DOFs on the boundary
        DofsBoundary = []
        for nodeTag in tags_nodes_boundary:
            if nodeTag is not None:
                nodeIndices = self.BeamModel.facets.find(nodeTag)
                if nodeIndices is not None and nodeIndices.size > 0:
                    nodesLocatedDofs = fem.locate_dofs_topological(
                        self._V, self.domain.topology.dim - 1, nodeIndices)
                    DofsBoundary.extend(nodesLocatedDofs)
        DofsBoundary = np.array(DofsBoundary, dtype=np.int32)

        if len(DofsBoundary) == 0:
            raise ValueError("No boundary DOFs found; check your boundary node tags.")

        # Create IS objects for boundary and interior DOFs
        is_boundary = PETSc.IS().createGeneral(DofsBoundary)
        num_dofs = self._K.getSize()[0]
        all_dofs = np.arange(num_dofs, dtype=np.int32)
        DofsInterior = np.setdiff1d(all_dofs, DofsBoundary)
        is_interior = PETSc.IS().createGeneral(DofsInterior)

        if len(DofsInterior) == 0:
            raise ValueError("No interior DOFs found; the Schur complement cannot be computed.")

        # Extract the submatrices
        K_II = self._K.createSubMatrix(is_interior, is_interior)
        K_IB = self._K.createSubMatrix(is_interior, is_boundary)
        K_BI = self._K.createSubMatrix(is_boundary, is_interior)
        K_BB = self._K.createSubMatrix(is_boundary, is_boundary)

        # Create a KSP object for solving K_II * U = K_IB
        ksp = PETSc.KSP().create()
        ksp.setOperators(K_II)
        ksp.setType('preonly')
        ksp.getPC().setType('lu')
        ksp.setFromOptions()

        # Solve K_II * U = K_IB
        # Create a dense matrix to store the solution
        num_interior = len(DofsInterior)
        num_boundary = len(DofsBoundary)
        U = PETSc.Mat().createDense([num_interior, num_boundary])
        U.setUp()

        # Solve the linear system for each column of U
        K_IB_dense = K_IB.convert("dense")
        ksp.matSolve(K_IB_dense, U)

        # Compute Schur complement : S = K_BB - K_BI * U
        temp = PETSc.Mat().createDense([num_boundary, num_boundary])
        K_BI.matMult(U, temp)
        S = K_BB.copy()
        S.axpy(-1.0, temp)

        S_dense = S.convert(PETSc.Mat.Type.DENSE)
        SchurNumpy = S_dense.getDenseArray().copy()

        for M in (K_II, K_IB, K_BI, K_BB, U, temp, S, S_dense, K_IB_dense):
            try:
                M.destroy()
            except Exception:
                pass

        return SchurNumpy, DofsBoundary