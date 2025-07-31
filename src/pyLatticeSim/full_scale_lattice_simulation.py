"""
Class for full-scale lattice simulation.
"""
import numpy as np


from dolfinx import fem
from ufl import dot

from .simulation_base import SimulationBase

class FullScaleLatticeSimulation(SimulationBase):
    """
    All type of lattice simulation
    """

    def __init__(self, BeamModel):
        super().__init__(BeamModel)
        self.prepare_simulation()



    def apply_displacement_all_nodes_with_lattice_data(self):
        """
        Applying displacement at all nodes with lattice data

        Parameters:
        ------------
        latticeData: lattice object
            Lattice object
        """
        alreadyDone = []
        nodePosition = self.domain.geometry.x
        nodePosition = np.round(nodePosition, 3)
        triplet_tuples = [tuple(row) for row in nodePosition]
        dictNode = {triplet: idx for idx, triplet in enumerate(triplet_tuples)}
        for cell in self.BeamModel.lattice.cells:
            for beam in cell.beams:
                for node in [beam.point1, beam.point2]:
                    if 1 in node.fixed_DOF:
                        nodeIndices = dictNode.get(tuple(np.array([node.x,node.y,node.z])), None)
                        # Find the degrees of freedom associated with these nodes
                        nodesLocatedDofs = fem.locate_dofs_topological(self._V, self.domain.topology.dim - 1,
                                                                       np.array([nodeIndices], dtype=np.int32))
                        # Filter the values to apply
                        nodesLocatedDofs_filtered = [val for i, val in enumerate(nodesLocatedDofs)
                                                     if node.fixed_DOF[i]  == 1]
                        displacement_filtered = [val for i, val in enumerate(node.displacement_vector)
                                                    if node.fixed_DOF[i]  == 1]

                        # Define the displacement function and set the values
                        u_bc = fem.Function(self._V)
                        u_bc.x.array[nodesLocatedDofs_filtered] = displacement_filtered
                        # Apply the boundary condition
                        self._bcs.append(fem.dirichletbc(u_bc, np.array(nodesLocatedDofs_filtered, dtype=np.int32)))
                        alreadyDone.append(node.index_boundary)

    def set_result_diplacement_on_lattice_object(self):
        """
        Set the result displacement on the lattice object.
        """
        nodePosition = self.domain.geometry.x
        nodePosition = np.round(nodePosition, 3)
        triplet_tuples = [tuple(row) for row in nodePosition]
        dictNode = {triplet: idx for idx, triplet in enumerate(triplet_tuples)}
        for cell in self.BeamModel.lattice.cells:
            for beam in cell.beams:
                for node in [beam.point1, beam.point2]:
                    nodeCoord = (node.x,node.y,node.z)
                    if nodeCoord in dictNode:
                        index_vertex = dictNode[nodeCoord]  # indice du sommet

                        # --- Récupérer les DOFs de déplacement (sub(0)) ---
                        dof_indices_disp = fem.locate_dofs_topological(
                            self._V.sub(0), 0, [index_vertex])
                        uvals_disp = self.u.vector.array[dof_indices_disp]

                        # --- Récupérer les DOFs de rotation (sub(1)) ---
                        dof_indices_rot = fem.locate_dofs_topological(
                            self._V.sub(1), 0, [index_vertex])
                        uvals_rot = self.u.vector.array[dof_indices_rot]

                        # Concaténer pour avoir [ux, uy, uz, rx, ry, rz]
                        dof_values = np.concatenate((uvals_disp, uvals_rot))

                        # Stocker dans le Node de ta structure
                        node.displacementValue = dof_values

    def apply_force_on_all_nodes_with_lattice_data(self):
        """
        Applying force at all nodes with lattice data.

        This function applies forces stored in the lattice structure onto the corresponding nodes in the finite element model.
        """
        nodesLocatedDofs = []
        nodePosition = self.domain.geometry.x
        nodePosition = np.round(nodePosition, 3)
        triplet_tuples = [tuple(row) for row in nodePosition]
        dictNode = {triplet: idx for idx, triplet in enumerate(triplet_tuples)}

        for cell in self.BeamModel.lattice.cells:
            for beam in cell.beams:
                for node in [beam.point1, beam.point2]:
                    if np.any(node.applied_force):  # Check if any force is applied
                        nodeIndices = dictNode.get((node.x, node.y, node.z), None)
                        if nodeIndices is not None:
                            # Dirac au noeud
                            delta_func = fem.Function(self._V)
                            with delta_func.vector.localForm() as local_vec:
                                local_vec.set(0.0)
                                dofs_node = fem.locate_dofs_topological(self._V.sub(0),
                                                                        self.domain.topology.dim - 1,
                                                                        [nodeIndices])
                                for index in dofs_node:
                                    local_vec[index] += 1.0
                            print(local_vec)
                            # Force à appliquer
                            forceValue = np.array(node.appliedForce)
                            print(forceValue)
                            force_expr = fem.Constant(self.domain, forceValue)
                            print(force_expr)

                            # Ajout de la contribution dans $self._l_form$
                            if self._l_form is None:
                                self._l_form = dot(force_expr, self._V.sub(0)) * delta_func * self._dx
                            else:
                                self._l_form += dot(force_expr, self._V.sub(0)) * delta_func * self._dx

            # forceValue = np.array(node.appliedForce)
                            #
                            # # Create a delta function to apply nodal force
                            # delta_func = fem.Function(self._V)
                            # with delta_func.vector.localForm() as local_vec:
                            #     local_vec.set(0.0)
                            #     for index in nodesLocatedDofs:
                            #         local_vec[index] += 1.0  # Dirac function at node
                            #
                            # force_expr = fem.Constant(self.domain, forceValue)
                            #
                            # # Add nodal force contribution to the L-form
                            # if self._l_form is None:
                            #     self._l_form = dot(force_expr, delta_func) * self._dx
                            # else:
                            #     self._l_form += dot(force_expr, delta_func) * self._dx

    def print_number_DOFs(self):
        """
        Print the total number of degrees of freedom (DOFs) in the simulation.
        """
        num_dofs = self._V.dofmap.index_map.size_global * self._V.dofmap.index_map_bs
        print(f"Total number of DOFs: {num_dofs}")
