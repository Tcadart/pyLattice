"""
Class for full-scale lattice simulation.
"""
from basix.ufl import element
import numpy as np
import ufl
from dolfinx import fem

from .simulation_base import SimulationBase

class FullScaleLatticeSimulation(SimulationBase):
    """
    A class to handle full-scale lattice simulation using FenicsX.
    """

    def __init__(self, BeamModel):
        super().__init__(BeamModel)
        self.prepare_simulation()

    def apply_displacement_all_nodes_with_lattice_data(self):
        """
        Applying displacement at all nodes with lattice data.
        """
        alreadyDone = []
        nodePosition = np.round(self.domain.geometry.x, 3)
        triplet_tuples = [tuple(row) for row in nodePosition]
        dictNode = {triplet: idx for idx, triplet in enumerate(triplet_tuples)}
        for cell in self.BeamModel.lattice.cells:
            for beam in cell.beams_cell:
                for node in [beam.point1, beam.point2]:
                    if 1 in node.fixed_DOF:
                        key = tuple(np.round([node.x, node.y, node.z], 3))
                        nodeIndices = dictNode.get(key, None)
                        # Find the degrees of freedom associated with these nodes
                        nodesLocatedDofs = fem.locate_dofs_topological(
                            self._V, self.domain.topology.dim - 1,
                            np.array([nodeIndices], dtype=np.int32)
                        )
                        # Filter the values to apply
                        nodesLocatedDofs_filtered = [
                            val for i, val in enumerate(nodesLocatedDofs)
                            if node.fixed_DOF[i] == 1
                        ]
                        displacement_filtered = [
                            val for i, val in enumerate(node.displacement_vector)
                            if node.fixed_DOF[i] == 1
                        ]
                        # Define the displacement function and set the values
                        u_bc = fem.Function(self._V)
                        u_bc.x.array[nodesLocatedDofs_filtered] = displacement_filtered
                        # Apply the boundary condition
                        self._bcs.append(
                            fem.dirichletbc(u_bc, np.array(nodesLocatedDofs_filtered, dtype=np.int32))
                        )
                        alreadyDone.append(node.index_boundary)

    def set_result_diplacement_on_lattice_object(self):
        """
        Assigns the displacement and rotation values from the simulation to the lattice nodes.
        """
        # Displacement
        displacement_fem = self.u.sub(0).collapse()
        coords_disp = np.round(displacement_fem.function_space.tabulate_dof_coordinates(), 5)
        values_disp = displacement_fem.x.array.reshape((-1, 3))

        # Rotations
        rotation_fem = self.u.sub(1).collapse()
        coords_rot = np.round(rotation_fem.function_space.tabulate_dof_coordinates(), 5)
        values_rot = rotation_fem.x.array.reshape((-1, 3))

        # Mapping dictionaries
        pos_to_disp = {tuple(coord): disp for coord, disp in zip(coords_disp, values_disp)}
        pos_to_rot = {tuple(coord): rot for coord, rot in zip(coords_rot, values_rot)}

        # Node assignment
        for cell in self.BeamModel.lattice.cells:
            for beam in cell.beams_cell:
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

    def set_reaction_force_on_lattice_with_FEM_results(self):
        """
        Set reaction force on boundary condition nodes with FEM results.
        """
        for cell in self.BeamModel.lattice.cells:
            for beam in cell.beams_cell:
                for node in [beam.point1, beam.point2]:
                    if 1 in node.fixed_DOF:
                        RF = self.calculate_reaction_force_and_moment_at_position(
                            np.array([node.x, node.y, node.z]))
                        node.set_reaction_force(RF)

    # def apply_force_on_all_nodes_with_lattice_data(self):
    #     """
    #     Applying force at all nodes with lattice data.
    #
    #     This function applies forces stored in the lattice structure onto the corresponding nodes in the finite element model.
    #     """
    #     nodePosition = np.round(self.domain.geometry.x, 3)
    #     triplet_tuples = [tuple(row) for row in nodePosition]
    #     dictNode = {triplet: idx for idx, triplet in enumerate(triplet_tuples)}
    #
    #     for cell in self.BeamModel.lattice.cells:
    #         for beam in cell.beams_cell:
    #             for node in [beam.point1, beam.point2]:
    #                 if np.any(node.applied_force):  # Check if any force is applied
    #                     key = tuple(np.round([node.x, node.y, node.z], 3))
    #                     node_idx = dictNode.get(key, None)
    #                     if node_idx is None:
    #                         continue
    #
    #                     # DOFs at this vertex for the displacement subspace
    #                     entities = np.array([node_idx], dtype=np.int32)
    #                     dofs_node = fem.locate_dofs_topological(
    #                         self._V.sub(0), self.domain.topology.dim - 1, entities
    #                     )
    #
    #                     # Build a Dirac-like indicator on the *mixed* space so it can appear in the UFL form
    #                     indicator = fem.Function(self._V)
    #                     with indicator.x.petsc_vec.localForm() as local_vec:
    #                         local_vec.set(0.0)
    #                         for idx in np.atleast_1d(dofs_node):
    #                             local_vec[idx] += 1.0
    #
    #                     # Split test function and indicator
    #                     (w_, theta_) = ufl.split(self._u_)
    #                     (ind_w, ind_theta) = ufl.split(indicator)
    #
    #                     f = np.asarray(node.applied_force, dtype=float)
    #                     for i in range(3):
    #                         if not np.isclose(f[i], 0.0):
    #                             term = f[i] * ind_w[i] * w_[i] * self._dx
    #                             self._l_form = term if self._l_form is None else self._l_form + term
    def apply_force_on_all_nodes_with_lattice_data(self):
        """
        Register TRUE nodal point loads by collecting the global DOF indices
        and values (no smearing with dx). They will be injected in the RHS
        vector just before solving.
        """
        self._point_loads.clear()
        node_pos = np.round(self.domain.geometry.x, 3)
        map_pos2idx = {tuple(p): i for i, p in enumerate(node_pos)}

        for cell in self.BeamModel.lattice.cells:
            for beam in cell.beams_cell:
                for node in (beam.point1, beam.point2):
                    if not np.any(node.applied_force):
                        continue
                    key = tuple(np.round([node.x, node.y, node.z], 3))
                    vertex = map_pos2idx.get(key, None)
                    if vertex is None:
                        continue
                    entities = np.array([vertex], dtype=np.int32)
                    # displacement subspace (0), component i = 0,1,2
                    for i in range(3):
                        fi = float(node.applied_force[i])
                        if fi == 0.0:
                            continue
                        dofs = fem.locate_dofs_topological(
                            self._V.sub(0).sub(i), self.domain.topology.dim - 1, entities
                        )
                        # dofs are indices into the PARENT mixed space -> usable directly on RHS Vec
                        for di in np.atleast_1d(dofs):
                            self._point_loads.append((int(di), fi))

    def print_number_DOFs(self):
        """
        Print the total number of degrees of freedom (DOFs) in the simulation.
        """
        num_dofs = self._V.dofmap.index_map.size_global * self._V.dofmap.index_map_bs
        print(f"Total number of DOFs: {num_dofs}")
