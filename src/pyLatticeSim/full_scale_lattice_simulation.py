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

    def apply_force_on_all_nodes_with_lattice_data(self):
        """
        Applying force at all nodes with lattice data.

        This function applies forces stored in the lattice structure onto the corresponding nodes in the finite element model.
        """
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
                            # Dirac at node
                            # delta_func = fem.Function(self._V)
                            scalar_element = element("Lagrange", self.domain.basix_cell(), 1)
                            scalar_space = fem.functionspace(self.domain, scalar_element)
                            delta_func = fem.Function(scalar_space)
                            with delta_func.x.petsc_vec.localForm() as local_vec:
                                local_vec.set(0.0)
                                entities = np.array([nodeIndices], dtype=np.int32)
                                dofs_node = fem.locate_dofs_topological(self._V.sub(0),
                                                                        self.domain.topology.dim - 1,
                                                                        entities)
                                for index in dofs_node:
                                    local_vec[index] += 1.0
                            # Force to apply
                            forceValue = np.array(node.applied_force)
                            v = ufl.TestFunction(self._V)
                            for i in range(len(forceValue)):
                                if not np.isclose(forceValue[i], 0.0):
                                    if self._l_form is None:
                                        self._l_form = forceValue[i] * v[i] * delta_func * self._dx
                                    else:
                                        self._l_form += forceValue[i] * v[i] * delta_func * self._dx

    def print_number_DOFs(self):
        """
        Print the total number of degrees of freedom (DOFs) in the simulation.
        """
        num_dofs = self._V.dofmap.index_map.size_global * self._V.dofmap.index_map_bs
        print(f"Total number of DOFs: {num_dofs}")
