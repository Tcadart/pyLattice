"""
Superclass of lattice for optimization purposes
"""
from statistics import mean
from typing import TYPE_CHECKING

import numpy as np

from pyLatticeSim.lattice_sim import LatticeSim

if TYPE_CHECKING:
    from mesh_file.mesh_trimmer import MeshTrimmer

class latticeOpti(LatticeSim):
    """
    Superclass of lattice for optimization purposes
    """

    def __init__(self, name_file: str, mesh_trimmer: "MeshTrimmer" = None, verbose: int = 0):
        super().__init__(name_file, mesh_trimmer, verbose)
        self._simulation_flag = True
        self._optimization_flag = True

        self.load_relative_density_model()

        self.objectif_data = None
        self.initial_relative_density_constraint = None
        self.initial_continuity_constraint = None
        self.relative_density_poly = []
        self.relative_density_poly_deriv = []
        self.parameter_optimization = []
        self.kriging_model_relative_density = None
        self.initial_value_objective = None


    def get_relative_density_constraint(self, relativeDensityMax, geomScheme) -> float:
        """
        Get relative density of the lattice
        """
        relativeDensity = self.get_relative_density(geomScheme)
        print("Relative density: ", relativeDensity)
        error = relativeDensity - relativeDensityMax
        print("Relative density maximum: ", relativeDensityMax)
        print("Relative density error: ", error)
        return error

    def get_relative_density(self, geom_scheme=None) -> float:
        """
        Get mean relative density of all cells in lattice

        Returns:
        --------
        meanRelDens: float
            Mean relative density of the lattice
        """
        cellRelDens = []
        for cell in self.cells:
            if self._simulation_flag:
                if self.kriging_model_relative_density is not None:
                    cellRelDens.append(cell.get_relative_density_kriging(self.kriging_model_relative_density, geom_scheme))
            else:
                cellRelDens.append(cell.relative_density)
        meanRelDens = mean(cellRelDens)
        return meanRelDens

    def get_relative_density_gradient(self) -> list[float]:
        """
        Get relative density gradient of the lattice

        Returns:
        --------
        grad: list of float
            Gradient of relative density
        """
        if len(self.relative_density_poly) == 0:
            self.define_relative_density_function()
        if len(self.cells[0].radii) != len(self.relative_density_poly):
            raise ValueError("Invalid radii data.")

        grad = []
        for cell in self.cells:
            grad.append(cell.get_relative_density_gradient())
        return grad



    def get_relative_density_gradient_kriging(self, geom_scheme=None) -> list[float]:
        """
        Get relative density gradient of the lattice using kriging model

        Returns:
        --------
        grad: list of float
            Gradient of relative density
        """
        grad = []
        numberOfCells = len(self.cells)
        if geom_scheme is None or len(geom_scheme) != 3:
            geom_scheme = [i < len(self.radii) for i in range(3)]

        for cell in self.cells:
            gradient3Geom = cell.get_relative_density_gradient_kriging(self.kriging_model_relative_density,
                                                                       geom_scheme) / numberOfCells
            grad.extend(gradient3Geom[geom_scheme])
        return grad

    def get_radius_continuity_difference(self, delta: float = 0.01) -> list[float]:
        """
        Get the difference in radii between connected beams in the lattice

        Parameters:
        -----------
        delta: float
            Minimum difference in radii between connected cells
        """
        radiusContinuityDifference = []
        for cell in self.cells:
            radiusCell = cell.radii
            for neighbours in cell.neighbour_cells:
                for rad in range(len(radiusCell)):
                    radiusContinuityDifference.append((radiusCell[rad] - neighbours.radii[rad]) ** 2 - delta ** 2)
        return radiusContinuityDifference

    def get_radius_continuity_jacobian(self) -> np.ndarray:
        """
        Compute the Jacobian of the radii continuity constraint.

        Returns:
        --------
        np.ndarray
            Jacobian matrix of shape (num_constraints, num_radii)
        """
        rows = []
        cols = []
        values = []
        constraint_index = 0

        for cell in self.cells:
            radiusCell = cell.radii
            for neighbour in cell.neighbour_cells:
                radiusNeighbour = neighbour.radii
                for rad in range(len(radiusCell)):
                    i = cell.index * len(radiusCell) + rad
                    j = neighbour.index * len(radiusCell) + rad
                    diff = radiusCell[rad] - radiusNeighbour[rad]

                    rows.append(constraint_index)
                    cols.append(i)
                    values.append(2 * diff)

                    rows.append(constraint_index)
                    cols.append(j)
                    values.append(-2 * diff)

                    constraint_index += 1

        jacobian = np.zeros((constraint_index, self.get_number_parameters_optimization()))
        for r, c, v in zip(rows, cols, values):
            jacobian[r, c] = v

        return jacobian

    def define_relative_density_function(self, degree: int = 3) -> None:
        """
        Define relative density function
        Possible to define a more complex function with dependency on hybrid cells

        Parameters:
        -----------
        degree: int
            Degree of the polynomial function
        """
        if len(self.relative_density_poly) == 0:
            fictiveCell = Cell([0, 0, 0], [self.cell_size_x, self.cell_size_y, self.cell_size_z], [0, 0, 0],
                               self.geom_types, self.radii, self.grad_radius, self.grad_dim, self.grad_mat,
                               self.uncertainty_node, self._verbose)
            domainRadius = np.linspace(0.01, 0.1, 10)
            for idxRad in range(len(self.radii)):
                radius = np.zeros(len(self.radii))
                relativeDensity = []
                for domainIdx in domainRadius:
                    radius[idxRad] = domainIdx
                    fictiveCell.change_beam_radius([radius], self.grad_radius)
                    relativeDensity.append(fictiveCell.relative_density())
                poly_coeffs = np.polyfit(domainRadius, relativeDensity, degree).flatten()
                poly = np.poly1d(poly_coeffs)
                self.relative_density_poly.append(poly)
                self.relative_density_poly_deriv.append(poly.deriv())

    def set_optimization_parameters(self, optimizationParameters: list[float], geomScheme: list[bool]) -> None:
        """
        Set optimization parameters for the lattice

        Parameters:
        -----------
        optimizationParameters: list of float
            List of optimization parameters
        geomScheme: list of bool
            List of N boolean values indicating the scheme of geometry to optimize
        """
        if len(optimizationParameters) != self.get_number_parameters_optimization(geomScheme):
            raise ValueError("Invalid number of optimization parameters.")

        if geomScheme is None:
            numberOfParametersPerCell = len(self.geom_types)
        else:
            numberOfParametersPerCell = sum(geomScheme)

        for cell in self.cells:
            startIdx = cell.index * numberOfParametersPerCell
            endIdx = (cell.index + 1) * numberOfParametersPerCell
            radius = optimizationParameters[startIdx:endIdx]

            if len(radius) != len(cell.radii):
                # Reconstruct the full radii vector based on geomScheme
                full_radius = []
                i = 0  # index for radii (optimization vector)
                for keep, old in zip(geomScheme, cell.radii):
                    if keep:
                        full_radius.append(radius[i])
                        i += 1
                    else:
                        full_radius.append(old)
                radius = full_radius

            cell.change_beam_radius(radius, self.grad_radius)

    def calculate_objective(self, typeObjective: str) -> float:
        """
        Calculate objective function for the lattice optimization

        Parameters
        ----------
        typeObjective: str
            Type of objective function to calculate (Compliance...)

        Returns
        -------
        objectiveValue: float
            Objective function value
        """
        if typeObjective == "Compliance":
            reactionForce = self.get_global_reaction_force(appliedForceAdded=True)
            reaction_force_array = np.array(list(reactionForce.values())).flatten()
            displacement = np.array(self.get_global_displacement(OnlyImposed=True)[0])
            objective = 0.5 * np.dot(reaction_force_array, displacement)
            if self._verbose > 2:
                np.set_printoptions(threshold=np.inf)
                print("Reaction force: ", reaction_force_array[displacement != 0])
                print("Displacement: ", displacement[displacement != 0])
                print("Compliance: ", objective)
        elif typeObjective == "Displacement":
            setNode = self.find_point_on_lattice_surface(surfaceNames=self.objectif_data["surface"])
            displacements = []
            for node in setNode:
                for dof in self.objectif_data["DOF"]:
                    if dof < 0 or dof > 5:
                        raise ValueError("Invalid degree of freedom index.")
                    displacements.append(node.displacement_vector[dof])
            displacements = np.array(displacements)
            objective = sum(abs(displacements)) / len(displacements)

        elif typeObjective == "Stiffness":
            pass
        return objective

    def get_number_parameters_optimization(self, geomScheme) -> int:
        """
        Get number of parameters for optimization

        Returns:
        --------
        numParameters: int
            Number of parameters for optimization
        """
        numParameters = 0
        for cell in self.cells:
            if geomScheme is None:
                numParameters += len(cell.radii)
            else:
                numParameters += sum(geomScheme)
        return numParameters

    @timing.timeit
    def load_relative_density_model(self, model_path="Lattice/saved_lattice_file/RelativeDensityKrigingModel.pkl"):
        """
        Load the relative density model from a file

        Parameters:
        -----------
        model_path: str
            Path to the model file

        Returns:
        --------
        model: Kriging
            The loaded model
        """
        if not os.path.exists(model_path):
            print(f"Model file not found: {model_path}")
        else:
            gpr = joblib.load(model_path)
            self.kriging_model_relative_density = gpr


    def set_objective_data(self, objectifData: dict) -> None:
        """
        Add objective data to the lattice

        Parameters:
        -----------
        objectif_data: dict
            Dictionary containing objective data
        """
        self.objectif_data = objectifData