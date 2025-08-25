"""
Superclass of lattice for optimization purposes
"""
import os
from statistics import mean
from typing import TYPE_CHECKING

import joblib
import numpy as np
from scipy.optimize import NonlinearConstraint, Bounds, minimize
from colorama import Fore, Style
from matplotlib import pyplot as plt
import matplotlib

from pyLattice.cell import Cell
from pyLattice.utils import open_lattice_parameters

matplotlib.use('TkAgg')


from pyLatticeSim.lattice_sim import LatticeSim

if TYPE_CHECKING:
    from mesh_file.mesh_trimmer import MeshTrimmer

from pyLattice.timing import *
timingOpti = Timing()

class LatticeOpti(LatticeSim):
    """
    Superclass of lattice for optimization purposes
    """

    def __init__(self, name_file: str, mesh_trimmer: "MeshTrimmer" = None, verbose: int = 0):
        super().__init__(name_file, mesh_trimmer, verbose)
        self._simulation_flag = True
        self._optimization_flag = True

        self.solution = None
        self.objective = None
        self.initial_parameters = None
        self.bounds = None
        self.constraints = []
        self.callback_function = None
        self.max_iteration = 1000
        self.objectif_data = None
        self.objective_function = None
        self.position_objective = None
        self.objective_type = None
        self.number_parameters = None
        self.enable_normalization = False
        self.optimization_parameters = None
        self.constraints_dict = {}
        self.kriging_model_relative_density = None
        self.min_radius = 0.01
        self.max_radius = 0.1



        self.initial_relative_density_constraint = None
        self.initial_continuity_constraint = None
        self.relative_density_poly = []
        self.relative_density_poly_deriv = []
        self.parameter_optimization = []
        self.initial_value_objective = None

        self._get_optimization_parameters(name_file)
        self._set_number_parameters_optimization()

        self.load_relative_density_model() # TODO modifier le comportement



    def optimize_lattice(self):
        """
        Runs the optimization process using SLSQP.
        """
        self.solution = minimize(
            fun=self.objective,
            # jac=self.gradient,
            x0=self.initial_parameters,
            method='SLSQP',
            bounds=self.bounds,
            constraints=self.constraints,
            callback=self.callback_function,
            options={
                'maxiter': self.max_iteration,
                'ftol': 1e-3,
                # 'gtol': 1e-2,
                'disp': True,
                'eps': 1e-2
            }
        )
        plt.show()
        if self.solution.success:
            print("\n✅ Optimization succeeded!")
            print("Optimal parameters:", self.solution.x)
        else:
            print("\n⚠️ Optimization failed!")
            print(self.solution.message)

    def _get_optimization_parameters(self, name_file: str) -> None:
        """
        Define optimization parameters from the input file

        Parameters:
        -----------
        name_file: str
            Name of the input file
        """
        lattice_parameters = open_lattice_parameters(name_file)

        optimization_informations = lattice_parameters.get("optimization_informations", {})
        self.objective_function = optimization_informations.get("objective_function", None)
        self.objective_type = optimization_informations.get("objective_type", None)
        self.position_objective = optimization_informations.get("position_objective", None)
        self.constraints_dict = optimization_informations.get("constraints", {})
        self.optimization_parameters = optimization_informations.get("optimization_parameters", None)
        if self.optimization_parameters is None:
            raise ValueError("No optimization parameters defined.")

    def define_optimization_parameters(self, enable_normalization: bool = False):
        """
        Define optimization parameters
        """
        self.enable_normalization = enable_normalization

    def _initialize_optimization_solver(self, Radius):
        """
        Initialize the solver.
        """
        if self.enable_normalization:
            borneMin = 0
            borneMax = 1
        else:
            borneMin = 0.01
            borneMax = 0.1
        self.bounds = Bounds(lb=[borneMin] * self.number_parameters, ub=[borneMax] * self.number_parameters)

        if self.optimization_parameters["type"] == "unit_cell":
            self.initial_parameters = np.tile(Radius, self.number_parameters // len(Radius))
        elif self.optimization_parameters["type"] == "linear":
            self.initial_parameters = [1.0] * (len(self.optimization_parameters.get("direction", [])) + 1)
        else:
            raise ValueError("Invalid optimization parameters type.")

    def _add_constraint_density(self):
        """
        Add the density constraint to the list of constraints.
        """
        if "relative_density" not in self.constraints_dict:
            return
        self.relative_density_objective = self.constraints_dict.get("relative_density", {}).get("value", None)
        density_nl_constraint = NonlinearConstraint(
            fun=self.density_constraint,
            lb=-np.inf,
            ub=0,
            jac=self.density_constraint_gradient
        )
        self.constraints.append(density_nl_constraint)

    def density_constraint(self, r):
        """
        Density constraint function

        Parameters:
        r: list of float
            List of optimization parameters

        """
        self.set_optimization_parameters(r)
        densConstraint = self.get_relative_density_constraint()
        # if self.densConstraintInitial is None:
        #     self.densConstraintInitial = densConstraint
        # densConstraint = densConstraint/self.densConstraintInitial
        return densConstraint

    def density_constraint_gradient(self, r):
        self.set_optimization_parameters(r)
        gradDensConstraint = self.get_relative_density_gradient_kriging()
        # gradDensConstraint = gradDensConstraint/self.densConstraintInitial
        return gradDensConstraint


    def get_relative_density_constraint(self) -> float:
        """
        Get relative density of the lattice
        """
        relativeDensity = self.get_relative_density()
        error = relativeDensity - self.relative_density_objective
        if self._verbose > 0:
            print("Relative density: ", relativeDensity)
            print("Relative density maximum: ", self.relative_density_objective)
            print("Relative density error: ", error)
        return error

    def get_relative_density(self) -> float:
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
                    cellRelDens.append(cell.get_relative_density_kriging(self.kriging_model_relative_density))
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

    def set_optimization_parameters(self, optimization_parameters_actual: list[float]) -> None:
        """
        Set optimization parameters for the lattice

        Parameters:
        -----------
        optimizationParameters: list of float
            List of optimization parameters
        geomScheme: list of bool
            List of N boolean values indicating the scheme of geometry to optimize
        """
        if len(optimization_parameters_actual) != self.number_parameters:
            raise ValueError("Invalid number of optimization parameters.")


        if self.optimization_parameters["type"] == "unit_cell":
            number_parameters_per_cell = len(self.geom_types)
            for cell in self.cells:
                startIdx = cell.index * number_parameters_per_cell
                endIdx = (cell.index + 1) * number_parameters_per_cell
                radius = optimization_parameters_actual[startIdx:endIdx]
                cell.change_beam_radius(radius, self.grad_radius)
        elif self.optimization_parameters["type"] == "linear":
            dirs = self.optimization_parameters.get("direction", [])
            if not dirs:
                raise ValueError("No directions provided for linear optimization.")

            valid_dirs = {"x", "y", "z"}
            if any(d not in valid_dirs for d in dirs):
                raise ValueError(f"Invalid direction in {dirs}; valid are 'x', 'y', 'z'.")

            expected = len(dirs) + 1  # one coeff per listed direction + intercept d
            if len(optimization_parameters_actual) != expected:
                raise ValueError(
                    f"Invalid number of optimization parameters for linear optimization. "
                    f"Expected {expected} (one per {dirs} + intercept)."
                )

            # Build coefficients a, b, c mapped to x, y, z (missing ones = 0), and intercept d
            coeffs = {"x": 0.0, "y": 0.0, "z": 0.0}
            for i, dkey in enumerate(dirs):
                coeffs[dkey] = float(optimization_parameters_actual[i])
            d_intercept = float(optimization_parameters_actual[-1])

            # Evaluate f(x,y,z) = a*x + b*y + c*z + d at each cell center and set radius
            for cell in self.cells:
                cx, cy, cz = cell.center_point
                value = coeffs["x"] * cx + coeffs["y"] * cy + coeffs["z"] * cz + d_intercept
                value = max(self.min_radius, min(self.max_radius, value))  # Clamp to min max
                # Apply same scalar to all geom types for this cell
                radius_vec = [float(value)] * len(self.geom_types)
                cell.change_beam_radius(radius_vec, self.grad_radius)

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

    def _set_number_parameters_optimization(self):
        """
        Set number of parameters for optimization
        """
        if self.optimization_parameters["type"] == "unit_cell":
            numParameters = 0
            for cell in self.cells:
                numParameters += len(cell.radii)
            self.number_parameters = numParameters
        elif self.optimization_parameters["type"] == "linear":
            self.number_parameters = 2 * len(self.optimization_parameters["direction"])
        else:
            raise ValueError("Invalid optimization parameters type.")


    @timingOpti.timeit
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