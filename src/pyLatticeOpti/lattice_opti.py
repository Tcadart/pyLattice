"""
Superclass of lattice for optimization purposes
"""
import os
from pathlib import Path
from statistics import mean
from typing import TYPE_CHECKING

import joblib
import numpy as np
from scipy.optimize import NonlinearConstraint, Bounds, minimize
from colorama import Fore, Style
from matplotlib import pyplot as plt

from pyLattice.cell import Cell
from pyLattice.utils import open_lattice_parameters
from pyLatticeSim.utils_simulation import solve_FEM_FenicsX
from pyLatticeOpti.plotting_lattice_optim import OptimizationPlotter
from pyLatticeSim.lattice_sim import LatticeSim

if TYPE_CHECKING:
    from data.inputs.mesh_file.mesh_trimmer import MeshTrimmer

from pyLattice.timing import *
timingOpti = Timing()

class LatticeOpti(LatticeSim):
    """
    Superclass of lattice for optimization purposes
    """

    def __init__(self, name_file: str, mesh_trimmer: "MeshTrimmer" = None, verbose: int = 0,
                 convergence_plotting : bool = False):
        super().__init__(name_file, mesh_trimmer, verbose)
        self._convergence_plotting = convergence_plotting
        self._simulation_flag = True
        self._optimization_flag = True

        self.solution = None
        self.actual_objective = None
        self.denorm_objective = None
        self.initial_value_objective = None
        self.initial_parameters = None
        self.bounds = None
        self.constraints = []
        self.iteration = 0
        self.optim_max_iteration = 1000
        self.optim_ftol = 1e-3
        self.optim_disp = True
        self.optim_eps = 1e-2
        self.objectif_data = None
        self.objective_function = None
        self.position_objective = None
        self.objective_type = None
        self.actual_optimization_parameters = []
        self.number_parameters = None
        self.enable_normalization = False
        self.optimization_parameters = None
        self.constraints_dict = {}
        self.kriging_model_relative_density = None
        self.kriging_model_geometries_types = None
        self.min_radius = 0.01
        self.max_radius = 0.1
        self.radius_field = None  # callable: (x,y,z) -> r
        self.radius_field_info = {}  # metadata (e.g., dirs, centers, widths)
        self._flag_plotting_initialized = False
        self.plotting_densities = []
        self.plotting_objective = []
        self.plotter = None




        self.initial_relative_density_constraint = None
        self.initial_continuity_constraint = None
        self.relative_density_poly = []
        self.relative_density_poly_deriv = []
        self.parameter_optimization = []

        self._get_optimization_parameters(name_file)
        self._set_number_parameters_optimization()

        self.load_relative_density_model() # TODO modifier le comportement



    def optimize_lattice(self):
        """
        Runs the optimization process using SLSQP.
        """
        self._initialize_optimization_solver()
        self._add_constraint_density()
        self.solution = minimize(
            fun=self.objective,
            # jac=self.gradient,
            x0=self.initial_parameters,
            method='SLSQP',
            bounds=self.bounds,
            constraints=self.constraints,
            callback=self.callback_function,
            options={
                'maxiter': self.optim_max_iteration,
                'ftol': self.optim_ftol,
                # 'gtol': 1e-2,
                'disp': self.optim_disp,
                'eps': self.optim_eps
            }
        )
        if self._convergence_plotting and self.plotter is not None:
            try:
                self.plotter.update(self.actual_objective, self.solution.x)
            except Exception:
                pass
            self.plotter.finalize(block=True)

        if self.solution.success:
            print("\n✅ Optimization succeeded!")
            print("Optimal parameters:", self.solution.x)
        else:
            print("\n⚠️ Optimization failed!")
            print(self.solution.message)

    def redefine_optim_parameters(self, max_iteration: int = None, ftol: float = None, disp: bool = None,
                                  eps: float = None) -> None:
        """
        Redefine optimization parameters

        Parameters:
        -----------
        max_iteration: int
            Maximum number of iterations for the optimizer
        ftol: float
            Tolerance for termination by the optimizer
        disp: bool
            Whether to display optimization messages
        eps: float
            Step size for numerical approximation of the Jacobian
        """
        if max_iteration is not None:
            self.optim_max_iteration = max_iteration
        if ftol is not None:
            self.optim_ftol = ftol
        if disp is not None:
            self.optim_disp = disp
        if eps is not None:
            self.optim_eps = eps

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
        self._simulation_type = optimization_informations.get("simulation_type", None)
        if self._simulation_type not in {"FEM", "DDM"}:
            print("Simulation type for optimization:", self._simulation_type)
            raise ValueError("Invalid simulation type for optimization. Choose 'FEM' or 'DDM'.")

        self.enable_normalization = optimization_informations.get("enable_parameter_normalization", False)

    def _clamp_radius(self, v: float) -> float:
        return max(self.min_radius, min(self.max_radius, float(v)))

    def _initialize_optimization_solver(self):
        """
        Initialize the solver.
        """
        if self.enable_normalization:
            borneMin = 0
            borneMax = 1
        else:
            borneMin = self.min_radius
            borneMax = self.max_radius
        self.bounds = Bounds(lb=[borneMin] * self.number_parameters, ub=[borneMax] * self.number_parameters)

        if self.optimization_parameters["type"] == "unit_cell":
            raise NotImplementedError("Unit cell optimization not implemented yet.")
            self.initial_parameters = np.tile(borneMax, self.number_parameters // len(borneMax))
        elif self.optimization_parameters["type"] == "linear":
            if self.radius_field is None:
                self._build_radius_field()
            self.initial_parameters = [0.0] * (self.number_parameters - 1)
            self.initial_parameters.append(1)
        elif self.optimization_parameters["type"] == "constant":
            self.initial_parameters = [0.5] * self.number_parameters
        else:
            raise ValueError("Invalid optimization parameters type.")

    def _build_radius_field(self):
        """
        Prepare a parametric radius field r(x,y,z) based on self.optimization_parameters.
        Supported fields:
          - field='linear' with 'direction' subset of ['x','y','z'] -> f=[a, b,c,d] (restricted to listed dirs +
          intercept)
          - field='poly2' with 'terms' subset of ['x','y','z','x2','y2','z2','xy','xz','yz'] + intercept
        """
        opt = self.optimization_parameters or {}
        field_type = opt.get("type", "linear")

        if field_type == "linear":
            dirs = opt.get("direction", ["x", "y", "z"])
            valid = {"x", "y", "z"}
            if any(d not in valid for d in dirs):
                raise ValueError(f"Invalid direction in {dirs}; valid are {valid}.")
            # parameter vector f = [coeff for each dir in order given] + [intercept]
            n = len(dirs) + 1

            def f(x, y, z, theta):
                """Evaluate linear field at (x,y,z) with parameters theta."""
                coeff = dict(zip(dirs, theta[:len(dirs)]))
                d0 = theta[-1]
                return (coeff.get("x", 0.0) * x
                        + coeff.get("y", 0.0) * y
                        + coeff.get("z", 0.0) * z
                        + d0)

            self.radius_field = f
            self.number_parameters = n
            self.radius_field_info = {"type": "linear", "dirs": dirs}

        elif field_type == "poly2":
            # polynomial up to quadratic with chosen terms + intercept
            terms = opt.get("terms", ["x", "y", "z"])  # add 'x2','y2','z2','xy','xz','yz' as needed
            valid = {"x", "y", "z", "x2", "y2", "z2", "xy", "xz", "yz"}
            if any(t not in valid for t in terms):
                raise ValueError(f"Invalid term in {terms}; valid are {valid}.")
            n = len(terms) + 1  # + intercept

            def f(x, y, z, theta):
                """Evaluate poly2 field at (x,y,z) with parameters theta."""
                coeffs = dict(zip(terms, theta[:len(terms)]))
                d0 = theta[-1]
                v = d0
                v += coeffs.get("x", 0.0) * x
                v += coeffs.get("y", 0.0) * y
                v += coeffs.get("z", 0.0) * z
                v += coeffs.get("x2", 0.0) * (x * x)
                v += coeffs.get("y2", 0.0) * (y * y)
                v += coeffs.get("z2", 0.0) * (z * z)
                v += coeffs.get("xy", 0.0) * (x * y)
                v += coeffs.get("xz", 0.0) * (x * z)
                v += coeffs.get("yz", 0.0) * (y * z)
                return v

            self.radius_field = f
            self.number_parameters = n
            self.radius_field_info = {"type": "poly2", "terms": terms}
        else:
            raise ValueError(f"Unknown field type '{field_type}'.")

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
        print(Fore.LIGHTGREEN_EX + "Density constraint function" + Fore.RESET)
        self.set_optimization_parameters(r)
        densConstraint = self.get_relative_density_constraint()
        # if self.densConstraintInitial is None:
        #     self.densConstraintInitial = densConstraint
        # densConstraint = densConstraint/self.densConstraintInitial
        if self._verbose > 0:
            print("Density constraint: ", densConstraint)
        return densConstraint

    def density_constraint_gradient(self, r):
        """
        Density constraint gradient function

        Parameters:
        r: list of float
            List of optimization parameters
        """
        print(Fore.LIGHTBLUE_EX + "Density constraint gradient function" + Fore.RESET)
        self.set_optimization_parameters(r)
        gradDensConstraint = self.get_relative_density_gradient_kriging()
        print("Density constraint gradient: ", gradDensConstraint)
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
                    cellRelDens.append(cell.get_relative_density_kriging(self.kriging_model_relative_density,
                                                                         self.kriging_model_geometries_types))
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


    def get_relative_density_gradient_kriging(self) -> np.array:
        """
        Get relative density gradient of the lattice using kriging model

        Returns:
        --------
        grad: list of float
            Gradient of relative density
        """
        if self.optimization_parameters["type"] == "unit_cell":
            number_cells = len(self.cells)
            grad = np.zeros(self.number_parameters)

            for cell in self.cells:
                gradient3Geom = cell.get_relative_density_gradient_kriging(
                    self.kriging_model_relative_density, self.kriging_model_geometries_types) / number_cells
                grad[cell.index] = gradient3Geom

            return grad
        elif self.optimization_parameters["type"] == "linear":
            raise NotImplementedError("Gradient for linear optimization not implemented yet.")
            numberOfCells = len(self.cells)
            grad = np.zeros(self.number_parameters)
            dirs = self.optimization_parameters.get("direction", [])
            if not dirs:
                raise ValueError("No directions provided for linear optimization.")
            valid_dirs = {"x", "y", "z"}
            if any(d not in valid_dirs for d in dirs):
                raise ValueError(f"Invalid direction in {dirs}; valid are 'x', 'y', 'z'.")
            coeffs = {"x": 0.0, "y": 0.0, "z": 0.0}
            for i, dkey in enumerate(dirs):
                coeffs[dkey] = float(self.initial_parameters[i])
            d_intercept = float(self.initial_parameters[-1])

            for cell in self.cells:
                cx, cy, cz = cell.center_point
                value = coeffs["x"] * cx + coeffs["y"] * cy + coeffs["z"] * cz + d_intercept
                value = max(self.min_radius, min(self.max_radius, value))
                gradient3Geom = cell.get_relative_density_gradient_kriging(
                    self.kriging_model_relative_density, self.kriging_model_geometries_types) / numberOfCells
                for i, dkey in enumerate(dirs):
                    grad[i] += gradient3Geom[0] * cx if dkey == "x" else gradient3Geom[1] * cy if dkey == "y" else gradient3Geom[2] * cz
                grad[-1] += sum(gradient3Geom)
            return grad
        elif self.optimization_parameters["type"] == "constant":
            n_cells = len(self.cells)
            hybrid = bool(self.optimization_parameters.get("hybrid", False))
            n_geom = len(self.geom_types)
            scale = (self.max_radius - self.min_radius) if self.enable_normalization else 1.0

            if hybrid:
                grad = np.zeros(n_geom, dtype=float)
                for cell in self.cells:
                    g_cell = np.asarray(
                        cell.get_relative_density_gradient_kriging(
                            self.kriging_model_relative_density,
                            self.kriging_model_geometries_types
                        ),
                        dtype=float
                    )
                    grad += g_cell
                grad /= n_cells
                grad *= scale
                return grad
            else:
                g_total = 0.0
                for cell in self.cells:
                    g_cell = np.asarray(
                        cell.get_relative_density_gradient_kriging(
                            self.kriging_model_relative_density,
                            self.kriging_model_geometries_types
                        ),
                        dtype=float
                    )
                    g_total += float(np.sum(g_cell))
                g_total /= n_cells
                g_total *= scale
                return np.array([g_total], dtype=float)
        else:
            raise ValueError("Invalid optimization parameters type.")

    def objective(self, r) -> float:
        """
        Objective function for the optimization

        Parameters:
        -----------
        r: list of float
            List of optimization parameters

        Returns:
        --------
        objectiveValue: float
            Value of the objective function
        """
        if self._verbose >= 1:
            print(Fore.GREEN + "Objective function" + Fore.RESET)
        print("Parameters:", r)
        self.set_optimization_parameters(r)
        self._simulate_lattice_equilibrium()
        objective = self.calculate_objective()
        print("objective", objective)
        self.denorm_objective = objective
        objectiveNorm = self.normalize_objective_data(objective, objective_type = True)

        if self.objective_function == "max":
            objectiveNorm = -objectiveNorm
        elif self.objective_function == "min":
            pass
        else:
            raise ValueError("objective_function must be 'min' or 'max'")

        print("Normalized objective", objectiveNorm)
        self.actual_objective = objectiveNorm
        print("Actual objective", self.actual_objective)
        return self.actual_objective

    def normalize_objective_data(self, value, objective_type: bool = True):
        """
        Normalize the objective value if normalization is enabled.

        Parameters:
        -----------
        value: float
            The objective value to normalize.
        objectif : bool
            True if normalizing objective, False if normalizing gradient.

        Returns:
        --------
        float
            Normalized objective value.
        """
        if self.enable_normalization:
            if self.initial_value_objective is None and objective_type:
                # Use a positive scale to avoid sign flips when normalizing
                scale = abs(value)
                if scale == 0.0:
                    scale = 1.0  # fallback to avoid division by zero
                self.initial_value_objective = scale
                print("Initial objective value (scale): ", self.initial_value_objective)

            value_normalized = value / self.initial_value_objective
            if self._verbose >= 1:
                if objective_type:
                    print("Objective normalized: ", value_normalized)
                else:
                    print("Gradient normalized: ", value_normalized)
            return value_normalized
        else:
            return value


    def _simulate_lattice_equilibrium(self):
        """
        Simulate the lattice equilibrium using the internal solver.
        """
        self._initialize_simulation_parameters()
        if self._simulation_type == "FEM":
            solve_FEM_FenicsX(self)
        elif self._simulation_type == "DDM":
            if self._parameters_define is False:
                self.define_parameters(True, 100)
                self.solve_DDM()
        else:
            raise ValueError("Invalid simulation type for optimization. Choose 'FEM' or 'DDM'.")

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
            for neighbours in cell.get_all_cell_neighbours():
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
            for neighbour in cell.get_all_cell_neighbours():
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
                    fictiveCell.change_beam_radius([radius])
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
        if (    self.actual_optimization_parameters is not None
                and len(self.actual_optimization_parameters) == len(optimization_parameters_actual)
                and np.allclose(optimization_parameters_actual, self.actual_optimization_parameters,
                                rtol=1e-10, atol=1e-10)):
            return

        self.actual_optimization_parameters = optimization_parameters_actual.copy()

        if self._verbose >= 2:
            print(Style.DIM + "Optimization parameters: ", self.actual_optimization_parameters, Style.RESET_ALL)

        if self.optimization_parameters["type"] == "constant":
            for cell in self.cells:
                radius = self.denormalize_optimization_parameters([float(optimization_parameters_actual[0])])
                cell.change_beam_radius(radius * len(self.geom_types))
        elif self.optimization_parameters["type"] == "unit_cell":
            number_parameters_per_cell = len(self.geom_types)
            for cell in self.cells:
                startIdx = cell.index * number_parameters_per_cell
                endIdx = (cell.index + 1) * number_parameters_per_cell
                radius = optimization_parameters_actual[startIdx:endIdx]
                cell.change_beam_radius(radius)
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
                coeffs[dkey] = self.denormalize_optimization_parameters([float(optimization_parameters_actual[i])])[0]
            d_intercept = self.denormalize_optimization_parameters([float(optimization_parameters_actual[-1])])[0]


            # Evaluate f(x,y,z) = a*x + b*y + c*z + d at each cell center and set radius
            for cell in self.cells:
                cx, cy, cz = cell.center_point
                value = coeffs["x"] * cx + coeffs["y"] * cy + coeffs["z"] * cz + d_intercept
                value = max(self.min_radius, min(self.max_radius, value))  # Clamp to min max
                # Apply same scalar to all geom types for this cell
                radius_vec = [float(value)] * len(self.geom_types)
                cell.change_beam_radius(radius_vec)
        else:
            raise ValueError("Invalid optimization parameters type.")


    def denormalize_optimization_parameters(self, r_norm: list[float]) -> list[float]:
        """
        Denormalize optimization parameters

        Parameters:
        -----------
        r_norm: list of float
            List of normalized optimization parameters

        Returns:
        --------
        r: list of float
            List of denormalized optimization parameters
        """
        if not self.enable_normalization:
            return r_norm
        r = []
        for val in r_norm:
            denorm_val = self._clamp_radius(val * (self.max_radius - self.min_radius) + self.min_radius)
            r.append(denorm_val)
        return r

    def calculate_objective(self) -> float:
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
        if self.objective_type == "compliance":
            reactionForce = self.get_global_reaction_force(appliedForceAdded=True)
            reaction_force_array = np.array(list(reactionForce.values())).flatten()
            displacement = np.array(self.get_global_displacement_DDM(OnlyImposed=True)[0])
            objective = - 0.5 * np.dot(reaction_force_array, displacement)
            if self._verbose > 2:
                np.set_printoptions(threshold=np.inf)
                print("Reaction force: ", reaction_force_array[displacement != 0])
                print("Displacement: ", displacement[displacement != 0])
                print("Compliance: ", objective)
        elif self.objective_type == "displacement":
            setNode = self.find_point_on_lattice_surface(surfaceNames=self.objectif_data["surface"])
            displacements = []
            for node in setNode:
                for dof in self.objectif_data["DOF"]:
                    if dof < 0 or dof > 5:
                        raise ValueError("Invalid degree of freedom index.")
                    displacements.append(node.displacement_vector[dof])
            displacements = np.array(displacements)
            objective = sum(abs(displacements)) / len(displacements)

        elif self.objective_type == "stiffness":
            raise NotImplementedError("Stiffness objective not implemented yet.")
        else:
            raise ValueError("Invalid objective function type.")
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
            self._build_radius_field()
        elif self.optimization_parameters["type"] == "constant":
            if self.optimization_parameters.get("hybrid", False):
                self.number_parameters = len(self.geom_types)
            else:
                self.number_parameters = 1
        else:
            raise ValueError("Invalid optimization parameters type.")

    def gradient(self, r):
        """
        Gradient function for the optimization
        """
        if self._verbose >= 1:
            print(Fore.BLUE + "Gradient function"+ Fore.RESET)
        self.set_optimization_parameters(r)
        print("Calculating gradient...") # TODO
        # grad, _ = self.conjugateGradient.calculateGradient(self.GeomScheme)
        # print("Gradient", grad)
        # gradNorm = self.normalizeObjective(grad, objectif=False)
        # self.actualGradient = gradNorm
        # print("Actual gradient", self.actualGradient)
        # return gradNorm

    @timingOpti.timeit
    def load_relative_density_model(self, model_name="RelativeDensityKrigingModel"):
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
        path_model = Path(__file__).parents[2] / "src" / "pyLatticeOpti" / model_name
        if path_model.suffix != ".pkl":
            path_model = str(path_model) + ".pkl"

        if not os.path.exists(path_model):
            print(f"Model file not found: {path_model}")
        else:
            gpr = joblib.load(path_model)
            self.kriging_model_relative_density = gpr
            self.kriging_model_geometries_types = ["BCC", "Hybrid1", "Hybrid4"] # TODO a modifier

    def callback_function(self, r):
        """
        Callback function for the optimization (Printing and plotting)
        """
        self.iteration += 1
        print(f"[Itération {self.iteration}] Objective : {self.actual_objective}")

        if "relative_density" in self.constraints_dict:
            print("Relative density = ", self.get_relative_density())

        if self._convergence_plotting:
            if self.plotter is None:
                self.plotter = OptimizationPlotter(self)
            self.plotter.update(self.actual_objective, r)





