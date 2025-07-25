import math
import random

from Timing import *

timing = Timing()


@timing.timeit
def getGradSettings(numCellsX, numCellsY, numCellsZ, gradProperties: list) -> list[list[float]]:
    """
    Generate gradient settings based on the provided rule, direction, and parameters.

    Parameters:
    -----------
    gradProperties: list[Rule, Direction, Parameters]
        All types of properties for gradient definition.

    Return:
    ---------
    gradientData: list[list[float]]
        Generated gradient settings (list of lists).
    """

    def apply_rule(i: int, total_cells: int, param_value: float, rule: str) -> float:
        """
        Apply a specific gradient rule to calculate the gradient factor.

        Parameters:
        -----------
        i : int
            Current cell index in active direction.
        total_cells : int
            Total number of cells in the active direction.
        param_value : float
            Gradient parameter value.
        rule : str
            The gradient rule to apply ('constant', 'linear', etc.).

        Returns:
        --------
        float
            Calculated gradient factor.
        """
        mid = total_cells / 2
        match rule:
            case 'constant':
                return 1.0
            case 'linear':
                return 1.0 + i * param_value
            case 'parabolic':
                if i < mid:
                    return 1.0 + (i / mid) * param_value
                else:
                    return 1.0 + ((total_cells - i - 1) / mid) * param_value
            case 'sinusoide':
                return 1.0 + param_value * math.sin((i / total_cells) * math.pi)
            case 'exponential':
                return 1.0 + math.exp(i * param_value)
            case _:
                raise ValueError(f"Unknown gradient rule: {rule}")

    # Extract gradient properties
    rule, direction, parameters = gradProperties

    # Determine the number of cells in each direction
    number_cells = [numCellsX, numCellsY, numCellsZ]

    indices = [0, 0, 0]

    gradientData = []

    for _ in range(max(number_cells)):
        gradientData.append([
            apply_rule(indices[dim], number_cells[dim], parameters[dim], rule) if direction[dim] == 1 else 1.0
            for dim in range(3)
        ])

        for dim in range(3):
            if direction[dim] == 1 and indices[dim] < number_cells[dim] - 1:
                indices[dim] += 1
    return gradientData


def grad_settings_constant(num_cells_x: int, num_cells_y: int, num_cells_z: int, material_gradient: bool = False) -> (
        list)[list]:
    """
    Generate constant gradient settings (i.e., all values = 1.0).

    Parameters:
    -----------
    num_cells_x : int
    num_cells_y : int
    num_cells_z : int

    Returns:
    --------
    list[list[float]]:
        A list of [1.0, 1.0, 1.0] repeated for the total number of cells.
    """
    if material_gradient:
        return [[[1 for _ in range(num_cells_x)] for _ in range(num_cells_y)] for _ in range(num_cells_z)]
    else:
        total_cells = num_cells_x * num_cells_y * num_cells_z
        return [[1.0, 1.0, 1.0] for _ in range(total_cells)]


@timing.timeit
def gradMaterialSetting(numCellsX, numCellsY, numCellsZ, gradMatProperty: list) -> list:
    """
    Define gradient material settings.

    Parameters:
    ------------
    gradMatProperty: list[Multimat, GradMaterialDirection]
        Set of properties for material gradient.

    Returns:
    --------
    gradMat: list
        3D list representing the material type in the structure.
    """
    multimat, direction = gradMatProperty

    # Initialize gradMat based on `multimat` value
    if multimat == -1:  # Random materials
        return [[[random.randint(1, 3) for _ in range(numCellsX)] for _ in range(numCellsY)] for _ in
                range(numCellsZ)]

    if multimat == 0:  # Single material
        return [[[1 for _ in range(numCellsX)] for _ in range(numCellsY)] for _ in range(numCellsZ)]

    if multimat == 1:  # Graded materials
        # Generate gradient based on the direction
        return [
            [
                [
                    X if direction == 1 else Y if direction == 2 else Z
                    for X in range(numCellsX)
                ]
                for Y in range(numCellsY)
            ]
            for Z in range(numCellsZ)
        ]

    # Default case: return an empty gradMat if no valid `multimat` is provided
    return []
