import random

from Lattice import *
import matplotlib.pyplot as plt
import math

# Function to visualize the chosen cell
def visu_Cellule(ax):
    ax.set_title("Cellule choisie")
    lattice.cells[0].visualize_3d(ax)

# Function to visualize the generated lattice
def visu_lattice(ax):
    ax.set_title("Lattice généré")
    lattice.visualize_3d(ax)

# Function to display lattice points and beams
def display_lattice_points_beams():
    fig = plt.figure()
    ax1 = fig.add_subplot(1, 2, 1, projection='3d')
    visu_Cellule(ax1)
    ax2 = fig.add_subplot(1, 2, 2, projection='3d')
    visu_lattice(ax2)
    plt.show()

# Function to display only the lattice
def display_only_lattice():
    fig = plt.figure()
    ax = fig.add_subplot(111, projection='3d')
    visu_lattice(ax)
    plt.show()

def gradSettings(rule,direction,parameters,number_cell_X,number_cell_Y, number_cell_Z):
    """
    Generate gradient settings based on the provided rule, direction, and parameters.

    :param rule: Gradient rule ('constant', 'linear', 'parabolic', 'sinusoide', 'exponential')
    :param direction: Direction of gradient for each axis (list of 3 integers [X,Y,Z])
    :param parameters: Gradient parameters for each axis (list of 3 floats [X,Y,Z])
    :param number_cell_X: Number of cells along X axis (int)
    :param number_cell_Y: Number of cells along Y axis (int)
    :param number_cell_Z: Number of cells along Z axis (int)
    :return: Generated gradient settings (list of lists)
    """
    def apply_rule(i, number_cell, dir_value, param_value, rule):
        if i < number_cell:
            if i >= 1 and dir_value == 1:
                if rule == 'constant':
                    return 1.0
                elif rule == 'linear':
                    return i * param_value
                elif rule == 'parabolic':
                    return i * param_value if i < number_cell / 2 else (number_cell - i - 1) * param_value
                elif rule == 'sinusoide':
                    if i < number_cell / 4:
                        return i * parameters
                    elif i < number_cell / 2:
                        return (number_cell / 2 - i) * parameters
                    elif i == number_cell / 2:
                        return 1.0
                    elif i < 3 / 4 * number_cell:
                        return (3 / 4 * number_cell - i) * (1 / parameters)
                elif rule == 'exponential':
                    return math.exp(i * param_value)
            return 1.0
        return 0.0  

    # Initialization matrix
    max_cells = max(number_cell_X, number_cell_Y, number_cell_Z)
    gradient = [[0.0, 0.0, 0.0] for _ in range(max_cells)]

    # Processing multiple rules
    for i in range(max_cells):
        number_cells = [number_cell_X, number_cell_Y, number_cell_Z]
        for dim_index in range(3):
            gradient[i][dim_index] = apply_rule(i, number_cells[dim_index], direction[dim_index], parameters[dim_index], rule)


    return gradient

def gradMaterialSetting(Multimat,direction,number_cell_X,number_cell_Y, number_cell_Z):

    # def calculatePositionPointOnPlane(CoeffPlan, PositionPlan , x, y, z):
    #     value = CoeffPlan[0]*x + CoeffPlan[1]*y + CoeffPlan[2]*z + PositionPlan
    #     if value > 0:
    #         return 0
    #     else:
    #         return 1
    if Multimat == -1: #Random
        gradMat = [[[random.randint(1,3) for X in range(number_cell_X)] for Y in range(number_cell_Y)] for Z in range(number_cell_Z)]
    if Multimat == 0: # Mono material
        gradMat = [[[1 for X in range(number_cell_X)] for Y in range(number_cell_Y)] for Z in range(number_cell_Z)]
    elif Multimat == 1: # Graded material
        if direction == 1:
            gradMat = [[[X for X in range(number_cell_X)] for Y in range(number_cell_Y)] for Z in range(number_cell_Z)]
        if direction == 2:
            gradMat = [[[Y for X in range(number_cell_X)] for Y in range(number_cell_Y)] for Z in range(number_cell_Z)]
        if direction == 3:
            gradMat = [[[Z for X in range(number_cell_X)] for Y in range(number_cell_Y)] for Z in range(number_cell_Z)]
    # elif Multimat == 2:
    #     gradMat = [[[0 for X in range(number_cell_X)] for Y in range(number_cell_Y)] for Z in range(number_cell_Z)]
    #     centerCell = []
    #     for x in range(number_cell_X):
    #         for y in range(number_cell_Y):
    #             for z in range(number_cell_Z):
    #                 centerCell.append([0.5+x,0.5+y,0.5+z])
    #     for PositionPlan in range(number_cell_Z):
    #         IdxMat = PositionPlan%2
    #         for IdxCenterCells in range(len(lattice.centerCell)):
    #             if calculatePositionPointOnPlane([1,0,0], PositionPlan , lattice.centerCell[IdxCenterCells,0], lattice.centerCell[IdxCenterCells,1], lattice.centerCell[IdxCenterCells,2]):
    #                 gradMat[lattice.centerCell[IdxCenterCells,3]][lattice.centerCell[IdxCenterCells,4]][lattice.centerCell[IdxCenterCells,5]] = IdxMat
    return gradMat

#*******************************************************************************************************************
#*******************************************************************************************************************

                                    #Variables

#*******************************************************************************************************************
#*******************************************************************************************************************

name_model = 'Lattice_cube'
name_Job = 'Job_1'
name_Part = 'Lattice_Part'
name_Assembly = 'Lattice_assembly'
VectorOrientation = [0,0,-1]
Radius = 0.5
cell_size = 1
cell_size_X = cell_size
cell_size_Y = cell_size
cell_size_Z = cell_size
number_cell = 3
number_cell_X = number_cell
number_cell_Y = number_cell
number_cell_Z = number_cell

Lattice_Type = -1
# -1 => Full random
# 0 => BCC
# 1 => Octet
# 2 => OctetExt
# 3 => OctetInt
# 4 => BCCZ
# 5 => Cubic
# 6 => OctahedronZ
# 7 => OctahedronZcross
# 8 => Kelvin
# 9 => Cubic formulation 2 (centered)
# 10 => Cubic V3
# 11 => Cubic V4
# 12 => New lattice (non connu) GPT generated
# 13 => Diamond
# 14 => Auxetic

# Gradient on cell dimensions
GradDimRule = 'constant'
GradDimDirection = [1,0,0]
GradDimParameters = [1.1,0.0,2.0] #Float
# Gradient on radius of beams
GradRadRule = 'constant'
GradRadDirection = [0,0,1]
GradRadParameters = [1.0,0.0,2.0]
#Gradient Rule
# - constant
# - linear
# - parabolic
# - sinusoide
# - exponential

Multimat = -1
# -1 => Full random
# 0 -> materiaux
# 1 -> multimat par couche
GradMaterialDirection = 3 # 1:X / 2:Y / 3:Z

AnalysisType = 0
# 0 Modelisation lattice only
# 1 Compression Z

MethodSim = 0
# 0 No modification
# 1 Node Modification

#*******************************************************************************************************************
#*******************************************************************************************************************

                                    #Main

#*******************************************************************************************************************
#*******************************************************************************************************************
# Gradient properties
gradDim = gradSettings(GradDimRule,GradDimDirection,GradDimParameters,number_cell_X,number_cell_Y, number_cell_Z) 
gradRadius = gradSettings(GradRadRule,GradRadDirection,GradRadParameters,number_cell_X,number_cell_Y, number_cell_Z) 
gradMat = gradMaterialSetting(Multimat,GradMaterialDirection,number_cell_X,number_cell_Y, number_cell_Z)


#Generate data from lattice
lattice = Lattice(cell_size_X,cell_size_Y,cell_size_Z, number_cell_X,number_cell_Y,number_cell_Z,Lattice_Type, Radius,gradRadius,gradDim,gradMat,MethodSim)

display_only_lattice()
