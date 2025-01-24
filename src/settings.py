# settings.py

# *******************************************************************************************************************
# *******************************************************************************************************************

# Variables

# *******************************************************************************************************************
# *******************************************************************************************************************
# Lattice properties
Radius = 0.01
materialName = 'VeroClear'
cell_size = 1
cell_size_X = 1
cell_size_Y = cell_size
cell_size_Z = cell_size
number_cell = 1
number_cell_X = 10
number_cell_Y = 10
number_cell_Z = 2

Lattice_Type = 1000
# -2 => Method random cell
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
# 1000 => Hybrid Lattice with data in hybridLatticeData
hybridLatticeData = [0.1, 0.0, 0.0]

# Gradient on cell dimensions
GradDimRule = 'constant'
GradDimDirection = [1, 0, 0]
GradDimParameters = [2.5, 0.0, 1.5]  # Float
# Gradient on radius of beams
GradRadRule = 'linear'
GradRadDirection = [1, 0, 0]
GradRadParameters = [0.01, 0.0, 1]
# Gradient Rule
# - constant
# - linear
# - parabolic
# - sinusoide
# - exponential

Multimat = 0
# -1 => Full random
# 0 -> materiaux
# 1 -> multimat par couche
GradMaterialDirection = 3  # 1:X / 2:Y / 3:Z

MethodSim = 0
# 0 No modification
# 1 Node Modification

uncertaintyNodeSD = 0.0
# Value of the standard deviation for the uncertainty on the nodes

# erasedParts = [(30.0, 0.0, 0.0, 19.0, 50.0, 19.0)]
# List of erased parts in the lattice
# [(xStart, yStart, zStart, xDim, yDim, zDim), ...] of the erased region


# Gradient properties
gradDimProperty = [GradDimRule, GradDimDirection, GradDimParameters]
gradRadiusProperty = [GradRadRule, GradRadDirection, GradRadParameters]
gradMatProperty = [Multimat, GradMaterialDirection]
