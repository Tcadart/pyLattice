from src.settings import *
from src.Lattice import *

#Generate data from lattice
lattice = Lattice(cell_size_X, cell_size_Y, cell_size_Z, number_cell_X, number_cell_Y, number_cell_Z, Lattice_Type,
                  Radius, gradRadiusProperty, gradDimProperty, gradMatProperty, MethodSim, uncertaintyNodeSD)

fig = lattice.visualizeLattice3D("Type", deformedForm=True, plotCellIndex=False, voxelViz=False)
# fig.show()

