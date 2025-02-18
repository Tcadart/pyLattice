from src.Lattice import *
from src.Utils import *

from Preset_Lattice.Helmet import *

lattice = Lattice(cell_size_X, cell_size_Y, cell_size_Z, number_cell_X, number_cell_Y, number_cell_Z,
                  Lattice_Type, Radius, materialName, gradRadiusProperty, gradDimProperty, gradMatProperty,
                  MethodSim, uncertaintyNodeSD, periodicity=True, randomHybrid=False)

# print(lattice.getRelativeDensity())
# saveJSONToGrasshopper(lattice, "LatticeTest")
print("Start vizualisation")
visualizeLattice3D(lattice.cells, lattice.latticeDimensionsDict, "Radius", voxelViz=False)
# fig.show()
