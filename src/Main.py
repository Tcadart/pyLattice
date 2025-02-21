from src.Lattice import *
from src.Utils import *
# from Preset_Lattice.Helmet import *
from src.settings import *

meshObject = mesh("Lion")
meshObject.scaleMesh(5)
meshObject.moveMeshToOrigin()


vizualizer = LatticeUtils()
lattice = Lattice(cell_size_X, cell_size_Y, cell_size_Z, number_cell_X, number_cell_Y, number_cell_Z,
                  Lattice_Type, Radius, materialName, gradRadiusProperty, gradDimProperty, gradMatProperty,
                  MethodSim, uncertaintyNodeSD, periodicity=False, randomHybrid=False, meshObject=meshObject)

# print(lattice.getRelativeDensity())
vizualizer.saveJSONToGrasshopper(lattice, "BelfortLion")
vizualizer.visualizeLattice3D(lattice.cells, lattice.latticeDimensionsDict, "Radius", voxelViz=False, explodeVoxel=0.1,
                   plotting = False)

# vizualizer.visualizeMesh(meshObject)

vizualizer.show()
# fig.show()
