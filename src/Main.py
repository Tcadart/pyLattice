from src.Lattice import *
from src.Utils import *
from Mesh.Mesh import *
# from Preset_Lattice.Helmet import *
from src.settings import *

# meshObject = mesh("CutedBone2.stl")
# meshObject.scaleMesh(4)
# meshObject.saveMesh("CutedBone2_scaled")

vizualizer = LatticeUtils()
lattice = Lattice(cell_size_X, cell_size_Y, cell_size_Z, number_cell_X, number_cell_Y, number_cell_Z,
                  Lattice_Type, Radius, materialName, gradRadiusProperty, gradDimProperty, gradMatProperty,
                  MethodSim, uncertaintyNodeSD, periodicity=False, randomHybrid=False, erasedParts=erasedParts,
                  meshObject=None)

# vizualizer.saveLatticeObject(lattice, "Beam3PointFlexion")
# lattice.cutBeamsAtMeshIntersection()
# lattice.printStatistics()

# print(lattice.getRelativeDensity())
# vizualizer.saveJSONToGrasshopper(lattice, "Beam3PointFlexion", multipleParts=1)
# vizualizer.visualizeLattice3D(lattice.cells, lattice.latticeDimensionsDict, "Raduis", voxelViz=False, plotting = False)

lattice.generateMeshLattice(15, cutMeshAtBoundary=True)
# lattice.cutMeshLatticeAtBoundary()

# vizualizer.visualizeMesh(meshObject)
vizualizer.visualizeMesh(lattice.meshLattice)
vizualizer.saveMeshLattice("testLatticeMesh", lattice.meshLattice)

vizualizer.show()
# fig.show()
