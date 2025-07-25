from Lattice import *
from Utils import *
from LatticePlotting import *
from Mesh import *
from preset_lattice.settings import *


# meshObject = mesh("CutedBone2.stl")
# meshObject = mesh("Lattice.stl")
# meshObject.scaleMesh(4)
# meshObject.saveMesh("CutedBone2_scaled")

vizualizer = LatticePlotting()
lattice = Lattice(cell_size_X, cell_size_Y, cell_size_Z, number_cell_X, number_cell_Y, number_cell_Z, Lattice_Type,
                  Radius, materialName, gradRadiusProperty, gradDimProperty, gradMatProperty, MethodSim,
                  uncertaintyNodeSD, eraser_blocks=erasedParts, meshObject=None)

# lattice.are_cells_identical()

# vizualizer.saveLatticeObject(lattice, "Beam3PointFlexion")
# lattice.cutBeamsAtMeshIntersection()
# lattice.printStatistics()

# print(lattice.getRelativeDensity())
# vizualizer.saveJSONToGrasshopper(lattice, "Beam3PointFlexion", multipleParts=1)
vizualizer.visualizeLattice3D(lattice.cells, lattice.latticeDimensionsDict, "radii", plotNodeIndex=False,
                              plotting=False)

# lattice.generateMeshLattice(15, cutMeshAtBoundary=True, remeshLattice=False)
# lattice.cutMeshLatticeAtBoundary()
# lattice.generateMeshLatticeGmsh(cutMeshAtBoundary=True, meshSize=0.1, saveMesh=True, saveSTL=False, runGmshApp=False)

# meshObject = mesh("Lattice.stl")
#
# vizualizer.visualizeMesh(meshObject)
# vizualizer.visualizeMesh(lattice.meshLattice)
# vizualizer.saveMeshLattice("testLatticeMesh", lattice.meshLattice)

vizualizer.show()
# fig.show()

# Timing
timing.summary()
