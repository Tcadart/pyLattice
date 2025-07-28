from lattice import *
from utils import *
from plotting_lattice import *
from Mesh import *
from preset_lattice.settings import *


# mesh_object = mesh("CutedBone2.stl")
# mesh_object = mesh("Lattice.stl")
# mesh_object.scaleMesh(4)
# mesh_object.saveMesh("CutedBone2_scaled")

vizualizer = LatticePlotting()
lattice = Lattice(cell_size_X, cell_size_Y, cell_size_Z, number_cell_X, number_cell_Y, number_cell_Z, Lattice_Type,
                  Radius, materialName, gradRadiusProperty, gradDimProperty, gradMatProperty, MethodSim,
                  uncertaintyNodeSD, eraser_blocks=erasedParts, mesh_object=None)

# lattice.are_cells_identical()

# vizualizer.saveLatticeObject(lattice, "Beam3PointFlexion")
# lattice.cutBeamsAtMeshIntersection()
# lattice.printStatistics()

# print(lattice.getRelativeDensity())
# vizualizer.saveJSONToGrasshopper(lattice, "Beam3PointFlexion", multipleParts=1)
vizualizer.visualize_lattice_3D(lattice.cells, lattice.lattice_dimension_dict, "radii", plotNodeIndex=False,
                                plotting=False)

# lattice.generateMeshLattice(15, cutMeshAtBoundary=True, remeshLattice=False)
# lattice.cutMeshLatticeAtBoundary()
# lattice.generateMeshLatticeGmsh(cutMeshAtBoundary=True, meshSize=0.1, saveMesh=True, saveSTL=False, runGmshApp=False)

# mesh_object = mesh("Lattice.stl")
#
# vizualizer.visualizeMesh(mesh_object)
# vizualizer.visualizeMesh(lattice.mesh_lattice)
# vizualizer.saveMeshLattice("testLatticeMesh", lattice.mesh_lattice)

vizualizer.show()
# fig.show()

# Timing
timing.summary()
