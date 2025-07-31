from lattice import *
from utils import *
from plotting_lattice import *
from mesh_file import *
from preset_lattice.settings import *

name_mesh = "bike-helmet_0_5"
mesh_trimmer = MeshTrimmer(name_mesh)
# mesh_trimmer.plot_mesh()

vizualizer = LatticePlotting()
lattice = Lattice(cell_size_X, cell_size_Y, cell_size_Z, number_cell_X, number_cell_Y, number_cell_Z, Lattice_Type,
                  Radius, materialName, gradRadiusProperty, gradDimProperty, gradMatProperty, mesh_trimmer=mesh_trimmer)

lattice.cut_beam_with_mesh_trimmer()

save_lattice_object(lattice, "Kelvin_helmet")

# lattice.are_cells_identical()

# vizualizer.saveLatticeObject(lattice, "Beam3PointFlexion")
# lattice.cutBeamsAtMeshIntersection()
# lattice.printStatistics()

# print(lattice.getRelativeDensity())
# vizualizer.saveJSONToGrasshopper(lattice, "Beam3PointFlexion", multipleParts=1)
vizualizer.visualize_lattice(lattice, "radii", plotNodeIndex=False, plotting=False)

# lattice.generateMeshLattice(15, cutMeshAtBoundary=True, remeshLattice=False)
# lattice.cutMeshLatticeAtBoundary()
# lattice.generateMeshLatticeGmsh(cutMeshAtBoundary=True, meshSize=0.1, saveMesh=True, saveSTL=False, runGmshApp=False)

# mesh_trimmer = mesh("Lattice.stl")
#
# vizualizer.visualizeMesh(mesh_trimmer)
# vizualizer.visualizeMesh(lattice.mesh_lattice)
# vizualizer.saveMeshLattice("testLatticeMesh", lattice.mesh_lattice)

vizualizer.show()
# fig.show()

# Timing
timing.summary()
