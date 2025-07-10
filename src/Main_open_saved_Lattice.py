from src.Utils import *
from src.Lattice import *

nameLattice = "BeamFlexionOptimization"
# nameLattice = "LatticeTest"

lattice = Lattice.loadLatticeObject(nameLattice)

vizualizer = LatticeUtils()
vizualizer.visualizeLattice3D(lattice.cells, lattice.latticeDimensionsDict, "Radius", voxelViz=False,
                              explodeVoxel=0.1, plotting = False)
vizualizer.show()
