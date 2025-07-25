from src.Utils import *
from src.Lattice import *

nameLattice = "BeamFlexionOptimization"
# nameLattice = "LatticeTest"

lattice = Lattice.loadLatticeObject(nameLattice)

vizualizer = LatticeUtils()
vizualizer.visualizeLattice3D(lattice.cells, lattice.latticeDimensionsDict, "radii", explode_voxel=0.1, plotting=False)
vizualizer.show()
