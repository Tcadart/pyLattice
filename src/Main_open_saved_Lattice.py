from Lattice.src.Utils import *
from Lattice.src.Lattice import *

nameLattice = "BeamFlexionOptimization"

lattice = Lattice.loadLatticeObject(nameLattice)

vizualizer = LatticeUtils()
vizualizer.visualizeLattice3D(lattice.cells, lattice.latticeDimensionsDict, "Radius", voxelViz=False,
                              explodeVoxel=0.1, plotting = False)
vizualizer.show()
