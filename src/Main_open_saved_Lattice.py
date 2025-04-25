from Utils import LatticeUtils
from src.Lattice import *

nameLattice = "BeamFlexionOptimization.json"

lattice = Lattice.loadLatticeObject(nameLattice)

vizualizer = LatticeUtils()
vizualizer.visualizeLattice3D(lattice.cells, lattice.latticeDimensionsDict, "Radius", voxelViz=False,
                              explodeVoxel=0.1, plotting = False)
vizualizer.show()
