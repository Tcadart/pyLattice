from src.Utils import *
from src.Lattice import *

nameLattice = "BeamFlexionOptimization"
# nameLattice = "LatticeTest"

lattice = Lattice.pickle_lattice(nameLattice)

vizualizer = LatticeUtils()
vizualizer.visualizeLattice3D(lattice.cells, lattice.lattice_dimension_dict, "radii", explode_voxel=0.1, plotting=False)
vizualizer.show()
