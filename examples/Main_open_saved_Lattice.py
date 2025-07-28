from src.utils import *
from src.lattice import *

nameLattice = "BeamFlexionOptimization"
# nameLattice = "LatticeTest"

lattice = Lattice.pickle_lattice(nameLattice)

vizualizer = LatticeUtils()
vizualizer.visualize_lattice_3D(lattice.cells, lattice.lattice_dimension_dict, "radii", explode_voxel=0.1,
                                plotting=False)
vizualizer.show()

#TODO : A refaire

