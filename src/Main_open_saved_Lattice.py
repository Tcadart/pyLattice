from src.Lattice import *

nameLattice = "HybridRandom"

lattice = Lattice.loadLatticeObject(nameLattice)

fig = lattice.visualizeLattice3D("Radius", deformedForm=True, plotCellIndex=False, voxelViz=False)

