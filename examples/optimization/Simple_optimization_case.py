"""
Examples of a simple optimization case.
"""

from pyLatticeOpti.lattice_opti import LatticeOpti

name_file = "optimization_beam_flexion"

lattice_object = LatticeOpti(name_file, verbose=1, convergence_plotting = True)

lattice_object.optimize_lattice()
