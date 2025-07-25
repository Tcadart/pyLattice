"""
Example of lattice structure design with all entries from the package.
"""
import os
import sys
from src.Lattice import Lattice

from LatticePlotting import LatticePlotting

project_root = os.path.abspath(os.path.join(os.path.dirname(__file__), ".."))
if project_root not in sys.path:
    sys.path.insert(0, project_root)

json_path = os.path.join(project_root, "preset_lattice", "full_entry_lattice.json")
lattice_object = Lattice.from_json(json_path)

vizualizer = LatticePlotting()
vizualizer.visualizeLattice3D(lattice_object, beam_color_type="radii")
