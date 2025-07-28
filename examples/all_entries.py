"""
Example of lattice structure design with all entries from the package.
"""
from pathlib import Path

from src.lattice import Lattice
from plotting_lattice import LatticePlotting

project_root = Path(__file__).resolve().parent.parent
json_path = project_root / "preset_lattice" / "full_entry_lattice.json"

lattice_object = Lattice.from_json(json_path)

vizualizer = LatticePlotting()
vizualizer.visualize_lattice_3D(lattice_object, beam_color_type="radii")
