"""
Simple example of how to plot a simple BCC lattice using Matplotlib.
"""
from pathlib import Path

from src.lattice import Lattice
from plotting_lattice import LatticePlotting


project_root = Path(__file__).resolve().parent.parent
json_path = project_root / "preset_lattice" / "simple_BCC.json"

lattice_object = Lattice.from_json(json_path)

vizualizer = LatticePlotting()
vizualizer.visualize_lattice_3D(lattice_object, beam_color_type="radii")