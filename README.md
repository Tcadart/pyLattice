![pyLattice_logo](/docs/pyLattice_logo.png)

# pyLattice
Design and simulation of truss lattice structures

## Overview

This Python project generates 3D lattice structures with customizable geometry, material properties, gradient, and 
specific modification for simualtions (More information on An optimal penalty method for the joint stiffening in 
beam models of additively manufactured lattice structures: [Paper Link](https://doi.org/10.1016/j.ijsolstr.2024.113107)). 
The structures are defined by cells, beams, and nodes. The code supports various lattice types, including Body-Centered Cubic (BCC), Octet, Diamond, and more (You can also add your own to personalized).

## Key Features
- **Lattice Types**: Generate lattice structures of different geometric configurations, including:
  - BCC (Body-Centered Cubic)
  - Octet
  - Kelvin
  - Cubic (multiple formulations)
  - Hybrid and custom types
- **Gradient Settings**: Apply gradients for cell size, beam radius, and material properties across the lattice.
- **Simulation Preparation**: Supports structural modification for finite element simulations (e.g., node modifications, periodicity).
- **Visualization**: Use matplotlib or Plotly for 3D visualization of the generated structures.

## Project Structure
- **src/pyLattice/lattice.py**: Main class for generating and managing lattice structures.
- **src/pyLattice/beam.py**: Defines the properties and methods of beams connecting nodes in the lattice.
- **src/pyLattice/cell.py**: Defines cells in the lattice structure.
- **src/pyLattice/point.py**: Defines points (nodes) in 3D space for the lattice.
- **src/pyLattice/materials.py**: Handles the material properties.
- **src/pyLattice/geometries/**: Contains geometry definitions for various lattice types.
- **src/pyLatticeSim/**: Simulation backend using FEniCSx for finite element analysis.
- **src/pyLatticeOpti/**: Optimization tools for lattice design.
- **examples/**: Example scripts demonstrating various use cases.
- **main.py**: Main example script to generate and visualize lattice structures.

## Installation

1. Clone the repository:
   ```bash
   git clone https://github.com/Tcadart/pyLattice.git
   cd pyLattice
   ```
2. Set up the project:
   ```bash
   pip install -e .
   ```

For detailed installation instructions including optional dependencies for simulation and mesh operations, see the [Installation Guide](docs/source/Installation_tutorial.md).

## Usage

### Basic Lattice Generation
To generate a simple lattice structure:

```python
from pyLattice.lattice import Lattice
from pyLattice.plotting_lattice import LatticePlotting

# Load lattice from JSON parameter file
name_file = "design/simple_BCC"
lattice = Lattice(name_file)

# Visualize the lattice
visualizer = LatticePlotting()
visualizer.visualize_lattice(lattice, beam_color_type="radii")
```

### Advanced Lattice with Custom Parameters
For more complex lattice configurations, create a JSON parameter file with custom settings:

```json
{
  "cell_size": [1.0, 1.0, 1.0],
  "num_cells": [5, 5, 5],
  "geom_types": ["BCC"],
  "beam_radius": 0.1,
  "materials": {"E": 200000, "nu": 0.3, "rho": 7800}
}
```

Then load and visualize:
```python
lattice = Lattice("path/to/your/config")
visualizer = LatticePlotting()
visualizer.visualize_lattice(lattice)
```

## Configuration

Lattice parameters are configured using JSON files stored in the `data/inputs/` directory. These files define:
- **Cell dimensions**: Size of unit cells
- **Lattice geometry**: Type of lattice structure (BCC, Octet, Kelvin, etc.)
- **Material properties**: Elastic modulus, Poisson's ratio, density
- **Beam properties**: Radius and gradient settings
- **Simulation settings**: Boundary conditions and solver parameters

Example parameter files can be found in the `examples/` directory.

## Dependencies

### Core Dependencies
- Python 3.12+
- NumPy >= 1.26.4
- Matplotlib >= 3.8.0
- Gmsh >= 4.14.0
- SymPy >= 1.12

### Optional Dependencies
- **Simulation**: FEniCSx, UFL, Basix (install via conda-forge)
- **Mesh Operations**: Trimesh, RTree, PyEmbree
- **Optimization**: Scikit-learn

See [Installation Guide](docs/source/Installation_tutorial.md) for detailed setup instructions.

## License

This project is licensed under the MIT License.
