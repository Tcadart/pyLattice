# Lattice Structure Generator

## Overview

This Python project generates 3D lattice structures with customizable geometry, material properties, gradient, and specific modification for simualtions (More information on An optimal penalty method for the joint stiffening in beam models of additively manufactured lattice structures). The structures are defined by cells, beams, and nodes. The code supports various lattice types, including Body-Centered Cubic (BCC), Octet, Diamond, and more (You can also add your own to personalized).

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
- **Lattice.py**: Main class for generating and managing lattice structures.
- **Beam.py**: Defines the properties and methods of beams connecting nodes in the lattice.
- **Cell.py**: Defines cells in the lattice structure.
- **Point.py**: Defines points (nodes) in 3D space for the lattice.
- **Materials.py**: Handles the material properties.
- **Geometry_Lattice.py**: Definition of geometries and contains methods for specific geometric manipulations (e.g., cylindrical transformation).
- **Main.py**: Example script to generate and visualize lattice structures.
- **settings.py**: Configuration file to set lattice settings.

## Installation

1. Clone the repository:
   ```bash
   git clone https://github.com/yourusername/lattice-generator.git
   ```
2. Navigate to the project directory and install dependencies:
   ```bash
   pip install -r requirements.txt
   ```
3. Set up the project:
   ```bash
   python setup.py install
   ```

## Usage

### Basic Lattice Generation
To generate a simple lattice structure:
```python
from Lattice import Lattice

lattice = Lattice.simpleLattice(
    cell_size_x=1.0, 
    cell_size_y=1.0, 
    cell_size_z=1.0, 
    num_cells_x=5, 
    num_cells_y=5, 
    num_cells_z=5, 
    Lattice_Type=0, 
    Radius=0.1
)
```

### Custom Lattice with Gradients
To generate a lattice with custom gradient settings:
```python
gradDimProperty = ['linear', [1, 0, 1], [0.2, 0.0, 0.5]]
gradRadiusProperty = ['parabolic', [0, 0, 1], [0.0, 0.0, 0.1]]

lattice = Lattice(
    cell_size_x=1.0, 
    cell_size_y=1.0, 
    cell_size_z=1.0, 
    num_cells_x=10, 
    num_cells_y=10, 
    num_cells_z=10,
    Lattice_Type=5, 
    Radius=0.1,
    gradRadiusProperty=gradRadiusProperty, 
    gradDimProperty=gradDimProperty,
    gradMatProperty=[0, 3]
)
```

### Visualizing the Lattice
To visualize the generated lattice structure in 3D:
```python
lattice.visualizeLattice3D(beamColor="Material")
```

### Simulation Preparation
Prepare the lattice for a simulation with node uncertainty and boundary conditions:
```python
lattice = Lattice.hybridgeometry(
    cell_size_x=1.0, 
    cell_size_y=1.0, 
    cell_size_z=1.0, 
    simMethod=1, 
    uncertaintyNode=1,
    hybridLatticeData=[0.1, 0.2, 0.15]
)
```

## Configuration

The configuration file `settings.py` allows setting default parameters such as:
- `Radius`: Initial beam radius.
- `Lattice_Type`: Type of lattice geometry to generate.
- `GradDimRule`: Rule for gradient of cell dimensions.
- `GradRadRule`: Rule for gradient of beam radii.
- `MethodSim`: Simulation method (e.g., node modification).
- `uncertaintyNodeSD`: Standard deviation for node uncertainty.

## Dependencies

- Python 3.x
- Matplotlib
- Plotly
- NumPy
- SciPy

## License

This project is licensed under the MIT License.
