"""
Utility functions for plotting and saving lattice structures in 3D, validating inputs, and saving data in JSON
format for Grasshopper compatibility.
"""
import json
import math
import os
import pickle
from pathlib import Path
from typing import Tuple

import numpy as np
import matplotlib.colors as mcolors

def open_lattice_parameters(file_name: str):
    """
    Open a JSON file containing lattice parameters.

    Parameters:
    -----------
    file_name: str
        Name of the JSON file containing lattice parameters.
    """
    project_root = Path(__file__).resolve().parents[2]
    json_path = project_root / "preset_lattice" / file_name
    if json_path.suffix != ".json":
        json_path = json_path.with_suffix('.json')

    try:
        with open(json_path, 'r') as file:
            lattice_parameters = json.load(file)
    except FileNotFoundError:
        raise FileNotFoundError(f"The file {json_path} does not exist.")
    return lattice_parameters

def _validate_inputs_lattice(cell_size_x, cell_size_y, cell_size_z,
                     num_cells_x, num_cells_y, num_cells_z,
                     geom_types, radii, grad_radius_property, grad_dim_property, grad_mat_property,
                     uncertainty_node, eraser_blocks):
    # Check cell sizes
    assert isinstance(cell_size_x, (int, float)) and cell_size_x > 0, "cell_size_x must be a positive number"
    assert isinstance(cell_size_y, (int, float)) and cell_size_y > 0, "cell_size_y must be a positive number"
    assert isinstance(cell_size_z, (int, float)) and cell_size_z > 0, "cell_size_z must be a positive number"

    # Check number of cells
    assert isinstance(num_cells_x, int) and num_cells_x > 0, "num_cells_x must be a positive integer"
    assert isinstance(num_cells_y, int) and num_cells_y > 0, "num_cells_y must be a positive integer"
    assert isinstance(num_cells_z, int) and num_cells_z > 0, "num_cells_z must be a positive integer"

    # Check lattice type_beam
    assert isinstance(geom_types, list), "Lattice_Type must be a list"
    assert all(isinstance(lt, str) for lt in geom_types), "All elements of Lattice_Type must be strings"

    # Check radii
    assert isinstance(radii, list), "radii must be a list"
    assert all(isinstance(r, float) for r in radii), "All radii values must be floats"
    assert len(radii) == len(geom_types), "The number of radii must be equal to the number of lattice types"

    # Check gradient properties
    if grad_radius_property is not None:
        assert isinstance(grad_radius_property, list), "gradRadiusProperty must be a list"
        assert len(grad_radius_property) == 3, "gradRadiusProperty must be a list of 3 elements"
    if grad_dim_property is not None:
        assert isinstance(grad_dim_property, list), "gradDimProperty must be a list"
        assert len(grad_dim_property) == 3, "gradDimProperty must be a list of 3 elements"
    if grad_mat_property is not None:
        assert len(grad_mat_property) == 2, "gradMatProperty must be a list of 2 elements"
        assert isinstance(grad_mat_property[0], int), "gradMatProperty[0] must be an integer"
        assert isinstance(grad_mat_property[1], int), "gradMatProperty[1] must be an integer"

    # Check optional parameters
    assert isinstance(uncertainty_node, float), "uncertainty_node must be a float"

    if eraser_blocks is not None:
        for erasedPart in eraser_blocks:
            assert len(erasedPart) == 6 and all(
                isinstance(x, float) for x in erasedPart), "eraser_blocks must be a list of 6 floats"

def _validate_inputs_cell(
        pos: list,
        initial_size: list,
        coordinate: list,
        geom_types: list[str],
        radii: list[float],
        grad_radius: list,
        grad_dim: list,
        grad_mat: list,
        uncertainty_node: float,
        _verbose: int,
):
    """Validate inputs for the class constructor."""

    if not isinstance(pos, list) or len(pos) != 3:
        raise TypeError(f"'pos' must be a list of length 3, got {pos}")

    if not isinstance(initial_size, list) or len(initial_size) != 3:
        raise TypeError(f"'initial_size' must be a list of length 3, got {initial_size}")

    if not isinstance(coordinate, list) or len(coordinate) != 3:
        raise TypeError(f"'coordinate' must be a list of length 3, got {coordinate}")

    if not isinstance(geom_types, list) or not all(isinstance(x, str) for x in geom_types):
        raise TypeError(f"'geom_types' must be a list of str, got {geom_types}")

    if not isinstance(radii, list) or not all(isinstance(x, (float, int)) for x in radii):
        raise TypeError(f"'radii' must be a list of float, got {radii}")

    if grad_radius is not None and not isinstance(grad_radius, list):
        raise TypeError(f"'grad_radius' must be a list or None, got {grad_radius}")

    if grad_dim is not None and not isinstance(grad_dim, list):
        raise TypeError(f"'grad_dim' must be a list or None, got {grad_dim}")

    if grad_mat is not None and not isinstance(grad_mat, list):
        raise TypeError(f"'grad_mat' must be a list or None, got {grad_mat}")

    if not isinstance(uncertainty_node, (float, int)):
        raise TypeError(f"'uncertainty_node' must be a float, got {uncertainty_node}")

    if not isinstance(_verbose, int):
        raise TypeError(f"'_verbose' must be an int, got {_verbose}")


def function_penalization_Lzone(radiusAngleData: Tuple[float, float]) -> float:
    """
    Calculate the penalization length based on radii and angle data.

    Args:
        radiusAngleData (Tuple[float, float]): (radii, angle).

    Returns:
        float: Length of the penalization zone.
    """
    if radiusAngleData is None or len(radiusAngleData) != 2:
        raise ValueError("radiusAngleData must be a tuple of (radii, angle).")
    radius, angle = radiusAngleData
    # Case beam quasi-aligned, avoid division by zero
    if angle > 170:
        return 0.0000001
    return radius / math.tan(math.radians(angle) / 2)


def save_lattice_object(lattice, file_name: str = "LatticeObject") -> None:
    """
    Save the lattice object to a file.

    Parameters:
    -----------
    file_name: str
        Name of the file to save (with or without the '.pkl' extension).
    """
    project_root = Path(__file__).resolve().parent.parent
    path = project_root / "saved_lattice_file" / file_name
    if path.suffix != ".pkl":
        path = path.with_suffix('.pkl')

    try:
        with open(path, "wb") as file:
            pickle.dump(lattice, file)
    except Exception as e:
        raise IOError(f"Failed to save lattice pickle: {e}")

    print(f"Lattice pickle saved successfully to {path}")


def _prepare_lattice_plot_data(beam, deformedForm: bool = False):
    beamDraw = set()
    lines = []
    index = []
    nodes = []

    if beam.radius != 0.0 and beam not in beamDraw:
        node1 = beam.point1.deformed_coordinates if deformedForm else (beam.point1.x, beam.point1.y, beam.point1.z)
        node2 = beam.point2.deformed_coordinates if deformedForm else (beam.point2.x, beam.point2.y, beam.point2.z)
        lines.append([node1, node2])
        nodes.append(beam.point1)
        nodes.append(beam.point2)
        index.append(beam.point1.index)
        index.append(beam.point2.index)

    return lines, nodes, index


def _get_beam_color(beam, color_palette, beamColor, idxColor, cells, nbRadiusBins):
    beamColor = beamColor.lower()

    def _to_scalar_radius(r):
        arr = np.atleast_1d(r)
        return float(arr[0])

    if beamColor == "material":
        mat = int(getattr(beam, "material", 0))
        colorBeam = color_palette[mat % len(color_palette)]

    elif beamColor == "type":
        t = int(getattr(beam, "type_beam", getattr(beam, "geom_types", 0)))
        colorBeam = color_palette[t % len(color_palette)]

    elif beamColor == "radii":
        r = _to_scalar_radius(getattr(beam, "radius", 0.0))
        if r not in idxColor:
            idxColor.append(r)
        colorBeam = color_palette[idxColor.index(r) % len(color_palette)]

    elif beamColor == "radiusbin":
        # Construire les bords des classes (bin edges) paresseusement dans idxColor
        if not idxColor:
            all_radii = [
                _to_scalar_radius(getattr(b, "radius", 0.0))
                for c in cells for b in c.beams
                if _to_scalar_radius(getattr(b, "radius", 0.0)) > 0.0
            ]
            if not all_radii:
                idxColor = [0.0, 1.0]
            else:
                min_r, max_r = min(all_radii), max(all_radii)
                # nbRadiusBins classes => nbRadiusBins+1 bornes
                idxColor = list(np.linspace(min_r, max_r, nbRadiusBins + 1))

        r = _to_scalar_radius(getattr(beam, "radius", 0.0))
        bin_idx = np.digitize([r], idxColor, right=False)[0] - 1
        bin_idx = max(0, min(len(idxColor) - 2, bin_idx))
        colorBeam = color_palette[bin_idx % len(color_palette)]

    else:
        colorBeam = "blue"

    return colorBeam, idxColor


def get_boundary_condition_color(fixed_DOF: list[bool]) -> str:
    """
    Generate a color based on the fixed DOFs using a bitmask approach.
    """
    # Convert fixed_DOF to a bitmask integer
    bitmask = sum(2**i for i, val in enumerate(fixed_DOF) if val)

    # Create a color palette (reproducible)
    base_colors = list(mcolors.TABLEAU_COLORS.values()) + list(mcolors.CSS4_COLORS.values())
    color_index = bitmask % len(base_colors)

    return base_colors[color_index]


def visualize_lattice_3D_interactive(lattice, beamColor: str = "Material", voxelViz: bool = False,
                                     deformedForm: bool = False, plotCellIndex: bool = False) -> "go.Figure":
    """
    Visualizes the lattice in 3D using Plotly.

    Parameters:
    -----------
    beamColor: string (default: "Material")
        "Material" -> color by material
        "Type" -> color by type_beam
    voxelViz: boolean (default: False)
        True -> voxel visualization
        False -> beam visualization
    deformedForm: boolean (default: False)
        True -> deformed form
    plotCellIndex: boolean (default: False)
        True -> plot the index of each cell
    """

    color_list = ['blue', 'green', 'red', 'yellow', 'orange', 'purple', 'cyan', 'magenta']
    fig = go.Figure()

    if not voxelViz:
        beamDraw = set()
        nodeDraw = set()
        node_coords = []
        node_colors = []
        lines_x = []
        lines_y = []
        lines_z = []
        line_colors = []
        node1 = None
        node2 = None

        for cell in lattice.cells:
            for beam in cell.beams:
                if beam not in beamDraw:
                    if deformedForm:
                        node1 = beam.point1.deformed_coordinates
                        node2 = beam.point2.deformed_coordinates
                    else:
                        node1 = (beam.point1.x, beam.point1.y, beam.point1.z)
                        node2 = (beam.point2.x, beam.point2.y, beam.point2.z)

                    # Add the beam to the figure
                    lines_x.extend([node1[0], node2[0], None])
                    lines_y.extend([node1[1], node2[1], None])
                    lines_z.extend([node1[2], node2[2], None])

                    # Determine the color of the beam
                    if beamColor == "Material":
                        colorBeam = color_list[beam.material % len(color_list)]
                    elif beamColor == "Type":
                        colorBeam = color_list[beam.type_beam % len(color_list)]
                    else:
                        colorBeam = 'grey'

                    line_colors.extend([colorBeam, colorBeam, colorBeam])

                    beamDraw.add(beam)

                # Add the nodes to the figure
                for node in [node1, node2]:
                    if node not in nodeDraw:
                        node_coords.append(node)
                        nodeDraw.add(node)
                        # Determine the color of the node
                        node_colors.append('black')

            if plotCellIndex:
                cell_center = cell.center_point
                fig.add_trace(go.Scatter3d(
                    x=[cell_center[0]],
                    y=[cell_center[1]],
                    z=[cell_center[2]],
                    mode='text',
                    text=str(cell.index),
                    textposition="top center",
                    showlegend=False
                ))

        # Add the beams to the figure
        fig.add_trace(go.Scatter3d(
            x=lines_x,
            y=lines_y,
            z=lines_z,
            mode='lines',
            line=dict(color=line_colors, width=5),
            hoverinfo='none',
            showlegend=False
        ))

        # Add the nodes to the figure
        if node_coords:
            node_x, node_y, node_z = zip(*node_coords)
            fig.add_trace(go.Scatter3d(
                x=node_x,
                y=node_y,
                z=node_z,
                mode='markers',
                marker=dict(size=4, color=node_colors),
                hoverinfo='none',
                showlegend=False
            ))

    else:
        # Vizualize the lattice as a voxel grid
        for cell in lattice.cells:
            x, y, z = cell.coordinate
            dx, dy, dz = cell.size

            if beamColor == "Material":
                colorCell = color_list[cell.beams[0].material % len(color_list)]
            elif beamColor == "Type":
                colorCell = color_list[int(str(cell.geom_types)[0]) % len(color_list)]
            else:
                colorCell = 'grey'

            # Create the voxel
            fig.add_trace(go.Mesh3d(
                x=[x, x + dx, x + dx, x, x, x + dx, x + dx, x],
                y=[y, y, y + dy, y + dy, y, y, y + dy, y + dy],
                z=[z, z, z, z, z + dz, z + dz, z + dz, z + dz],
                color=colorCell,
                opacity=0.5,
                showlegend=False
            ))

    # Configure the layout
    limMin = min(lattice.x_min, lattice.y_min, lattice.z_min)
    limMax = max(lattice.x_max, lattice.y_max, lattice.z_max)
    fig.update_layout(
        scene=dict(
            xaxis=dict(title='X', range=[limMin, limMax], backgroundcolor='white', showgrid=True, zeroline=True),
            yaxis=dict(title='Y', range=[limMin, limMax], backgroundcolor='white', showgrid=True, zeroline=True),
            zaxis=dict(title='Z', range=[limMin, limMax], backgroundcolor='white', showgrid=True, zeroline=True),
            aspectmode='cube'
        ),
        margin=dict(l=0, r=0, b=0, t=0),
        showlegend=False
    )

    return fig  # Return the figure


def save_JSON_to_Grasshopper(lattice, nameLattice: str = "LatticeObject", multipleParts: int = 1) -> None:
    """
    Save the current lattice object to JSON files for Grasshopper compatibility, separating by cells.

    Parameters:
    -----------
    lattice: Lattice
        Lattice object to save.
    nameLattice: str
        Name of the lattice file to save.
    multipleParts: int, optional (default: 1)
        Number of parts to save.
    """
    folder = "saved_lattice_file"
    if not os.path.exists(folder):
        os.makedirs(folder)

    numberCell = len(lattice.cells)
    cellsPerPart = max(1, numberCell // multipleParts)

    for partIdx in range(multipleParts):
        partName = f"{nameLattice}_part{partIdx + 1}.json" if multipleParts > 1 else f"{nameLattice}.json"
        file_pathJSON = os.path.join(folder, partName)

        partNodesX = []
        partNodesY = []
        partNodesZ = []
        partRadius = []

        startIdx = partIdx * cellsPerPart
        endIdx = min((partIdx + 1) * cellsPerPart, numberCell)

        for cell in lattice.cells[startIdx:endIdx]:
            for beam in cell.beams:
                partNodesX.append(beam.point1.x)
                partNodesX.append(beam.point2.x)
                partNodesY.append(beam.point1.y)
                partNodesY.append(beam.point2.y)
                partNodesZ.append(beam.point1.z)
                partNodesZ.append(beam.point2.z)
                partRadius.append(beam.radius)

        obj = {
            "nodesX": partNodesX,
            "nodesY": partNodesY,
            "nodesZ": partNodesZ,
            "radii": partRadius,
            "maxX": lattice.x_max,
            "minX": lattice.x_min,
            "maxY": lattice.y_max,
            "minY": lattice.y_min,
            "maxZ": lattice.z_max,
            "minZ": lattice.z_min,
            "relativeDensity": lattice.get_relative_density()
        }

        with open(file_pathJSON, 'w') as f:
            json.dump(obj, f)

        print(f"Saved lattice part {partIdx + 1} to {file_pathJSON}")


def plot_coordinate_system(ax):
    """
    Plot a 3D coordinate system with arrows representing the X, Y, and Z axes.
    """
    origin = [0, 0, 0]
    axis_length = 0.8

    ax.quiver(*origin, axis_length, 0, 0, color='r', arrow_length_ratio=0.1)
    ax.quiver(*origin, 0, axis_length, 0, color='g', arrow_length_ratio=0.1)
    ax.quiver(*origin, 0, 0, axis_length, color='b', arrow_length_ratio=0.1)

    ax.text(origin[0] + axis_length, origin[1], origin[2], "X", color='r')
    ax.text(origin[0], origin[1] + axis_length, origin[2], "Y", color='g')
    ax.text(origin[0], origin[1], origin[2] + axis_length, "Z", color='b')
