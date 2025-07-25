"""
Utility functions for plotting and saving lattice structures in 3D, validating inputs, and saving data in JSON
format for Grasshopper compatibility.
"""
import json
import math
import os
import pickle
from typing import Tuple

import numpy as np
import plotly.graph_objects as go
import trimesh


def _validate_inputs(cell_size_x, cell_size_y, cell_size_z,
                     num_cells_x, num_cells_y, num_cells_z,
                     Lattice_Type, Radius, materialName, gradRadiusProperty, gradDimProperty, gradMatProperty,
                     uncertaintyNode, periodicity, erasedParts):
    # Check cell sizes
    assert isinstance(cell_size_x, (int, float)) and cell_size_x > 0, "cell_size_x must be a positive number"
    assert isinstance(cell_size_y, (int, float)) and cell_size_y > 0, "cell_size_y must be a positive number"
    assert isinstance(cell_size_z, (int, float)) and cell_size_z > 0, "cell_size_z must be a positive number"

    # Check number of cells
    assert isinstance(num_cells_x, int) and num_cells_x > 0, "num_cells_x must be a positive integer"
    assert isinstance(num_cells_y, int) and num_cells_y > 0, "num_cells_y must be a positive integer"
    assert isinstance(num_cells_z, int) and num_cells_z > 0, "num_cells_z must be a positive integer"

    # Check lattice type
    assert isinstance(Lattice_Type, list), "Lattice_Type must be a list"
    assert all(isinstance(lt, int) for lt in Lattice_Type), "All elements of Lattice_Type must be integers"

    # Check radius
    assert isinstance(Radius, list), "radii must be a list"
    assert all(isinstance(r, float) for r in Radius), "All radius values must be floats"
    assert len(Radius) == len(Lattice_Type), "The number of radius must be equal to the number of lattice types"

    # Check material name
    assert isinstance(materialName, str), "material_name must be a string"

    # Check gradient properties
    if gradRadiusProperty is not None:
        assert isinstance(gradRadiusProperty, list), "gradRadiusProperty must be a list"
        assert len(gradRadiusProperty) == 3, "gradRadiusProperty must be a list of 3 elements"
    if gradDimProperty is not None:
        assert isinstance(gradDimProperty, list), "gradDimProperty must be a list"
        assert len(gradDimProperty) == 3, "gradDimProperty must be a list of 3 elements"
    if gradMatProperty is not None:
        assert len(gradMatProperty) == 2, "gradMatProperty must be a list of 2 elements"
        assert isinstance(gradMatProperty[0], int), "gradMatProperty[0] must be an integer"
        assert isinstance(gradMatProperty[1], int), "gradMatProperty[1] must be an integer"

    # Check optional parameters
    assert isinstance(uncertaintyNode, float), "uncertainty_node must be a float"

    assert isinstance(periodicity, int), "enable_periodicity must be an integer"

    if erasedParts is not None:
        for erasedPart in erasedParts:
            assert len(erasedPart) == 6 and all(
                isinstance(x, float) for x in erasedPart), "eraser_blocks must be a list of 6 floats"


def functionPenalizationLzone(radiusAngleData: Tuple[float, float]) -> float:
    """
    Calculate the penalization length based on radius and angle data.

    Args:
        radiusAngleData (Tuple[float, float]): (radius, angle).

    Returns:
        float: Length of the penalization zone.
    """
    if radiusAngleData is None or len(radiusAngleData) != 2:
        raise ValueError("radiusAngleData must be a tuple of (radius, angle).")
    radius, angle = radiusAngleData
    # Case beam quasi-aligned, avoid division by zero
    if angle > 170:
        return 0.0000001
    return radius / math.tan(math.radians(angle) / 2)


def saveLatticeObject(lattice, file_name: str = "LatticeObject") -> None:
    """
    Save the lattice object to a file.

    Parameters:
    -----------
    file_name: str
        Name of the file to save (with or without the '.pkl' extension).
    """
    folder = "Saved_Lattice"
    os.makedirs(folder, exist_ok=True)

    if not file_name.endswith(".pkl"):
        file_name += ".pkl"

    file_path = os.path.join(folder, file_name)

    with open(file_path, "wb") as file:
        pickle.dump(lattice, file)

    print(f"Lattice saved successfully to {file_path}")


def _prepareLatticePlotData(beam, deformedForm: bool = False):
    """Prepare lines and node positions for lattice plotting."""
    beamDraw = set()
    lines = []
    index = []
    nodes = set()

    if beam.radius != 0.0 and beam not in beamDraw:
        node1 = beam.point1.getDeformedPos() if deformedForm else (
            beam.point1.x, beam.point1.y, beam.point1.z)
        node2 = beam.point2.getDeformedPos() if deformedForm else (
            beam.point2.x, beam.point2.y, beam.point2.z)
        lines.append([node1, node2])
        nodes.update([node1, node2])
        index.append(beam.point1.index)
        index.append(beam.point2.index)

    return lines, nodes, index


def _getBeamColor(beam, color_palette, beamColor, idxColor, cells, nbRadiusBins):
    # Assign colors based on the selected scheme
    if beamColor == "Material":
        colorBeam = color_palette[beam.material % len(color_palette)]
    elif beamColor == "Type":
        colorBeam = color_palette[beam.type % len(color_palette)]
    elif beamColor == "radii":
        if beam.radius not in idxColor:
            idxColor.append(beam.radius)
        colorBeam = color_palette[idxColor.index(beam.radius) % len(color_palette)]
    elif beamColor == "RadiusBin":
        dimRadius = len(cells[0].radius)
        # Ensure radii are extracted properly depending on the dimension
        if not idxColor:
            all_radii = sorted({
                b.radius[dimRadius] if hasattr(b.radius, '__len__') else b.radius
                for c in cells for b in c.beams if
                (b.radius[dimRadius] if hasattr(b.radius, '__len__') else b.radius) > 0
            })
            if not all_radii:
                idxColor = [0.0]
            else:
                min_r, max_r = min(all_radii), max(all_radii)
                bin_edges = np.linspace(min_r, max_r, nbRadiusBins + 1)
                idxColor = bin_edges

        # Get the beam radius at the specified dimension
        radius_value = beam.radius[dimRadius] if hasattr(beam.radius,
                                                         '__len__') else beam.radius
        bin_idx = np.digitize(radius_value, idxColor, right=True) - 1
        colorBeam = color_palette[bin_idx % len(color_palette)]
    else:
        colorBeam = "blue"
    return colorBeam, idxColor


def visualizeLattice3D_interactive(lattice, beamColor: str = "Material", voxelViz: bool = False,
                                   deformedForm: bool = False, plotCellIndex: bool = False) -> "go.Figure":
    """
    Visualizes the lattice in 3D using Plotly.

    Parameters:
    -----------
    beamColor: string (default: "Material")
        "Material" -> color by material
        "Type" -> color by type
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
                        node1 = beam.point1.getDeformedPos()
                        node2 = beam.point2.getDeformedPos()
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
                        colorBeam = color_list[beam.type % len(color_list)]
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
                cell_center = cell.centerPoint
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
            x, y, z = cell.coordinateCell
            dx, dy, dz = cell.cellSize

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
    limMin = min(lattice.xMin, lattice.yMin, lattice.zMin)
    limMax = max(lattice.xMax, lattice.yMax, lattice.zMax)
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


def saveJSONToGrasshopper(lattice, nameLattice: str = "LatticeObject", multipleParts: int = 1) -> None:
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
    folder = "Saved_Lattice"
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
                partRadius.append(max(beam.radius, 0.015))

        obj = {
            "nodesX": partNodesX,
            "nodesY": partNodesY,
            "nodesZ": partNodesZ,
            "radius": partRadius,
            "maxX": lattice.xMax,
            "minX": lattice.xMin,
            "maxY": lattice.yMax,
            "minY": lattice.yMin,
            "maxZ": lattice.zMax,
            "minZ": lattice.zMin,
            "relativeDensity": lattice.getRelativeDensity()
        }

        with open(file_pathJSON, 'w') as f:
            json.dump(obj, f)

        print(f"Saved lattice part {partIdx + 1} to {file_pathJSON}")


def saveLatticeObject(lattice, file_name: str = "LatticeObject") -> None:
    """
    Save the lattice object to a file.

    Parameters:
    -----------
    file_name: str
        Name of the file to save (with or without the '.pkl' extension).
    """
    folder = "Saved_Lattice"
    os.makedirs(folder, exist_ok=True)

    if not file_name.endswith(".pkl"):
        file_name += ".pkl"

    file_path = os.path.join(folder, file_name)

    with open(file_path, "wb") as file:
        pickle.dump(lattice, file)

    print(f"Lattice saved successfully to {file_path}")


def saveMeshLattice(outputPath: str, meshObject: trimesh.Trimesh = None):
    """
    Save the mesh to a file.

    Parameters
    ----------
    outputPath : str
        Path where the mesh should be saved.
    meshObject : mesh, optional
        Mesh object to save. If None, uses the default mesh object.
    """
    if not outputPath.endswith('.stl'):
        outputPath += '.stl'
    if not outputPath.startswith('Mesh/'):
        outputPath = 'Mesh/' + outputPath

    meshObject.export(outputPath)
    print(f"Mesh Lattice saved to {outputPath}")
