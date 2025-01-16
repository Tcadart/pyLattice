"""
Visualization and saving of lattice structures from lattice objects.

Created in 2025-01-16 by Cadart Thomas, University of technology Belfort Montbéliard.
"""
import json
import os
import pickle

import matplotlib.pyplot as plt
from mpl_toolkits.mplot3d.art3d import Line3DCollection, Poly3DCollection
import matplotlib.colors as mcolors
import plotly.graph_objects as go

from .Cell import Cell


def visualizeLattice3D(cells: list["Cell"], latticeDimDict: dict, beamColor: str = "Material",
                       voxelViz: bool = False, deformedForm: bool = False, nameSave: str = None,
                       plotCellIndex: bool = False) -> None:
    """
    Visualizes the lattice in 3D using matplotlib.

    Parameters:
    -----------
    cells: list of Cell
        List of cells to visualize.
    latticeDimDict: dict
        Dictionary containing lattice dimension information (xMin, xMax, yMin, yMax, zMin, zMax).
    beamColor: str, optional (default: "Material")
        Color scheme for beams. Options:
        - "Material": Color by material.
        - "Type": Color by type.
        - "Radius": Color by radius.
    voxelViz: bool, optional (default: False)
        If True, visualize as voxels; otherwise, use beam visualization.
    deformedForm: bool, optional (default: False)
        If True, use deformed node positions.
    nameSave: str, optional
        If provided, save the plot with this name.
    plotCellIndex: bool, optional (default: False)
        If True, plot cell indices.
    """

    def generate_colors(n: int) -> list:
        """Generate a list of `n` distinct colors."""
        base_colors = list(mcolors.TABLEAU_COLORS.values())
        if n <= len(base_colors):
            return base_colors[:n]
        return base_colors + list(mcolors.CSS4_COLORS.values())[:n - len(base_colors)]

    fig = plt.figure()
    ax = fig.add_subplot(111, projection='3d')
    ax.set_title("Lattice generated")

    # Generate a large color palette to avoid missing colors
    max_elements = max(len(cells), 20)  # Dynamically decide the number of colors
    color_palette = generate_colors(max_elements)
    idxColor = []

    if not voxelViz:
        beamDraw = set()
        lines = []
        colors = []
        nodeX, nodeY, nodeZ = [], [], []
        nodeDraw = set()

        for cell in cells:
            for beam in cell.beams:
                node1 = beam.point1.getDeformedPos() if deformedForm else (beam.point1.x, beam.point1.y, beam.point1.z)
                node2 = beam.point2.getDeformedPos() if deformedForm else (beam.point2.x, beam.point2.y, beam.point2.z)

                if beam not in beamDraw:
                    # Assign colors based on the selected scheme
                    if beamColor == "Material":
                        colorBeam = color_palette[beam.material % len(color_palette)]
                    elif beamColor == "Type":
                        colorBeam = color_palette[beam.type % len(color_palette)]
                    elif beamColor == "Radius":
                        if beam.radius not in idxColor:
                            idxColor.append(beam.radius)
                        colorBeam = color_palette[idxColor.index(beam.radius) % len(color_palette)]
                    else:
                        colorBeam = "blue"  # Default color

                    # Add line data
                    lines.append([(node1[0], node1[1], node1[2]), (node2[0], node2[1], node2[2])])
                    colors.append(colorBeam)
                    beamDraw.add(beam)

                # Add node data
                for node in [node1, node2]:
                    if node not in nodeDraw:
                        nodeDraw.add(node)
                        nodeX.append(node[0])
                        nodeY.append(node[1])
                        nodeZ.append(node[2])

            if plotCellIndex:
                ax.text(cell.centerPoint[0], cell.centerPoint[1], cell.centerPoint[2], str(cell.index),
                        color='black', fontsize=10)

        # Plot lines and nodes
        line_collection = Line3DCollection(lines, colors=colors, linewidths=2)
        ax.add_collection3d(line_collection)
        ax.scatter(nodeX, nodeY, nodeZ, c='black', s=5)

    else:  # Voxel visualization
        for cell in cells:
            x, y, z = cell.coordinateCell
            dx, dy, dz = cell.cellSize

            if beamColor == "Material":
                colorCell = color_palette[cell.beams[0].material % len(color_palette)]
            elif beamColor == "Type":
                colorCell = color_palette[cell.latticeType % len(color_palette)]
            else:
                colorCell = "blue"  # Default color

            ax.bar3d(x, y, z, dx, dy, dz, color=colorCell, alpha=1, shade=True, edgecolor='k')

    # Set axis labels and limits
    ax.set_xlabel('X')
    ax.set_ylabel('Y')
    ax.set_zlabel('Z')
    limMin = min(latticeDimDict["xMin"], latticeDimDict["yMin"], latticeDimDict["zMin"])
    limMax = max(latticeDimDict["xMax"], latticeDimDict["yMax"], latticeDimDict["zMax"])
    ax.set_xlim3d(limMin, limMax)
    ax.set_ylim3d(limMin, limMax)
    ax.set_zlim3d(limMin, limMax)

    # Save or show the plot
    if nameSave is not None:
        plt.savefig(nameSave)
    else:
        plt.show()


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
                colorCell = color_list[int(str(cell.latticeType)[0]) % len(color_list)]
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


def visualCellZoneBlocker(lattice, erasedParts: list[tuple]) -> None:
    """
    Visualize the lattice with erased parts

    Parameters:
    -----------
    erasedParts: list of tuple
        List of erased parts with (x_start, y_start, z_start, x_dim, y_dim
    """
    fig = plt.figure()
    ax = fig.add_subplot(111, projection='3d')

    # Plot global lattice cube
    x_max = lattice.xMax
    y_max = lattice.yMax
    z_max = lattice.zMax
    vertices_global = [[0, 0, 0], [x_max, 0, 0], [x_max, y_max, 0], [0, y_max, 0],
                       [0, 0, z_max], [x_max, 0, z_max], [x_max, y_max, z_max], [0, y_max, z_max]]
    ax.add_collection3d(
        Poly3DCollection([vertices_global], facecolors='grey', linewidths=1, edgecolors='black', alpha=0.3))

    # Plot erased region cube
    for erased in erasedParts:
        x_start, y_start, z_start, x_dim, y_dim, z_dim = erased
        vertices_erased = [[x_start, y_start, z_start], [x_start + x_dim, y_start, z_start],
                           [x_start + x_dim, y_start + y_dim, z_start], [x_start, y_start + y_dim, z_start],
                           [x_start, y_start, z_start + z_dim], [x_start + x_dim, y_start, z_start + z_dim],
                           [x_start + x_dim, y_start + y_dim, z_start + z_dim],
                           [x_start, y_start + y_dim, z_start + z_dim]]
        ax.add_collection3d(
            Poly3DCollection([vertices_erased], facecolors='red', linewidths=1, edgecolors='black', alpha=0.6))

    # Set labels and limits
    ax.set_xlabel('X')
    ax.set_ylabel('Y')
    ax.set_zlabel('Z')
    ax.set_xlim([0, x_max])
    ax.set_ylim([0, y_max])
    ax.set_zlim([0, z_max])

    plt.show()



def saveJSONToGrasshopper(lattice, nameLattice: str = "LatticeObject") -> None:
    """
    Save the current lattice object to a JSON file for Grasshopper compatibility.

    Parameters:
    -----------
    nameLattice: str
        Name of the lattice file to save.
    """
    folder = "Saved_Lattice"
    if not os.path.exists(folder):
        os.makedirs(folder)

    file_pathJSON = os.path.join(folder, nameLattice + ".json")

    outNodesX = []
    outNodesY = []
    outNodesZ = []
    outRadius = []
    c = 0
    beamAdded = []
    for cell in lattice.cells:
        for beam in cell.beams:
            if beam not in beamAdded:
                beamAdded.append(beam)
                outNodesX.append(beam.point1.x)
                outNodesX.append(beam.point2.x)
                outNodesY.append(beam.point1.y)
                outNodesY.append(beam.point2.y)
                outNodesZ.append(beam.point1.z)
                outNodesZ.append(beam.point2.z)
                outRadius.append(beam.radius)
                if outRadius[-1] < 0.015:
                    c += 1
                    outRadius[-1] = 0.015

    outMaxX = lattice.xMax
    outMinX = lattice.xMin
    outMaxY = lattice.yMax
    outMinY = lattice.yMin
    outMaxZ = lattice.zMax
    outMinZ = lattice.zMin
    relativeDensity = lattice.getRelativeDensity()

    obj = {
        "nodesX": outNodesX,
        "nodesY": outNodesY,
        "nodesZ": outNodesZ,
        "radius": outRadius,
        "maxX": outMaxX,
        "minX": outMinX,
        "maxY": outMaxY,
        "minY": outMinY,
        "maxZ": outMaxZ,
        "minZ": outMinZ,
        "relativeDensity": relativeDensity
    }
    # Sauvegarder les données au format JSON
    with open(file_pathJSON, 'w') as f:
        json.dump(obj, f)
