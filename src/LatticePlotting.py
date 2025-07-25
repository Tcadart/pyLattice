"""
Visualization and saving of lattice structures from lattice objects.

Created in 2025-01-16 by Cadart Thomas, University of technology Belfort-Montbéliard.
"""
from mpl_toolkits.mplot3d.art3d import Line3DCollection, Poly3DCollection

from src.Cell import Cell
from src.Utils import *
from src.Utils import _getBeamColor, _prepareLatticePlotData

import matplotlib
import matplotlib.colors as mcolors
import matplotlib.pyplot as plt

matplotlib.use('TkAgg')  # Or 'Qt5Agg' if you prefer Qt backend


class LatticePlotting:
    """
    Class for visualizing lattice structures in 3D.
    """

    def __init__(self, initFig: bool = False):
        if initFig:
            self.initFigure()
        self.fig = None
        self.ax = None
        self.minAxis = None
        self.maxAxis = None
        self.initFig = initFig
        self.axisSet = False

    def initFigure(self):
        """Initialize the 3D figure for plotting."""
        self.fig = plt.figure()
        self.ax = self.fig.add_subplot(111, projection='3d')
        self.initFig = True
        self.ax.set_axis_off()

    def _setMinMaxAxis(self, latticeDimDict: dict) -> None:
        limMin = min(latticeDimDict["xMin"], latticeDimDict["yMin"], latticeDimDict["zMin"])
        limMax = max(latticeDimDict["xMax"], latticeDimDict["yMax"], latticeDimDict["zMax"])
        self.minAxis = min(limMin, self.minAxis) if self.minAxis is not None else limMin
        self.maxAxis = max(limMax, self.maxAxis) if self.maxAxis is not None else limMax
        self.axisSet = True

    def visualizeLattice3D(self, lattice_object, beam_color_type: str = "Material",
                           voxelViz: bool = False, deformedForm: bool = False, file_save_path: str = None,
                           plotCellIndex: bool = False, plotNodeIndex: bool = False, explode_voxel: float = 0.0,
                           plotting: bool = True, nbRadiusBins: int = 5) -> None:

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
            - "radii": Color by radius.
        voxelViz: bool, optional (default: False)
            If True, visualize as voxels; otherwise, use beam visualization.
        deformedForm: bool, optional (default: False)
            If True, use deformed node positions.
        nameSave: str, optional
            If provided, save the plot with this name.
        plotCellIndex: bool, optional (default: False)
            If True, plot cell indices.
        """
        if self.initFig is False:
            self.initFigure()

        def generate_colors(n: int) -> list:
            """Generate a list of `n` distinct colors."""
            base_colors = list(mcolors.TABLEAU_COLORS.values())
            if n <= len(base_colors):
                return base_colors[:n]
            return base_colors + list(mcolors.CSS4_COLORS.values())[:n - len(base_colors)]

        cells = lattice_object.cells
        latticeDimDict = lattice_object.lattice_dimension_dict

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
                    if beam.radius != 0.0 and beam not in beamDraw:
                        colorBeam, idxColor = _getBeamColor(beam, color_palette, beam_color_type, idxColor, cells,
                                                                 nbRadiusBins)

                        # Add line and node data
                        beam_lines, beam_nodes, beam_indices = _prepareLatticePlotData(beam, deformedForm)
                        lines.extend(beam_lines)
                        colors.extend([colorBeam] * len(beam_lines))  # One color per line

                        for i, node in enumerate(beam_nodes):
                            if node not in nodeDraw:
                                nodeDraw.add(node)
                                nodeX.append(node[0])
                                nodeY.append(node[1])
                                nodeZ.append(node[2])
                                if plotNodeIndex:
                                    self.ax.text(node[0], node[1], node[2], str(beam_indices[i]), fontsize=6,
                                                 color='gray')

                        beamDraw.add(beam)

                if plotCellIndex:
                    self.ax.text(cell.centerPoint[0], cell.centerPoint[1], cell.centerPoint[2], str(cell.index),
                                 color='black', fontsize=10)

            # Plot lines and nodes
            line_collection = Line3DCollection(lines, colors=colors, linewidths=2)
            self.ax.add_collection3d(line_collection)
            self.ax.scatter(nodeX, nodeY, nodeZ, c='black', s=5)

        else:  # Voxel visualization
            for cell in cells:
                x, y, z = cell.coordinateCell
                dx, dy, dz = cell.cellSize

                if beam_color_type == "Material":
                    colorCell = color_palette[cell.beams[0].material % len(color_palette)]
                elif beam_color_type == "Type":
                    colorCell = color_palette[cell.geom_types % len(color_palette)]
                elif beam_color_type == "radii":
                    colorCell = cell.getRGBcolorDependingOfRadius()
                else:
                    colorCell = "green"  # Default color

                x_offset = explode_voxel * (x - latticeDimDict["xMin"]) / dx
                y_offset = explode_voxel * (y - latticeDimDict["yMin"]) / dy
                z_offset = explode_voxel * (z - latticeDimDict["zMin"]) / dz
                self.ax.bar3d(x + x_offset, y + y_offset, z + z_offset,
                              dx, dy, dz, color=colorCell, alpha=0.5, shade=True, edgecolor='k')

        if self.axisSet is False:
            self._setMinMaxAxis(latticeDimDict)

        # Save or show the plot
        if plotting:
            self.show()
        if file_save_path is not None:
            plt.savefig(file_save_path)


    def visualCellZoneBlocker(self, lattice, erasedParts: list[tuple]) -> None:
        """
        Visualize the lattice with erased parts

        Parameters:
        -----------
        eraser_blocks: list of tuple
            List of erased parts with (x_start, y_start, z_start, x_dim, y_dim
        """

        # Plot global lattice cube
        x_max = lattice.xMax
        y_max = lattice.yMax
        z_max = lattice.zMax
        vertices_global = [[0, 0, 0], [x_max, 0, 0], [x_max, y_max, 0], [0, y_max, 0],
                           [0, 0, z_max], [x_max, 0, z_max], [x_max, y_max, z_max], [0, y_max, z_max]]
        self.ax.add_collection3d(
            Poly3DCollection([vertices_global], facecolors='grey', linewidths=1, edgecolors='black', alpha=0.3))

        # Plot erased region cube
        for erased in erasedParts:
            x_start, y_start, z_start, x_dim, y_dim, z_dim = erased
            vertices_erased = [[x_start, y_start, z_start], [x_start + x_dim, y_start, z_start],
                               [x_start + x_dim, y_start + y_dim, z_start], [x_start, y_start + y_dim, z_start],
                               [x_start, y_start, z_start + z_dim], [x_start + x_dim, y_start, z_start + z_dim],
                               [x_start + x_dim, y_start + y_dim, z_start + z_dim],
                               [x_start, y_start + y_dim, z_start + z_dim]]
            self.ax.add_collection3d(
                Poly3DCollection([vertices_erased], facecolors='red', linewidths=1, edgecolors='black', alpha=0.6))

        # Set labels and limits
        self.ax.set_xlabel('X')
        self.ax.set_ylabel('Y')
        self.ax.set_zlabel('Z')
        self.ax.set_xlim([0, x_max])
        self.ax.set_ylim([0, y_max])
        self.ax.set_zlim([0, z_max])

        plt.show()

    def visualizeMesh(self, meshObject):
        """
        Visualize a mesh object in 3D.
        """
        if self.initFig is False:
            self.initFigure()

        if hasattr(meshObject, "mesh"):
            mesh = meshObject.mesh
        else:
            mesh = meshObject

        faces = mesh.vertices[mesh.faces]
        self.ax.add_collection3d(Poly3DCollection(faces, facecolors='cyan', linewidths=0.1, edgecolors='k'))

        minAxis = min(mesh.bounds[0])
        maxAxis = max(mesh.bounds[1])
        self.minAxis = min(minAxis, self.minAxis) if self.minAxis is not None else minAxis
        self.maxAxis = max(maxAxis, self.maxAxis) if self.maxAxis is not None else maxAxis

        # Ajouter un point rouge à l'origine
        self.ax.scatter([0], [0], [0], color='red', s=50, label="Origin (0,0,0)")

    def plotRadiusDistribution(self, cells: list["Cell"], nbRadiusBins: int = 5):
        """
        Plot the radius distribution of beams in the lattice.

        Parameters:
        -----------
        cells: list of Cell
            List of cells to visualize.
        latticeDimDict: dict
            Dictionary containing lattice dimension information (xMin, xMax, yMin, zMin, zMax).
        nbRadiusBins: int
            Number of bins for the histogram.
        """
        all_radii = []
        all_volumes = []

        for cell in cells:
            radius = cell.radius
            if hasattr(radius, '__len__'):
                all_radii.append(radius)
            else:
                all_radii.append([radius])
            all_volumes.append(cell.getVolumeGeomSeparated())

        all_radii = np.array(all_radii)
        dimRadius = all_radii.shape[1]
        all_volumes = np.array(all_volumes)
        sumVolume = np.sum(all_volumes, axis=0)
        ratio_volume = sumVolume / np.sum(sumVolume) * 100

        colors = plt.cm.tab10.colors  # Up to 10 colors predefined

        plt.figure(figsize=(7, 5))
        bins = np.linspace(np.min(all_radii), np.max(all_radii), nbRadiusBins + 1)
        bin_width = (bins[1] - bins[0]) / (dimRadius + 1)

        for i in range(dimRadius):
            shifted_bins = bins[:-1] + i * bin_width
            plt.bar(shifted_bins, np.histogram(all_radii[:, i], bins=bins)[0],
                    width=bin_width, align='edge', color=colors[i % len(colors)], edgecolor='black',
                    label=f'Geometry {i}, Ratio Volume: {ratio_volume[i]:.2f}%')

        plt.title('radii Distribution')
        plt.xlabel('radii')
        plt.ylabel('Count')
        plt.legend()
        plt.show()

    def subplotLatticeGeometries(self, cells: list["Cell"], latticeDimDict: dict, nbRadiusBins: int = 5,
                                 explodeVoxel: float = 0.0):
        """
        Create subplots:
        - One subplot per geometry (radius index) with voxel visualization.
        """
        rmin = 0
        rmax = 0.1

        # Determine number of geometries
        dimRadius = len(cells[0].radius) if hasattr(cells[0].radius, '__len__') else 1
        fig, axs = plt.subplots(1, dimRadius, figsize=(5 * dimRadius, 5), subplot_kw={'projection': '3d'})
        axs = [axs] if dimRadius == 1 else axs  # Ensure axs is always iterable
        for ax in axs:
            ax.set_axis_off()

        for rad in range(dimRadius):
            ax = axs[rad]
            for cell in cells:
                x, y, z = cell.coordinateCell
                dx, dy, dz = cell.cellSize

                # Get color based on the radius value for current geometry
                radius = cell.radius
                radius_value = radius[rad] if hasattr(radius, '__len__') else radius
                import matplotlib.cm as cm  # ajouter en haut si pas encore fait

                # Define the colormap from blue to red
                colormap = cm.get_cmap('coolwarm')

                # Normalize radius between 0 and 1
                radius_norm = (radius_value - rmin) / (rmax - rmin)
                radius_norm = np.clip(radius_norm, 0.0, 1.0)
                colorCell = colormap(radius_norm)

                x_offset = explodeVoxel * (x - latticeDimDict["xMin"]) / dx
                y_offset = explodeVoxel * (y - latticeDimDict["yMin"]) / dy
                z_offset = explodeVoxel * (z - latticeDimDict["zMin"]) / dz

                ax.bar3d(x + x_offset, y + y_offset, z + z_offset, dx, dy, dz,
                         color=colorCell, alpha=0.5, shade=True, edgecolor='k')

            if self.axisSet is False:
                self._setMinMaxAxis(latticeDimDict)

            ax.set_title(f'Geometry {rad}')
            ax.set_xlim3d(self.minAxis, self.maxAxis)
            ax.set_ylim3d(self.minAxis, self.maxAxis)
            ax.set_zlim3d(self.minAxis, self.maxAxis)
            ax.set_box_aspect([1, 1, 1])

        plt.tight_layout()
        plt.show()

    def show(self):
        """
        Show the 3D plot with axis labels and limits.
        """
        self.ax.set_xlabel('X')
        self.ax.set_ylabel('Y')
        self.ax.set_zlabel('Z')
        self.ax.set_xlim3d(self.minAxis, self.maxAxis)
        self.ax.set_ylim3d(self.minAxis, self.maxAxis)
        self.ax.set_zlim3d(self.minAxis, self.maxAxis)
        self.ax.set_box_aspect([1, 1, 1])
        plt.show()
