"""
Visualization and saving of lattice structures from lattice objects.

Created in 2025-01-16 by Cadart Thomas, University of technology Belfort-Montbéliard.
"""
from typing import Tuple

import numpy as np
from matplotlib.lines import Line2D
from matplotlib.patches import Patch
from mpl_toolkits.mplot3d.art3d import Line3DCollection, Poly3DCollection
import matplotlib
import matplotlib.colors as mcolors
import matplotlib.pyplot as plt

from .cell import Cell

from .utils import _get_beam_color, _prepare_lattice_plot_data, plot_coordinate_system, get_boundary_condition_color

matplotlib.use('TkAgg')  # Or 'Qt5Agg' if you prefer Qt backend


class LatticePlotting:
    """
    Class for visualizing lattice structures in 3D.
    """

    def __init__(self, initFig: bool = False):
        if initFig:
            self.init_figure()
        self.fig = None
        self.ax = None
        self.minAxis = None
        self.maxAxis = None
        self.initFig = initFig
        self.axisSet = False

    def init_figure(self):
        """Initialize the 3D figure for plotting."""
        self.fig = plt.figure()
        self.ax = self.fig.add_subplot(111, projection='3d')
        self.initFig = True
        self.ax.set_axis_off()

    def _set_min_max_axis(self, latticeDimDict: dict) -> None:
        limMin = min(latticeDimDict["x_min"], latticeDimDict["y_min"], latticeDimDict["z_min"])
        limMax = max(latticeDimDict["x_max"], latticeDimDict["y_max"], latticeDimDict["z_max"])
        self.minAxis = min(limMin, self.minAxis) if self.minAxis is not None else limMin
        self.maxAxis = max(limMax, self.maxAxis) if self.maxAxis is not None else limMax
        self.axisSet = True

    def visualize_lattice(self, lattice_object, beam_color_type: str = "Material",
                             voxelViz: bool = False, deformedForm: bool = False, file_save_path: str = None,
                             plotCellIndex: bool = False, plotNodeIndex: bool = False, explode_voxel: float = 0.0,
                             plotting: bool = True, nbRadiusBins: int = 5,
                             enable_system_coordinates: bool = True, enable_boundary_conditions: bool = False,
                             camera_position: Tuple[float, float] = None) -> None:

        """
        Visualizes the lattice in 3D using matplotlib.

        Parameters:
        -----------
        cells: list of Cell
            List of cells to visualize.
        latticeDimDict: dict
            Dictionary containing lattice dimension information (x_min, x_max, y_min, y_max, z_min, z_max).
        beamColor: str, optional (default: "Material")
            Color scheme for beams. Options:
            - "Material": Color by material.
            - "Type": Color by type_beam.
            - "radii": Color by radii.
        voxelViz: bool, optional (default: False)
            If True, visualize as voxels; otherwise, use beam visualization.
        deformedForm: bool, optional (default: False)
            If True, use deformed node positions.
        nameSave: str, optional
            If provided, save the plot with this name_lattice.
        plotCellIndex: bool, optional (default: False)
            If True, plot cell indices.
        plotNodeIndex: bool, optional (default: False)
            If True, plot node indices.
        explode_voxel: float, optional (default: 0.0)
            If greater than 0, apply an explosion effect to the voxel visualization.
        plotting: bool, optional (default: True)
            If True, display the plot after creation.
        nbRadiusBins: int, optional (default: 5)
            Number of bins for the histogram of radii distribution.
        enable_system_coordinates: bool, optional (default: True)
            If True, plot the coordinate system axes.
        enable_boundary_conditions: bool, optional (default: False)
            If True, visualize boundary conditions on nodes.
        camera_position: Tuple[float, float], optional
            If provided, set the camera position for the 3D plot as (elevation, azimuth).
        """
        if self.initFig is False:
            self.init_figure()

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
        legend_map = {}  # {nb_fixed: color}

        if not voxelViz:
            beamDraw = set()
            lines = []
            colors = []
            nodeX, nodeY, nodeZ = [], [], []
            nodeDraw = set()

            for cell in cells:
                for beam in cell.beams:
                    if beam.radius != 0.0 and beam not in beamDraw:
                        colorBeam, idxColor = _get_beam_color(beam, color_palette, beam_color_type, idxColor, cells,
                                                              nbRadiusBins)

                        # Add line and node data
                        beam_lines, beam_nodes, beam_indices, = _prepare_lattice_plot_data(beam, deformedForm)
                        lines.extend(beam_lines)
                        colors.extend([colorBeam] * len(beam_lines))  # One color per line

                        for i, node in enumerate(beam_nodes):
                            if node not in nodeDraw:
                                nodeDraw.add(node)
                                nodeX.append(node.x if not deformedForm else node.deformed_coordinates[0])
                                nodeY.append(node.y if not deformedForm else node.deformed_coordinates[1])
                                nodeZ.append(node.z if not deformedForm else node.deformed_coordinates[2])
                                if plotNodeIndex:
                                    self.ax.text(nodeX[-1], nodeY[-1], nodeZ[-1], str(beam_indices[i]), fontsize=5,
                                                 color='gray')
                                if enable_boundary_conditions:
                                    if any(node.fixed_DOF):
                                        bc_color = get_boundary_condition_color(node.fixed_DOF)
                                        nb_fixed = sum(node.fixed_DOF)
                                        legend_map.setdefault(nb_fixed, bc_color)

                                        self.ax.scatter(nodeX[-1], nodeY[-1], nodeZ[-1], c=bc_color, s=70)
                                        # Apply boundary condition labels
                                        for i, (is_fixed, d_val) in enumerate(zip(node.fixed_DOF, node.displacement_vector)):
                                            if is_fixed and abs(d_val) > 1e-10:
                                                self.ax.text(nodeX[-1] + cell.size[0]/4, nodeY[-1], nodeZ[-1] + cell.size[2]/4,
                                                             f"u{i}={d_val:.2e}", fontsize=10, color=bc_color)
                                        if deformedForm:
                                            # plot initial position of nodes of boundary conditions
                                            self.ax.scatter(node.x, node.y, node.z, facecolor='none', edgecolor = 'k',
                                                            s=70, label="Initial Position")



                        beamDraw.add(beam)

                if plotCellIndex:
                    self.ax.text(cell.center_point[0], cell.center_point[1], cell.center_point[2], str(cell.index),
                                 color='black', fontsize=10)

            # Plot lines and nodes
            line_collection = Line3DCollection(lines, colors=colors, linewidths=2)
            self.ax.add_collection3d(line_collection)
            self.ax.scatter(nodeX, nodeY, nodeZ, c='black', s=5)

        else:  # Voxel visualization
            for cell in cells:
                x, y, z = cell.coordinate
                dx, dy, dz = cell.size

                beam_color_type = beam_color_type.lower()
                if beam_color_type == "material":
                    colorCell = color_palette[cell.beams[0].material % len(color_palette)]
                elif beam_color_type == "type":
                    colorCell = color_palette[cell.geom_types % len(color_palette)]
                elif beam_color_type == "radii":
                    colorCell = cell.get_RGBcolor_depending_of_radius()
                elif beam_color_type == "random":
                    rng = np.random.default_rng()
                    colorCell = rng.random(3,)
                else:
                    colorCell = "green"  # Default color

                x_offset = explode_voxel * (x - latticeDimDict["x_min"]) / dx
                y_offset = explode_voxel * (y - latticeDimDict["y_min"]) / dy
                z_offset = explode_voxel * (z - latticeDimDict["z_min"]) / dz
                self.ax.bar3d(x + x_offset, y + y_offset, z + z_offset,
                              dx, dy, dz, color=colorCell, alpha=0.5, shade=True, edgecolor='k')

        if self.axisSet is False:
            self._set_min_max_axis(latticeDimDict)



        if enable_system_coordinates:
            plot_coordinate_system(self.ax)

        if enable_boundary_conditions and legend_map:
            legend_elements = [Patch(facecolor=color, label=f"{n} DOF")
                               for n, color in sorted(legend_map.items())]

            legend_elements.append(
                Line2D([0], [0], marker='o', color='k', label="Initial Position",
                       markerfacecolor='none', markersize=8, linestyle='None')
            )

            self.ax.legend(handles=legend_elements, title="Boundary Conditions", loc='upper right')

        if camera_position is not None:
            self.ax.view_init(elev=camera_position[0], azim=camera_position[1])

        # Save or show the plot
        if plotting:
            self.show()
        if file_save_path is not None:
            plt.savefig(file_save_path)



    def visual_cell_zone_blocker(self, lattice, erasedParts: list[tuple]) -> None:
        """
        Visualize the lattice with erased parts

        Parameters:
        -----------
        eraser_blocks: list of tuple
            List of erased parts with (x_start, y_start, z_start, x_dim, y_dim
        """

        # Plot global lattice cube
        x_max = lattice.x_max
        y_max = lattice.y_max
        z_max = lattice.z_max
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

    def plot_radius_distribution(self, cells: list["Cell"], nbRadiusBins: int = 5):
        """
        Plot the radii distribution of beams in the lattice.

        Parameters:
        -----------
        cells: list of Cell
            List of cells to visualize.
        latticeDimDict: dict
            Dictionary containing lattice dimension information (x_min, x_max, y_min, z_min, z_max).
        nbRadiusBins: int
            Number of bins for the histogram.
        """
        all_radii = []
        all_volumes = []

        for cell in cells:
            radius = cell.radii
            if hasattr(radius, '__len__'):
                all_radii.append(radius)
            else:
                all_radii.append([radius])
            all_volumes.append(cell.volume_each_geom)

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

    def subplot_lattice_geometries(self, cells: list["Cell"], latticeDimDict: dict, nbRadiusBins: int = 5,
                                   explodeVoxel: float = 0.0):
        """
        Create subplots:
        - One subplot per geometry (radii index) with voxel visualization.
        """
        rmin = 0
        rmax = 0.1

        # Determine number of geometries
        dimRadius = len(cells[0].radii) if hasattr(cells[0].radii, '__len__') else 1
        fig, axs = plt.subplots(1, dimRadius, figsize=(5 * dimRadius, 5), subplot_kw={'projection': '3d'})
        axs = [axs] if dimRadius == 1 else axs  # Ensure axs is always iterable
        for ax in axs:
            ax.set_axis_off()

        for rad in range(dimRadius):
            ax = axs[rad]
            for cell in cells:
                x, y, z = cell.coordinate
                dx, dy, dz = cell.size

                # Get color based on the radii value for current geometry
                radius = cell.radii
                radius_value = radius[rad] if hasattr(radius, '__len__') else radius
                import matplotlib.cm as cm  # ajouter en haut si pas encore fait

                # Define the colormap from blue to red
                colormap = cm.get_cmap('coolwarm')

                # Normalize radii between 0 and 1
                radius_norm = (radius_value - rmin) / (rmax - rmin)
                radius_norm = np.clip(radius_norm, 0.0, 1.0)
                colorCell = colormap(radius_norm)

                x_offset = explodeVoxel * (x - latticeDimDict["x_min"]) / dx
                y_offset = explodeVoxel * (y - latticeDimDict["y_min"]) / dy
                z_offset = explodeVoxel * (z - latticeDimDict["z_min"]) / dz

                ax.bar3d(x + x_offset, y + y_offset, z + z_offset, dx, dy, dz,
                         color=colorCell, alpha=0.5, shade=True, edgecolor='k')

            if self.axisSet is False:
                self._set_min_max_axis(latticeDimDict)

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
