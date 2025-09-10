# `MeshTrimmer` Class

The `MeshTrimmer` class provides functionalities for manipulating 3D meshes (typically in STL format) and trimming 
lattice structures using these meshes. It can be used to load, save, scale, move, and perform geometric checks on 
meshes. Additionally, it allows for cutting lattice beams that intersect with a mesh surface.

---

## 📦 Features

* Load and save `.stl` mesh files
* Automatically reposition mesh to origin
* Check whether a point or a lattice cell is inside the mesh
* Trim beams that intersect the mesh surface
* Use with a lattice generation process for geometrically conforming trimming

---

## 🧰 Dependencies

* [`trimesh`](https://trimsh.org/)
* [`numpy`](https://numpy.org/)
* [`matplotlib`](https://matplotlib.org/)
* Internal dependencies: `Beam`, `Point`, `Cell` class

---

## 🧱 Class Usage

Typically used in conjunction with the `Lattice` class:

```python
from data.inputs.mesh_file.mesh_trimmer import MeshTrimmer
from pyLattice.lattice import Lattice

trimmer = MeshTrimmer("example_mesh.stl")
lattice = Lattice(..., mesh_trimmer=trimmer)
```

---

## 🧩 Methods

### `__init__(mesh_name: str = None)`

Initializes a `MeshTrimmer` object. Optionally loads a mesh at initialization.

---

### `load_mesh(mesh_name: str)`

Loads a mesh from the `Mesh/` directory (STL format). Automatically handles extension and path formatting.

---

### `save_mesh(output_name: str)`

Saves the current mesh to the `Mesh/` directory, with automatic `.stl` extension handling.

---

### `scale_mesh(scale: float)`

Scales the mesh uniformly and moves its minimum corner to the origin.

---

### `move_mesh_to_origin()`

Translates the mesh so that its minimum bounding box corner is at the origin.

---

### `is_inside_mesh(point: array-like) -> bool`

Checks if a given 3D point is located inside the mesh volume.

---

### `is_cell_in_mesh(cell: Cell) -> bool`

Returns `True` if at least one corner of the given `Cell` object is inside the mesh. Used to filter cells during lattice generation.

---

### `cut_beams_at_mesh_intersection(cells: list[Cell])`

Cuts or removes beams in the provided `Cell` list depending on their intersection with the mesh surface.

* Fully outside: removed
* Partially intersecting: shortened
* Inside: unchanged

---

### `find_intersection_with_mesh(beam: Beam) -> tuple[float, float, float] | None`

Uses `pyembree` acceleration to find the first intersection between a beam and the mesh surface. Returns intersection point or `None` if no hit.

---

### `plot_mesh()`

Visualizes the current mesh using `matplotlib` in a 3D view.

* The mesh is displayed with translucent cyan faces and black edges.
* The camera is centered around the mesh centroid, with equal zoom in all directions.
* The 3D axes and background are hidden for a cleaner visualization.

---



## 🧬 Integration with `Lattice`

When passed to the `Lattice` class via the `mesh_trimmer` argument, `MeshTrimmer` is automatically used during cell generation to:

* Filter out cells not intersecting with the mesh
* Optionally trim beams based on the geometry

```python
from data.inputs.mesh_file.mesh_trimmer import MeshTrimmer
from pyLattice.lattice import Lattice

trimmer = MeshTrimmer("example_mesh.stl")
lattice = Lattice(..., mesh_trimmer=trimmer)
lattice.cut_beam_with_mesh_trimmer()  # Use the trimmer to cut beams at mesh intersection
```
---

## 📁 File Structure Convention

* Mesh files are expected to be stored under a `Mesh/` folder.
* STL format is enforced.

---

## 🛠️ Limitations

* Only beams that haven't already been modified (beam_mod == False) can be trimmed, due to the lack of 
  implementation for simulating trimmed lattice structures.