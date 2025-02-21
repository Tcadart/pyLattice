import trimesh
import matplotlib.pyplot as plt

mesh = trimesh.load_mesh("Mesh/bike-helmet.stl")
mesh.apply_scale(0.1)

from mpl_toolkits.mplot3d.art3d import Poly3DCollection

# Extraire les faces
fig = plt.figure()
ax = fig.add_subplot(111, projection="3d")

# Obtenir les sommets des faces
faces = mesh.vertices[mesh.faces]
ax.add_collection3d(Poly3DCollection(faces, facecolors='cyan', linewidths=0.2, edgecolors='k', alpha=0.5))

# Définir les limites des axes
minAxis = min(mesh.bounds[0])
maxAxis = max(mesh.bounds[1])
ax.set_xlim([minAxis, maxAxis])
ax.set_ylim([minAxis, maxAxis])
ax.set_zlim([minAxis, maxAxis])
ax.set_box_aspect([1, 1, 1])  # Égaliser les axes X, Y et Z


plt.show()
