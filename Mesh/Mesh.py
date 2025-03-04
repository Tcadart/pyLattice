import trimesh

class mesh:
    def __init__(self, meshPath: str = None):
        """
        Initialize a mesh object.
        """
        self.mesh = None
        if meshPath is not None:
            self.loadMesh(meshPath)

    def loadMesh(self, meshPath: str):
        """
        Load a mesh from a file.

        Parameters
        ----------
        meshPath : str
            Path to the mesh file.
        """
        if not meshPath.endswith('.stl'):
            meshPath += '.stl'
        if not meshPath.startswith('Mesh/'):
            meshPath = 'Mesh/' + meshPath
        self.mesh = trimesh.load_mesh(meshPath)

    def saveMesh(self, outputPath: str):
        """
        Save the mesh to a file.

        Parameters
        ----------
        outputPath : str
            Path where the mesh should be saved.
        """
        if not outputPath.endswith('.stl'):
            outputPath += '.stl'
        if not outputPath.startswith('Mesh/'):
            outputPath = 'Mesh/' + outputPath

        self.mesh.export(outputPath)
        print(f"Mesh saved to {outputPath}")


    def scaleMesh(self, scale: float):
        """
        Scale the mesh.

        Parameters
        ----------
        scale : float
            Scale factor.
        """
        self.mesh.apply_scale(scale)
        self.moveMeshToOrigin()

    def moveMeshToOrigin(self):
        translation = -self.mesh.bounds[0]

        # Appliquer la translation
        self.mesh.apply_translation(translation)

    def is_inside_mesh(self, point):
        """
        Check if a point is inside the mesh.
        """
        return self.mesh.contains([point])[0]
