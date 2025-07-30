import json
import os


class MatProperties:
    """
    A class to represent the properties of a material loaded from a file.
    """

    MATERIALS_DIR = "src/materials/"  # Directory containing material files

    def __init__(self, material_file):
        """
        Initialize the MatProperties object by loading data from a file.

        :param material_file: Filename of the material JSON (str)
        """
        self.file_path = os.path.join(self.MATERIALS_DIR, material_file)
        self.name, self.density, self.elastic, self.plastic = self.load_material()

    def load_material(self):
        """
        Loads material properties from a JSON file.

        :return: Material name_lattice, density, elastic properties, and plastic properties
        """
        if not os.path.exists(self.file_path):
            raise FileNotFoundError(f"Material file not found: {self.file_path}")

        with open(self.file_path, "r") as file:
            data = json.load(file)

        return (
            data["name_lattice"],
            data["density"],
            data["elastic"],
            data["plastic"]
        )
