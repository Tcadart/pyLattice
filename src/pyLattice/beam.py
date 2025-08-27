"""
# beam.py
"""

from typing import List, Tuple, Optional
import math

from .point import Point
from .utils import function_penalization_Lzone

class Beam(object):
    """
    Class Beam represents a beam element in lattice structures.
    """

    def __init__(self, point1: 'Point', point2: 'Point', radius: float, material: int, type_beam: int) -> None:
        """
        Initialize a Beam object representing a beam element.

        Args:
            point1 (Point): First endpoint of the beam.
            point2 (Point): Second endpoint of the beam.
            radius (float): radii of the beam.
            material (int): Material identifier of the beam.
            type_beam (int): Type of the beam (0: normal, 1: modified, 2: boundary beam).
        """
        self.point1: 'Point' = point1
        self.point2: 'Point' = point2
        self.radius: float = radius
        self.material: int = material
        self.type_beam: int = type_beam
        self.index: Optional[int] = None
        self.angle_point_1: dict = {"radius": None, "angle": None, "L_zone": None}
        self.angle_point_2: dict = {"radius": None, "angle": None, "L_zone": None}
        self.length: float = self.get_length()
        self.volume: float = self.get_volume(sectionType="circular")
        self.beam_mod: bool = False
        self.penalization_coefficient: float = 1.5  # Fixed with previous optimization
        self.initial_radius: Optional[float] = None

    def __repr__(self) -> str:
        return f"Beam({self.point1}, {self.point2}, radii:{self.radius}, Type:{self.type_beam}, Index:{self.index})"

    def __eq__(self, other: object) -> bool:
        return isinstance(other, Beam) and self.point1 == other.point1 and self.point2 == other.point2

    def __hash__(self) -> int:
        return hash((self.point1, self.point2))

    @property
    def data(self) -> List[int]:
        """
        Property to retrieve beam data for exporting.

        Returns:
            List[int]: [beam_index, point1_index, point2_index, beam_type].
        """
        return [self.index, self.point1.index, self.point2.index, self.type_beam]

    def get_length(self) -> float:
        """
        Calculate the length of the beam.

        Returns:
            float: Length of the beam.
        """
        x1, y1, z1 = self.point1.x, self.point1.y, self.point1.z
        x2, y2, z2 = self.point2.x, self.point2.y, self.point2.z
        length = round(math.sqrt((x2 - x1) ** 2 + (y2 - y1) ** 2 + (z2 - z1) ** 2), 4)
        return length

    def get_volume(self, sectionType: str = "circular") -> float:
        """
        Calculate the volume of the beam in case of a circular section.

        Args:
            sectionType (str): Type of the beam section. Currently only "circular" is supported.
        Returns:
            float: Volume of the beam.
        """
        if sectionType.lower() != "circular":
            raise ValueError("Currently only circular section type_beam is supported.")
        return math.pi * (self.radius ** 2) * self.length

    def get_angle_between_beams(self, other: 'Beam', periodicity: bool) -> float:
        """
        Calculates the angle between two beams

        Parameters:
        ----------
        other : Beam
            The other beam to calculate the angle with.
        periodicity : bool
            If True, considers periodic boundary conditions.

        Return:
        --------
        Angle: float
            angle in degrees
        """
        p1, p2 = None, None
        if periodicity:
            # Tag for corners (1000-1007)
            for idx1, p1_candidate in enumerate([self.point1, self.point2]):
                if p1_candidate.tag and 1000 <= p1_candidate.tag[0] <= 1007:
                    p1 = idx1
                    for idx2, p2_candidate in enumerate([other.point1, other.point2]):
                        if p2_candidate.tag and 1000 <= p2_candidate.tag[0] <= 1007:
                            p2 = idx2  # Enregistrer l'indice pour beam2
                            break

            # Tags for edges (100-111)
            if p1 is None and p2 is None:
                list_tag_edge = [[102, 104, 106, 107], [100, 108, 105, 111], [101, 109, 103, 110]]
                for tag_list in list_tag_edge:
                    for idx1, p1_candidate in enumerate([self.point1, self.point2]):
                        if p1_candidate.tag and p1_candidate.tag[0] in tag_list:
                            p1 = idx1
                            for idx2, p2_candidate in enumerate([other.point1, other.point2]):
                                if p2_candidate.tag and p2_candidate.tag[0] in tag_list:
                                    p2 = idx2
                                    break
            # Tags for faces (10-15)
            if p1 is None and p2 is None:
                list_face_tag = [[10, 15], [11, 14], [12, 13]]
                for face_tag in list_face_tag:
                    for idx1, p1_candidate in enumerate([self.point1, self.point2]):
                        if p1_candidate.tag and p1_candidate.tag[0] in face_tag:
                            p1 = idx1
                            for idx2, p2_candidate in enumerate([other.point1, other.point2]):
                                if p2_candidate.tag and p2_candidate.tag[0] in face_tag:
                                    p2 = idx2
                                    break

        if self.point1 == other.point1 or (p1 == 0 and p2 == 0):
            u = self.point2 - self.point1
            v = other.point2 - other.point1
        elif self.point1 == other.point2 or (p1 == 0 and p2 == 1):
            u = self.point2 - self.point1
            v = other.point1 - other.point2
        elif self.point2 == other.point1 or (p1 == 1 and p2 == 0):
            u = self.point1 - self.point2
            v = other.point2 - other.point1
        elif self.point2 == other.point2 or (p1 == 1 and p2 == 1):
            u = self.point1 - self.point2
            v = other.point1 - other.point2
        else:
            raise ValueError("beams are not connected at any point")

        dot_product = sum(a * b for a, b in zip(u, v))
        u_norm = math.sqrt(sum(a * a for a in u))
        v_norm = math.sqrt(sum(b * b for b in v))
        cos_theta = dot_product / (u_norm * v_norm)
        cos_theta = max(min(cos_theta, 1.0), -1.0)
        angle_rad = math.acos(cos_theta)
        angle_deg = math.degrees(angle_rad)
        return angle_deg

    def get_point_on_beam_at_distance(self, distance: float, start_point: int) -> "point":
        """
        Calculate the coordinates of a point on the beam at a specific distance from an endpoint.

        Args:
            distance (float): Distance from the specified endpoint.
            start_point (int): Index of the endpoint (1 for point1, 2 for point2).

        Returns:
            List[float]: Coordinates [x, y, z] of the calculated point.

        Raises:
            ValueError: If the point index is not 1 or 2.
        """
        if start_point == 1:
            start_point = self.point1
            end_point = self.point2
        elif start_point == 2:
            start_point = self.point2
            end_point = self.point1
        else:
            raise ValueError("point must be 1 or 2.")

        direction_ratios = [
            (end_point.x - start_point.x) / self.length,
            (end_point.y - start_point.y) / self.length,
            (end_point.z - start_point.z) / self.length,
        ]

        factors = [dr * distance for dr in direction_ratios]

        point_mod = [
            start_point.x + factors[0],
            start_point.y + factors[1],
            start_point.z + factors[2]
        ]

        point_mod = Point(*point_mod)

        return point_mod

    def is_point_on_beam(self, node: 'point') -> bool:
        """
        Check if a given node lies on the beam.

        Args:
            node (Point): The point to check.

        Returns:
            bool: True if the node lies on the beam, False otherwise.
        """
        vector1 = (self.point2.x - self.point1.x, self.point2.y - self.point1.y, self.point2.z - self.point1.z)
        vector2 = (node.x - self.point1.x, node.y - self.point1.y, node.z - self.point1.z)

        if (node.x == self.point1.x and node.y == self.point1.y and node.z == self.point1.z) or (
                node.x == self.point2.x and node.y == self.point2.y and node.z == self.point2.z):
            return False
        cross_product = (
            vector1[1] * vector2[2] - vector1[2] * vector2[1],
            vector1[2] * vector2[0] - vector1[0] * vector2[2],
            vector1[0] * vector2[1] - vector1[1] * vector2[0]
        )

        if cross_product == (0, 0, 0):
            dot_product = (vector2[0] * vector1[0] + vector2[1] * vector1[1] + vector2[2] * vector1[2])
            vector1_length_squared = (vector1[0] ** 2 + vector1[1] ** 2 + vector1[2] ** 2)
            return 0 <= dot_product <= vector1_length_squared
        else:
            return False

    def set_angle(self, radius: float, angle: float, point: "Point") -> None:
        """
        Assign angle and radius data to one of the beam's endpoints.

        Parameters:
        ----------
        radius : float
            Radius at the point.
        angle : float
            Angle at the point in degrees.
        point : Point
            The point (endpoint) of the beam to which the data is assigned.
        """
        if point == self.point1:
            self.angle_point_1["radius"] = radius
            self.angle_point_1["angle"] = angle
            self.angle_point_1["L_zone"] = function_penalization_Lzone(radius, angle)
        elif point == self.point2:
            self.angle_point_2["radius"] = radius
            self.angle_point_2["angle"] = angle
            self.angle_point_2["L_zone"] = function_penalization_Lzone(radius, angle)
        else:
            raise ValueError("The specified point is not an endpoint of the beam.")

    def get_length_mod(self) -> Tuple[float, float]:
        """
        Calculate the modification length for the penalization method.

        Returns:
            Tuple[float, float]: Length modifications for point1 and point2.
        """
        L1 = self.angle_point_1["L_zone"]
        L2 = self.angle_point_2["L_zone"]
        return L1, L2

    def set_beam_mod(self):
        """
        Set the beam as modified.
        """
        self.beam_mod = True
        self.initial_radius = self.radius
        self.radius *= self.penalization_coefficient

    def change_beam_radius(self, new_radius: float):
        """
        Change the radius of the beam.

        Args:
            new_radius (float): The new radius to set for the beam.
        """
        if new_radius <= 0:
            raise ValueError("Radius must be a positive value.")
        if self.beam_mod:
            self.radius = new_radius * self.penalization_coefficient
        else:
            self.radius = new_radius

    def is_identical_to(self, other: "Beam", cell_size: list) -> bool:
        """
        Check if this beam is identical to another beam.

        Parameters
        ----------
        other : Beam
            The other beam to compare with.
        cell_size : list
            cell size [cell_size_X, cell_size_Y, cell_size_Z].
        """
        lengthtest = math.isclose(self.length, other.length, rel_tol=1e-5)
        radiustest = math.isclose(self.radius, other.radius, rel_tol=1e-5)
        point1test = self.point1.is_identical_to(other.point1, cell_size)
        point2test = self.point2.is_identical_to(other.point2, cell_size)
        materialtest = self.material == other.material
        typetest = self.type_beam == other.type_beam
        return (
                lengthtest
                and radiustest
                and point1test
                and point2test
                and materialtest
                and typetest
        )
