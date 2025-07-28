import sys
import os

sys.path.append(os.path.abspath(os.path.join(os.path.dirname(__file__), "../")))

import pytest

from src.point import Point
from src.beam import Beam


def test_beam_length():
    point1 = Point(0, 0, 0)
    point2 = Point(1, 0, 0)
    beam1 = Beam(point1, point2, 0.05, 1, 0)

    assert beam1.length == 1.0


def test_beam_equality():
    point1 = Point(0, 0, 0)
    point2 = Point(1, 0, 0)
    beam1 = Beam(point1, point2, 0.05, 1, 0)

    point3 = Point(0, 0, 0)
    point4 = Point(1, 0, 0)
    beam2 = Beam(point3, point4, 0.05, 1, 0)

    assert beam1.is_identical_to(beam2)

