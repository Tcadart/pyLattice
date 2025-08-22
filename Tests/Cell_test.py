import numpy as np
import pytest
from unittest.mock import MagicMock

from src.pyLattice.cell import Cell
from src.pyLattice.beam import Beam
from src.pyLattice.point import Point


@pytest.fixture
def basic_cell():
    """
    Fixture to create a basic Cell object for testing.
    """
    pos = [0, 0, 0]
    initial_size = [1.0, 1.0, 1.0]
    start_pos = [0.0, 0.0, 0.0]
    lattice_type = ["BCC"]
    radii = [0.05]
    grad_radius = None
    grad_dim = None
    grad_mat = None
    return Cell(pos, initial_size, start_pos, lattice_type, radii, grad_radius, grad_dim, grad_mat)


def test_cell_initialization(basic_cell):
    assert basic_cell.pos == [0, 0, 0]
    assert basic_cell.coordinate == [0.0, 0.0, 0.0]
    assert isinstance(basic_cell.beams, list)
    assert basic_cell.size is not None


def test_get_cell_volume(basic_cell):
    expected_volume = 1.0
    assert basic_cell.volume == pytest.approx(expected_volume)


def test_add_remove_beam(basic_cell):
    p1 = Point(0, 0, 0)
    p2 = Point(1, 0, 0)
    b = Beam(p1, p2, 0.05, 0, 0)
    initial_count = len(basic_cell.beams)
    basic_cell.add_beam(b)
    assert len(basic_cell.beams) == initial_count + 1
    basic_cell.remove_beam(b)
    assert len(basic_cell.beams) == initial_count


def test_cell_center(basic_cell):
    expected_center = [0.5, 0.5, 0.5]
    assert all(np.isclose(a, b) for a, b in zip(basic_cell.center_point, expected_center))


def test_get_cell_corner_coordinates(basic_cell):
    corners = basic_cell.corner_coordinates
    assert len(corners) == 8
    assert all(len(c) == 3 for c in corners)


def test_beam_radius_modification_warning(basic_cell, capsys):
    basic_cell.change_beam_radius([0.05])
    captured = capsys.readouterr()
    assert "WARNING" in captured.out


def test_get_rgb_color(basic_cell):
    rgb = basic_cell.get_RGBcolor_depending_of_radius()
    assert isinstance(rgb, tuple)
    assert all(isinstance(v, float) for v in rgb)
