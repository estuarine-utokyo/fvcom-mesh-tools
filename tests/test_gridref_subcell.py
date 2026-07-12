"""Sub-cell ("F12-c3") addressing of the atlas grid."""

import pytest

from fvcom_mesh_tools.gridref import TOKYO_BAY_GRID as G


def test_subcell_roundtrip():
    lon, lat = 139.8043, 35.5259          # Haneda site
    ref = G.point_to_subcell(lon, lat)
    x0, y0, x1, y1 = G.subcell_bounds(ref)
    assert x0 <= lon <= x1 and y0 <= lat <= y1
    # parent consistency
    parent = ref.split("-")[0]
    px0, py0, px1, py1 = G.cell_bounds(parent)
    assert px0 <= x0 and x1 <= px1 and py0 <= y0 and y1 <= py1


def test_subcell_size_is_cell_over_subn():
    x0, y0, x1, y1 = G.subcell_bounds("H8-a1")
    assert abs((x1 - x0) - G.dlon / G.SUB_N) < 1e-12
    assert abs((y1 - y0) - G.dlat / G.SUB_N) < 1e-12
    # a1 = NW corner of the parent
    px0, py0, px1, py1 = G.cell_bounds("H8")
    assert abs(x0 - px0) < 1e-12 and abs(y1 - py1) < 1e-12


def test_cell_bounds_dispatches_sub_refs():
    assert G.cell_bounds("H8-c3") == G.subcell_bounds("H8-c3")


def test_bad_subcell_raises():
    with pytest.raises(ValueError):
        G.subcell_bounds("H8-z9")
    with pytest.raises(ValueError):
        G.subcell_bounds("H8-c0")
