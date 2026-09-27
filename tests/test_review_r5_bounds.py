"""Coastline rules review, round 5: bounds that hold, and inputs refused."""

import hashlib

import numpy as np
import pytest
import shapely
from shapely.geometry import LineString, Polygon, box

from fvcom_mesh_tools import provenance
from fvcom_mesh_tools.patch import (
    base_size_field,
    filter_shoreline_local,
    land_an_element_fits,
    patch_sizing,
    size_lower_bound,
    size_upper_bound,
)
from fvcom_mesh_tools.walls import close_wall_pockets


def _const(v):
    return lambda q: np.full(len(np.atleast_2d(q)), float(v))


def test_lower_bound_covers_the_edges_of_a_thin_box():
    c = np.array([50.0, 0.5])

    def size(q):
        return 30.0 + np.linalg.norm(np.atleast_2d(q)[:, :2] - c, axis=1)

    lo = size_lower_bound(box(0, 0, 100, 1), size, 20.0)
    assert lo <= 30.0
    hi = size_upper_bound(box(0, 0, 100, 1), size, 20.0)
    assert hi >= 30.0 + np.hypot(50.0, 0.5)


def test_pocket_over_a_fine_spot_beside_its_edge_stays_water():
    wall = LineString([(0, 0), (0, 250), (100, 250), (100, 150), (20, 150)])
    c = np.array([48.75, 150.1])

    def size(q):
        return 59.9 + 0.8 * np.linalg.norm(np.atleast_2d(q)[:, :2] - c, axis=1)

    added, _, rep = close_wall_pockets([wall], box(-2000, -2000, 2000, 0), size,
                                       min_h=60.000001, size_floor=30.0)
    assert rep["n_pockets_closed"] == 0 and added.is_empty


@pytest.mark.parametrize("width,y", [(90, 150), (100, 48), (90, 48)])
def test_pocket_area_and_clearance_are_judged_at_the_coarsest_element(width, y):
    wall = LineString([(0, 0), (0, y + 100), (width, y + 100), (width, y), (20, y)])
    _, _, rep = close_wall_pockets([wall], box(-2000, -2000, 2000, 0), _const(100.0),
                                   min_h=60.000001, size_floor=30.0)
    assert rep["n_pockets_closed"] == 0


def test_a_pocket_that_passes_every_test_is_still_closed():
    wall = LineString([(0, 0), (0, 400), (110, 400), (110, 290), (20, 290)])
    _, _, rep = close_wall_pockets([wall], box(-2000, -2000, 2000, 0), _const(100.0),
                                   min_h=60.000001, size_floor=30.0)
    assert rep["n_pockets_closed"] == 1


def _voronoi_mesh():
    xy = np.array([(0, 0), (100, 0), (100, 100), (0, 100),
                   (300, 0), (300, 100), (700, 0), (700, 100)], dtype=float)
    tri = np.array([(0, 1, 3), (1, 2, 3), (1, 4, 2), (4, 5, 2), (4, 6, 5), (6, 7, 5)])
    return xy, tri


def test_nearest_extension_bounds_hold_across_voronoi_jumps():
    xy, tri = _voronoi_mesh()
    f = base_size_field(xy, tri, outside="nearest")
    a, b = f([[49.999999, -10], [50.000001, -10]])
    assert abs(a - b) > 1.0                           # the premise the slope bound needed
    for geom in (box(40, -30, 60, -5), box(-50, -50, 750, -1), box(150, 20, 250, 80),
                 box(250, -40, 400, 30)):
        lo = size_lower_bound(geom, f, 1.0)
        hi = size_upper_bound(geom, f, 1.0)
        x0, y0, x1, y1 = geom.bounds
        gx, gy = np.meshgrid(np.linspace(x0, x1, 301), np.linspace(y0, y1, 301))
        v = f(np.column_stack([gx.ravel(), gy.ravel()]))
        assert lo <= v.min() + 1e-9 and hi >= v.max() - 1e-9


def test_patch_sizing_bounds_hold_with_a_region():
    xy, tri = _voronoi_mesh()
    region = (box(120, 20, 160, 60), 5.0, 80.0)
    h = patch_sizing(xy, tri, [region], distmesh_scale=1.0, outside="nearest")
    for geom in (box(140, -40, 260, -1), box(100, 0, 300, 100), box(500, -60, 650, -5)):
        lo = size_lower_bound(geom, h, 1.0)
        hi = size_upper_bound(geom, h, 1.0)
        x0, y0, x1, y1 = geom.bounds
        gx, gy = np.meshgrid(np.linspace(x0, x1, 301), np.linspace(y0, y1, 301))
        v = h(np.column_stack([gx.ravel(), gy.ravel()]))
        assert lo <= v.min() + 1e-9 and hi >= v.max() - 1e-9


def test_land_guard_reaches_polygons_nested_in_collections():
    raw = Polygon([(0, 0), (100, 0), (100, 100), (0, 100), (0, 0),
                   (-100, 0), (-100, -100), (-200, -100), (-200, 0), (-100, 0), (0, 0)])
    land = shapely.make_valid(raw)
    assert land.geom_type == "GeometryCollection"
    assert len(land_an_element_fits(land, _const(30.0), 30.0)) == 2


@pytest.mark.parametrize("bad", [-10.0, 0.0, np.nan, np.inf])
def test_local_filter_refuses_a_bad_size_field(bad):
    with pytest.raises(ValueError, match="size field"):
        filter_shoreline_local(box(0, 0, 100, 100), _const(bad), 30.0,
                               box(-50, -50, 150, 150))


def test_local_filter_still_filters_with_a_good_field():
    out, _ = filter_shoreline_local(box(0, 0, 100, 100), _const(30.0), 30.0,
                                    box(-50, -50, 150, 150))
    assert out.area == pytest.approx(10000.0, rel=0.05)


def _no_git(monkeypatch):
    monkeypatch.setattr(provenance, "git_state", lambda p: None)


def test_code_state_marks_an_unreadable_source(tmp_path, monkeypatch):
    _no_git(monkeypatch)
    (tmp_path / "a.py").write_text("x = 1\n")
    (tmp_path / "b.py").write_text("y = 2\n")
    real = provenance.file_sha256
    calls = []

    def flaky(f):
        calls.append(f)
        return None if f.name == "b.py" else real(f)

    monkeypatch.setattr(provenance, "file_sha256", flaky)
    st = provenance.code_state("nope-not-installed", tmp_path)
    assert st["source_sha256"] is None
    assert st["source_unreadable"] == [str(tmp_path / "b.py")]
    assert len(calls) == 2                             # read once each


def test_code_state_hashes_readable_sources(tmp_path, monkeypatch):
    _no_git(monkeypatch)
    (tmp_path / "a.py").write_text("x = 1\n")
    st = provenance.code_state("nope-not-installed", tmp_path)
    want = hashlib.sha256()
    want.update(b"a.py")
    want.update(provenance.file_sha256(tmp_path / "a.py").encode())
    assert st["source_sha256"] == want.hexdigest() and st["source_unreadable"] == []


def test_patch_sizing_bound_follows_the_ramp_not_the_target():
    xy, tri = _voronoi_mesh()
    region = (box(120, 20, 160, 60), 5.0, 80.0)
    h = patch_sizing(xy, tri, [region], distmesh_scale=1.0, outside="nearest")
    base = base_size_field(xy, tri, outside="nearest")
    geom = box(200, -60, 260, -40)                    # 60 m from the region, u >= 0.75
    lo_b = size_lower_bound(geom, base, 1.0)
    assert size_lower_bound(geom, h, 1.0) >= 5.0 + (lo_b - 5.0) * 0.75 - 1e-9
