"""Coastline rules review, round 12."""

import numpy as np
from shapely.geometry import Point, Polygon, box

from fvcom_mesh_tools.patch import hole_polygon, island_rings, rim_repair


def _const(v):
    return lambda q: np.full(len(np.atleast_2d(q)), float(v))


def _with_rings(water, rings):
    p = np.asarray(water.exterior.coords, dtype=float)[:-1]
    e = [np.array([(k, (k + 1) % len(p)) for k in range(len(p))])]
    for r in rings:
        k0 = len(p)
        p = np.vstack([p, r])
        e.append(np.array([(k0 + k, k0 + (k + 1) % len(r)) for k in range(len(r))]))
    return hole_polygon(p, np.vstack(e))


def test_an_island_keeps_its_lake():
    water = box(0, 0, 1000, 1000)
    land = Polygon(box(100, 100, 900, 900).exterior.coords,
                   [box(300, 300, 700, 700).exterior.coords])
    rings, rep = island_rings(land, water, _const(30.0))
    assert rep["n_islands_added"] == 1 and rep["n_lakes_added"] == 1
    after = _with_rings(water, rings)
    assert after.contains(Point(500, 500)) and not after.contains(Point(200, 200))
    assert abs(after.area - 520000.0) < 2000.0


def test_a_lake_too_close_to_its_island_coast_is_reported():
    land = Polygon(box(100, 100, 900, 900).exterior.coords,
                   [box(105, 300, 300, 700).exterior.coords])
    rings, rep = island_rings(land, box(0, 0, 1000, 1000), _const(30.0))
    assert rep["n_lakes_added"] == 0
    assert any("lake" in s["why"] for s in rep["skipped"])


def test_an_edge_refused_before_is_tried_again_after_the_rim_changes():
    p = np.array([[36.1, 23.5], [24.0, 46.8], [5.1, 76.2], [-7.4, 19.5],
                  [-9.2, 23.5], [-74.3, -12.9], [4.9, -35.2], [58.2, -5.4]])
    e = np.array([(k, (k + 1) % 8) for k in range(8)])
    out, eg, _b, _m, rep = rim_repair(p, e, np.full(8, -1), Polygon(p), _const(30.0))
    L = np.linalg.norm(out[eg[:, 0]] - out[eg[:, 1]], axis=1)
    assert (L >= 15.0).all()


def test_short_edges_left_are_read_off_the_returned_rim():
    p = np.array([(0.0, 0.0), (10.0, 0.0), (100.0, 20.0)])
    e = np.array([(0, 1), (1, 2), (2, 0)])
    out, eg, _b, _m, rep = rim_repair(p, e, np.full(3, -1), Polygon(p), _const(30.0))
    L = np.linalg.norm(out[eg[:, 0]] - out[eg[:, 1]], axis=1)
    assert rep["n_short_edges_left"] == int((L < 15.0).sum())
    for s in rep["short_edges_left"]:
        assert any(abs(s["length_m"] - x) < 0.01 for x in L)
