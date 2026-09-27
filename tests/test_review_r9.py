"""Coastline rules review, round 9."""

import numpy as np
from shapely.geometry import GeometryCollection, LineString, Point, Polygon, box

from fvcom_mesh_tools.patch import blunt_acute_corners, hole_polygon, rim_repair
from fvcom_mesh_tools.walls import close_wall_pockets
from tests.test_local_refine_driver import DRIVER, driver_function


def _const(v):
    return lambda q: np.full(len(np.atleast_2d(q)), float(v))


def _cycle(n, k0=0):
    return np.array([(k0 + k, k0 + (k + 1) % n) for k in range(n)])


def test_the_first_blunting_pins_what_it_changed():
    p = np.array([(0, 0), (0, 600), (-600, 900), (-600, 600)], dtype=float)
    e = _cycle(4)
    b = np.array([-1, 1, 2, 3])
    p2, e2, _b2, rep = blunt_acute_corners(p, e, b, Polygon(p), _const(30.0))
    assert rep["n_corners_blunted"] == 1
    env = {"np": np}
    driver_function("_nbr_xy", env)
    changed = driver_function("_rim_changed", env)(p, e, p2, e2)
    got = {tuple(np.round(p2[k], 3)) for k in changed}
    assert (0.0, 30.0) in got and (-21.213, 21.213) in got
    src = DRIVER.read_text()
    assert "for k in _rim_changed(rc[\"pfix\"], rc[\"egfix\"], _p, _e)}" in src
    assert "in PINNED_XY]" in src


def _outer_and_island():
    outer = np.array([(0, 0), (200, 600), (0, 900), (-200, 600)], dtype=float)
    island = np.array([(-1, 10), (1, 10), (1, 12), (-1, 12)], dtype=float)
    p = np.vstack([outer, island])
    e = np.vstack([_cycle(4), _cycle(4, 4)])
    b = np.array([-1, -1, -1, -1, 4, 5, 6, 7])
    return p, e, b


def test_blunting_does_not_strand_an_island_outside_its_shell():
    p, e, b = _outer_and_island()
    water = hole_polygon(p, e)
    assert not water.contains(Point(0, 11))
    for kw in ({}, {"operations": ("angles",)}):
        out, eg, _b, _m, _rep = rim_repair(p, e, b, water, _const(30.0), **kw)
        after = hole_polygon(out, eg)
        assert not after.contains(Point(0, 11))


def test_blunting_still_cuts_a_corner_with_no_island_near():
    p = np.array([(0, 0), (200, 600), (0, 900), (-200, 600)], dtype=float)
    out, eg, _b, rep = blunt_acute_corners(p, _cycle(4), np.full(4, -1), Polygon(p),
                                           _const(30.0))
    assert rep["n_corners_blunted"] >= 1


def test_pockets_take_the_polygons_of_a_mixed_land_collection():
    wall = LineString([(0, 0), (0, 400), (110, 400), (110, 290), (20, 290)])
    land = GeometryCollection([box(-2000, -2000, 2000, 0),
                               LineString([(1000, 1000), (1100, 1000)])])
    _, _, rep = close_wall_pockets([wall], land, _const(100.0), min_h=60.0, size_floor=30.0)
    _, _, ref = close_wall_pockets([wall], box(-2000, -2000, 2000, 0), _const(100.0),
                                   min_h=60.0, size_floor=30.0)
    assert rep["n_pockets_closed"] == ref["n_pockets_closed"] == 1
