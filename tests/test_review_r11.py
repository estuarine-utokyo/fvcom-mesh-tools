"""Coastline rules review, round 11: one rim-edit guard for every path."""

import numpy as np
from shapely.geometry import Point

from fvcom_mesh_tools.patch import _rim_edit_ok, hole_polygon
from tests.test_review_r10 import _blunt_env, _nested


def _cycle(n, k0=0):
    return np.array([(k0 + k, k0 + (k + 1) % n) for k in range(n)])


def _rim(*rings):
    p = np.vstack([np.asarray(r, dtype=float) for r in rings])
    e, k = [], 0
    for r in rings:
        e.append(_cycle(len(r), k))
        k += len(r)
    return p, np.vstack(e)


def _root_env(p, e, base):
    env = _blunt_env(p, e, base)
    env.update(h_achieved=lambda q: np.full(len(np.atleast_2d(q)), 30.0), _moved=[])
    exec(compile(_nested("_root"), "<driver>", "exec"), env)
    return env


def test_a_root_move_that_strands_an_island_is_refused():
    outer = [(100, -5), (300, 0), (300, 300), (-100, 300), (-100, 0)]
    island = [(99, -4.95), (101, -4.95), (101, -4.9), (99, -4.9)]
    p, e = _rim(outer, island)
    base = np.arange(len(p))
    base[0] = -1
    before = hole_polygon(p, e)
    env = _root_env(p, e, base)
    env["_root"](np.array([110.0, -4.75]))
    after = hole_polygon(env["rim_xy"], env["rim_eg"])
    assert before.contains(Point(100, -4.925)) == after.contains(Point(100, -4.925))


def test_a_fold_that_crosses_an_island_is_refused():
    p, e = _rim([(0, 0), (10, 0), (300, 0), (300, 300), (0, 300)],
                [(6, 10), (6, 20), (1, 20), (1, 10)])
    env = _blunt_env(p, e, np.full(len(p), -1))
    env["_blunt"](0, 30.0)
    hole_polygon(env["rim_xy"], env["rim_eg"])       # raises if the rings touch


def test_a_fold_that_swaps_an_island_and_a_lake_is_refused():
    p, e = _rim([(0, 0), (10, 10), (300, 0), (300, 300), (0, 300)],
                [(1, 20), (2, 20), (2, 21), (1, 21)],
                [(8, 6), (9, 6), (9, 7), (8, 7)])
    before = hole_polygon(p, e)
    env = _blunt_env(p, e, np.full(len(p), -1))
    env["_blunt"](0, 30.0)
    after = hole_polygon(env["rim_xy"], env["rim_eg"])
    for q in (Point(1.5, 20.5), Point(8.5, 6.5)):
        assert before.contains(q) == after.contains(q)


def test_a_fold_pins_the_root_it_makes():
    p, e = _rim([(0, 0), (10, 0), (300, 0), (300, 300), (0, 100)])
    env = _blunt_env(p, e, np.array([-1, -1, 2, 3, 4]))
    env["_blunt"](0, 30.0)
    assert env["_folded"] == [(1, 0)]
    assert (5.0, 0.0) in env["PINNED_XY"] and (300.0, 0.0) in env["PINNED_XY"]


def test_the_guard_accepts_a_plain_edit_and_refuses_a_swap():
    p, e = _rim([(0, 0), (300, 0), (300, 300), (0, 300)], [(100, 100), (110, 100), (110, 110)])
    moved = p.copy()
    moved[1] = (310, 0)
    assert _rim_edit_ok(p, e, moved, e)
    out = p.copy()
    out[4:] += (500, 0)                               # the island leaves its shell
    assert not _rim_edit_ok(p, e, out, e)
