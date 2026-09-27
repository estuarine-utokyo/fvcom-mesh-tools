"""Coastline rules review, round 10."""

import ast

import numpy as np
from shapely.geometry import Point

from fvcom_mesh_tools.patch import _rim_edit_ok, hole_polygon, rim_repair
from tests.test_local_refine_driver import DRIVER


def _const(v):
    return lambda q: np.full(len(np.atleast_2d(q)), float(v))


def _cycle(n, k0=0):
    return np.array([(k0 + k, k0 + (k + 1) % n) for k in range(n)])


def test_a_tip_does_not_step_back_over_a_lake():
    outer = [(0, 0), (190, 0), (190, 297), (200, 297), (210, 297),
             (210, 0), (600, 0), (600, 300), (0, 300)]
    lake = [(191, 280), (193, 280), (193, 282), (191, 282)]
    p = np.array(outer + lake, dtype=float)
    e = np.vstack([_cycle(9), _cycle(4, 9)])
    b = np.arange(len(p))
    b[2] = -1
    water = hole_polygon(p, e)
    assert water.contains(Point(192, 281))
    out, eg, _b, _m, rep = rim_repair(p, e, b, water, _const(30.0), operations=("slits",),
                                      rounds=1, focus=[[190, 297]], focus_factor=0.001)
    assert hole_polygon(out, eg).contains(Point(192, 281))
    assert _rim_edit_ok(p, e, out, eg)


def _nested(name):
    tree = ast.parse(DRIVER.read_text())
    fn = next(n for n in ast.walk(tree) if isinstance(n, ast.FunctionDef) and n.name == name)
    return ast.Module(body=[fn], type_ignores=[])


def _blunt_env(xy, eg, base, pinned=()):
    env = {"np": np, "rim_xy": np.asarray(xy, dtype=float),
           "rim_eg": np.asarray(eg, dtype=np.int64), "rim_base": np.asarray(base),
           "_segs": [], "_folded": [], "PINNED_XY": set(pinned),
           "_rim_edit_ok": _rim_edit_ok}
    for name in ("_blunt", "_pin"):
        exec(compile(_nested(name), "<driver>", "exec"), env)
    return env


def test_a_wall_root_fold_does_not_collapse_a_triangle():
    env = _blunt_env([(0, 0), (10, 0), (0, 10)], _cycle(3), np.full(3, -1))
    env["_blunt"](0, 30.0)
    hole = hole_polygon(env["rim_xy"], env["rim_eg"])
    assert hole.area > 0


def test_a_fold_still_merges_a_short_edge_and_keeps_its_pin():
    xy = [(0, 0), (10, 0), (300, 0), (300, 300), (0, 300)]
    env = _blunt_env(xy, _cycle(5), np.full(5, -1), pinned={(10.0, 0.0)})
    env["_blunt"](0, 30.0)
    assert env["_folded"] == [(1, 0)]
    assert (5.0, 0.0) in env["PINNED_XY"]
