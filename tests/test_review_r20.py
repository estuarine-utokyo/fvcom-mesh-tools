"""Coastline rules review, round 20."""

import numpy as np
from shapely.affinity import translate
from shapely.geometry import Polygon, box

from fvcom_mesh_tools.patch import island_rings
from tests.test_review_r12 import _with_rings


def test_lakes_needing_different_outlines_all_survive_in_either_order():
    a = box(100, 100, 900, 900).difference(box(100, 100, 110, 500))
    b = box(-100, 200, 94, 450)
    c = translate(box(100, 100, 900, 900).union(
        Polygon([(100, 485), (90, 500), (100, 515)])), yoff=1500)
    d = translate(box(-100, 475, 80, 525), yoff=1500)
    shell = box(-300, -100, 1100, 2600)
    water = box(-1000, -1000, 2000, 3500)
    for lakes in ([a, b, c, d], [b, a, d, c]):
        land = Polygon(shell.exterior.coords, [lk.exterior.coords for lk in lakes])
        rings, rep = island_rings(land, water, 30.0)
        assert rep["n_lakes_added"] == 4, rep["skipped"]
        after = _with_rings(water, rings)
        assert after.intersection(b).area > 48000.0


def test_a_small_rim_island_is_not_swallowed_by_a_large_shell():
    outer = box(0, 0, 22000, 22000)
    water = outer.difference(box(10995, 10995, 11005, 11005))
    land = Polygon(box(1000, 1000, 21000, 21000).exterior.coords,
                   [box(10994, 10994, 11006, 11006).exterior.coords])
    rings, rep = island_rings(land, water, 200.0)
    existing = [np.asarray(water.interiors[0].coords)[:-1]]
    after = _with_rings(outer, existing + list(rings))
    assert after.difference(water).area < 1e-6
