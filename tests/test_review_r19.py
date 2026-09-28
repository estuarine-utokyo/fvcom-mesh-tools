"""Coastline rules review, round 19."""

import numpy as np
from shapely.geometry import Polygon, box

from fvcom_mesh_tools.patch import island_rings
from tests.test_review_r12 import _with_rings


def test_a_filled_lake_may_not_cover_a_rim_island():
    outer = box(0, 0, 1000, 1000)
    water = outer.difference(box(450, 450, 550, 550))
    land = Polygon(box(100, 100, 900, 900).exterior.coords,
                   [box(440, 440, 560, 560).exterior.coords])
    rings, rep = island_rings(land, water, 30.0)
    existing = [np.asarray(water.interiors[0].coords)[:-1]]
    after = _with_rings(outer, existing + list(rings))
    assert after.difference(water).area < 1.0


def test_sibling_lakes_survive_in_either_order():
    a = box(100, 100, 900, 900).difference(box(100, 100, 110, 500))
    b = box(-100, 200, 94, 450)
    shell = box(-300, -100, 1100, 1100)
    water = box(-1000, -1000, 2000, 2000)
    for lakes in ([a, b], [b, a]):
        land = Polygon(shell.exterior.coords, [lk.exterior.coords for lk in lakes])
        rings, rep = island_rings(land, water, 30.0)
        assert rep["n_lakes_added"] == 2, rep["skipped"]
