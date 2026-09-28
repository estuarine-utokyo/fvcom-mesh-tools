"""Coastline rules review, round 18."""

import numpy as np
from shapely.geometry import Polygon, box

from fvcom_mesh_tools.patch import island_rings
from tests.test_review_r12 import _with_rings


def _h(q):
    return np.full(len(np.atleast_2d(q)), 30.0)


_h.size_min = _h.size_max = 30.0


def test_a_lake_the_source_shell_keeps_is_not_filled():
    water = box(0, 0, 1000, 1000)
    shell = box(100, 100, 900, 900).union(box(90, 100, 100, 500))
    lake = box(106, 200, 700, 450)
    land = Polygon(shell.exterior.coords, [lake.exterior.coords])
    rings, rep = island_rings(land, water, _h)
    assert rep["n_islands_added"] == 1 and rep["n_lakes_added"] == 1


def test_an_island_may_surround_a_rim_island_in_its_lake():
    water = box(0, 0, 1000, 1000).difference(box(450, 450, 550, 550))
    land = Polygon(box(100, 100, 900, 900).exterior.coords,
                   [box(300, 300, 700, 700).exterior.coords])
    assert land.within(water)
    rings, rep = island_rings(land, water, _h)
    assert rep["n_islands_added"] == 1 and rep["n_lakes_added"] == 1
    hole = [np.asarray(water.interiors[0].coords)[:-1]]
    after = _with_rings(box(0, 0, 1000, 1000), hole + list(rings))
    assert after.symmetric_difference(water.difference(land)).area < 50.0
