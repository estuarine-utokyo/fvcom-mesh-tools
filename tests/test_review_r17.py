"""Coastline rules review, round 17."""

import numpy as np
from shapely.geometry import LinearRing, MultiPolygon, Polygon, box

from fvcom_mesh_tools.patch import _clear_of, island_rings


def test_a_ring_placed_first_is_judged_against_its_later_neighbour():
    def h(q):
        return np.clip(590 - np.atleast_2d(q)[:, 0], 30, 200)

    h.size_min, h.size_max = 30, 200
    land = MultiPolygon([box(100, 100, 500, 900), box(540, 200, 740, 800)])
    rings, rep = island_rings(land, box(-500, -500, 1500, 1500), h)
    for i, a in enumerate(rings):
        for j, b in enumerate(rings):
            if i != j:
                assert _clear_of(a, LinearRing(b), h, 0.5)[3]       # possible
    assert rep["n_islands_added"] + len(rep["skipped"]) == 2


def test_a_lake_falls_back_to_its_source_outline():
    lake = Polygon([(100, 100), (900, 100), (900, 900), (600, 900),
                    (600, 550), (400, 550), (400, 900), (100, 900)])
    shell = box(0, 0, 1000, 1000).difference(box(440, 570, 560, 1100))
    land = Polygon(shell.exterior.coords, [lake.exterior.coords])

    def h(q):
        return np.clip(30 + 0.6 * (np.atleast_2d(q)[:, 1] - 100), 30, 300)

    rings, rep = island_rings(land, box(-1000, -1000, 2000, 2000), h, fine_h=60.000001)
    assert rep["n_lakes_added"] == 1
