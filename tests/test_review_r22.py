"""Coastline rules review, round 22."""

from shapely.geometry import MultiPolygon, Polygon, box

from fvcom_mesh_tools.patch import _preference_search, island_rings
from tests.test_review_r12 import _with_rings

A = box(100, 100, 900, 900).difference(box(100, 100, 110, 500))
B = box(-100, 200, 94, 450)


def test_a_sibling_island_is_not_lost_to_its_neighbours_outline():
    water = box(-1000, -1000, 2000, 2000)
    rings, rep = island_rings(MultiPolygon([A, B]), water, 30.0)
    assert rep["n_islands_added"] == 2, rep["skipped"]
    assert _with_rings(water, rings).intersection(B).area < 1.0


def test_an_island_in_a_lake_is_not_lost_to_the_lake_outline():
    lake = box(100, 100, 900, 900).union(box(90, 100, 100, 500))
    child = box(106, 200, 700, 450)
    parent = Polygon(box(0, 0, 1000, 1000).exterior.coords, [lake.exterior.coords])
    water = box(-1000, -1000, 2000, 2000)
    rings, rep = island_rings(MultiPolygon([parent, child]), water, 30.0)
    assert rep["n_islands_added"] == 2 and rep["n_lakes_added"] == 1, rep["skipped"]


def test_the_large_search_is_not_spent_on_duplicates():
    target = (True,) * 6 + (False,) * 3

    def run(prefs):
        return ((0 if prefs == target else 1), prefs), None

    best, _why = _preference_search(9, [(True,) * 9, (False,) * 9], run)
    assert best[0] == 0


def test_invalid_starts_still_expand():
    target = (False,) + (True,) * 8

    def run(prefs):
        if prefs == target:
            return (0, prefs), None
        return None, "invalid"

    best, _why = _preference_search(9, [(True,) * 9, (False,) * 9], run)
    assert best is not None and best[0] == 0
