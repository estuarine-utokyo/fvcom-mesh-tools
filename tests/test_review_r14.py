"""Coastline rules review, round 14."""

import os

import numpy as np
from shapely.geometry import MultiPolygon, Polygon, box

from fvcom_mesh_tools import provenance
from fvcom_mesh_tools.patch import island_rings
from tests.test_review_r12 import _with_rings


def _const(v):
    return lambda q: np.full(len(np.atleast_2d(q)), float(v))


def test_an_island_in_a_large_lake_comes_after_the_island_around_it():
    water = box(0, 0, 1000, 1000)
    outer = Polygon(box(100, 100, 900, 900).exterior.coords,
                    [box(105, 105, 895, 895).exterior.coords])
    inner = box(200, 200, 800, 800)
    rings, rep = island_rings(MultiPolygon([outer, inner]), water, _const(30.0))
    after = _with_rings(water, rings)
    assert after.intersection(inner).area < 1.0


def test_sibling_lakes_stay_siblings():
    water = box(0, 0, 1000, 1000)
    lake1 = Polygon([(200, 200), (700, 200), (700, 440), (400, 440),
                     (400, 450), (700, 450), (700, 700), (200, 700)])
    lake2 = box(420, 447, 425, 449)
    land = Polygon(box(100, 100, 900, 900).exterior.coords,
                   [lake1.exterior.coords, lake2.exterior.coords])
    rings, rep = island_rings(land, water, _const(30.0))
    polys = [Polygon(r) for r in rings]
    for i, a in enumerate(polys):
        for j, b in enumerate(polys):
            if i != j and a.contains(b):
                assert i == 0                          # only the island holds lakes
    after = _with_rings(water, rings)
    assert after.contains(lake1.representative_point())


def test_a_parallel_gap_is_judged_along_its_length():
    def h(q):
        q = np.atleast_2d(q)
        return 30.0 + 170.0 * np.maximum(0.0, 1.0 - np.abs(q[:, 1] - 500.0) / 400.0)

    rings, rep = island_rings(box(50, 100, 950, 900), box(0, 0, 1000, 1000), h)
    assert rep["n_islands_added"] == 1 and rep["tight"]


def test_constant_and_column_sizes_still_work():
    for size in (30.0, lambda q: np.full((len(q), 1), 30.0)):
        rings, rep = island_rings(box(100, 100, 900, 900), box(0, 0, 1000, 1000), size)
        assert rep["n_islands_added"] == 1


def _pkg_with_aliases(tmp_path, names):
    pkg, shared = tmp_path / "pkg", tmp_path / "shared"
    pkg.mkdir(parents=True)
    shared.mkdir(parents=True)
    (pkg / "__init__.py").write_text("")
    (shared / "mod.py").write_text("x = 1\n")
    for n in names:
        os.symlink(shared, pkg / n)
    return pkg


def test_aliases_all_count_and_the_digest_is_stable(tmp_path, monkeypatch):
    monkeypatch.setattr(provenance, "git_state", lambda p: None)
    both = provenance.code_state("nope", _pkg_with_aliases(tmp_path / "1", ["a", "b"]))
    again = provenance.code_state("nope", _pkg_with_aliases(tmp_path / "2", ["b", "a"]))
    one = provenance.code_state("nope", _pkg_with_aliases(tmp_path / "3", ["b"]))
    assert both["source_sha256"] == again["source_sha256"] != one["source_sha256"]


def test_a_link_to_an_ancestor_ends_the_walk(tmp_path, monkeypatch):
    monkeypatch.setattr(provenance, "git_state", lambda p: None)
    pkg = tmp_path / "pkg"
    (pkg / "sub").mkdir(parents=True)
    (pkg / "__init__.py").write_text("")
    os.symlink(pkg, pkg / "sub" / "up")
    st = provenance.code_state("nope", pkg)
    assert st["source_sha256"] and st["source_unreadable"] == []
