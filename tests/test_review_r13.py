"""Coastline rules review, round 13."""

import os

import numpy as np
from shapely.geometry import MultiPolygon, Point, Polygon, box

from fvcom_mesh_tools import provenance
from fvcom_mesh_tools.patch import island_rings
from tests.test_review_r12 import _with_rings


def _const(v):
    return lambda q: np.full(len(np.atleast_2d(q)), float(v))


def test_two_lakes_resampled_apart_do_not_overlap():
    water = box(0, 0, 1000, 1000)
    lake1 = Polygon([(200, 200), (700, 200), (700, 440), (400, 440),
                     (400, 450), (700, 450), (700, 700), (200, 700)])
    lake2 = box(580, 441, 620, 442)
    land = Polygon(box(100, 100, 900, 900).exterior.coords,
                   [lake1.exterior.coords, lake2.exterior.coords])
    assert land.is_valid
    rings, rep = island_rings(land, water, _const(30.0))
    _with_rings(water, rings)                         # raises if rings touch


def test_an_island_in_a_filled_lake_stays_land():
    water = box(0, 0, 1000, 1000)
    outer = Polygon(box(100, 100, 900, 900).exterior.coords,
                    [box(105, 300, 125, 700).exterior.coords])    # no element fits
    inner = box(110, 400, 120, 410)
    rings, rep = island_rings(MultiPolygon([outer, inner]), water, _const(30.0))
    after = _with_rings(water, rings)
    assert not after.contains(inner.centroid)
    assert rep["n_lakes_added"] == 0
    assert any(s["why"] == "inside a lake meshed as land" for s in rep["skipped"])


def test_an_island_in_a_kept_lake_is_added():
    water = box(0, 0, 1000, 1000)
    outer = Polygon(box(100, 100, 900, 900).exterior.coords,
                    [box(300, 300, 700, 700).exterior.coords])
    inner = box(450, 450, 550, 550)
    rings, rep = island_rings(MultiPolygon([outer, inner]), water, _const(30.0))
    after = _with_rings(water, rings)
    assert after.contains(Point(350, 350)) and not after.contains(Point(500, 500))
    assert rep["n_islands_added"] == 2 and rep["n_lakes_added"] == 1


def test_lake_clearance_is_judged_at_the_local_size():
    def h(q):
        q = np.atleast_2d(q)
        return np.maximum(30.0, 30.0 + 0.1 * ((q[:, 0] - 150) + (q[:, 1] - 150)))

    land = Polygon(box(100, 100, 900, 900).exterior.coords,
                   [box(150, 150, 850, 850).exterior.coords])
    rings, rep = island_rings(land, box(-1000, -1000, 2000, 2000), h)
    # 50 m from its coast where the element is 170 m: named, not passed
    # silently; an element fits in it, so it stays water
    assert rep["n_lakes_added"] == 1
    tight = [s for s in rep["tight"] if "lake" in s["why"]]
    assert tight and float(tight[0]["why"].split(" where the element is ")[1][:-2]) > 150


def test_code_state_hashes_a_symlinked_subpackage(tmp_path, monkeypatch):
    monkeypatch.setattr(provenance, "git_state", lambda p: None)
    pkg, ext = tmp_path / "pkg", tmp_path / "ext"
    pkg.mkdir()
    ext.mkdir()
    (pkg / "__init__.py").write_text("")
    (ext / "__init__.py").write_text("a = 1\n")
    os.symlink(ext, pkg / "subpkg")
    os.symlink(pkg, ext / "loop")                    # a cycle ends the walk
    first = provenance.code_state("nope", pkg)["source_sha256"]
    (ext / "__init__.py").write_text("a = 2\n")
    second = provenance.code_state("nope", pkg)["source_sha256"]
    assert first and second and first != second
