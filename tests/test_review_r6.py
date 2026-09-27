"""Coastline rules review, round 6."""

import os

import numpy as np
import pytest
import shapely
from shapely.geometry import Polygon, box

from fvcom_mesh_tools import provenance
from fvcom_mesh_tools.patch import filter_shoreline_local, island_rings, unresolvable_water


def _const(v):
    return lambda q: np.full(len(np.atleast_2d(q)), float(v))


def test_continuous_width_pad_holds_discs_centred_outside_the_footprint():
    def size(q):
        q = np.atleast_2d(q)
        return np.minimum(200.0, 100.0 + 0.165 * np.maximum(q[:, 0], 0.0))

    added, rep = unresolvable_water(box(-1000, -4000, 0, 4000), size,
                                    box(0, -1000, 5, 1000),
                                    radius_factor=0.75, min_h=60.0, spacing=1.0)
    assert added.is_empty and not rep["unbounded_pad"]


def test_continuous_width_still_closes_a_dead_end():
    land = shapely.difference(box(0, 0, 400, 400), box(150, 200, 170, 400))
    added, rep = unresolvable_water(land, _const(100.0), box(0, 0, 400, 400),
                                    radius_factor=0.75, min_h=60.0, spacing=5.0)
    assert added.area > 0.5 * 20 * 200 and not rep["unbounded_pad"]


def test_continuous_width_closes_nothing_when_the_pad_cannot_settle():
    def size(q):                                   # grows faster than any disc
        q = np.atleast_2d(q)
        return 100.0 + 5.0 * np.abs(q[:, 0])

    added, rep = unresolvable_water(box(-1000, -4000, 0, 4000), size,
                                    box(0, -1000, 5, 1000),
                                    radius_factor=0.75, min_h=60.0, spacing=5.0)
    assert added.is_empty and rep["unbounded_pad"]


def test_island_rings_reach_polygons_nested_in_collections():
    raw = Polygon([(0, 0), (100, 0), (100, 100), (0, 100), (0, 0),
                   (-100, 0), (-100, -100), (-200, -100), (-200, 0), (-100, 0), (0, 0)])
    land = shapely.make_valid(raw)
    assert land.geom_type == "GeometryCollection"
    rings, rep = island_rings(land, box(-1000, -1000, 1000, 1000), _const(30.0))
    assert rep["n_islands_added"] + len(rep["skipped"]) == 2


def test_local_filter_refuses_an_empty_footprint():
    with pytest.raises(ValueError, match="footprint is empty"):
        filter_shoreline_local(box(0, 0, 100, 100), _const(30.0), 30.0, Polygon())


def test_code_state_reports_a_directory_it_cannot_scan(tmp_path, monkeypatch):
    monkeypatch.setattr(provenance, "git_state", lambda p: None)
    (tmp_path / "a.py").write_text("x = 1\n")

    def refuse(*a, **k):
        raise PermissionError(13, "injected", str(tmp_path))

    monkeypatch.setattr(os, "scandir", refuse)
    st = provenance.code_state("nope-not-installed", tmp_path)
    assert st["source_sha256"] is None and st["source_unreadable"]


def test_a_larger_pad_only_adds_cells_around_the_same_raster(monkeypatch):
    import fvcom_mesh_tools.patch as patch_mod

    land = shapely.difference(box(0, 0, 400, 400), box(150, 200, 173, 400))
    kw = dict(radius_factor=0.75, min_h=60.0, spacing=7.0)
    base, _ = unresolvable_water(land, _const(100.0), box(0, 0, 400, 400), **kw)
    real = patch_mod.size_upper_bound
    monkeypatch.setattr(patch_mod, "size_upper_bound", lambda *a, **k: 1.9 * real(*a, **k))
    wide, rep = unresolvable_water(land, _const(100.0), box(0, 0, 400, 400), **kw)
    assert not rep["unbounded_pad"] and not base.is_empty
    assert shapely.equals_exact(base.normalize(), wide.normalize(), 1e-6)
