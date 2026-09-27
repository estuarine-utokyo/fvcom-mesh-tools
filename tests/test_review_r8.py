"""Coastline rules review, round 8."""

import numpy as np
import pytest
from shapely.geometry import LineString, Polygon

from fvcom_mesh_tools import provenance
from fvcom_mesh_tools.patch import _source_substring, hole_polygon, rim_repair
from fvcom_mesh_tools.walls import close_wall_pockets


def _const(v):
    return lambda q: np.full(len(np.atleast_2d(q)), float(v))


def test_a_snapped_slit_chord_that_crosses_the_rim_is_refused():
    p = np.array([(0, 0), (0, 5), (20, 5), (20, 2), (10, 2),
                  (10, -2), (8, -1), (8, -20), (-100, -20), (-100, 0)], dtype=float)
    e = np.column_stack([np.arange(10), np.roll(np.arange(10), -1)])
    out, eg, _b, _m, rep = rim_repair(p, e, np.full(10, -1), Polygon(p), _const(30.0),
                                      protect=(5, 6), operations=("slits",), rounds=1,
                                      retreat_tips=False, focus=[[0, 0]], focus_factor=0.001)
    hole = hole_polygon(out, eg)                     # raises on an invalid rim
    assert hole.is_valid and hole.area > 0


def test_hole_polygon_refuses_an_invalid_rim():
    p = np.array([(0, 0), (10, 0), (0, 10), (10, 10)], dtype=float)   # a bow tie
    e = np.array([(0, 1), (1, 2), (2, 3), (3, 0)])
    with pytest.raises(ValueError, match="valid polygon"):
        hole_polygon(p, e)


def test_equal_arcs_do_not_depend_on_the_ring_start_or_direction():
    pts = np.array([(0, 0), (60, 40), (100, 100)], dtype=float)
    rings = [LineString([(0, 0), (100, 0), (100, 100), (0, 100), (0, 0)]),
             LineString([(100, 0), (100, 100), (0, 100), (0, 0), (100, 0)]),
             LineString([(0, 0), (0, 100), (100, 100), (100, 0), (0, 0)])]
    got = [LineString(_source_substring(pts, r)) for r in rings]
    want = LineString([(0, 0), (100, 0), (100, 100)])     # the arc that fits
    for g in got:
        assert g.equals(want)


def test_pockets_close_with_no_land_at_all():
    wall = LineString([(0, 0), (0, 400), (110, 400), (110, 290), (20, 290)])
    added, _, rep = close_wall_pockets([wall], Polygon(), _const(100.0),
                                       min_h=60.000001, size_floor=30.0)
    assert rep["n_pockets_closed"] == 1 and added.area > 10000.0


def test_unlisted_sidecars_null_the_dataset_digest(tmp_path, monkeypatch):
    from pathlib import Path

    shp = tmp_path / "land.shp"
    shp.write_bytes(b"x")
    (tmp_path / "land.dbf").write_bytes(b"y")

    def refuse(self):
        raise PermissionError(13, "injected")

    monkeypatch.setattr(Path, "iterdir", refuse)
    files = provenance.dataset_files(shp)
    assert files[0] == shp and len(files) == 2 and provenance.UNLISTED in str(files[1])
    rec = provenance.collect(files={"land": files}, libraries=())
    assert rec["files"]["land"]["sha256"] is None


def test_listed_sidecars_still_hash(tmp_path):
    shp = tmp_path / "land.shp"
    shp.write_bytes(b"x")
    (tmp_path / "land.dbf").write_bytes(b"y")
    rec = provenance.collect(files={"land": provenance.dataset_files(shp)}, libraries=())
    assert rec["files"]["land"]["sha256"] is not None


def test_a_filtered_away_island_is_legitimate_water():
    # the driver's rule: an empty source is refused only for anchored stretches
    from fvcom_mesh_tools.patch import rim_constraints
    from tests.test_review_r7 import _island_mesh

    nodes, elements, sel = _island_mesh()
    frozen = set(sel.frozen_nodes.tolist())
    anchored = any(any(int(v) in frozen for v in r) and any(int(v) not in frozen for v in r)
                   for r in sel.rings)
    assert not anchored
    rc = rim_constraints(nodes, sel, size=_const(30.0), coastline="resolve", shoreline=[])
    assert rc["n_islands_left_to_source"] == 1
    assert hole_polygon(rc["pfix"], rc["egfix"]).area == pytest.approx(800.0 * 800.0)


def test_each_search_pass_names_its_rejected_meshes(tmp_path):
    from tests.test_local_refine_driver import DRIVER, driver_function, seed_search_env

    passes = []
    env = seed_search_env(tmp_path, {3: (np.array([[3.0]]),)}, {3: 0})

    def attempt(seed):
        passes.append(env["SEARCH_PASS"])
        raise ValueError("stub: no candidate")

    env.update(attempt=attempt, SEARCH_PASS=1, save_report=lambda: None)
    env["reports"]["attempts"] = []
    search = driver_function("seed_search", env)
    search([3])
    search([3], search_pass=2)
    assert passes == [1, 2]
    src = DRIVER.read_text()
    assert 'f"rejected_pass{SEARCH_PASS}_seed{seed}.npz"' in src
    assert "seed_search(seeds, search_pass=2)" in src
