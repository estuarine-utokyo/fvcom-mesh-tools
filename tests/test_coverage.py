"""Synthetic geometry checks for missing-water coverage, in metres."""
import json

import numpy as np
import pytest

pytest.importorskip("shapely")
pytest.importorskip("skimage")
from shapely.geometry import Polygon, box
from shapely.ops import unary_union

from fvcom_mesh_tools.coverage import coverage_gaps
from fvcom_mesh_tools.io import Fort14Mesh


@pytest.fixture
def mesh():
    nodes = np.array([(x, y) for y in range(0, 101, 10) for x in range(0, 41, 10)], float)
    triangles = []
    for j in range(10):
        for i in range(4):
            n = j * 5 + i
            triangles.extend([(n, n + 1, n + 5), (n + 1, n + 6, n + 5)])
    return Fort14Mesh("synthetic", nodes, np.ones(len(nodes)), np.array(triangles), [], [])


def run(mesh, gap, **kwargs):
    domain = box(0, 0, 200, 100)
    land = domain.difference(unary_union([box(0, 0, 40, 100), gap]))
    return coverage_gaps(mesh, land, domain, pixel_size_m=1, min_area_m2=0, **kwargs)


def test_known_connected_gap(mesh):
    result = run(mesh, box(40, 30, 140, 70))
    r, = result["regions"]
    assert r["area_m2"] == pytest.approx(4000)
    assert r["h_m"] == 10
    assert r["width_p90_m"] == pytest.approx(40, abs=2)
    assert 25 < r["width_p50_m"] <= 40
    assert r["rows_p50"] >= 2
    assert r["classification"] == "clear-defect"
    assert r["touches_mesh"]
    assert r["inland_reach_m"] == pytest.approx(100, abs=2)
    # Rectangle medial axis: 60 m trunk plus four 20*sqrt(2) m branches.
    assert r["length_m"] == pytest.approx(60 + 80 * np.sqrt(2), abs=20)
    assert all(c["feasible"] for c in r["corridors"].values())
    assert result["clear_defects"] == 1
    json.dumps(result, allow_nan=False)


def test_neck_blocks_corridor_despite_wide_basin(mesh):
    gap = unary_union([box(40, 47, 75, 53), box(70, 30, 170, 70)])
    r, = run(mesh, gap)["regions"]
    assert r["classification"] == "clear-defect"
    assert r["w_h_p50"] > 2
    assert not r["corridors"]["1"]["feasible"]
    assert not r["corridors"]["2"]["feasible"]


@pytest.mark.parametrize("width,classification,rows", [(6, "subtarget", 0), (16, "marginal", 1)])
def test_narrow_gap(mesh, width, classification, rows):
    r, = run(mesh, box(40, 50 - width / 2, 140, 50 + width / 2))["regions"]
    assert r["classification"] == classification
    assert r["rows_p50"] == rows
    assert r["width_p50_m"] == pytest.approx(width, abs=2)
    assert r["corridors"]["1"]["feasible"] == (rows == 1)
    assert not r["corridors"]["2"]["feasible"]


def test_isolated_and_policy_separation(mesh):
    gap = unary_union([box(40, 30, 140, 70), box(160, 30, 190, 70)])
    result = run(mesh, gap, policy_land=box(40, 30, 140, 70))
    r, = result["regions"]
    assert not r["touches_mesh"]
    assert r["inland_reach_m"] is None
    assert r["distance_to_mesh_m"] == 120
    assert result["policy_filled_area_m2"] == 4000
    p, = result["policy_filled"]
    assert p["classification"] == "policy-filled"
    assert p["area_m2"] == 4000
    assert result["omitted_water_area_m2"] == 5200


def test_domain_excludes_outside_and_covered_water(mesh):
    result = coverage_gaps(mesh, Polygon(), box(0, 0, 40, 100))
    assert result["regions"] == []
    assert result["omitted_water_area_m2"] == 0


def test_partial_policy_fill_and_area_filter(mesh):
    result = run(mesh, box(40, 30, 140, 70), policy_land=box(90, 30, 140, 70))
    assert sum(r["area_m2"] for r in result["regions"]) == 2000
    assert result["policy_filled_area_m2"] == 2000
    result = coverage_gaps(mesh, Polygon(), box(0, 0, 40.01, 100), min_area_m2=10)
    assert result["regions"] == []
    assert result["below_area_threshold"]["area_m2"] == pytest.approx(1)


def test_budget_is_explicitly_unmeasured(mesh):
    r, = run(mesh, box(40, 30, 140, 70), max_raster_cells=10)["regions"]
    assert r["classification"] == "unmeasured"
    assert "limit" in r["measurement_note"]


@pytest.mark.parametrize("kwargs", [{"pixel_size_m": 0}, {"pixel_size_m": float("nan")},
                                   {"min_area_m2": -1}, {"corridor_fraction": 0}])
def test_bad_parameters(mesh, kwargs):
    with pytest.raises(ValueError):
        coverage_gaps(mesh, Polygon(), box(0, 0, 200, 100), **kwargs)


def test_invalid_geometry_and_empty_mesh(mesh):
    with pytest.raises(ValueError, match="valid polygonal"):
        coverage_gaps(mesh, Polygon([(0, 0), (2, 2), (0, 2), (2, 0)]), box(0, 0, 200, 100))
    mesh.elements = np.empty((0, 3), int)
    with pytest.raises(ValueError, match="nonempty"):
        coverage_gaps(mesh, Polygon(), box(0, 0, 200, 100))


def test_certified_runner_domain_cut_retains_unmeshed_bays():
    """The explicit OBC cut removes the sea side without shrinking to a mesh hull."""
    import runpy
    from pathlib import Path

    pytest.importorskip("geopandas")
    pytest.importorskip("matplotlib")
    pyproj = pytest.importorskip("pyproj")
    from shapely.geometry import Point

    runner = runpy.run_path(str(Path(__file__).resolve().parents[1]
                               / "notebooks/389_coverage_gaps.py"))
    project = pyproj.Transformer.from_crs(4326, 32654, always_xy=True)
    lonlat = np.array([[139.70, 35.0], [139.70, 35.75],
                       [139.80, 35.4], [139.90, 35.4], [139.85, 35.5]])
    nodes = np.column_stack(project.transform(lonlat[:, 0], lonlat[:, 1]))
    mesh = Fort14Mesh("OBC", nodes, np.ones(5), np.array([[2, 3, 4]]),
                      [np.array([0, 1])], [])
    domain = runner["domain_of_interest"](mesh)
    assert domain.contains(Point(*project.transform(139.73, 35.60)))
    assert not domain.contains(Point(*project.transform(139.65, 35.60)))
    assert not domain.contains(Point(*project.transform(140.20, 35.60)))
