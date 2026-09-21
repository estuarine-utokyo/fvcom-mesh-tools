"""Tests for the local-refinement specification and its pre-flight checks."""

from __future__ import annotations

import numpy as np
import pytest
from shapely.geometry import Polygon

from fvcom_mesh_tools.refine import (
    RefineRegion,
    frozen_changes,
    load_refine,
    preflight,
    transition_width_m,
)

AMBIENT = 350.0
GRAD = 0.165


def _region(**over):
    spec = {
        "name": "fishery",
        "geometry": {"circle": {"center": [139.7881, 35.3228], "radius_m": 300}},
        "target_h_m": 30,
    }
    spec.update(over)
    return RefineRegion(spec)


def _flat_depth(value):
    return lambda lon, lat: np.full(np.size(lon), float(value))


def _land_north_of(lat):
    return Polygon([(139.0, lat), (140.5, lat), (140.5, lat + 1.0), (139.0, lat + 1.0)])


# --- the width a gradation needs -------------------------------------------

def test_transition_width_is_the_size_gap_over_the_gradation():
    assert transition_width_m(30, 350, 0.165) == pytest.approx(320 / 0.165)
    # A 30 m target in a 350 m field needs nearly 2 km, six times a 300 m core.
    assert transition_width_m(30, 350, 0.165) == pytest.approx(1939.4, abs=0.1)


def test_transition_width_is_zero_when_no_climb_is_needed():
    assert transition_width_m(350, 350, 0.165) == 0.0
    assert transition_width_m(400, 350, 0.165) == 0.0


def test_a_steeper_gradation_buys_a_shorter_transition():
    assert transition_width_m(30, 350, 0.33) < transition_width_m(30, 350, 0.165)


# --- the specification ------------------------------------------------------

def test_circle_geometry_has_the_requested_radius():
    r = _region()
    w, s, e, n = r.geometry.bounds
    half_ns = (n - s) / 2 * 111000.0
    half_ew = (e - w) / 2 * 111000.0 * np.cos(np.radians(35.3228))
    assert half_ns == pytest.approx(300, abs=1)
    assert half_ew == pytest.approx(300, abs=1)


@pytest.mark.parametrize("bad", [
    {"name": ""},
    {"target_h_m": 0},
    {"target_h_m": -30},
    {"priority": True},
    {"touch_coast": "yes"},
    {"transition_m": -1},
])
def test_invalid_region_fields_are_rejected(bad):
    with pytest.raises(ValueError):
        _region(**bad)


def test_unknown_region_keys_are_rejected():
    with pytest.raises(ValueError):
        RefineRegion({"name": "x", "geometry": {"bbox": [139.7, 35.3, 139.8, 35.4]},
                      "target_h_m": 30, "surprise": 1})


def test_recipe_round_trip(tmp_path):
    mesh = tmp_path / "base.14"
    mesh.write_text("stub\n")
    p = tmp_path / "r.yaml"
    p.write_text(
        "base_mesh: base.14\ndt_floor_s: 4.5\ngradation: 0.165\n"
        "refine:\n  - name: a\n"
        "    geometry: {circle: {center: [139.79, 35.32], radius_m: 300}}\n"
        "    target_h_m: 30\n"
    )
    cfg = load_refine(p)
    assert cfg["base_mesh"] == mesh.resolve()
    assert cfg["dt_floor_s"] == 4.5
    assert [r.name for r in cfg["refine"]] == ["a"]


def test_recipe_rejects_a_missing_base_mesh(tmp_path):
    p = tmp_path / "r.yaml"
    p.write_text(
        "base_mesh: nowhere.14\ndt_floor_s: 4.5\ngradation: 0.165\n"
        "refine:\n  - name: a\n    geometry: {bbox: [139.7, 35.3, 139.8, 35.4]}\n"
        "    target_h_m: 30\n"
    )
    with pytest.raises(ValueError, match="base_mesh"):
        load_refine(p)


def test_recipe_rejects_duplicate_region_names(tmp_path):
    mesh = tmp_path / "base.14"
    mesh.write_text("stub\n")
    p = tmp_path / "r.yaml"
    p.write_text(
        "base_mesh: base.14\ndt_floor_s: 4.5\ngradation: 0.165\nrefine:\n"
        + "".join(
            "  - name: a\n    geometry: {bbox: [139.7, 35.3, 139.8, 35.4]}\n"
            "    target_h_m: 30\n" for _ in range(2))
    )
    with pytest.raises(ValueError):
        load_refine(p)


# --- pre-flight -------------------------------------------------------------

def test_preflight_passes_and_reports_the_derived_width():
    # 2 m of water is shallow enough that a 30 m target clears a 4.5 s floor.
    rep = preflight(_region(touch_coast=True), gradation=GRAD, dt_floor_s=4.5,
                    ambient_h_m=AMBIENT, depth_of=_flat_depth(2.0), land=None)
    assert rep["transition_m"] == pytest.approx(1939.4, abs=0.1)
    assert rep["elements_core"] > 700
    assert rep["core_on_land"] == 0


def test_preflight_reports_the_altitude_time_step_not_the_edge_one():
    # The reported dt uses the minimum ALTITUDE (notebook 392), which for an
    # equilateral cell is sqrt(3)/2 of the edge. Quoting the edge figure
    # overstates the step by 15 % -- enough to clear a floor it does not meet.
    rep = preflight(_region(touch_coast=True), gradation=GRAD, dt_floor_s=4.0,
                    ambient_h_m=AMBIENT, depth_of=_flat_depth(2.0), land=None)
    edge = 30 / np.sqrt(9.81 * 2.0)
    assert rep["dt_by_shortest_edge_s"] == pytest.approx(edge, rel=1e-9)
    assert rep["dt_s"] == pytest.approx(np.sqrt(3) / 2 * edge, rel=1e-9)
    assert rep["dt_s"] < rep["dt_by_shortest_edge_s"]


def test_preflight_refuses_a_target_that_breaks_the_time_step():
    # 30 m over 20 m of water allows 1.9 s by altitude, far below a 4.5 s floor.
    with pytest.raises(ValueError, match="allows dt"):
        preflight(_region(touch_coast=True), gradation=GRAD, dt_floor_s=4.5,
                  ambient_h_m=AMBIENT, depth_of=_flat_depth(20.0), land=None)


def test_preflight_refuses_the_futtsu_target_on_its_own_floor():
    # 30 m over 4.15 m of water: 4.70 s by edge, 4.07 s by altitude. The first
    # Futtsu recipe was accepted on the edge figure and fails on the real one.
    with pytest.raises(ValueError, match="4.07 s by minimum altitude"):
        preflight(_region(touch_coast=True), gradation=GRAD, dt_floor_s=4.5,
                  ambient_h_m=AMBIENT, depth_of=_flat_depth(4.15), land=None)


def test_preflight_names_the_largest_target_the_water_permits():
    with pytest.raises(ValueError) as exc:
        preflight(_region(touch_coast=True), gradation=GRAD, dt_floor_s=4.5,
                  ambient_h_m=AMBIENT, depth_of=_flat_depth(20.0), land=None)
    permitted = 4.5 * np.sqrt(9.81 * 20.0) / (np.sqrt(3) / 2)
    assert f"{permitted:.0f} m" in str(exc.value)


def test_preflight_refuses_a_transition_too_short_for_the_gradation():
    with pytest.raises(ValueError, match="shorter than"):
        preflight(_region(transition_m=500, touch_coast=True), gradation=GRAD,
                  dt_floor_s=4.5, ambient_h_m=AMBIENT, depth_of=_flat_depth(2.0), land=None)


def test_preflight_accepts_a_transition_wider_than_required():
    rep = preflight(_region(transition_m=3000, touch_coast=True), gradation=GRAD,
                    dt_floor_s=4.5, ambient_h_m=AMBIENT, depth_of=_flat_depth(2.0),
                    land=None)
    assert rep["transition_m"] == 3000
    assert rep["transition_required_m"] == pytest.approx(1939.4, abs=0.1)


def test_preflight_refuses_a_core_on_land_unless_asked():
    land = _land_north_of(35.3228)          # covers the northern half of the core
    with pytest.raises(ValueError, match="touch_coast"):
        preflight(_region(), gradation=GRAD, dt_floor_s=4.5, ambient_h_m=AMBIENT,
                  depth_of=_flat_depth(2.0), land=land)
    rep = preflight(_region(touch_coast=True), gradation=GRAD, dt_floor_s=4.5,
                    ambient_h_m=AMBIENT, depth_of=_flat_depth(2.0), land=land)
    assert rep["core_on_land"] > 0


def test_preflight_requires_a_land_polygon_when_the_core_must_avoid_it():
    # Passing land=None used to be a silent pass, which is the wrong default
    # for a check whose whole purpose is to refuse.
    with pytest.raises(ValueError, match="land polygon is required"):
        preflight(_region(), gradation=GRAD, dt_floor_s=4.5, ambient_h_m=AMBIENT,
                  depth_of=_flat_depth(2.0), land=None)


def test_preflight_refuses_a_core_entirely_on_land():
    land = _land_north_of(35.0)
    with pytest.raises(ValueError, match="entirely on land"):
        preflight(_region(touch_coast=True), gradation=GRAD, dt_floor_s=4.5,
                  ambient_h_m=AMBIENT, depth_of=_flat_depth(2.0), land=land)


# --- the frozen-region contract --------------------------------------------

def test_frozen_changes_passes_when_only_affected_nodes_moved():
    base = np.array([[0.0, 0.0], [1.0, 0.0], [2.0, 0.0]])
    new = base.copy()
    new[0] += 5.0
    out = frozen_changes(base, new, np.array([True, False, False]))
    assert out["ok"] and out["n_moved_in_frozen"] == 0 and out["n_frozen"] == 2


def test_frozen_changes_catches_a_node_that_moved_outside_the_region():
    base = np.array([[0.0, 0.0], [1.0, 0.0], [2.0, 0.0]])
    new = base.copy()
    new[2, 0] += 0.01          # x only: the move is exactly 0.01 m
    out = frozen_changes(base, new, np.array([True, False, False]))
    assert not out["ok"]
    assert out["n_moved_in_frozen"] == 1
    assert out["max_move_in_frozen_m"] == pytest.approx(0.01)


def test_frozen_changes_requires_matching_arrays():
    with pytest.raises(ValueError):
        frozen_changes(np.zeros((3, 2)), np.zeros((4, 2)), np.zeros(3, bool))


# --- what the hole actually cuts -------------------------------------------

def _strip(nx=12, ny=4, h=100.0):
    """A regular right-triangle strip: a stand-in for a patch of real mesh."""
    xs, ys = np.meshgrid(np.arange(nx + 1) * h, np.arange(ny + 1) * h, indexing="ij")
    nodes = np.column_stack([xs.ravel(), ys.ravel()])

    def nid(i, j):
        return i * (ny + 1) + j

    tri = []
    for i in range(nx):
        for j in range(ny):
            tri.append([nid(i, j), nid(i + 1, j), nid(i + 1, j + 1)])
            tri.append([nid(i, j), nid(i + 1, j + 1), nid(i, j + 1)])
    return nodes, np.asarray(tri, dtype=np.int64)


def test_hole_clearance_sees_a_footprint_that_reaches_the_boundary():
    from shapely.geometry import Point

    from fvcom_mesh_tools.refine import hole_clearance

    nodes, tri = _strip()
    # Centred on the strip, a 400 m footprint reaches the top and bottom edges.
    rep = hole_clearance(nodes, tri, Point(600.0, 200.0).buffer(1.0),
                         transition_m=400.0)
    assert rep["n_selected"] > 0
    assert rep["reaches_boundary"] is True
    assert rep["n_physical_boundary_edges"] > 0


def test_hole_clearance_reports_an_interior_footprint_as_clear():
    from shapely.geometry import Point

    from fvcom_mesh_tools.refine import hole_clearance

    nodes, tri = _strip(nx=20, ny=20)
    rep = hole_clearance(nodes, tri, Point(1000.0, 1000.0).buffer(1.0),
                         transition_m=250.0)
    assert rep["n_selected"] > 0
    assert rep["reaches_boundary"] is False
    assert rep["n_interface_edges"] == rep["n_rim_edges"]


def test_hole_clearance_reports_the_reach_beyond_the_requested_envelope():
    from shapely.geometry import Point

    from fvcom_mesh_tools.refine import hole_clearance

    nodes, tri = _strip(nx=20, ny=20)
    rep = hole_clearance(nodes, tri, Point(1000.0, 1000.0).buffer(1.0),
                         transition_m=250.0)
    # Whole triangles are taken, so the selection always reaches past the
    # analytic envelope; the caller needs the number, not a promise.
    assert rep["selection_reach_m"] > rep["requested_reach_m"]


def test_hole_clearance_flags_the_open_boundary():
    from shapely.geometry import Point

    from fvcom_mesh_tools.refine import hole_clearance

    nodes, tri = _strip(nx=20, ny=20)
    rep = hole_clearance(nodes, tri, Point(1000.0, 1000.0).buffer(1.0),
                         transition_m=250.0, open_boundaries=[np.arange(0, 21)])
    assert rep["reaches_open_boundary"] is False
    near = hole_clearance(nodes, tri, Point(100.0, 100.0).buffer(1.0),
                          transition_m=250.0, open_boundaries=[np.arange(0, 21)])
    assert near["reaches_open_boundary"] is True
    assert near["n_open_boundary_nodes"] > 0


def test_frozen_changes_rejects_a_destroyed_coordinate():
    # NaN fails every comparison, so a bare `moved > tol` passed this.
    base = np.array([[0.0, 0.0], [1.0, 0.0]])
    new = base.copy()
    new[1, 0] = np.nan
    with pytest.raises(ValueError, match="non-finite"):
        frozen_changes(base, new, np.array([True, False]))
