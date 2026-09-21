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
    rep = preflight(_region(), gradation=GRAD, dt_floor_s=4.5, ambient_h_m=AMBIENT,
                    depth_of=_flat_depth(4.0), land=None)
    assert rep["transition_m"] == pytest.approx(1939.4, abs=0.1)
    assert rep["dt_s"] == pytest.approx(30 / np.sqrt(9.81 * 4.0), rel=1e-9)
    assert rep["elements_core"] > 700
    assert rep["core_on_land"] == 0


def test_preflight_refuses_a_target_that_breaks_the_time_step():
    # 30 m over 20 m of water allows 2.1 s, far below a 4.5 s floor.
    with pytest.raises(ValueError, match="allows dt"):
        preflight(_region(), gradation=GRAD, dt_floor_s=4.5, ambient_h_m=AMBIENT,
                  depth_of=_flat_depth(20.0), land=None)


def test_preflight_names_the_largest_target_the_water_permits():
    with pytest.raises(ValueError) as exc:
        preflight(_region(), gradation=GRAD, dt_floor_s=4.5, ambient_h_m=AMBIENT,
                  depth_of=_flat_depth(20.0), land=None)
    permitted = 4.5 * np.sqrt(9.81 * 20.0)
    assert f"{permitted:.0f} m" in str(exc.value)


def test_preflight_refuses_a_transition_too_short_for_the_gradation():
    with pytest.raises(ValueError, match="shorter than"):
        preflight(_region(transition_m=500), gradation=GRAD, dt_floor_s=4.5,
                  ambient_h_m=AMBIENT, depth_of=_flat_depth(4.0), land=None)


def test_preflight_accepts_a_transition_wider_than_required():
    rep = preflight(_region(transition_m=3000), gradation=GRAD, dt_floor_s=4.5,
                    ambient_h_m=AMBIENT, depth_of=_flat_depth(4.0), land=None)
    assert rep["transition_m"] == 3000
    assert rep["transition_required_m"] == pytest.approx(1939.4, abs=0.1)


def test_preflight_refuses_a_core_on_land_unless_asked():
    land = _land_north_of(35.3228)          # covers the northern half of the core
    with pytest.raises(ValueError, match="touch_coast"):
        preflight(_region(), gradation=GRAD, dt_floor_s=4.5, ambient_h_m=AMBIENT,
                  depth_of=_flat_depth(4.0), land=land)
    rep = preflight(_region(touch_coast=True), gradation=GRAD, dt_floor_s=4.5,
                    ambient_h_m=AMBIENT, depth_of=_flat_depth(4.0), land=land)
    assert rep["core_on_land"] > 0


def test_preflight_refuses_a_core_entirely_on_land():
    land = _land_north_of(35.0)
    with pytest.raises(ValueError, match="entirely on land"):
        preflight(_region(touch_coast=True), gradation=GRAD, dt_floor_s=4.5,
                  ambient_h_m=AMBIENT, depth_of=_flat_depth(4.0), land=land)


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
