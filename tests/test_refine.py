"""Tests for the local-refinement specification and its pre-flight checks."""

from __future__ import annotations

import numpy as np
import pytest
import shapely
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
        "base_mesh: base.14\ndt_expected_s: 4.5\ngradation: 0.165\n"
        "refine:\n  - name: a\n"
        "    geometry: {circle: {center: [139.79, 35.32], radius_m: 300}}\n"
        "    target_h_m: 30\n"
    )
    cfg = load_refine(p)
    assert cfg["base_mesh"] == mesh.resolve()
    assert cfg["dt_expected_s"] == 4.5
    assert [r.name for r in cfg["refine"]] == ["a"]


def test_recipe_rejects_a_missing_base_mesh(tmp_path):
    p = tmp_path / "r.yaml"
    p.write_text(
        "base_mesh: nowhere.14\ndt_expected_s: 4.5\ngradation: 0.165\n"
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
        "base_mesh: base.14\ndt_expected_s: 4.5\ngradation: 0.165\nrefine:\n"
        + "".join(
            "  - name: a\n    geometry: {bbox: [139.7, 35.3, 139.8, 35.4]}\n"
            "    target_h_m: 30\n" for _ in range(2))
    )
    with pytest.raises(ValueError):
        load_refine(p)


# --- pre-flight -------------------------------------------------------------

def test_preflight_passes_and_reports_the_derived_width():
    # 2 m of water is shallow enough that a 30 m target clears a 4.5 s floor.
    rep = preflight(_region(), gradation=GRAD, dt_expected_s=4.5,
                    ambient_h_m=AMBIENT, depth_of=_flat_depth(2.0), land=None)
    assert rep["transition_m"] == pytest.approx(1939.4, abs=0.1)
    assert rep["elements_core"] > 700
    assert rep["core_on_land"] == 0


def test_preflight_reports_the_altitude_time_step_not_the_edge_one():
    # The reported dt uses the minimum ALTITUDE (notebook 392), which for an
    # equilateral cell is sqrt(3)/2 of the edge. Quoting the edge figure
    # overstates the step by 15 % -- enough to clear a floor it does not meet.
    rep = preflight(_region(), gradation=GRAD, dt_expected_s=4.0,
                    ambient_h_m=AMBIENT, depth_of=_flat_depth(2.0), land=None)
    edge = 30 / np.sqrt(9.81 * 2.0)
    assert rep["dt_by_shortest_edge_s"] == pytest.approx(edge, rel=1e-9)
    assert rep["dt_s"] == pytest.approx(np.sqrt(3) / 2 * edge, rel=1e-9)
    assert rep["dt_s"] < rep["dt_by_shortest_edge_s"]


def test_a_short_time_step_alerts_but_does_not_veto():
    # The region is given: its position and resolution are inputs, so a time
    # step the caller did not expect is news, not grounds for refusal.
    rep = preflight(_region(), gradation=GRAD, dt_expected_s=4.5,
                    ambient_h_m=AMBIENT, depth_of=_flat_depth(20.0), land=None)
    assert rep["dt_alert"] is not None
    assert "cost" in rep["dt_alert"]
    assert rep["dt_step_cost_factor"] > 2.0


def test_the_alert_names_the_target_that_would_keep_the_expected_step():
    rep = preflight(_region(), gradation=GRAD, dt_expected_s=4.5,
                    ambient_h_m=AMBIENT, depth_of=_flat_depth(20.0), land=None)
    needed = 4.5 * np.sqrt(9.81 * 20.0) / (np.sqrt(3) / 2)
    assert f"{needed:.0f} m" in rep["dt_alert"]


def test_the_futtsu_target_alerts_at_its_expected_step():
    # 30 m over 4.15 m of water: 4.70 s by edge, 4.07 s by altitude.
    rep = preflight(_region(), gradation=GRAD, dt_expected_s=4.5,
                    ambient_h_m=AMBIENT, depth_of=_flat_depth(4.15), land=None)
    assert rep["dt_alert"] is not None
    assert "4.07 s by minimum altitude" in rep["dt_alert"]
    assert rep["dt_s"] == pytest.approx(4.07, abs=0.01)


def test_no_alert_when_the_expected_step_is_met():
    rep = preflight(_region(), gradation=GRAD, dt_expected_s=4.0,
                    ambient_h_m=AMBIENT, depth_of=_flat_depth(2.0), land=None)
    assert rep["dt_alert"] is None


def test_preflight_refuses_a_transition_too_short_for_the_gradation():
    with pytest.raises(ValueError, match="shorter than"):
        preflight(_region(transition_m=500), gradation=GRAD,
                  dt_expected_s=4.5, ambient_h_m=AMBIENT, depth_of=_flat_depth(2.0), land=None)


def test_preflight_accepts_a_transition_wider_than_required():
    rep = preflight(_region(transition_m=3000), gradation=GRAD,
                    dt_expected_s=4.5, ambient_h_m=AMBIENT, depth_of=_flat_depth(2.0),
                    land=None)
    assert rep["transition_m"] == 3000
    assert rep["transition_required_m"] == pytest.approx(1939.4, abs=0.1)


def test_a_partly_dry_core_alerts_but_does_not_veto():
    # The region is given; a dry patch of it is news, not grounds for refusal.
    land = _land_north_of(35.3228)          # covers the northern half of the core
    rep = preflight(_region(), gradation=GRAD, dt_expected_s=4.5,
                    ambient_h_m=AMBIENT, depth_of=_flat_depth(2.0), land=land)
    assert rep["core_on_land"] > 0
    assert rep["land_alert"] is not None


def test_no_land_alert_for_a_wet_core():
    rep = preflight(_region(), gradation=GRAD, dt_expected_s=4.5,
                    ambient_h_m=AMBIENT, depth_of=_flat_depth(2.0),
                    land=_land_north_of(35.40))
    assert rep["core_on_land"] == 0
    assert rep["land_alert"] is None


def test_preflight_refuses_a_core_entirely_on_land():
    land = _land_north_of(35.0)
    with pytest.raises(ValueError, match="entirely on land"):
        preflight(_region(), gradation=GRAD, dt_expected_s=4.5,
                  ambient_h_m=AMBIENT, depth_of=_flat_depth(2.0), land=land)


def test_coastline_mode_defaults_to_preserve(tmp_path):
    mesh = tmp_path / "base.14"
    mesh.write_text("stub\n")
    p = tmp_path / "r.yaml"
    p.write_text(
        "base_mesh: base.14\ndt_expected_s: 4.5\ngradation: 0.165\n"
        "refine:\n  - name: a\n"
        "    geometry: {circle: {center: [139.79, 35.32], radius_m: 300}}\n"
        "    target_h_m: 30\n"
    )
    cfg = load_refine(p)
    assert cfg["coastline"] == "preserve"
    assert cfg["coastline_tolerance_m"] == 100.0


def test_an_unknown_coastline_mode_is_rejected(tmp_path):
    mesh = tmp_path / "base.14"
    mesh.write_text("stub\n")
    p = tmp_path / "r.yaml"
    p.write_text(
        "base_mesh: base.14\ndt_expected_s: 4.5\ngradation: 0.165\n"
        "coastline: redraw\n"
        "refine:\n  - name: a\n"
        "    geometry: {bbox: [139.7, 35.3, 139.8, 35.4]}\n    target_h_m: 30\n"
    )
    with pytest.raises(ValueError, match="coastline must be one of"):
        load_refine(p)


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


# --- depths come from the base mesh, never re-sampled ------------------------

def _one_element():
    nodes = np.array([[0.0, 0.0], [100.0, 0.0], [0.0, 100.0]])
    tri = np.array([[0, 1, 2]], dtype=np.int64)
    return nodes, tri


def test_retained_nodes_keep_their_exact_depth():
    from fvcom_mesh_tools.refine import depths_from_base

    nodes, tri = _one_element()
    dep = np.array([5.0, 40.0, 12.0])
    out, outside = depths_from_base(nodes, tri, dep, nodes)
    assert outside == 0
    assert out == pytest.approx(dep)


def test_a_new_node_gets_the_interpolated_base_depth():
    from fvcom_mesh_tools.refine import depths_from_base

    nodes, tri = _one_element()
    dep = np.array([0.0, 30.0, 60.0])
    mid = np.array([[50.0, 0.0], [0.0, 50.0], [100.0 / 3, 100.0 / 3]])
    out, _ = depths_from_base(nodes, tri, dep, mid)
    assert out == pytest.approx([15.0, 30.0, 30.0])


def test_refinement_cannot_worsen_the_r_factor():
    # A value interpolated inside an element lies between that element's own
    # vertex depths, so any new pair's r is bounded by the element's worst
    # edge r. Refining the mesh can only leave the seabed slope alone or ease
    # it -- which is why the bathymetry needs no constrained re-smoothing.
    from fvcom_mesh_tools.refine import depths_from_base

    nodes, tri = _one_element()
    rng = np.random.default_rng(0)
    for _ in range(200):
        dep = rng.uniform(2.0, 300.0, 3)
        r_base = max(abs(dep[i] - dep[j]) / (dep[i] + dep[j])
                     for i, j in ((0, 1), (1, 2), (2, 0)))
        w = rng.dirichlet(np.ones(3), size=2)
        pts = w @ nodes
        out, _ = depths_from_base(nodes, tri, dep, pts)
        r_new = abs(out[0] - out[1]) / (out[0] + out[1])
        assert r_new <= r_base + 1e-12


def test_a_point_outside_the_base_mesh_falls_back_to_the_nearest_node():
    from fvcom_mesh_tools.refine import depths_from_base

    nodes, tri = _one_element()
    dep = np.array([5.0, 40.0, 12.0])
    out, outside = depths_from_base(nodes, tri, dep, np.array([[-10.0, -10.0]]))
    assert outside == 1
    assert out[0] == pytest.approx(5.0)


def test_base_depths_must_match_the_base_nodes():
    from fvcom_mesh_tools.refine import depths_from_base

    nodes, tri = _one_element()
    with pytest.raises(ValueError, match="one per base node"):
        depths_from_base(nodes, tri, np.array([1.0, 2.0]), nodes)


def test_base_depths_must_be_finite():
    from fvcom_mesh_tools.refine import depths_from_base

    nodes, tri = _one_element()
    with pytest.raises(ValueError, match="finite"):
        depths_from_base(nodes, tri, np.array([1.0, np.nan, 3.0]), nodes)


def test_the_base_r_factor_property_is_inherited_not_just_its_values():
    """A product's guarantee is part of what is inherited.

    `m7001tp_rfac0p2_cap300` means every base edge satisfies r <= 0.2.
    Interpolation does not carry that across a new edge joining different
    base elements: on the Futtsu patch ten new edges came out above it, the
    worst at 0.3075, while every wholly frozen edge stayed at 0.2.
    """
    import numpy as np

    from fvcom_mesh_tools.refine import limit_rfactor

    # a 3-node chain: the two ends are frozen at 3 and 6 m (r = 1/3), the
    # middle is new and interpolation put it at 3.2
    elements = np.array([[0, 1, 2], [1, 3, 2]])
    depths = np.array([3.0, 3.2, 6.0, 5.5])
    movable = np.array([False, True, False, False])
    out, info = limit_rfactor(elements, depths, movable, 0.2,
                              depth_min=3.0, depth_max=300.0)
    assert np.array_equal(out[~movable], depths[~movable]), "a base depth moved"
    assert info["n_depths_changed"] == 1
    e = np.array([[0, 1], [1, 2], [1, 3]])
    r = np.abs(out[e[:, 0]] - out[e[:, 1]]) / (out[e[:, 0]] + out[e[:, 1]])
    # the edge between the two FROZEN nodes cannot be fixed and is not counted
    assert info["max_r_on_new_edges"] <= 0.2 + 1e-9
    assert r.max() <= 0.2 + 1e-9


def test_limit_rfactor_never_moves_a_frozen_depth():
    import numpy as np

    from fvcom_mesh_tools.refine import limit_rfactor

    rng = np.random.default_rng(0)
    n = 30
    elements = np.array([[i, (i + 1) % n, (i + 2) % n] for i in range(n)])
    depths = 3.0 + 40.0 * rng.random(n)
    movable = rng.random(n) < 0.5
    out, info = limit_rfactor(elements, depths, movable, 0.2, depth_min=3.0)
    assert np.array_equal(out[~movable], depths[~movable])
    assert info["max_frozen_depth_change_m"] == 0.0


def test_an_impossible_r_factor_request_is_reported_not_hidden():
    """A movable node between two frozen depths too far apart has no answer."""
    import numpy as np

    from fvcom_mesh_tools.refine import limit_rfactor

    elements = np.array([[0, 1, 2]])
    depths = np.array([3.0, 6.0, 30.0])
    movable = np.array([False, True, False])
    _, info = limit_rfactor(elements, depths, movable, 0.2, depth_min=3.0,
                            rounds=50)
    assert not info["converged"]
    assert info["n_edges_over_rmax"] > 0


# --------------------------------------------------- a polygon from a file
#
# A real fishery boundary arrives as a shapefile or GeoJSON. Typing its
# coordinates into a recipe is tedious and is a place to put a typo nobody
# will ever find.


def _fishery(tmp_path, rows, crs="EPSG:4326", name="fishery.geojson"):
    import geopandas as gpd
    from shapely.geometry import Polygon

    gdf = gpd.GeoDataFrame(
        {k: [r[k] for r in rows] for k in rows[0] if k != "xy"},
        geometry=[Polygon(r["xy"]) for r in rows], crs=crs)
    path = tmp_path / name
    gdf.to_file(path)
    return path


def _square(lon, lat, half=0.002):
    return [(lon - half, lat - half), (lon + half, lat - half),
            (lon + half, lat + half), (lon - half, lat + half),
            (lon - half, lat - half)]


def test_a_polygon_is_read_from_a_file(tmp_path):
    from fvcom_mesh_tools.refine import RefineRegion

    path = _fishery(tmp_path, [{"NAME": "nori", "xy": _square(139.788, 35.323)}])
    r = RefineRegion({"name": "futtsu", "target_h_m": 30,
                      "geometry": {"file": str(path)}})
    assert r.kind == "file"
    assert r.geometry.geom_type == "Polygon"
    assert r.geometry.contains(shapely.Point(139.788, 35.323))
    assert r.source["file"] == str(path)
    assert r.source["n_vertices"] == 5


def test_a_file_polygon_is_reprojected_to_lon_lat(tmp_path):
    """A boundary in UTM is still a boundary; the vocabulary is lon/lat."""
    from pyproj import Transformer

    from fvcom_mesh_tools.refine import RefineRegion

    to_m = Transformer.from_crs(4326, 32654, always_xy=True)
    xy = [to_m.transform(x, y) for x, y in _square(139.788, 35.323)]
    path = _fishery(tmp_path, [{"NAME": "nori", "xy": xy}], crs="EPSG:32654",
                    name="utm.geojson")
    r = RefineRegion({"name": "futtsu", "target_h_m": 30,
                      "geometry": {"file": str(path)}})
    assert r.geometry.contains(shapely.Point(139.788, 35.323))


def test_an_attribute_filter_picks_the_feature(tmp_path):
    from fvcom_mesh_tools.refine import RefineRegion

    path = _fishery(tmp_path, [
        {"NAME": "nori", "xy": _square(139.788, 35.323)},
        {"NAME": "kaki", "xy": _square(139.900, 35.500)}])
    with pytest.raises(ValueError, match="selection is 2 features"):
        RefineRegion({"name": "x", "target_h_m": 30,
                      "geometry": {"file": str(path)}})
    r = RefineRegion({"name": "x", "target_h_m": 30,
                      "geometry": {"file": str(path), "where": {"NAME": "kaki"}}})
    assert r.geometry.contains(shapely.Point(139.900, 35.500))
    assert r.source["where"] == {"NAME": "kaki"}


def test_a_row_index_picks_the_feature(tmp_path):
    from fvcom_mesh_tools.refine import RefineRegion

    path = _fishery(tmp_path, [
        {"NAME": "a", "xy": _square(139.788, 35.323)},
        {"NAME": "b", "xy": _square(139.900, 35.500)}])
    r = RefineRegion({"name": "x", "target_h_m": 30,
                      "geometry": {"file": str(path), "index": 1}})
    assert r.geometry.contains(shapely.Point(139.900, 35.500))
    with pytest.raises(ValueError, match="outside the 2 selected"):
        RefineRegion({"name": "x", "target_h_m": 30,
                      "geometry": {"file": str(path), "index": 5}})


def test_an_unknown_column_says_which_columns_there_are(tmp_path):
    from fvcom_mesh_tools.refine import RefineRegion

    path = _fishery(tmp_path, [{"NAME": "nori", "xy": _square(139.788, 35.323)}])
    with pytest.raises(ValueError, match="no column 'FISHERY'"):
        RefineRegion({"name": "x", "target_h_m": 30,
                      "geometry": {"file": str(path), "where": {"FISHERY": "n"}}})


def test_a_shapefile_without_a_prj_is_refused(tmp_path):
    """GeoJSON is WGS84 by definition; a shapefile is whatever its .prj says.

    A boundary whose projection nobody recorded cannot be placed, and
    guessing lon/lat would put a UTM polygon in the Gulf of Guinea.
    """
    from fvcom_mesh_tools.refine import RefineRegion

    path = _fishery(tmp_path, [{"NAME": "n", "xy": _square(139.788, 35.323)}],
                    name="nocrs.shp")
    path.with_suffix(".prj").unlink()
    with pytest.raises(ValueError, match="no CRS"):
        RefineRegion({"name": "x", "target_h_m": 30,
                      "geometry": {"file": str(path)}})


def test_several_disjoint_parts_are_several_regions(tmp_path):
    """Taking the largest silently would be worse than refusing."""
    import geopandas as gpd
    from shapely.geometry import MultiPolygon, Polygon

    from fvcom_mesh_tools.refine import RefineRegion

    path = tmp_path / "multi.geojson"
    gpd.GeoDataFrame(
        {"NAME": ["both"]},
        geometry=[MultiPolygon([Polygon(_square(139.788, 35.323)),
                                Polygon(_square(139.900, 35.500))])],
        crs="EPSG:4326").to_file(path)
    with pytest.raises(ValueError, match="MultiPolygon of 2 parts"):
        RefineRegion({"name": "x", "target_h_m": 30,
                      "geometry": {"file": str(path)}})


def test_a_buffer_grows_the_polygon_by_metres(tmp_path):
    from fvcom_mesh_tools.refine import RefineRegion

    path = _fishery(tmp_path, [{"NAME": "n", "xy": _square(139.788, 35.323)}])
    plain = RefineRegion({"name": "x", "target_h_m": 30,
                          "geometry": {"file": str(path)}}).geometry
    grown = RefineRegion({"name": "x", "target_h_m": 30,
                          "geometry": {"file": str(path), "buffer_m": 100}}).geometry
    assert grown.contains(plain)
    # a 100 m buffer on a ~360 x 440 m square roughly doubles its area
    assert 1.8 < grown.area / plain.area < 2.6


def test_a_missing_geometry_file_says_so(tmp_path):
    from fvcom_mesh_tools.refine import RefineRegion

    with pytest.raises(ValueError, match="geometry file not found"):
        RefineRegion({"name": "x", "target_h_m": 30,
                      "geometry": {"file": str(tmp_path / "nope.geojson")}})
