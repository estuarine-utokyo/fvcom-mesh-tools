"""Fifth-review reproductions. Run from the repository root with pytest.

Tests named ``test_defect_*`` assert the desired contract and intentionally
fail on 95538f5. Other tests are controls and read-only artifact measurements.
No mesher, FVCOM run, or production-output write is performed.
"""

import ast
import json
from dataclasses import replace
from pathlib import Path
from types import SimpleNamespace

import numpy as np
import pytest
import shapely

from fvcom_mesh_tools import patch
from fvcom_mesh_tools.io.fort14 import Fort14Mesh, read_fort14
from fvcom_mesh_tools.io.fvcom_native import (
    apply_obc_depth_control,
    export_fvcom_case,
    read_fvcom_case,
    read_obc_types,
)
from fvcom_mesh_tools.mesh_clean import compact_nodes
from fvcom_mesh_tools.qa import QACheck, run_qa
from fvcom_mesh_tools.refine import limit_rfactor

ROOT = Path(__file__).resolve().parents[1]


def square():
    return Fort14Mesh(
        "review", np.array([[0., 0.], [100., 0.], [100., 100.], [0., 100.]]),
        np.full(4, 8.), np.array([[0, 1, 2], [0, 2, 3]]),
        [np.array([0, 1])], [(20, np.array([1, 2, 3, 0]))], obc_type=3,
    )


def gate(regions):
    """Execute the actual driver function and tolerance, without its top level."""
    path = ROOT / "notebooks/420_local_refine.py"
    tree = ast.parse(path.read_text())
    body = [n for n in tree.body if
            (isinstance(n, ast.FunctionDef) and n.name == "achieved_per_region") or
            (isinstance(n, ast.Assign) and any(
                isinstance(t, ast.Name) and t.id == "RESOLUTION_TOLERANCE"
                for t in n.targets))]
    env = dict(np=np, shapely=shapely, regions_m=regions)
    exec(compile(ast.Module(body=body, type_ignores=[]), str(path), "exec"), env)
    return env["achieved_per_region"]


def region(geom, target=30.):
    return geom, SimpleNamespace(name="requested", target_h_m=target)


@pytest.mark.parametrize("cap", [0, 1, 100])
def test_real_qa_attribution_is_cap_independent(cap):
    mesh = square()
    # Force node, element, and edge checks to fail, including optional dt gate.
    mesh.nodes[3, 1] = 20
    qa = run_qa(mesh, coords="metric", max_offenders=cap, min_angle_deg=80,
                max_angle_deg=85, max_area_change=.1, max_valence=1,
                min_depth_m=9, min_dt_s=100)
    for c in qa.checks:
        if c.status == "fail":
            assert len(c.offender_ids) == c.n_violations, c.check_id
    actual = patch.introduced_violations(qa.checks, 1, mesh.elements)
    reference = run_qa(mesh, coords="metric", max_offenders=100, min_angle_deg=80,
                       max_angle_deg=85, max_area_change=.1, max_valence=1,
                       min_depth_m=9, min_dt_s=100)
    expected = patch.introduced_violations(reference.checks, 1, mesh.elements)
    def key(v):
        return v["check"], patch._offender_key(v)
    assert sorted(map(key, actual)) == sorted(map(key, expected))


def test_duplicate_obc_is_counted_but_not_identified_and_still_blocks():
    mesh = square()
    mesh.open_boundaries = [np.array([0, 1, 0])]
    c = next(c for c in run_qa(mesh, coords="metric").checks
             if c.check_id == "obc_ordering")
    assert c.n_violations == 1 and c.offender_ids == []
    out = patch.introduced_violations([c], 2, mesh.elements)
    assert out[0]["kind"] == "unattributed" and out[0]["n_unattributed"] == 1


def test_defect_fallback_keeps_incident_elements():
    c = QACheck("edge_check", "quality", True, False, "test", "test", 1,
                offenders=[{"kind": "edge", "id": [0, 2], "elements": [0, 1]}])
    print("fallback identities", c.offender_ids)
    assert patch.introduced_violations([c], 2, square().elements) == []


def test_obc_key_collision_changes_decoration_not_gate():
    mesh = square()
    mesh.open_boundaries = [np.array([1, 3]), np.array([1, 3])]
    c = next(c for c in run_qa(mesh, coords="metric").checks
             if c.check_id == "obc_ordering")
    assert [o["segment"] for o in c.offender_ids] == [0, 1]
    out = patch.introduced_violations([c], 2, mesh.elements)
    assert len(out) == 2
    assert [o["segment"] for o in out] == [0, 0]


@pytest.mark.parametrize("ident", [False, np.bool_(True), -.5, -1, 2, "0", None])
def test_malformed_element_and_incident_ids_block(ident):
    for off in [{"kind": "element", "id": ident},
                {"kind": "edge", "id": [0, 2], "elements": [0, ident]}]:
        c = QACheck("test", "quality", True, False, "test", "test", 1,
                    offenders=[], offender_ids=[off])
        assert patch.introduced_violations([c], 2, square().elements)


def test_defect_unhashable_malformed_id_is_conservative():
    off = {"kind": "element", "id": {"index": 0}}
    c = QACheck("test", "quality", True, False, "test", "test", 1,
                offenders=[off])
    assert patch.introduced_violations([c], 2, square().elements)


def test_type_survives_control_and_native_roundtrip(tmp_path):
    mesh, _ = apply_obc_depth_control(square())
    files = export_fvcom_case(mesh, tmp_path, "kept", obc_depth_control=False)
    back = read_fvcom_case(files["grd"], files["dep"], files["obc"])
    assert mesh.obc_type == back.obc_type == 3


def test_defect_compaction_preserves_type(tmp_path):
    mesh = square()
    mesh = replace(mesh, nodes=np.vstack([mesh.nodes, [200., 200.]]),
                   depths=np.append(mesh.depths, 8.))
    compact, info = compact_nodes(mesh)
    assert info["n_orphans_removed"] == 1
    assert np.array_equal(compact.open_boundaries[0], mesh.open_boundaries[0])
    files = export_fvcom_case(compact, tmp_path, "compacted", obc_depth_control=False)
    assert read_obc_types(files["obc"]) == [3]


@pytest.mark.parametrize("x", [0., 4950.])
def test_defect_thin_conflict_area(x):
    coarse = shapely.box(0, 0, 10000, 10)
    fine = shapely.box(x, 0, x + 100, 10)
    regions = [(coarse, 30., 0.), (fine, 5., 0.)]
    exact = patch.region_conflicts(regions, ["coarse", "fine"])
    sampled = patch.region_conflicts(
        regions, ["coarse", "fine"], base_size=lambda p: np.full(len(p), 300.))
    assert exact["finer_than_declared"]["coarse"]["fraction"] == .01
    print("sampled conflict", sampled)
    fraction = sampled["finer_than_declared"].get("coarse", {}).get("fraction", 0.)
    assert fraction == pytest.approx(.01, abs=.005)


def test_shared_field_matches_independent_linear_ramp():
    mesh = square()
    regions = [(shapely.box(20, 20, 30, 30), 5., 50.),
               (shapely.box(60, 20, 70, 30), 12., 30.)]
    p = np.array([[25., 25.], [45., 25.], [65., 25.], [80., 80.]])
    base = patch.base_size_field(mesh.nodes, mesh.elements)(p)
    expected = base.copy()
    for geom, target, width in regions:
        distance = shapely.distance(shapely.points(p), geom)
        expected = np.minimum(expected, target + (base-target)*np.clip(distance/width, 0, 1))
    assert patch.patch_sizing(mesh.nodes, mesh.elements, regions,
                              distmesh_scale=1.)(p) == pytest.approx(expected)


def test_defect_transverse_slope_not_reported_as_measured_flat():
    # A 40 m channel at the DEFAULT 25 m spacing: only one interior row.
    geom = shapely.box(0, 0, 1000, 40)
    def field(p):
        return 100. + p[:, 1]
    coarse = patch.field_gradation(field, geom)
    fine = patch.field_gradation(field, geom, spacing=5.)
    assert fine["max_slope"] == pytest.approx(1.)
    print("coarse transverse gradation", coarse)
    assert not coarse["measured"] or coarse["max_slope"] >= .99


def test_outside_mesh_values_do_not_enter_gradation():
    geom = shapely.box(0, 0, 100, 100)
    def field(p):
        return np.where(shapely.contains_xy(geom, p[:, 0], p[:, 1]),
                        100. + .2*p[:, 0], 10000.)
    measured = patch.field_gradation(field, geom, spacing=10.)
    assert measured["max_slope"] == pytest.approx(.2)
    assert measured["n_partial_samples"] == 0


@pytest.mark.parametrize("value", ["nan", "inf", "-inf"])
@pytest.mark.parametrize("axis", [1, 2])
def test_grid_nonfinite_either_coordinate_refused(tmp_path, value, axis):
    files = export_fvcom_case(square(), tmp_path, "invalid", obc_depth_control=False)
    rows = files["grd"].read_text().splitlines()
    fields = rows[4].split()
    fields[axis] = value
    rows[4] = " ".join(fields)
    files["grd"].write_text("\n".join(rows) + "\n")
    with pytest.raises(ValueError, match="non-finite"):
        read_fvcom_case(files["grd"], files["dep"], files["obc"])


@pytest.mark.parametrize("tol", [np.nan, np.inf, -1.])
def test_invalid_coordinate_tolerance_refused(tol):
    with pytest.raises(ValueError, match="coord_tol_m"):
        read_fvcom_case("unused", "unused", "unused", coord_tol_m=tol)


def test_limiter_unused_dry_node_and_positive_used_nodes():
    depths, info = limit_rfactor(np.array([[0, 1, 2]]), [3., 30., 3., 0.],
                                [False, True, False, False], .2)
    assert info["converged"] and depths[3] == 0.
    assert depths[0] == depths[2] == 3.
    assert max(abs(depths[i]-depths[j])/(depths[i]+depths[j])
               for i, j in [(0, 1), (1, 2), (0, 2)]) <= .2 + 1e-9


def test_defect_limiter_cannot_create_invalid_depth_and_certify_it():
    with np.errstate(invalid="ignore", divide="ignore"):
        try:
            depths, info = limit_rfactor(np.array([[0, 1, 2]]), [1., 10., 1.],
                                        [True, True, True], .2, depth_max=0.)
        except ValueError:
            return
    print("invalid limiter result", depths, info)
    assert not info["converged"]


def test_resolution_empty_and_uniformly_coarse_are_rejected():
    mesh = square()
    for geom in [shapely.box(1, 1, 2, 2), shapely.box(-1, -1, 101, 101)]:
        assert gate([region(geom)])(mesh)[1] == ["requested"]


@pytest.mark.parametrize("side, missed", [(31.49, False), (31.51, True)])
def test_resolution_tolerance_boundary(side, missed):
    mesh = Fort14Mesh("equilateral", np.array([
        [0., 0.], [side, 0.], [side/2, side*np.sqrt(3)/2]]),
        np.full(3, 8.), np.array([[0, 1, 2]]), [], [])
    assert bool(gate([region(shapely.box(-1, -1, 40, 40))])(mesh)[1]) == missed


def artifact(name):
    path = ROOT / f"outputs/v4_{name}/report.json"
    if not path.exists():
        pytest.skip("local production artifacts unavailable")
    report = json.loads(path.read_text())
    mesh_path = path.parent / Path(report["mesh"]).name
    if not mesh_path.exists():
        pytest.skip("local production mesh unavailable")
    return read_fort14(mesh_path), report


@pytest.mark.parametrize("name", ["futtsu_nori", "futtsu_nori_polygon", "futtsu_two_beds"])
def test_delivered_artifacts_measured_read_only(name):
    mesh, report = artifact(name)
    regions = [(shapely.Polygon(r["xy"]), SimpleNamespace(
        name=r["name"], target_h_m=r["target_h_m"])) for r in report["regions"]]
    stats, missed = gate(regions)(mesh)
    qa = run_qa(mesh, coords="metric", max_offenders=0)
    retained = report["selection"]["n_elements_retained"]
    introduced = patch.introduced_violations(qa.checks, retained, mesh.elements)
    tri = mesh.nodes[mesh.elements]
    sides = np.linalg.norm(tri - np.roll(tri, -1, axis=1), axis=2)
    u, v = tri[:, 1] - tri[:, 0], tri[:, 2] - tri[:, 0]
    double_area = abs(u[:, 0]*v[:, 1] - u[:, 1]*v[:, 0])
    dt = (double_area / sides.max(axis=1)) / np.sqrt(9.81*mesh.depths[mesh.elements].max(axis=1))
    print(name, dict(nodes=mesh.n_nodes, elements=mesh.n_elements,
                     gates=qa.n_gate_total, failed=qa.n_gate_failed,
                     introduced=len(introduced), dt=float(dt.min()), achieved=stats))
    assert not missed and not introduced
    assert float(dt.min()) == pytest.approx(report["achieved"]["dt_min_s"])


def test_defect_real_mesh_expanded_request_is_not_delivered():
    mesh, report = artifact("futtsu_nori")
    geom = shapely.Polygon(report["regions"][0]["xy"]).buffer(300.)
    stats, missed = gate([region(geom)])(mesh)
    # No mesh is edited. This deliberately changes only the requested geometry.
    assert stats["requested"]["p90_m"] > 60.
    assert stats["requested"]["max_m"] > 140.
    print("expanded request on actual mesh", stats)
    assert missed == ["requested"]


def test_defect_real_mesh_unsampled_lobe_is_not_hidden_by_fine_lobe():
    mesh, report = artifact("futtsu_nori")
    fine = shapely.Polygon(report["regions"][0]["xy"])
    c = fine.centroid
    remote = shapely.Point(c.x, c.y + 5000.).buffer(300., quad_segs=64)
    corridor = shapely.LineString([(c.x, c.y), (c.x, c.y + 5000.)]).buffer(5.)
    geom = shapely.union_all([fine, remote, corridor])
    triangles = shapely.polygons(mesh.nodes[mesh.elements])
    wet = shapely.union_all(triangles)
    assert geom.geom_type == "Polygon" and geom.is_valid
    assert geom.difference(wet).area == 0.
    assert remote.area / geom.area > .46
    alone, alone_missed = gate([region(remote)])(mesh)
    assert alone_missed and alone["requested"]["n_edges"] == 0
    stats, missed = gate([region(geom)])(mesh)
    print("unsampled lobe", dict(area_m2=geom.area, remote_area_m2=remote.area,
                                 outside_mesh_m2=geom.difference(wet).area,
                                 achieved=stats))
    assert missed == ["requested"]


def test_saved_futtsu_field_on_actual_patch_footprint():
    mesh, report = artifact("futtsu_nori")
    paths = [Path(report[k]) for k in ("base_mesh", "base_depth", "base_obc")]
    if not all(p.exists() for p in paths):
        pytest.skip("native base inputs unavailable")
    base = read_fvcom_case(*paths)
    assert np.isfinite(base.depths).all() and base.depths.min() > 0.
    retained = report["selection"]["n_elements_retained"]
    hole = shapely.union_all(shapely.polygons(mesh.nodes[mesh.elements[retained:]]))
    geom = shapely.Polygon(report["regions"][0]["xy"])
    regions = [(geom, 30., report["preflight"][0]["transition_m"])]
    field = patch.patch_sizing(base.nodes, base.elements, regions, distmesh_scale=1.)
    default = patch.field_gradation(field, hole)
    finer = patch.field_gradation(field, hole, spacing=default["spacing_m"]/2)
    print("real Futtsu field", dict(default=default, finer=finer,
                                    base_depth_range=[base.depths.min(), base.depths.max()]))
    assert default["measured"] and finer["measured"]
    assert .3 < default["max_slope"] < .5
    assert .3 < finer["max_slope"] < .5
