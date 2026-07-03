"""Tests for the unified QA gate (fvcom_mesh_tools.qa / fmesh-mesh-qa).

Fixtures are small metric-coordinate grids built in memory. The
pristine 4x4-node grid (1 km squares, right triangles) is constructed
to pass every gated check: the diagonal of the square touching the
upper open/land junction is flipped so that no element carries the
junction node together with two solid nodes (the {2,1,1} ISBCE=2
mis-classification), mirroring what a real mesh must do near an OBC
junction.
"""

from __future__ import annotations

import json

import numpy as np
import pytest

from fvcom_mesh_tools.cli import meshqa
from fvcom_mesh_tools.diagnostics import channel_width_metric
from fvcom_mesh_tools.io import Fort14Mesh, write_fort14
from fvcom_mesh_tools.qa import compute_isonb, detect_coords, format_report, run_qa

N = 4  # nodes per side
SPACING = 1000.0  # metres


def _nid(i: int, j: int) -> int:
    return j * N + i


def _grid_mesh(
    *,
    shear: float = 0.0,
    obc: list[int] | None = None,
    land: list[int] | None = None,
    flip_squares: set[tuple[int, int]] | None = None,
    depth: float = 5.0,
) -> Fort14Mesh:
    """Right-triangle grid over an (N-1) x (N-1) square arrangement."""
    flip_squares = flip_squares or set()
    nodes = np.array(
        [[i * SPACING + shear * j * SPACING, j * SPACING]
         for j in range(N) for i in range(N)],
        dtype=np.float64,
    )
    elements = []
    for j in range(N - 1):
        for i in range(N - 1):
            a, b = _nid(i, j), _nid(i + 1, j)
            c, d = _nid(i + 1, j + 1), _nid(i, j + 1)
            if (i, j) in flip_squares:
                elements.append([a, b, d])
                elements.append([b, c, d])
            else:
                elements.append([a, b, c])
                elements.append([a, c, d])
    open_boundaries = [np.asarray(obc, dtype=np.int64)] if obc else []
    if land is None:
        # Full outer loop from the top of the default left-side OBC
        # around to its bottom.
        land = (
            [_nid(0, 2), _nid(0, 3)]
            + [_nid(i, 3) for i in range(1, N)]
            + [_nid(N - 1, j) for j in range(N - 2, -1, -1)]
            + [_nid(i, 0) for i in range(N - 2, -1, -1)]
            + [_nid(0, 1)]
        )
    return Fort14Mesh(
        title="qa-test",
        nodes=nodes,
        depths=np.full(len(nodes), depth),
        elements=np.asarray(elements, dtype=np.int64),
        open_boundaries=open_boundaries,
        land_boundaries=[(0, np.asarray(land, dtype=np.int64))],
    )


def _pristine() -> Fort14Mesh:
    return _grid_mesh(
        obc=[_nid(0, 1), _nid(0, 2)],
        flip_squares={(0, 2)},
    )


def _check(report, check_id):
    matches = [c for c in report.checks if c.check_id == check_id]
    assert len(matches) == 1, f"check {check_id} missing"
    return matches[0]


# ---------------------------------------------------------------------------
# Happy path
# ---------------------------------------------------------------------------


def test_pristine_mesh_passes_all_gates():
    report = run_qa(_pristine())
    failed = [c.check_id for c in report.checks
              if c.gate and not c.skipped and not c.passed]
    assert failed == []
    assert report.passed
    assert report.coords == "metric"
    assert report.n_obc_nodes == 2


def test_pristine_informational_values():
    report = run_qa(_pristine())
    dt = _check(report, "implied_dt")
    # min edge 1000 m, depth 5 m -> 1000 / sqrt(9.81 * 5) ~ 142.9 s
    assert dt.data["dt_min_s"] == pytest.approx(142.9, abs=1.0)
    delaunay = _check(report, "delaunay_local")
    assert delaunay.data["n_non_delaunay"] == 0
    perp = _check(report, "obc_perpendicularity")
    assert perp.data["worst_deviation_deg"] == pytest.approx(0.0, abs=1e-9)


def test_detect_coords():
    assert detect_coords(np.array([[139.8, 35.5], [140.0, 35.6]])) == "lonlat"
    assert detect_coords(np.array([[380000.0, 3.9e6]])) == "metric"


def test_channel_width_metric_metric_coords():
    info = channel_width_metric(_pristine(), coords="metric")
    ratio = info["w_h_ratio"]
    finite = ratio[np.isfinite(ratio)]
    assert finite.size > 0
    # Domain is 3 km wide with ~1-1.4 km elements: widths must be of
    # kilometre order, not degree-projection garbage.
    widths = info["channel_width_m"][np.isfinite(info["channel_width_m"])]
    assert widths.min() > 100.0
    assert widths.max() < 20000.0


def test_compute_isonb_matches_tge_semantics():
    mesh = _pristine()
    from fvcom_mesh_tools.qa import _edge_topology

    topo = _edge_topology(mesh.elements, mesh.n_nodes)
    isonb = compute_isonb(
        mesh.n_nodes, topo.uv[topo.counts == 1], mesh.open_boundaries,
    )
    assert isonb[_nid(0, 1)] == 2 and isonb[_nid(0, 2)] == 2
    assert isonb[_nid(0, 0)] == 1 and isonb[_nid(0, 3)] == 1
    assert isonb[_nid(1, 1)] == 0 and isonb[_nid(2, 2)] == 0


# ---------------------------------------------------------------------------
# FVCOM-fatal failures
# ---------------------------------------------------------------------------


def test_r4_mixed_boundary_detected():
    # Standard diagonals + corner node in the OBC list: the top-left
    # corner triangle gets an open edge plus a solid third node
    # (ISONB sum = 5) -- FVCOM's tge.F PSTOP case.
    land = (
        [_nid(0, 3)]
        + [_nid(i, 3) for i in range(1, N)]
        + [_nid(N - 1, j) for j in range(N - 2, -1, -1)]
        + [_nid(i, 0) for i in range(N - 2, -1, -1)]
        + [_nid(0, 1)]
    )
    mesh = _grid_mesh(obc=[_nid(0, 1), _nid(0, 2), _nid(0, 3)], land=land)
    report = run_qa(mesh, channel_check=False)
    r4 = _check(report, "r4_mixed_boundary")
    assert not r4.passed
    assert r4.n_violations >= 1
    assert not report.passed


def test_flipped_element_detected():
    mesh = _pristine()
    mesh.elements[0] = mesh.elements[0][[0, 2, 1]]
    report = run_qa(mesh, channel_check=False)
    ccw = _check(report, "ccw_all_elements")
    assert not ccw.passed
    assert ccw.n_violations == 1
    assert ccw.offenders[0]["id"] == 0


def test_duplicate_and_degenerate_nodes_detected():
    mesh = _pristine()
    mesh.nodes[_nid(1, 1)] = mesh.nodes[_nid(2, 1)]
    report = run_qa(mesh, channel_check=False)
    assert not _check(report, "no_duplicate_nodes").passed
    assert not _check(report, "no_tiny_area").passed


def test_orphan_node_detected():
    base = _pristine()
    mesh = Fort14Mesh(
        title=base.title,
        nodes=np.vstack([base.nodes, [[10000.0, 10000.0]]]),
        depths=np.append(base.depths, 5.0),
        elements=base.elements,
        open_boundaries=base.open_boundaries,
        land_boundaries=base.land_boundaries,
    )
    report = run_qa(mesh, channel_check=False)
    orphan = _check(report, "no_orphan_nodes")
    assert not orphan.passed
    assert orphan.offenders[0]["id"] == mesh.n_nodes - 1


def test_bowtie_mesh_fails_manifold_and_isolation():
    nodes = np.array([
        [0.0, 0.0], [1000.0, 0.0], [500.0, 500.0],
        [0.0, 1000.0], [1000.0, 1000.0],
    ])
    elements = np.array([[0, 1, 2], [2, 4, 3]])
    mesh = Fort14Mesh(
        title="bowtie",
        nodes=nodes,
        depths=np.full(5, 5.0),
        elements=elements,
        open_boundaries=[],
        land_boundaries=[(0, np.array([0, 1, 2, 4, 3, 2, 0]))],
    )
    report = run_qa(mesh, channel_check=False)
    assert not _check(report, "manifold_boundary").passed
    assert not _check(report, "no_isolated_elements").passed


# ---------------------------------------------------------------------------
# OBC failures
# ---------------------------------------------------------------------------


def test_broken_obc_chain_detected():
    mesh = _grid_mesh(
        obc=[_nid(0, 1), _nid(0, 3)],  # not mesh-adjacent
        flip_squares={(0, 2)},
    )
    report = run_qa(mesh, channel_check=False)
    chain = _check(report, "obc_chain_adjacency")
    order = _check(report, "obc_ordering")
    assert not chain.passed
    assert chain.n_violations == 2
    assert not order.passed


def test_obc_perpendicularity_gate():
    # Sheared grid, OBC on the (unsheared) bottom row: the best
    # interior edge per OBC node deviates atan(0.4) ~ 21.8 deg.
    land = (
        [_nid(2, 0)]
        + [_nid(N - 1, j) for j in range(N)]
        + [_nid(i, 3) for i in range(N - 2, -1, -1)]
        + [_nid(0, j) for j in range(N - 2, -1, -1)]
        + [_nid(1, 0)]
    )
    mesh = _grid_mesh(shear=0.4, obc=[_nid(1, 0), _nid(2, 0)], land=land)
    report = run_qa(mesh, channel_check=False)
    perp = _check(report, "obc_perpendicularity")
    assert not perp.passed
    assert perp.data["worst_deviation_deg"] == pytest.approx(21.8, abs=0.1)

    relaxed = run_qa(mesh, channel_check=False, max_obc_perp_dev_deg=25.0)
    assert _check(relaxed, "obc_perpendicularity").passed


def test_obc_checks_skipped_without_open_boundary():
    land = (
        [_nid(0, 0)]
        + [_nid(0, j) for j in range(1, N)]
        + [_nid(i, 3) for i in range(1, N)]
        + [_nid(N - 1, j) for j in range(N - 2, -1, -1)]
        + [_nid(i, 0) for i in range(N - 2, 0, -1)]
        + [_nid(0, 0)]
    )
    mesh = _grid_mesh(obc=None, land=land)
    report = run_qa(mesh, channel_check=False)
    assert _check(report, "obc_perpendicularity").skipped
    assert _check(report, "obc_reachable").skipped
    assert report.passed  # skipped checks do not fail the gate


# ---------------------------------------------------------------------------
# Quality criteria failures
# ---------------------------------------------------------------------------


def test_c1_min_angle_detected():
    mesh = _pristine()
    mesh.nodes[_nid(1, 1)] = [1300.0, 700.0]
    report = run_qa(mesh, channel_check=False)
    c1 = _check(report, "c1_min_angle")
    assert not c1.passed
    assert c1.data["min_angle_deg"] < 30.0
    assert _check(report, "ccw_all_elements").passed


def test_c4_area_change_detected():
    nodes = np.array([
        [0.0, 0.0], [1000.0, 0.0], [500.0, 1000.0], [500.0, -3000.0],
    ])
    mesh = Fort14Mesh(
        title="c4",
        nodes=nodes,
        depths=np.full(4, 5.0),
        elements=np.array([[0, 1, 2], [0, 3, 1]]),
        open_boundaries=[],
        land_boundaries=[(0, np.array([0, 3, 1, 2, 0]))],
    )
    report = run_qa(mesh, channel_check=False)
    c4 = _check(report, "c4_area_change")
    assert not c4.passed
    assert c4.data["max_area_change"] == pytest.approx(2.0 / 3.0, abs=1e-6)


def test_c5_valence_detected():
    n_ring = 10
    ang = np.deg2rad(36.0 * np.arange(n_ring))
    nodes = np.vstack([
        [[0.0, 0.0]],
        np.column_stack([1000.0 * np.cos(ang), 1000.0 * np.sin(ang)]),
    ])
    elements = np.array([[0, i, i + 1] for i in range(1, n_ring)])
    mesh = Fort14Mesh(
        title="fan",
        nodes=nodes,
        depths=np.full(n_ring + 1, 5.0),
        elements=elements,
        open_boundaries=[],
        land_boundaries=[(0, np.arange(0, n_ring + 1))],
    )
    report = run_qa(mesh, channel_check=False)
    c5 = _check(report, "c5_valence")
    assert not c5.passed
    assert c5.data["max_valence"] == 9
    assert c5.offenders[0]["id"] == 0


def test_min_depth_detected():
    mesh = _pristine()
    mesh.depths[_nid(1, 1)] = 1.0
    report = run_qa(mesh, channel_check=False)
    depth = _check(report, "min_depth_clip")
    assert not depth.passed
    assert depth.data["min_depth_m"] == pytest.approx(1.0)


# ---------------------------------------------------------------------------
# Report rendering + CLI
# ---------------------------------------------------------------------------


def test_format_report_both_languages():
    report = run_qa(_pristine(), channel_check=False)
    ja = format_report(report, lang="ja")
    en = format_report(report, lang="en")
    assert "総合判定" in ja and "PASS" in ja
    assert "overall: PASS" in en
    with pytest.raises(ValueError):
        format_report(report, lang="de")


def test_report_json_roundtrip():
    report = run_qa(_pristine(), channel_check=False)
    payload = json.loads(json.dumps(report.to_dict()))
    assert payload["passed"] is True
    assert payload["mesh"]["n_elements"] == 18
    ids = {c["check_id"] for c in payload["checks"]}
    assert {"r4_mixed_boundary", "c1_min_angle", "implied_dt"} <= ids


def test_cli_pass_and_fail(tmp_path, capsys):
    good = tmp_path / "good.14"
    write_fort14(_pristine(), good)
    assert meshqa.main([str(good), "--no-channel"]) == 0
    assert (tmp_path / "good_qa.json").exists()
    out = capsys.readouterr().out
    assert "総合判定" in out

    bad_mesh = _pristine()
    bad_mesh.elements[0] = bad_mesh.elements[0][[0, 2, 1]]
    bad = tmp_path / "bad.14"
    write_fort14(bad_mesh, bad)
    assert meshqa.main([str(bad), "--no-channel", "--lang", "en"]) == 1
    out = capsys.readouterr().out
    assert "overall: FAIL" in out
    payload = json.loads((tmp_path / "bad_qa.json").read_text())
    assert payload["passed"] is False


def test_cli_missing_input(tmp_path):
    assert meshqa.main([str(tmp_path / "nope.14")]) == 2
