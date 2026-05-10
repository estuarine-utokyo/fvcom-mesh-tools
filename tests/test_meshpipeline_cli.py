"""Tests for ``fmesh-mesh-pipeline``."""

from __future__ import annotations

import json

import numpy as np

from fvcom_mesh_tools.cli import meshpipeline
from fvcom_mesh_tools.io import Fort14Mesh, write_fort14


def _square_mesh() -> Fort14Mesh:
    """A simple mesh that already passes loose thresholds."""
    nodes = np.array(
        [[0.0, 0.0], [1.0, 0.0], [1.0, 1.0], [0.0, 1.0], [0.5, 0.5]],
        dtype=float,
    )
    elements = np.array(
        [[0, 1, 4], [1, 2, 4], [2, 3, 4], [3, 0, 4]], dtype=np.int64,
    )
    return Fort14Mesh(
        title="t", nodes=nodes, depths=np.zeros(5, dtype=float),
        elements=elements,
        open_boundaries=[np.array([0, 1])],
        land_boundaries=[(0, np.array([1, 2, 3, 0]))],
    )


def _write(path) -> None:
    write_fort14(_square_mesh(), path)


def test_pipeline_no_thresholds_runs_all_rungs_and_exits_zero(tmp_path) -> None:
    src = tmp_path / "in.14"
    out = tmp_path / "out.14"
    summary = tmp_path / "p.json"
    _write(src)

    rc = meshpipeline.main(
        [str(src), str(out), "--summary", str(summary), "--quiet"],
    )
    assert rc == 0
    payload = json.loads(summary.read_text())
    assert payload["max_iters"] == 3
    # No thresholds => final.passed is None (gate not evaluated).
    assert payload["final"]["passed"] is None
    # All 3 rungs were attempted.
    assert len(payload["history"]) == 3
    assert [h["rung_label"] for h in payload["history"]] == [
        "rung0:A+B+C", "rung1:+D+F+G", "rung2:+E",
    ]


def test_pipeline_stops_at_first_passing_rung(tmp_path) -> None:
    """Loose threshold passes at rung 0; subsequent rungs must be
    skipped."""
    src = tmp_path / "in.14"
    out = tmp_path / "out.14"
    summary = tmp_path / "p.json"
    _write(src)

    rc = meshpipeline.main(
        [str(src), str(out),
         "--min-alpha", "0.5",          # easy to satisfy
         "--max-flipped", "0",          # trivially true
         "--summary", str(summary), "--quiet"],
    )
    assert rc == 0
    payload = json.loads(summary.read_text())
    # Should stop at rung 0.
    assert len(payload["history"]) == 1
    assert payload["history"][0]["rung_label"] == "rung0:A+B+C"
    assert payload["history"][0]["passed"] is True
    assert payload["final"]["rung_label"] == "rung0:A+B+C"
    assert payload["final"]["passed"] is True


def test_pipeline_threshold_failure_exits_one(tmp_path) -> None:
    """Demand alpha >= 0.95; the unit-square mesh has alpha ≈ 0.866,
    so no rung will satisfy it. Pipeline must exhaust all 3 rungs and
    exit 1."""
    src = tmp_path / "in.14"
    out = tmp_path / "out.14"
    summary = tmp_path / "p.json"
    _write(src)

    rc = meshpipeline.main(
        [str(src), str(out), "--min-alpha", "0.95",
         "--summary", str(summary), "--quiet"],
    )
    assert rc == 1
    payload = json.loads(summary.read_text())
    # All 3 rungs attempted (none passed).
    assert len(payload["history"]) == 3
    assert all(not h["passed"] for h in payload["history"])
    assert payload["final"]["passed"] is False


def test_pipeline_max_iters_caps_attempts(tmp_path) -> None:
    src = tmp_path / "in.14"
    out = tmp_path / "out.14"
    summary = tmp_path / "p.json"
    _write(src)

    rc = meshpipeline.main(
        [str(src), str(out), "--max-iters", "1",
         "--min-alpha", "0.95",
         "--summary", str(summary), "--quiet"],
    )
    # Failed at rung 0, --max-iters=1 means no further rungs.
    assert rc == 1
    payload = json.loads(summary.read_text())
    assert len(payload["history"]) == 1


def test_pipeline_writes_output_fort14(tmp_path) -> None:
    src = tmp_path / "in.14"
    out = tmp_path / "out.14"
    _write(src)

    rc = meshpipeline.main(
        [str(src), str(out), "--quiet"],
    )
    assert rc == 0
    assert out.exists()
    # The output is a fort.14; first non-empty line should match the
    # convention 'NE NP'.
    body = out.read_text().splitlines()
    assert len(body) > 1


def test_pipeline_invalid_max_iters_rejected(tmp_path) -> None:
    src = tmp_path / "in.14"
    _write(src)
    rc = meshpipeline.main(
        [str(src), str(tmp_path / "out.14"), "--max-iters", "5", "--quiet"],
    )
    assert rc == 2


def test_pipeline_missing_input_returns_2(tmp_path) -> None:
    rc = meshpipeline.main(
        [str(tmp_path / "no.14"), str(tmp_path / "out.14"), "--quiet"],
    )
    assert rc == 2


def test_pipeline_history_records_phases_per_rung(tmp_path) -> None:
    """The rung 0 history entry should list only A+B+C phases; rung 1
    adds D, F, G; rung 2 also adds E. Verifies that the rung overlay
    correctly toggles the underlying clean_mesh phase set."""
    src = tmp_path / "in.14"
    out = tmp_path / "out.14"
    summary = tmp_path / "p.json"
    _write(src)

    # Force exhaustion (impossible threshold) so all 3 rungs run.
    rc = meshpipeline.main(
        [str(src), str(out), "--min-alpha", "0.99",
         "--summary", str(summary), "--quiet"],
    )
    assert rc == 1
    payload = json.loads(summary.read_text())

    rung0_phases = set(payload["history"][0]["phases_run"])
    rung1_phases = set(payload["history"][1]["phases_run"])
    rung2_phases = set(payload["history"][2]["phases_run"])

    # Phase D / E / F / G all turn into named phases in clean_mesh's
    # info dict only when they actually run.
    assert "repair_overconnected_nodes" not in rung0_phases
    assert "repair_skewed_elements" not in rung0_phases
    assert "smooth_mesh_laplacian" not in rung0_phases
    assert "repair_under_resolved_channels" not in rung0_phases

    assert "repair_overconnected_nodes" in rung1_phases
    assert "repair_skewed_elements" in rung1_phases
    assert "smooth_mesh_laplacian" in rung1_phases
    assert "repair_under_resolved_channels" not in rung1_phases

    assert "repair_under_resolved_channels" in rung2_phases
