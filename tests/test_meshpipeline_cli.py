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


# ---------------------------------------------------------------------------
# --best-rung selection
# ---------------------------------------------------------------------------


def _fake_rung_results(entries: list[dict]) -> list[tuple[dict, object]]:
    """Build the (history-dict, mesh-stand-in) tuples ``_select_rung``
    expects from a list of compact entry dicts. The mesh stand-in is
    the entry's ``rung_label`` so the caller can identify which mesh
    would have been picked."""
    out = []
    for i, e in enumerate(entries):
        out.append((
            {
                "rung_index": i,
                "rung_label": e["label"],
                "metrics": {"alpha_mean": e["alpha"]},
                "passed": bool(e.get("passed", False)),
            },
            e["label"],   # stand in for the mesh
        ))
    return out


def test_select_rung_default_returns_last_attempted() -> None:
    """``best_rung=False`` must reproduce the legacy "last attempted"
    behaviour regardless of alpha."""
    from fvcom_mesh_tools.cli.meshpipeline import _select_rung

    rung_results = _fake_rung_results([
        {"label": "rung0", "alpha": 0.97, "passed": True},
        {"label": "rung1", "alpha": 0.95, "passed": True},
    ])
    idx, reason = _select_rung(rung_results, thresholds=True, best_rung=False)
    assert idx == 1
    assert "first-passing" in reason


def test_select_rung_best_picks_highest_alpha_among_passing() -> None:
    from fvcom_mesh_tools.cli.meshpipeline import _select_rung

    rung_results = _fake_rung_results([
        {"label": "rung0", "alpha": 0.95, "passed": True},
        {"label": "rung1", "alpha": 0.97, "passed": True},
        {"label": "rung2", "alpha": 0.96, "passed": False},
    ])
    idx, reason = _select_rung(rung_results, thresholds=True, best_rung=True)
    assert idx == 1
    assert "best-alpha-mean" in reason
    assert "gate-passing" in reason


def test_select_rung_best_ignores_failing_rungs_when_thresholds_set() -> None:
    """rung2 has the highest alpha but fails the gate; rung0 is the
    only passing rung and should win even though rung1 (failing) has
    higher alpha."""
    from fvcom_mesh_tools.cli.meshpipeline import _select_rung

    rung_results = _fake_rung_results([
        {"label": "rung0", "alpha": 0.91, "passed": True},
        {"label": "rung1", "alpha": 0.95, "passed": False},
        {"label": "rung2", "alpha": 0.98, "passed": False},
    ])
    idx, _ = _select_rung(rung_results, thresholds=True, best_rung=True)
    assert idx == 0


def test_select_rung_best_ties_prefer_lower_rung() -> None:
    """Two rungs with identical alpha; the lighter repair (lower
    index) must win."""
    from fvcom_mesh_tools.cli.meshpipeline import _select_rung

    rung_results = _fake_rung_results([
        {"label": "rung0", "alpha": 0.96, "passed": True},
        {"label": "rung1", "alpha": 0.96, "passed": True},
    ])
    idx, _ = _select_rung(rung_results, thresholds=True, best_rung=True)
    assert idx == 0


def test_select_rung_best_no_threshold_picks_highest_alpha_overall() -> None:
    from fvcom_mesh_tools.cli.meshpipeline import _select_rung

    rung_results = _fake_rung_results([
        {"label": "rung0", "alpha": 0.92, "passed": False},
        {"label": "rung1", "alpha": 0.95, "passed": False},
        {"label": "rung2", "alpha": 0.93, "passed": False},
    ])
    idx, reason = _select_rung(rung_results, thresholds=False, best_rung=True)
    assert idx == 1
    assert "best-alpha-mean" in reason


def test_select_rung_best_falls_back_when_no_rung_passes() -> None:
    """All rungs fail the gate; pick the highest-alpha one and report
    the fallback reason."""
    from fvcom_mesh_tools.cli.meshpipeline import _select_rung

    rung_results = _fake_rung_results([
        {"label": "rung0", "alpha": 0.91, "passed": False},
        {"label": "rung1", "alpha": 0.92, "passed": False},
        {"label": "rung2", "alpha": 0.93, "passed": False},
    ])
    idx, reason = _select_rung(rung_results, thresholds=True, best_rung=True)
    assert idx == 2
    assert "no rung passed" in reason


def test_select_rung_best_handles_nan_alpha() -> None:
    """Empty meshes give NaN alpha; they must never beat a finite
    alpha entry."""
    from fvcom_mesh_tools.cli.meshpipeline import _select_rung

    rung_results = _fake_rung_results([
        {"label": "rung0", "alpha": float("nan"), "passed": True},
        {"label": "rung1", "alpha": 0.50, "passed": True},
    ])
    idx, _ = _select_rung(rung_results, thresholds=True, best_rung=True)
    assert idx == 1


def test_pipeline_best_rung_runs_every_rung_and_records_selection(
    tmp_path,
) -> None:
    """Loose threshold passes from rung 0; ``--best-rung`` must still
    run rungs 1 and 2 and record the selection_reason in the JSON."""
    src = tmp_path / "in.14"
    out = tmp_path / "out.14"
    summary = tmp_path / "p.json"
    _write(src)

    rc = meshpipeline.main([
        str(src), str(out),
        "--min-alpha", "0.5",
        "--max-flipped", "0",
        "--best-rung",
        "--summary", str(summary), "--quiet",
    ])
    assert rc == 0
    payload = json.loads(summary.read_text())

    # All three rungs were attempted (no early-stop with --best-rung).
    assert len(payload["history"]) == 3
    assert payload["best_rung_mode"] is True
    assert "selection_reason" in payload["final"]
    # Final payload includes the rung_index of the chosen rung.
    assert "rung_index" in payload["final"]
    assert payload["final"]["selection_reason"].startswith("best-alpha-mean")


def test_pipeline_best_rung_default_off_keeps_first_pass_stop(tmp_path) -> None:
    """Without ``--best-rung``, the first-pass-stop behaviour holds."""
    src = tmp_path / "in.14"
    out = tmp_path / "out.14"
    summary = tmp_path / "p.json"
    _write(src)

    rc = meshpipeline.main([
        str(src), str(out),
        "--min-alpha", "0.5",
        "--max-flipped", "0",
        "--summary", str(summary), "--quiet",
    ])
    assert rc == 0
    payload = json.loads(summary.read_text())
    assert payload["best_rung_mode"] is False
    assert len(payload["history"]) == 1   # stopped at rung 0
