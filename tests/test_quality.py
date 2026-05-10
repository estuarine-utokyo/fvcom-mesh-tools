"""Tests for ``fvcom_mesh_tools.quality`` and the ``fmesh-mesh-quality`` CLI."""

from __future__ import annotations

import json

import numpy as np
import pytest

from fvcom_mesh_tools.cli import meshquality
from fvcom_mesh_tools.io import Fort14Mesh, write_fort14
from fvcom_mesh_tools.quality import (
    METRIC_KEYS,
    check_thresholds,
    compute_metrics,
    format_comparison_table,
    format_threshold_table,
)


def _square_mesh() -> Fort14Mesh:
    """Unit square split into 4 right triangles around a centre node."""
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


def _empty_mesh() -> Fort14Mesh:
    return Fort14Mesh(
        title="empty",
        nodes=np.empty((0, 2), dtype=float),
        depths=np.empty((0,), dtype=float),
        elements=np.empty((0, 3), dtype=np.int64),
        open_boundaries=[], land_boundaries=[],
    )


# ---------------------------------------------------------------------------
# compute_metrics
# ---------------------------------------------------------------------------


def test_compute_metrics_returns_all_documented_keys() -> None:
    metrics = compute_metrics(_square_mesh())
    assert set(metrics.keys()) == set(METRIC_KEYS)


def test_compute_metrics_square_mesh_values() -> None:
    metrics = compute_metrics(_square_mesh())
    assert metrics["n_nodes"] == 5
    assert metrics["n_elements"] == 4
    assert metrics["n_components"] == 1
    assert metrics["n_disjoint_elems"] == 0
    assert metrics["n_flipped"] == 0
    # The 4 right triangles have 90° / 45° / 45° interior angles.
    np.testing.assert_allclose(metrics["min_angle_p50_deg"], 45.0, atol=1e-6)
    # alpha for a 1-1-sqrt(2) right triangle is ~0.866.
    np.testing.assert_allclose(metrics["alpha_mean"], 0.8660254, atol=1e-6)
    # Centre node has valence 4; corners have valence 2.
    assert metrics["max_valence"] == 4


def test_compute_metrics_empty_mesh() -> None:
    metrics = compute_metrics(_empty_mesh())
    assert metrics["n_elements"] == 0
    assert metrics["n_nodes"] == 0
    # Shape metrics are nan, not 0 — they're undefined on an empty mesh.
    assert np.isnan(metrics["alpha_mean"])
    assert np.isnan(metrics["frac_lt_20deg"])


# ---------------------------------------------------------------------------
# check_thresholds
# ---------------------------------------------------------------------------


def test_check_thresholds_pass() -> None:
    metrics = compute_metrics(_square_mesh())
    passed, checks = check_thresholds(
        metrics, min_alpha_mean=0.5, max_valence=4,
    )
    assert passed
    assert len(checks) == 2
    assert all(c.passed for c in checks)


def test_check_thresholds_fail_on_alpha() -> None:
    metrics = compute_metrics(_square_mesh())
    passed, checks = check_thresholds(metrics, min_alpha_mean=0.95)
    assert not passed
    assert checks[0].metric == "alpha_mean"
    assert checks[0].op == "≥"
    assert not checks[0].passed


def test_check_thresholds_no_thresholds_means_pass() -> None:
    metrics = compute_metrics(_square_mesh())
    passed, checks = check_thresholds(metrics)
    assert passed
    assert checks == []


def test_check_thresholds_nan_actual_fails() -> None:
    metrics = compute_metrics(_empty_mesh())
    passed, _ = check_thresholds(metrics, min_alpha_mean=0.5)
    assert not passed


# ---------------------------------------------------------------------------
# Pretty printing
# ---------------------------------------------------------------------------


def test_format_comparison_table_single_mesh() -> None:
    metrics = compute_metrics(_square_mesh())
    text = format_comparison_table([("only", metrics)])
    assert "metric" in text
    assert "only" in text
    assert "delta" not in text  # delta only when 2 rows
    for k in METRIC_KEYS:
        assert k in text


def test_format_comparison_table_two_meshes_show_delta() -> None:
    m1 = compute_metrics(_square_mesh())
    m2 = dict(m1)
    m2["n_elements"] = m1["n_elements"] - 1
    text = format_comparison_table([("a", m1), ("b", m2)])
    assert "delta" in text
    # n_elements decreased by 1 → "-1" should appear.
    assert "-1" in text


def test_format_threshold_table_includes_pass_fail() -> None:
    metrics = compute_metrics(_square_mesh())
    _, checks = check_thresholds(metrics, min_alpha_mean=0.5, max_valence=2)
    text = format_threshold_table(checks)
    assert "PASS" in text
    assert "FAIL" in text


# ---------------------------------------------------------------------------
# CLI integration
# ---------------------------------------------------------------------------


def _write_fixture(path) -> None:
    write_fort14(_square_mesh(), path)


def test_cli_single_mesh_writes_summary_and_exits_zero(tmp_path) -> None:
    fort14 = tmp_path / "mesh.14"
    _write_fixture(fort14)
    summary = tmp_path / "out.json"

    rc = meshquality.main(
        [str(fort14), "--summary", str(summary), "--quiet"],
    )
    assert rc == 0
    assert summary.exists()
    payload = json.loads(summary.read_text())
    assert payload["passed"] is True
    assert len(payload["meshes"]) == 1
    assert set(payload["meshes"][0]["metrics"].keys()) == set(METRIC_KEYS)


def test_cli_threshold_gate_fails_with_exit_1(tmp_path) -> None:
    fort14 = tmp_path / "mesh.14"
    _write_fixture(fort14)
    summary = tmp_path / "out.json"

    # The square mesh has alpha_mean ≈ 0.866, so a 0.95 floor fails.
    rc = meshquality.main(
        [str(fort14), "--min-alpha", "0.95",
         "--summary", str(summary), "--quiet"],
    )
    assert rc == 1
    payload = json.loads(summary.read_text())
    assert payload["passed"] is False
    assert any(not c["passed"] for c in payload["checks"])


def test_cli_two_meshes_compare(tmp_path) -> None:
    """Two-input comparison mode — JSON has both meshes, exit 0 (no
    thresholds supplied)."""
    a = tmp_path / "before.14"
    b = tmp_path / "after.14"
    _write_fixture(a)
    _write_fixture(b)
    summary = tmp_path / "out.json"

    rc = meshquality.main(
        [str(a), str(b), "--labels", "before", "after",
         "--summary", str(summary), "--quiet"],
    )
    assert rc == 0
    payload = json.loads(summary.read_text())
    assert [m["label"] for m in payload["meshes"]] == ["before", "after"]


def test_cli_label_count_mismatch_rejected(tmp_path, capsys) -> None:
    fort14 = tmp_path / "mesh.14"
    _write_fixture(fort14)

    rc = meshquality.main(
        [str(fort14), "--labels", "a", "b", "--quiet"],
    )
    assert rc == 2
    err = capsys.readouterr().err
    assert "labels" in err.lower()


def test_cli_missing_input_returns_2(tmp_path) -> None:
    rc = meshquality.main([str(tmp_path / "nonexistent.14"), "--quiet"])
    assert rc == 2


def test_cli_invalid_max_nbr_elem_rejected(tmp_path) -> None:
    fort14 = tmp_path / "mesh.14"
    _write_fixture(fort14)
    rc = meshquality.main(
        [str(fort14), "--max-nbr-elem", "2", "--quiet"],
    )
    assert rc == 2


def test_cli_threshold_evaluated_against_last_input(tmp_path) -> None:
    """Build two meshes where the FIRST passes a threshold and the
    SECOND fails. The CLI should use the second (last) and exit 1."""
    a = tmp_path / "good.14"
    b = tmp_path / "bad.14"
    # Good: square mesh.
    write_fort14(_square_mesh(), a)
    # Bad: same but with one element flipped to violate n_flipped<=0.
    bad = _square_mesh()
    bad.elements[0] = bad.elements[0][::-1]  # invert orientation
    write_fort14(bad, b)
    summary = tmp_path / "q.json"
    rc = meshquality.main(
        [str(a), str(b), "--max-flipped", "0",
         "--summary", str(summary), "--quiet"],
    )
    # With the second mesh having a flipped triangle, exit 1.
    assert rc == 1


@pytest.mark.parametrize(
    "key", ["alpha_mean", "frac_lt_20deg", "max_valence", "n_flipped"],
)
def test_metric_key_present_in_table(key) -> None:
    metrics = compute_metrics(_square_mesh())
    text = format_comparison_table([("a", metrics)])
    assert key in text
