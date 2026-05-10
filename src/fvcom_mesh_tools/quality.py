"""Unified mesh-quality metrics for FVCOM fort.14 meshes.

This module is the single source of truth for the quality numbers
that ``fmesh-mesh-quality`` reports and that the existing PoCs
compute ad-hoc. It deliberately re-uses
:mod:`fvcom_mesh_tools.algorithms.quality` (alpha, min interior
angle, signed areas) and :mod:`fvcom_mesh_tools.diagnostics`
(face-face adjacency, valence, over-connected detector) so the
numbers shown by ``fmesh-mesh-quality`` agree with those reported
by ``fmesh-mesh-check`` and ``fmesh-mesh-clean``.

Public entry points:

* :func:`compute_metrics` returns a flat ``dict`` of float / int
  metrics suitable for JSON serialisation and for direct comparison
  against a threshold dict.
* :func:`check_thresholds` evaluates a metrics dict against
  user-supplied thresholds (e.g. ``min_alpha_mean=0.95``) and
  returns ``(passed, checks)`` where ``checks`` is a list of
  per-threshold :class:`ThresholdCheck` records.
* :func:`format_comparison_table` renders a side-by-side text table
  for one or more meshes, with an automatic "delta" column when two
  meshes are passed.
"""

from __future__ import annotations

from dataclasses import dataclass
from typing import Any

import numpy as np
from scipy.sparse.csgraph import connected_components

from fvcom_mesh_tools.algorithms import (
    alpha_quality,
    min_interior_angle,
    signed_areas,
)
from fvcom_mesh_tools.diagnostics import (
    DEFAULT_MAX_NBR_ELEM,
    face_face_adjacency,
    node_valence,
    overconnected_nodes_flag,
)
from fvcom_mesh_tools.io import Fort14Mesh

# Stable order for printing / JSON. Keeping this in one place lets
# fmesh-mesh-quality, the JSON dumper, and any test agree.
METRIC_KEYS: tuple[str, ...] = (
    "n_nodes",
    "n_elements",
    "n_open_boundaries",
    "n_land_boundaries",
    "alpha_mean",
    "alpha_p05",
    "alpha_p50",
    "min_angle_p05_deg",
    "min_angle_p50_deg",
    "frac_lt_20deg",
    "max_valence",
    "n_overconnected",
    "n_flipped",
    "n_components",
    "n_disjoint_elems",
)

# Per-metric formatting hints for the comparison table.
_INT_KEYS = {
    "n_nodes", "n_elements", "n_open_boundaries", "n_land_boundaries",
    "max_valence", "n_overconnected", "n_flipped",
    "n_components", "n_disjoint_elems",
}
_PCT_KEYS = {"frac_lt_20deg"}


def compute_metrics(
    mesh: Fort14Mesh,
    *,
    max_nbr_elem: int = DEFAULT_MAX_NBR_ELEM,
) -> dict[str, Any]:
    """Return the standard quality-metrics dict for ``mesh``.

    Empty meshes return zero counts and ``nan`` for shape-dependent
    metrics — by convention, shape metrics over no triangles are not
    well-defined and we'd rather show ``nan`` than mislead with 0.0.
    """
    NP = int(mesh.n_nodes)
    NE = int(mesh.n_elements)
    n_open = int(len(mesh.open_boundaries))
    n_land = int(len(mesh.land_boundaries))

    if NE == 0:
        return {
            "n_nodes": NP,
            "n_elements": 0,
            "n_open_boundaries": n_open,
            "n_land_boundaries": n_land,
            "alpha_mean": float("nan"),
            "alpha_p05": float("nan"),
            "alpha_p50": float("nan"),
            "min_angle_p05_deg": float("nan"),
            "min_angle_p50_deg": float("nan"),
            "frac_lt_20deg": float("nan"),
            "max_valence": 0,
            "n_overconnected": 0,
            "n_flipped": 0,
            "n_components": 0,
            "n_disjoint_elems": 0,
        }

    alpha = alpha_quality(mesh)
    angles = min_interior_angle(mesh)  # degrees
    sa = signed_areas(mesh)
    valence = node_valence(mesh.elements, n_nodes=NP)
    over_flag, _ = overconnected_nodes_flag(
        mesh.elements, n_nodes=NP, max_nbr=int(max_nbr_elem),
    )

    adj = face_face_adjacency(mesh.elements)
    n_comp, labels = connected_components(adj, directed=False, return_labels=True)
    comp_sizes = np.bincount(labels, minlength=int(n_comp))
    largest = int(comp_sizes.max()) if comp_sizes.size else 0
    n_disjoint = NE - largest

    return {
        "n_nodes": NP,
        "n_elements": NE,
        "n_open_boundaries": n_open,
        "n_land_boundaries": n_land,
        "alpha_mean": float(alpha.mean()),
        "alpha_p05": float(np.percentile(alpha, 5)),
        "alpha_p50": float(np.median(alpha)),
        "min_angle_p05_deg": float(np.percentile(angles, 5)),
        "min_angle_p50_deg": float(np.median(angles)),
        "frac_lt_20deg": float((angles < 20.0).sum() / NE),
        "max_valence": int(valence.max()) if valence.size else 0,
        "n_overconnected": int(over_flag.sum()),
        "n_flipped": int((sa < 0).sum()),
        "n_components": int(n_comp),
        "n_disjoint_elems": int(n_disjoint),
    }


# ---------------------------------------------------------------------------
# Threshold checks
# ---------------------------------------------------------------------------


@dataclass
class ThresholdCheck:
    """One pass/fail record for a threshold."""

    metric: str
    op: str            # "≥" or "≤"
    threshold: float
    actual: float
    passed: bool

    def to_dict(self) -> dict[str, Any]:
        return {
            "metric": self.metric,
            "op": self.op,
            "threshold": self.threshold,
            "actual": self.actual,
            "passed": bool(self.passed),
        }


def check_thresholds(
    metrics: dict[str, Any],
    *,
    min_alpha_mean: float | None = None,
    max_frac_lt_20deg: float | None = None,
    max_valence: int | None = None,
    max_overconnected: int | None = None,
    max_flipped: int | None = None,
    max_disjoint_elems: int | None = None,
) -> tuple[bool, list[ThresholdCheck]]:
    """Evaluate ``metrics`` against the supplied thresholds.

    Each unspecified threshold is skipped. Returns ``(passed, checks)``
    where ``passed`` is True iff every supplied check passed.
    ``nan`` actuals (e.g. shape metrics on an empty mesh) fail any
    threshold they're evaluated against.
    """
    checks: list[ThresholdCheck] = []

    def _ge(metric: str, threshold: float) -> None:
        actual = metrics.get(metric)
        if actual is None:
            return
        passed = (
            isinstance(actual, (int, float))
            and not (isinstance(actual, float) and np.isnan(actual))
            and float(actual) >= float(threshold)
        )
        checks.append(ThresholdCheck(metric, "≥", float(threshold),
                                     float(actual), passed))

    def _le(metric: str, threshold: float) -> None:
        actual = metrics.get(metric)
        if actual is None:
            return
        passed = (
            isinstance(actual, (int, float))
            and not (isinstance(actual, float) and np.isnan(actual))
            and float(actual) <= float(threshold)
        )
        checks.append(ThresholdCheck(metric, "≤", float(threshold),
                                     float(actual), passed))

    if min_alpha_mean is not None:
        _ge("alpha_mean", min_alpha_mean)
    if max_frac_lt_20deg is not None:
        _le("frac_lt_20deg", max_frac_lt_20deg)
    if max_valence is not None:
        _le("max_valence", int(max_valence))
    if max_overconnected is not None:
        _le("n_overconnected", int(max_overconnected))
    if max_flipped is not None:
        _le("n_flipped", int(max_flipped))
    if max_disjoint_elems is not None:
        _le("n_disjoint_elems", int(max_disjoint_elems))

    overall = all(c.passed for c in checks) if checks else True
    return overall, checks


# ---------------------------------------------------------------------------
# Pretty printing
# ---------------------------------------------------------------------------


def _fmt_value(metric: str, value: Any) -> str:
    if value is None:
        return "n/a"
    if isinstance(value, float) and np.isnan(value):
        return "n/a"
    if metric in _INT_KEYS:
        return f"{int(value):,}"
    if metric in _PCT_KEYS:
        return f"{float(value) * 100:.4f}%"
    if isinstance(value, float):
        return f"{value:.4f}"
    return str(value)


def _fmt_delta(metric: str, before: Any, after: Any) -> str:
    """Pretty signed delta string for a metric. Returns '' if either
    side is non-numeric or NaN."""
    if before is None or after is None:
        return ""
    try:
        b = float(before)
        a = float(after)
    except (TypeError, ValueError):
        return ""
    if np.isnan(b) or np.isnan(a):
        return ""
    diff = a - b
    if metric in _INT_KEYS:
        if diff == 0:
            return "0"
        return f"{int(diff):+,}"
    if metric in _PCT_KEYS:
        return f"{diff * 100:+.4f} pp"
    return f"{diff:+.4f}"


def format_comparison_table(
    rows: list[tuple[str, dict[str, Any]]],
    *,
    keys: tuple[str, ...] = METRIC_KEYS,
) -> str:
    """Render a side-by-side comparison table for one or more meshes.

    ``rows`` is ``[(label, metrics_dict), ...]``. With exactly two
    rows the table also gains a "delta" column (``after - before``).
    """
    if not rows:
        return ""
    metric_w = max(len(k) for k in keys)
    metric_w = max(metric_w, len("metric"))
    label_w = max(
        max((len(_fmt_value(k, m.get(k))) for m in (md for _, md in rows)
             for k in keys), default=0),
        max(len(label) for label, _ in rows),
        len("delta"),
        12,
    )

    show_delta = len(rows) == 2
    headers = ["metric"] + [label for label, _ in rows]
    if show_delta:
        headers.append("delta")
    widths = [metric_w] + [label_w] * (len(headers) - 1)
    header_line = "  ".join(
        h.ljust(widths[i]) if i == 0 else h.rjust(widths[i])
        for i, h in enumerate(headers)
    )
    sep = "  ".join("-" * w for w in widths)
    out_lines: list[str] = [header_line, sep]

    for k in keys:
        cells = [k.ljust(widths[0])]
        values = [md.get(k) for _, md in rows]
        for i, v in enumerate(values, start=1):
            cells.append(_fmt_value(k, v).rjust(widths[i]))
        if show_delta:
            cells.append(_fmt_delta(k, values[0], values[1]).rjust(widths[-1]))
        out_lines.append("  ".join(cells))

    return "\n".join(out_lines)


def format_threshold_table(checks: list[ThresholdCheck]) -> str:
    if not checks:
        return "(no thresholds supplied)"
    metric_w = max(max(len(c.metric) for c in checks), len("metric"))
    threshold_w = max(
        max(len(_fmt_value(c.metric, c.threshold)) for c in checks), 9,
    )
    actual_w = max(
        max(len(_fmt_value(c.metric, c.actual)) for c in checks), 6,
    )
    lines = [
        "  ".join([
            "metric".ljust(metric_w), "op".center(2),
            "threshold".rjust(threshold_w),
            "actual".rjust(actual_w), "result",
        ]),
        "  ".join([
            "-" * metric_w, "--", "-" * threshold_w,
            "-" * actual_w, "------",
        ]),
    ]
    for c in checks:
        lines.append("  ".join([
            c.metric.ljust(metric_w),
            c.op.center(2),
            _fmt_value(c.metric, c.threshold).rjust(threshold_w),
            _fmt_value(c.metric, c.actual).rjust(actual_w),
            "PASS" if c.passed else "FAIL",
        ]))
    return "\n".join(lines)


__all__ = [
    "METRIC_KEYS",
    "ThresholdCheck",
    "check_thresholds",
    "compute_metrics",
    "format_comparison_table",
    "format_threshold_table",
]
