"""PoC #48: FVCOM manual mesh requirements check.

User-provided FVCOM mesh requirements (from the FVCOM manual,
SMS-checkable):

  1. Minimum interior angle: >= 30.0°
  2. Maximum interior angle: <= 130.0°
  3. Maximum slope:          <= 0.1   (depth gradient between adjacent nodes)
  4. Element area change:    <= 0.5   (ratio between adjacent element areas)
  5. Connecting elements:    <= 8     (max node valence)

This script evaluates a fort.14 mesh against all five criteria,
reporting per-element / per-edge / per-node fail counts. Read-only
— no mesh is written.

Two interpretations are reported for criteria 3 and 4 (SMS does not
explicitly cite formulas in its UI; the convention below is the
most common one in coastal-mesh literature; the user can request a
different convention if SMS uses something else):

  slope (per edge)         = |h_a - h_b| / max(0.5*(|h_a|+|h_b|), eps)
  area_change (per edge)   = |A_i - A_j| / max(A_i, A_j)

For each fort.14 input we report:

  * total element / node count
  * per-criterion fail count + share
  * percentile distributions (p05 / p50 / p95 / p99 / max for
    elements; corresponding values for edge / node metrics)

Multiple inputs can be given on the command line. Default inputs
are the Phase G output and the Phase H v3 output for direct A/B.

Outputs:
    outputs/48_fvcom_requirements_<stem>.txt
    outputs/48_fvcom_requirements_<stem>.json
"""
from __future__ import annotations

import argparse
import json
import sys
from collections import defaultdict
from pathlib import Path

import numpy as np

from fvcom_mesh_tools.algorithms.quality import (
    edge_lengths_planar,
    min_interior_angle,
)
from fvcom_mesh_tools.diagnostics import (
    node_valence,
)
from fvcom_mesh_tools.io import Fort14Mesh, read_fort14

REPO = Path(__file__).resolve().parent.parent
OUT_DIR = REPO / "outputs"

# FVCOM manual thresholds (per user spec):
MIN_ANGLE_DEG = 30.0
MAX_ANGLE_DEG = 130.0
MAX_SLOPE = 0.1
MAX_AREA_CHANGE = 0.5
MAX_VALENCE = 8


def max_interior_angle(mesh: Fort14Mesh) -> np.ndarray:
    """Per-triangle largest interior angle, in degrees.

    Mirrors ``min_interior_angle`` but returns the max instead of
    the min over the three vertex angles.
    """
    ll = edge_lengths_planar(mesh)
    a = ll[:, 1]
    b = ll[:, 2]
    c = ll[:, 0]

    def _angle(opp: np.ndarray, e1: np.ndarray, e2: np.ndarray) -> np.ndarray:
        cos = (e1 ** 2 + e2 ** 2 - opp ** 2) / np.where(
            e1 * e2 == 0, 1.0, 2.0 * e1 * e2,
        )
        return np.arccos(np.clip(cos, -1.0, 1.0))

    angle0 = _angle(a, c, b)
    angle1 = _angle(b, a, c)
    angle2 = _angle(c, a, b)
    return np.maximum(np.maximum(angle0, angle1), angle2) * 180.0 / np.pi


def element_areas(mesh: Fort14Mesh) -> np.ndarray:
    """Per-element signed area (positive for CCW)."""
    p0 = mesh.nodes[mesh.elements[:, 0]]
    p1 = mesh.nodes[mesh.elements[:, 1]]
    p2 = mesh.nodes[mesh.elements[:, 2]]
    return 0.5 * (
        (p1[:, 0] - p0[:, 0]) * (p2[:, 1] - p0[:, 1])
        - (p1[:, 1] - p0[:, 1]) * (p2[:, 0] - p0[:, 0])
    )


def _edges_with_buddy(mesh: Fort14Mesh) -> tuple[np.ndarray, np.ndarray, np.ndarray]:
    """Return ``(edge_uv, elem_pair)`` — for every internal edge
    shared by exactly two elements, give (u, v) and (e_i, e_j).
    Boundary edges are skipped.
    Also returns the boundary edge count for sanity.
    """
    NE = mesh.elements.shape[0]
    edges = np.vstack([
        mesh.elements[:, [0, 1]],
        mesh.elements[:, [1, 2]],
        mesh.elements[:, [2, 0]],
    ])
    edges_sorted = np.sort(edges, axis=1)
    elem_of = np.tile(np.arange(NE), 3)
    keys = edges_sorted[:, 0].astype(np.int64) * (mesh.n_nodes + 1) \
        + edges_sorted[:, 1].astype(np.int64)
    order = np.argsort(keys, kind="stable")
    keys_sorted = keys[order]
    edges_s = edges_sorted[order]
    elem_s = elem_of[order]

    # Same key in consecutive entries ⇒ shared edge.
    same = (np.diff(keys_sorted) == 0)
    pair_idx = np.where(same)[0]
    edge_uv = edges_s[pair_idx]
    elem_pair = np.column_stack([elem_s[pair_idx], elem_s[pair_idx + 1]])

    # Boundary edges: unique keys.
    n_unique = int(np.unique(keys_sorted).size)
    n_boundary = n_unique - int(pair_idx.size)
    return edge_uv, elem_pair, np.array([n_boundary], dtype=np.int64)


def _summarize_dist(arr: np.ndarray, label: str) -> dict[str, float]:
    if arr.size == 0:
        return {f"{label}_count": 0}
    return {
        f"{label}_count": int(arr.size),
        f"{label}_min": float(arr.min()),
        f"{label}_p05": float(np.percentile(arr, 5)),
        f"{label}_p50": float(np.percentile(arr, 50)),
        f"{label}_p95": float(np.percentile(arr, 95)),
        f"{label}_p99": float(np.percentile(arr, 99)),
        f"{label}_max": float(arr.max()),
        f"{label}_mean": float(arr.mean()),
    }


def evaluate(mesh: Fort14Mesh) -> dict:
    """Return a results dict covering all five FVCOM criteria."""
    result: dict = {
        "n_nodes": int(mesh.n_nodes),
        "n_elements": int(mesh.n_elements),
    }

    # Criterion 1 + 2: interior angles.
    min_ang = min_interior_angle(mesh)
    max_ang = max_interior_angle(mesh)

    fail_min_angle = min_ang < MIN_ANGLE_DEG
    fail_max_angle = max_ang > MAX_ANGLE_DEG

    result["criterion_1_min_angle"] = {
        "threshold_deg": MIN_ANGLE_DEG,
        "n_fail": int(fail_min_angle.sum()),
        "share_fail": (
            float(fail_min_angle.sum()) / max(mesh.n_elements, 1)
        ),
        **_summarize_dist(min_ang, "min_angle_deg"),
    }
    result["criterion_2_max_angle"] = {
        "threshold_deg": MAX_ANGLE_DEG,
        "n_fail": int(fail_max_angle.sum()),
        "share_fail": (
            float(fail_max_angle.sum()) / max(mesh.n_elements, 1)
        ),
        **_summarize_dist(max_ang, "max_angle_deg"),
    }

    # Criterion 3: max slope (depth gradient across an edge).
    edge_uv, elem_pair, _n_bnd = _edges_with_buddy(mesh)
    if mesh.depths is not None and mesh.depths.size == mesh.n_nodes:
        h = mesh.depths
        h_a = h[edge_uv[:, 0]]
        h_b = h[edge_uv[:, 1]]
        denom = np.maximum(0.5 * (np.abs(h_a) + np.abs(h_b)), 1e-9)
        slope_per_edge = np.abs(h_a - h_b) / denom
        fail_slope = slope_per_edge > MAX_SLOPE
        result["criterion_3_max_slope"] = {
            "threshold": MAX_SLOPE,
            "formula": "|h_a - h_b| / max(0.5*(|h_a|+|h_b|), eps) per edge",
            "n_edges_checked": int(slope_per_edge.size),
            "n_fail": int(fail_slope.sum()),
            "share_fail": (
                float(fail_slope.sum()) / max(slope_per_edge.size, 1)
            ),
            **_summarize_dist(slope_per_edge, "slope"),
        }
    else:
        result["criterion_3_max_slope"] = {
            "note": "depths missing or wrong size; skipped",
        }

    # Criterion 4: element area change (ratio between adjacent elements).
    areas = element_areas(mesh)
    a_i = areas[elem_pair[:, 0]]
    a_j = areas[elem_pair[:, 1]]
    larger = np.maximum(np.abs(a_i), np.abs(a_j))
    smaller = np.minimum(np.abs(a_i), np.abs(a_j))
    area_change = (larger - smaller) / np.maximum(larger, 1e-30)
    fail_area = area_change > MAX_AREA_CHANGE
    result["criterion_4_area_change"] = {
        "threshold": MAX_AREA_CHANGE,
        "formula": "(max(A_i, A_j) - min(A_i, A_j)) / max(A_i, A_j)",
        "n_edges_checked": int(area_change.size),
        "n_fail": int(fail_area.sum()),
        "share_fail": (
            float(fail_area.sum()) / max(area_change.size, 1)
        ),
        **_summarize_dist(area_change, "area_change"),
    }

    # Criterion 5: connecting elements (max node valence).
    valence = node_valence(mesh.elements, mesh.n_nodes)
    fail_valence = valence > MAX_VALENCE
    result["criterion_5_valence"] = {
        "threshold": MAX_VALENCE,
        "n_fail": int(fail_valence.sum()),
        "share_fail": (
            float(fail_valence.sum()) / max(valence.size, 1)
        ),
        "max_valence": int(valence.max()) if valence.size else 0,
        **_summarize_dist(valence.astype(float), "valence"),
    }

    # Composite: count of elements failing any of criteria 1 + 2.
    fail_any_element = fail_min_angle | fail_max_angle
    result["composite_fail_any_element_criterion"] = {
        "n_fail": int(fail_any_element.sum()),
        "share_fail": (
            float(fail_any_element.sum()) / max(mesh.n_elements, 1)
        ),
    }

    return result


def _format_block(label: str, r: dict) -> list[str]:
    lines = [f"  {label}"]
    for key, val in r.items():
        if isinstance(val, float):
            lines.append(f"    {key:<32} {val:>14.6g}")
        else:
            lines.append(f"    {key:<32} {val!s:>14}")
    return lines


def _format_report(path: Path, r: dict) -> str:
    lines: list[str] = [
        f"PoC #48 — FVCOM manual requirements check on {path.name}",
        f"  NP={r['n_nodes']:,}  NE={r['n_elements']:,}",
        "",
        "FVCOM manual thresholds:",
        f"  min interior angle >= {MIN_ANGLE_DEG}°",
        f"  max interior angle <= {MAX_ANGLE_DEG}°",
        f"  max slope          <= {MAX_SLOPE}",
        f"  element area change<= {MAX_AREA_CHANGE}",
        f"  max valence        <= {MAX_VALENCE}",
        "",
        f"  {'criterion':<26}  {'fail':>10}  {'share':>10}  pass?",
        "  " + "-" * 60,
    ]
    rows = [
        ("1. min interior angle", r["criterion_1_min_angle"]),
        ("2. max interior angle", r["criterion_2_max_angle"]),
        ("3. max slope (depth)", r["criterion_3_max_slope"]),
        ("4. element area change", r["criterion_4_area_change"]),
        ("5. max valence", r["criterion_5_valence"]),
    ]
    for label, blk in rows:
        nf = blk.get("n_fail")
        sf = blk.get("share_fail")
        if nf is None:
            lines.append(
                f"  {label:<26}  {'-':>10}  {'-':>10}  (skipped: " +
                blk.get("note", "n/a") + ")"
            )
        else:
            ok = "PASS ✓" if nf == 0 else "FAIL ✗"
            lines.append(
                f"  {label:<26}  {nf:>10,}  "
                f"{sf:>10.4%}  {ok}"
            )

    composite = r["composite_fail_any_element_criterion"]
    lines.append("")
    lines.append(
        f"  composite (elements failing C1 or C2): "
        f"{composite['n_fail']:,} "
        f"({composite['share_fail']:.4%})"
    )

    lines.append("")
    lines.append("Distribution detail:")
    for label, blk in rows:
        lines.extend(_format_block(label, blk))
    return "\n".join(lines) + "\n"


def main() -> int:
    p = argparse.ArgumentParser()
    p.add_argument("inputs", nargs="*", type=Path,
                   help="fort.14 paths to evaluate")
    p.add_argument("--out-dir", type=Path, default=OUT_DIR)
    args = p.parse_args()

    if not args.inputs:
        # Default: Phase G output (rung-1 passing) + Phase H v3 output.
        args.inputs = [
            REPO / "outputs" / "33_pipeline_passing.14",
            REPO / "outputs" / "43_phase_h_v3_optimized.14",
        ]
    for inp in args.inputs:
        if not inp.exists():
            print(f"missing: {inp}", file=sys.stderr)
            return 2

    args.out_dir.mkdir(parents=True, exist_ok=True)
    combined: dict = defaultdict(dict)
    for inp in args.inputs:
        mesh = read_fort14(inp)
        r = evaluate(mesh)
        stem = inp.stem
        combined[stem] = r
        report = _format_report(inp, r)
        txt_path = args.out_dir / f"48_fvcom_requirements_{stem}.txt"
        json_path = args.out_dir / f"48_fvcom_requirements_{stem}.json"
        txt_path.write_text(report, encoding="utf-8")
        json_path.write_text(
            json.dumps(r, indent=2, default=str), encoding="utf-8",
        )
        print(report)
        print(f"wrote {txt_path}")
        print(f"wrote {json_path}")
        print()

    if len(args.inputs) >= 2:
        # Compact A/B summary.
        lines = [
            "A/B summary (criteria 1-5; PASS/FAIL only):",
            f"  {'input':<48}  C1 C2 C3 C4 C5",
            "  " + "-" * 75,
        ]
        for inp in args.inputs:
            stem = inp.stem
            r = combined[stem]
            verdict = []
            for key in (
                "criterion_1_min_angle",
                "criterion_2_max_angle",
                "criterion_3_max_slope",
                "criterion_4_area_change",
                "criterion_5_valence",
            ):
                blk = r.get(key, {})
                nf = blk.get("n_fail")
                if nf is None:
                    verdict.append("--")
                else:
                    verdict.append("✓" if nf == 0 else "✗")
            lines.append(
                f"  {inp.name:<48}  "
                + "  ".join(f"{v:>2}" for v in verdict)
            )
        print()
        print("\n".join(lines))
    return 0


if __name__ == "__main__":
    raise SystemExit(main())
