"""PoC #42: Phase H v2 end-to-end validation (boundary-aware operators).

PoC #41 ran Phase H v1 on the pipeline-rung-1 output of PoC #19 and
left 11,182 abandoned fail elements — of which PoC #40 measured 84 %
touch a boundary. v2 adds two boundary-aware operators:

  * Boundary-tangent smooth in Pass A — segment-interior boundary
    nodes move along the prev-next tangent line via 1-ring centroid
    projection (segment endpoints / corners stay pinned).
  * ``edge_split_boundary`` in Pass B — insert a midpoint on a
    topological boundary edge, split the single incident triangle
    into two, and thread the new node into the segment.

This PoC runs Phase H v2 with the same defaults as PoC #41 and
reports the delta against v1.

Outputs:
    outputs/42_phase_h_v2_optimized.14
    outputs/42_phase_h_v2_summary.txt
    outputs/42_phase_h_v2_summary.json
"""
from __future__ import annotations

import json
import time
from pathlib import Path

import numpy as np

from fvcom_mesh_tools.algorithms.quality import (
    alpha_quality,
    min_interior_angle,
)
from fvcom_mesh_tools.io import read_fort14, write_fort14
from fvcom_mesh_tools.mesh_clean_phase_h import phase_h_optimize
from fvcom_mesh_tools.quality import compute_metrics

REPO = Path(__file__).resolve().parent.parent
INPUT = REPO / "outputs" / "33_pipeline_passing.14"
OUT_DIR = REPO / "outputs"
OUTPUT = OUT_DIR / "42_phase_h_v2_optimized.14"
SUMMARY_TXT = OUT_DIR / "42_phase_h_v2_summary.txt"
SUMMARY_JSON = OUT_DIR / "42_phase_h_v2_summary.json"

ALPHA_TARGET = 0.95
MIN_ANGLE_TARGET = 20.0
MAX_OUTER_ROUNDS = 10
MAX_TOPOLOGY_PER_ROUND = 10_000


def _signed_area_negative_count(mesh) -> int:
    p0 = mesh.nodes[mesh.elements[:, 0]]
    p1 = mesh.nodes[mesh.elements[:, 1]]
    p2 = mesh.nodes[mesh.elements[:, 2]]
    cross = (
        (p1[:, 0] - p0[:, 0]) * (p2[:, 1] - p0[:, 1])
        - (p1[:, 1] - p0[:, 1]) * (p2[:, 0] - p0[:, 0])
    )
    return int((cross <= 0).sum())


def main() -> int:
    OUT_DIR.mkdir(parents=True, exist_ok=True)
    if not INPUT.exists():
        raise SystemExit(f"input missing: {INPUT}")

    mesh_in = read_fort14(INPUT)
    print(f"[42] input: NP={mesh_in.n_nodes:,}  NE={mesh_in.n_elements:,}")
    a_before = alpha_quality(mesh_in)
    m_before = min_interior_angle(mesh_in)
    fail_before = ((a_before < ALPHA_TARGET) | (m_before < MIN_ANGLE_TARGET))
    n_fail_before = int(fail_before.sum())
    print(
        f"[42] fail elements: {n_fail_before:,}  "
        f"({n_fail_before / mesh_in.n_elements:.4%}) at "
        f"alpha>={ALPHA_TARGET}  min_angle>={MIN_ANGLE_TARGET}°"
    )

    print(
        f"[42] running Phase H v2 "
        f"(max_outer_rounds={MAX_OUTER_ROUNDS}, "
        f"max_topology_per_round={MAX_TOPOLOGY_PER_ROUND}) ..."
    )
    t0 = time.perf_counter()
    mesh_out, info = phase_h_optimize(
        mesh_in,
        alpha_target=ALPHA_TARGET,
        min_angle_target=MIN_ANGLE_TARGET,
        max_outer_rounds=MAX_OUTER_ROUNDS,
        max_topology_per_round=MAX_TOPOLOGY_PER_ROUND,
    )
    wall = time.perf_counter() - t0
    print(
        f"[42] Phase H v2 wall: {wall:.1f} s  iters={info['n_iters']:,}  "
        f"sweeps={info['n_smooth_sweeps']:,}  "
        f"rounds={info['n_outer_rounds']:,}"
    )

    a_after = alpha_quality(mesh_out)
    m_after = min_interior_angle(mesh_out)
    fail_after = ((a_after < ALPHA_TARGET) | (m_after < MIN_ANGLE_TARGET))
    n_fail_after = int(fail_after.sum())

    write_fort14(mesh_out, OUTPUT)

    metrics_before = compute_metrics(mesh_in)
    metrics_after = compute_metrics(mesh_out)

    payload = {
        "config": {
            "alpha_target": ALPHA_TARGET,
            "min_angle_target": MIN_ANGLE_TARGET,
            "max_outer_rounds": MAX_OUTER_ROUNDS,
            "max_topology_per_round": MAX_TOPOLOGY_PER_ROUND,
        },
        "input": {
            "path": str(INPUT.resolve()),
            "n_nodes": int(mesh_in.n_nodes),
            "n_elements": int(mesh_in.n_elements),
            "n_fail": n_fail_before,
            "metrics": metrics_before,
        },
        "output": {
            "path": str(OUTPUT.resolve()),
            "n_nodes": int(mesh_out.n_nodes),
            "n_elements": int(mesh_out.n_elements),
            "n_fail": n_fail_after,
            "metrics": metrics_after,
            "wall_seconds": wall,
            "n_signed_area_negative": _signed_area_negative_count(mesh_out),
        },
        "phase_h_info": {
            k: (dict(v) if isinstance(v, dict)
                else list(v) if isinstance(v, (list, tuple))
                else v)
            for k, v in info.items()
            if not isinstance(v, np.ndarray)
        },
    }
    SUMMARY_JSON.write_text(json.dumps(payload, indent=2, default=str),
                            encoding="utf-8")

    lines = [
        "PoC #42 Phase H v2 end-to-end validation",
        f"input: {INPUT.name}",
        f"thresholds: alpha>={ALPHA_TARGET}  min_angle>={MIN_ANGLE_TARGET}°",
        "",
        f"  {'metric':<24}  {'input':>14}  {'output':>14}  {'delta':>10}",
        "  " + "-" * 70,
    ]
    for k in (
        "n_nodes", "n_elements",
        "alpha_mean", "alpha_p05", "min_angle_p05_deg", "frac_lt_20deg",
        "max_valence", "n_overconnected", "n_flipped",
    ):
        if k == "n_nodes":
            b, a = mesh_in.n_nodes, mesh_out.n_nodes
        elif k == "n_elements":
            b, a = mesh_in.n_elements, mesh_out.n_elements
        else:
            b = metrics_before.get(k)
            a = metrics_after.get(k)
        if isinstance(b, (int, float)) and isinstance(a, (int, float)):
            d = f"{a - b:+}" if isinstance(b, int) else f"{a - b:+.4f}"
        else:
            d = "-"
        lines.append(f"  {k:<24}  {b!s:>14}  {a!s:>14}  {d:>10}")
    lines.append(
        f"  {'fail elements':<24}  {n_fail_before!s:>14}  {n_fail_after!s:>14}  "
        f"{n_fail_after - n_fail_before:+}"
    )
    lines.append(
        f"  {'fail fraction':<24}  "
        f"{n_fail_before / max(mesh_in.n_elements, 1):>14.4%}  "
        f"{n_fail_after / max(mesh_out.n_elements, 1):>14.4%}  -"
    )
    lines.append(
        f"  {'wall_seconds':<24}  {'-':>14}  {wall:>14.1f}  -"
    )
    lines.append("")
    lines.append("  operators applied:")
    for op_name, n in info.get("operators_applied", {}).items():
        lines.append(f"    {op_name:<24}  applied={n:>6,}")
    lines.append(
        f"  iterations                     : {info['n_iters']:,}"
    )
    lines.append(
        f"  smooth sweeps                  : {info['n_smooth_sweeps']:,}"
    )
    lines.append(
        f"  outer rounds                   : {info['n_outer_rounds']:,}"
    )
    lines.append(
        f"  abandoned (no op helped)       : {info['n_abandoned']:,}"
    )

    SUMMARY_TXT.write_text("\n".join(lines) + "\n", encoding="utf-8")
    print()
    print("\n".join(lines))
    print(f"\nwrote {OUTPUT}")
    print(f"wrote {SUMMARY_TXT}")
    print(f"wrote {SUMMARY_JSON}")
    return 0


if __name__ == "__main__":
    raise SystemExit(main())
