"""PoC #49 (H1): Phase H with FVCOM-compliant gate (min_angle=30°).

The user-provided FVCOM manual lists *min interior angle >= 30°* as
a hard mesh requirement. Phase H has historically defaulted to a
20° gate (and a per-element alpha >= 0.95 gate that PoC #48
showed is unrelated to FVCOM's actual requirements). PoC #48
measured Phase H v3's output against the FVCOM manual:

   C1 min angle >= 30°:  10 / 46 665 fail   (0.02 %)
   C2 max angle <= 130°:  1 / 46 665 fail   (0.002 %)
   C4 area change <= 0.5: 99 / 66 483 fail  (0.15 %)
   C5 valence <= 8:       0 / 27 349 fail   (PASS ✓)

— so Phase H v3 already removes 98 % of C1 failures despite never
explicitly targeting 30°. This PoC re-runs Phase H from the same
Phase G output ``outputs/33_pipeline_passing.14`` but with:

   alpha_target     = 0.5    (FVCOM-permissive; relaxed from 0.95)
   min_angle_target = 30.0   (FVCOM hard requirement)

so the optimiser drives strictly toward the 30° gate without
wasting effort on alpha refinements that have no FVCOM impact.

Coastline projection (v3 feature) stays on so new boundary nodes
land on the MLIT C23 polyline.

Outputs:
   outputs/49_phase_h_30deg_optimized.14
   outputs/49_phase_h_30deg_summary.txt
   outputs/49_phase_h_30deg_summary.json

Wall-time budget: PoC #43 (v3 with alpha=0.95, min_angle=20°) ran
in 1 735 s on this input. The 30° gate adds candidate elements but
the alpha=0.5 relaxation drops the fail set massively (from ~12 440
under v3 settings to ~452 under FVCOM-only settings), so the run
should be substantially shorter. Elapse cap 1 h.
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
from fvcom_mesh_tools.mesh_clean_phase_h import (
    build_coastline_projector,
    phase_h_optimize,
)
from fvcom_mesh_tools.quality import compute_metrics

REPO = Path(__file__).resolve().parent.parent
INPUT = REPO / "outputs" / "33_pipeline_passing.14"
COASTLINE = (
    REPO / "data" / "coastline" / "tokyo_bay"
    / "MLIT_C23" / "C23-06_TOKYOBAY.shp"
)
OUT_DIR = REPO / "outputs"
OUTPUT = OUT_DIR / "49_phase_h_30deg_optimized.14"
SUMMARY_TXT = OUT_DIR / "49_phase_h_30deg_summary.txt"
SUMMARY_JSON = OUT_DIR / "49_phase_h_30deg_summary.json"

V3_SUMMARY = OUT_DIR / "43_phase_h_v3_summary.json"

ALPHA_TARGET = 0.5
MIN_ANGLE_TARGET = 30.0
MAX_OUTER_ROUNDS = 10
MAX_TOPOLOGY_PER_ROUND = 10_000
MAX_SNAP_M = 500.0


def _signed_area_negative_count(mesh) -> int:
    p0 = mesh.nodes[mesh.elements[:, 0]]
    p1 = mesh.nodes[mesh.elements[:, 1]]
    p2 = mesh.nodes[mesh.elements[:, 2]]
    cross = (
        (p1[:, 0] - p0[:, 0]) * (p2[:, 1] - p0[:, 1])
        - (p1[:, 1] - p0[:, 1]) * (p2[:, 0] - p0[:, 0])
    )
    return int((cross <= 0).sum())


def _format_metric(b, a):
    if isinstance(b, (int, float)) and isinstance(a, (int, float)):
        return f"{a - b:+}" if isinstance(b, int) else f"{a - b:+.4f}"
    return "-"


def main() -> int:
    OUT_DIR.mkdir(parents=True, exist_ok=True)
    if not INPUT.exists():
        raise SystemExit(f"input missing: {INPUT}")
    if not COASTLINE.exists():
        raise SystemExit(f"coastline missing: {COASTLINE}")

    mesh_in = read_fort14(INPUT)
    print(
        f"[49] input: NP={mesh_in.n_nodes:,}  NE={mesh_in.n_elements:,}",
        flush=True,
    )
    a_before = alpha_quality(mesh_in)
    m_before = min_interior_angle(mesh_in)
    # Fail mask under FVCOM-only thresholds (relaxed alpha, strict angle).
    fail_before = ((a_before < ALPHA_TARGET) | (m_before < MIN_ANGLE_TARGET))
    n_fail_before = int(fail_before.sum())
    print(
        f"[49] fail elements (α<{ALPHA_TARGET} ∨ "
        f"min_ang<{MIN_ANGLE_TARGET}°): {n_fail_before:,} "
        f"({n_fail_before / mesh_in.n_elements:.4%})",
        flush=True,
    )
    n_below_30 = int((m_before < 30.0).sum())
    print(
        f"[49] elements with min_angle < 30°: {n_below_30:,}",
        flush=True,
    )

    print(
        f"[49] building coastline projector from {COASTLINE.name} "
        f"(max_snap_m={MAX_SNAP_M:g}) ...",
        flush=True,
    )
    t0 = time.perf_counter()
    projector = build_coastline_projector(
        [COASTLINE],
        max_snap_distance_m=MAX_SNAP_M,
        mean_latitude_deg=float(mesh_in.nodes[:, 1].mean()),
    )
    proj_build_wall = time.perf_counter() - t0
    if projector is None:
        raise SystemExit("coastline projector built no polylines")
    print(f"[49] projector built in {proj_build_wall:.1f} s", flush=True)

    print(
        f"[49] running Phase H v3 with FVCOM gate "
        f"(alpha={ALPHA_TARGET}, min_angle={MIN_ANGLE_TARGET}°, "
        f"max_outer_rounds={MAX_OUTER_ROUNDS}) ...",
        flush=True,
    )
    t0 = time.perf_counter()
    mesh_out, info = phase_h_optimize(
        mesh_in,
        alpha_target=ALPHA_TARGET,
        min_angle_target=MIN_ANGLE_TARGET,
        max_outer_rounds=MAX_OUTER_ROUNDS,
        max_topology_per_round=MAX_TOPOLOGY_PER_ROUND,
        coastline_projector=projector,
        # Lookahead / Pass D stay off — H1 is the canonical
        # baseline; opt-ins can be measured separately.
        lookahead_enabled=False,
        patch_recdt_enabled=False,
    )
    wall = time.perf_counter() - t0
    print(
        f"[49] Phase H wall: {wall:.1f} s  iters={info['n_iters']:,} "
        f"sweeps={info['n_smooth_sweeps']:,} "
        f"rounds={info['n_outer_rounds']:,}",
        flush=True,
    )

    a_after = alpha_quality(mesh_out)
    m_after = min_interior_angle(mesh_out)
    fail_after = ((a_after < ALPHA_TARGET) | (m_after < MIN_ANGLE_TARGET))
    n_fail_after = int(fail_after.sum())
    n_below_30_after = int((m_after < 30.0).sum())
    write_fort14(mesh_out, OUTPUT)
    print(
        f"[49] output: NP={mesh_out.n_nodes:,}  NE={mesh_out.n_elements:,}  "
        f"fail={n_fail_after:,}  min_angle<30°={n_below_30_after:,}",
        flush=True,
    )

    metrics_before = compute_metrics(mesh_in)
    metrics_after = compute_metrics(mesh_out)

    v3_payload = None
    if V3_SUMMARY.exists():
        try:
            v3_payload = json.loads(V3_SUMMARY.read_text())
        except json.JSONDecodeError:
            v3_payload = None

    payload = {
        "config": {
            "alpha_target": ALPHA_TARGET,
            "min_angle_target": MIN_ANGLE_TARGET,
            "max_outer_rounds": MAX_OUTER_ROUNDS,
            "max_topology_per_round": MAX_TOPOLOGY_PER_ROUND,
            "max_snap_m": MAX_SNAP_M,
            "coastline": str(COASTLINE.resolve()),
            "lookahead_enabled": False,
            "patch_recdt_enabled": False,
        },
        "input": {
            "path": str(INPUT.resolve()),
            "n_nodes": int(mesh_in.n_nodes),
            "n_elements": int(mesh_in.n_elements),
            "n_fail": n_fail_before,
            "n_below_30deg": n_below_30,
            "metrics": metrics_before,
        },
        "output": {
            "path": str(OUTPUT.resolve()),
            "n_nodes": int(mesh_out.n_nodes),
            "n_elements": int(mesh_out.n_elements),
            "n_fail": n_fail_after,
            "n_below_30deg": n_below_30_after,
            "metrics": metrics_after,
            "wall_seconds": wall,
            "n_signed_area_negative": _signed_area_negative_count(mesh_out),
            "projector_build_wall_seconds": proj_build_wall,
        },
        "phase_h_info": {
            k: (dict(v) if isinstance(v, dict)
                else list(v) if isinstance(v, (list, tuple))
                else v)
            for k, v in info.items()
            if not isinstance(v, np.ndarray)
        },
        "v3_reference": v3_payload,
    }
    SUMMARY_JSON.write_text(
        json.dumps(payload, indent=2, default=str), encoding="utf-8",
    )

    lines = [
        "PoC #49 (H1) Phase H with FVCOM gate (alpha=0.5, min_angle=30°)",
        f"input: {INPUT.name}",
        f"coastline: {COASTLINE.name}  max_snap_m={MAX_SNAP_M:g}",
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
        lines.append(
            f"  {k:<24}  {b!s:>14}  {a!s:>14}  {_format_metric(b, a):>10}"
        )
    lines.append(
        f"  {'min_angle<30°':<24}  {n_below_30!s:>14}  {n_below_30_after!s:>14}  "
        f"{n_below_30_after - n_below_30:+}"
    )
    lines.append(
        f"  {'fail_in_gate':<24}  {n_fail_before!s:>14}  {n_fail_after!s:>14}  "
        f"{n_fail_after - n_fail_before:+}"
    )
    lines.append(
        f"  {'wall_seconds':<24}  {'-':>14}  {wall:>14.1f}  -"
    )
    lines.append("")
    lines.append("  operators applied (Pass A + Pass B):")
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

    SUMMARY_TXT.write_text("\n".join(lines) + "\n", encoding="utf-8")
    print()
    print("\n".join(lines))
    print()
    print(f"wrote {OUTPUT}")
    print(f"wrote {SUMMARY_TXT}")
    print(f"wrote {SUMMARY_JSON}")
    return 0


if __name__ == "__main__":
    raise SystemExit(main())
