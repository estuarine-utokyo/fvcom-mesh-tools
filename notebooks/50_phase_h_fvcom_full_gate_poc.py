"""PoC #50 (I1+I2): Phase H with the full FVCOM gate
(alpha>=0.95, min_angle>=30°, max_angle<=130°).

Context — what changed since PoC #49:

  * PoC #48 measured Phase H v3 default output against the FVCOM
    manual and showed v3 (alpha=0.95, min_angle=20°) already removes
    98 % of the C1 (min angle >= 30°) failures incidentally, while
    leaving 1 / 46 665 element violating C2 (max angle <= 130°).
  * PoC #49 (H1) relaxed alpha to 0.5 and raised min_angle to 30°,
    expecting the optimiser to crack the residual C1 fails. The
    relaxed alpha gate proved *counter-productive* — alpha quality
    and min-angle correlate strongly, so v3's strict alpha=0.95
    indirectly fixed more 30° fails than the explicit but slack
    PoC #49 configuration. C1 went 10 -> 13 instead of 10 -> 0.

I1+I2 keeps the strict alpha=0.95 gate (so the indirect lift
survives), tightens C1 to 30° (so the optimiser stops calling 21°
elements "passing"), and adds the C2 max-angle gate at 130° so
Phase H is finally aware of the obtuse-angle requirement. This is
the canonical FVCOM-aware Phase H configuration; PoC #48 is the
no-op baseline to compare against.

Configuration:

   alpha_target     = 0.95   (keep v3's strict gate)
   min_angle_target = 30.0   (FVCOM C1 hard requirement)
   max_angle_target = 130.0  (FVCOM C2 hard requirement, NEW in I2)

Coastline projection stays on so new boundary nodes land on the
MLIT C23 polyline. Pass C lookahead and Pass D patch re-CDT stay
*off* — this PoC is meant to baseline I1+I2 alone before stacking
optional opt-ins on top.

Outputs:
   outputs/50_phase_h_fvcom_full_gate_optimized.14
   outputs/50_phase_h_fvcom_full_gate_summary.txt
   outputs/50_phase_h_fvcom_full_gate_summary.json

Wall-time budget: PoC #43 (v3 with alpha=0.95, min_angle=20°) ran
in 1 735 s and PoC #49 (alpha=0.5, min_angle=30°) ran in 291 s.
With the strict alpha=0.95 gate restored the fail set widens vs.
PoC #49, so wall time is expected somewhere in between
(~30-50 min on the same Tokyo Bay input). Elapse cap 1 h 30 min.
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
OUTPUT = OUT_DIR / "50_phase_h_fvcom_full_gate_optimized.14"
SUMMARY_TXT = OUT_DIR / "50_phase_h_fvcom_full_gate_summary.txt"
SUMMARY_JSON = OUT_DIR / "50_phase_h_fvcom_full_gate_summary.json"

V3_SUMMARY = OUT_DIR / "43_phase_h_v3_summary.json"
H1_SUMMARY = OUT_DIR / "49_phase_h_30deg_summary.json"

ALPHA_TARGET = 0.95
MIN_ANGLE_TARGET = 30.0
MAX_ANGLE_TARGET = 130.0
MAX_OUTER_ROUNDS = 10
MAX_TOPOLOGY_PER_ROUND = 10_000
MAX_SNAP_M = 500.0


def _max_interior_angle(mesh) -> np.ndarray:
    """Per-element max interior angle in degrees."""
    p0 = mesh.nodes[mesh.elements[:, 0]]
    p1 = mesh.nodes[mesh.elements[:, 1]]
    p2 = mesh.nodes[mesh.elements[:, 2]]
    e0 = np.linalg.norm(p1 - p0, axis=1)  # opp v2
    e1 = np.linalg.norm(p2 - p1, axis=1)  # opp v0
    e2 = np.linalg.norm(p0 - p2, axis=1)  # opp v1

    def _ang(opp, ea, eb):
        denom = 2.0 * ea * eb
        safe = np.where(denom > 0, denom, 1.0)
        cos = np.where(
            denom > 0, (ea * ea + eb * eb - opp * opp) / safe, 0.0,
        )
        return np.arccos(np.clip(cos, -1.0, 1.0))

    angle_v0 = _ang(e1, e2, e0)
    angle_v1 = _ang(e2, e0, e1)
    angle_v2 = _ang(e0, e1, e2)
    return np.degrees(
        np.maximum(np.maximum(angle_v0, angle_v1), angle_v2),
    )


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
        f"[50] input: NP={mesh_in.n_nodes:,}  NE={mesh_in.n_elements:,}",
        flush=True,
    )
    a_before = alpha_quality(mesh_in)
    m_before = min_interior_angle(mesh_in)
    M_before = _max_interior_angle(mesh_in)
    fail_before = (
        (a_before < ALPHA_TARGET)
        | (m_before < MIN_ANGLE_TARGET)
        | (M_before > MAX_ANGLE_TARGET)
    )
    n_fail_before = int(fail_before.sum())
    n_c1_before = int((m_before < MIN_ANGLE_TARGET).sum())
    n_c2_before = int((M_before > MAX_ANGLE_TARGET).sum())
    print(
        f"[50] fail elements (α<{ALPHA_TARGET} ∨ min_ang<{MIN_ANGLE_TARGET}° "
        f"∨ max_ang>{MAX_ANGLE_TARGET}°): {n_fail_before:,} "
        f"({n_fail_before / mesh_in.n_elements:.4%})",
        flush=True,
    )
    print(
        f"[50] C1 violations (min_ang < 30°): {n_c1_before:,}",
        flush=True,
    )
    print(
        f"[50] C2 violations (max_ang > 130°): {n_c2_before:,}",
        flush=True,
    )

    print(
        f"[50] building coastline projector from {COASTLINE.name} "
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
    print(f"[50] projector built in {proj_build_wall:.1f} s", flush=True)

    print(
        f"[50] running Phase H with FVCOM full gate "
        f"(alpha={ALPHA_TARGET}, min_angle={MIN_ANGLE_TARGET}°, "
        f"max_angle={MAX_ANGLE_TARGET}°, "
        f"max_outer_rounds={MAX_OUTER_ROUNDS}) ...",
        flush=True,
    )
    t0 = time.perf_counter()
    mesh_out, info = phase_h_optimize(
        mesh_in,
        alpha_target=ALPHA_TARGET,
        min_angle_target=MIN_ANGLE_TARGET,
        max_angle_target=MAX_ANGLE_TARGET,
        max_outer_rounds=MAX_OUTER_ROUNDS,
        max_topology_per_round=MAX_TOPOLOGY_PER_ROUND,
        coastline_projector=projector,
        lookahead_enabled=False,
        patch_recdt_enabled=False,
    )
    wall = time.perf_counter() - t0
    print(
        f"[50] Phase H wall: {wall:.1f} s  iters={info['n_iters']:,} "
        f"sweeps={info['n_smooth_sweeps']:,} "
        f"rounds={info['n_outer_rounds']:,}",
        flush=True,
    )

    a_after = alpha_quality(mesh_out)
    m_after = min_interior_angle(mesh_out)
    M_after = _max_interior_angle(mesh_out)
    fail_after = (
        (a_after < ALPHA_TARGET)
        | (m_after < MIN_ANGLE_TARGET)
        | (M_after > MAX_ANGLE_TARGET)
    )
    n_fail_after = int(fail_after.sum())
    n_c1_after = int((m_after < MIN_ANGLE_TARGET).sum())
    n_c2_after = int((M_after > MAX_ANGLE_TARGET).sum())
    write_fort14(mesh_out, OUTPUT)
    print(
        f"[50] output: NP={mesh_out.n_nodes:,}  NE={mesh_out.n_elements:,}  "
        f"fail={n_fail_after:,}  C1={n_c1_after:,}  C2={n_c2_after:,}",
        flush=True,
    )

    metrics_before = compute_metrics(mesh_in)
    metrics_after = compute_metrics(mesh_out)

    def _load_json(p):
        if not p.exists():
            return None
        try:
            return json.loads(p.read_text())
        except json.JSONDecodeError:
            return None

    v3_payload = _load_json(V3_SUMMARY)
    h1_payload = _load_json(H1_SUMMARY)

    payload = {
        "config": {
            "alpha_target": ALPHA_TARGET,
            "min_angle_target": MIN_ANGLE_TARGET,
            "max_angle_target": MAX_ANGLE_TARGET,
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
            "n_c1_violations": n_c1_before,
            "n_c2_violations": n_c2_before,
            "metrics": metrics_before,
        },
        "output": {
            "path": str(OUTPUT.resolve()),
            "n_nodes": int(mesh_out.n_nodes),
            "n_elements": int(mesh_out.n_elements),
            "n_fail": n_fail_after,
            "n_c1_violations": n_c1_after,
            "n_c2_violations": n_c2_after,
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
        "h1_reference": h1_payload,
    }
    SUMMARY_JSON.write_text(
        json.dumps(payload, indent=2, default=str), encoding="utf-8",
    )

    lines = [
        "PoC #50 (I1+I2) Phase H with full FVCOM gate "
        f"(alpha={ALPHA_TARGET}, min_angle={MIN_ANGLE_TARGET}°, "
        f"max_angle={MAX_ANGLE_TARGET}°)",
        f"input: {INPUT.name}",
        f"coastline: {COASTLINE.name}  max_snap_m={MAX_SNAP_M:g}",
        (
            f"thresholds: alpha>={ALPHA_TARGET}  "
            f"min_angle>={MIN_ANGLE_TARGET}°  "
            f"max_angle<={MAX_ANGLE_TARGET}°"
        ),
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
        f"  {'C1: min_ang<30°':<24}  {n_c1_before!s:>14}  "
        f"{n_c1_after!s:>14}  {n_c1_after - n_c1_before:+}"
    )
    lines.append(
        f"  {'C2: max_ang>130°':<24}  {n_c2_before!s:>14}  "
        f"{n_c2_after!s:>14}  {n_c2_after - n_c2_before:+}"
    )
    lines.append(
        f"  {'fail_in_gate':<24}  {n_fail_before!s:>14}  "
        f"{n_fail_after!s:>14}  {n_fail_after - n_fail_before:+}"
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
