"""PoC #46: Phase H v4.1 end-to-end validation (target-exits-fail gate).

PoC #45 showed that the v4 ``union_penalty`` Pass C gate is too
permissive — accepts can drop the local union penalty without lifting
the failing target element out of fail status, leading to a thrash
where 20 000 accepts net to ≈ 0 fail-count change.

v4.1 replaces the gate with **``target_exits_fail``** — accept iff
the originally failing target element ``E`` exits fail status on the
post-op mesh (``alpha(E) >= alpha_target ∧ min_angle(E) >=
min_angle_target``) or was removed by the operator chain. The gate
is the SMS manual-edit standard: each Pass C accept is guaranteed to
fix at least one fail element.

v4.1 also enables an **op1-only** accept path under the new gate:
if op1 alone makes E exit fail, the search short-circuits without
trying op2 candidates. Faster and surfaces moves that Pass B's
strict gate rejected because they raised SOME local penalty
elsewhere even while fixing E.

This PoC re-runs PoC #43 / #45's exact pipeline-rung-1 input + MLIT
C23 Tokyo Bay coastline with ``lookahead_enabled=True`` and the new
``target_exits_fail`` gate, then writes a side-by-side A/B/C table
against v3 (#43) and v4 (#45).

Outputs:
    outputs/46_phase_h_v4_1_optimized.14
    outputs/46_phase_h_v4_1_summary.txt
    outputs/46_phase_h_v4_1_summary.json
"""
from __future__ import annotations

import json
import os
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
OUTPUT = OUT_DIR / "46_phase_h_v4_1_optimized.14"
SUMMARY_TXT = OUT_DIR / "46_phase_h_v4_1_summary.txt"
SUMMARY_JSON = OUT_DIR / "46_phase_h_v4_1_summary.json"

V3_SUMMARY = OUT_DIR / "43_phase_h_v3_summary.json"
V4_SUMMARY = OUT_DIR / "45_phase_h_v4_summary.json"

ALPHA_TARGET = 0.95
MIN_ANGLE_TARGET = 20.0
MAX_OUTER_ROUNDS = 10
MAX_TOPOLOGY_PER_ROUND = 10_000
MAX_LOOKAHEAD_PER_ROUND = 2_000
MAX_SNAP_M = 500.0


def _env_int(name: str, default: int) -> int:
    raw = os.environ.get(name)
    if raw is None or raw == "":
        return default
    return int(raw)


MAX_LOOKAHEAD_PER_ROUND = _env_int(
    "POC46_MAX_LOOKAHEAD", MAX_LOOKAHEAD_PER_ROUND,
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


def _load_summary(path: Path) -> dict | None:
    if not path.exists():
        return None
    try:
        return json.loads(path.read_text())
    except json.JSONDecodeError:
        return None


def main() -> int:
    OUT_DIR.mkdir(parents=True, exist_ok=True)
    if not INPUT.exists():
        raise SystemExit(f"input missing: {INPUT}")
    if not COASTLINE.exists():
        raise SystemExit(f"coastline missing: {COASTLINE}")

    mesh_in = read_fort14(INPUT)
    print(
        f"[46] input: NP={mesh_in.n_nodes:,}  NE={mesh_in.n_elements:,}",
        flush=True,
    )
    a_before = alpha_quality(mesh_in)
    m_before = min_interior_angle(mesh_in)
    fail_before = ((a_before < ALPHA_TARGET) | (m_before < MIN_ANGLE_TARGET))
    n_fail_before = int(fail_before.sum())
    print(
        f"[46] fail elements: {n_fail_before:,} "
        f"({n_fail_before / mesh_in.n_elements:.4%})",
        flush=True,
    )

    print(
        f"[46] building coastline projector from {COASTLINE.name} "
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
    print(f"[46] projector built in {proj_build_wall:.1f} s", flush=True)

    print(
        f"[46] running Phase H v4.1 "
        f"(lookahead_gate=target_exits_fail, "
        f"max_outer_rounds={MAX_OUTER_ROUNDS}, "
        f"max_lookahead_per_round={MAX_LOOKAHEAD_PER_ROUND}) ...",
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
        lookahead_enabled=True,
        max_lookahead_per_round=MAX_LOOKAHEAD_PER_ROUND,
        lookahead_gate="target_exits_fail",
    )
    wall = time.perf_counter() - t0
    print(
        f"[46] Phase H v4.1 wall: {wall:.1f} s  iters={info['n_iters']:,} "
        f"sweeps={info['n_smooth_sweeps']:,} "
        f"rounds={info['n_outer_rounds']:,}",
        flush=True,
    )

    a_after = alpha_quality(mesh_out)
    m_after = min_interior_angle(mesh_out)
    fail_after = ((a_after < ALPHA_TARGET) | (m_after < MIN_ANGLE_TARGET))
    n_fail_after = int(fail_after.sum())

    write_fort14(mesh_out, OUTPUT)

    metrics_before = compute_metrics(mesh_in)
    metrics_after = compute_metrics(mesh_out)

    v3_payload = _load_summary(V3_SUMMARY)
    v4_payload = _load_summary(V4_SUMMARY)

    payload = {
        "config": {
            "alpha_target": ALPHA_TARGET,
            "min_angle_target": MIN_ANGLE_TARGET,
            "max_outer_rounds": MAX_OUTER_ROUNDS,
            "max_topology_per_round": MAX_TOPOLOGY_PER_ROUND,
            "max_lookahead_per_round": MAX_LOOKAHEAD_PER_ROUND,
            "max_snap_m": MAX_SNAP_M,
            "coastline": str(COASTLINE.resolve()),
            "lookahead_enabled": True,
            "lookahead_gate": "target_exits_fail",
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
        "v4_reference": v4_payload,
    }
    SUMMARY_JSON.write_text(
        json.dumps(payload, indent=2, default=str), encoding="utf-8",
    )

    lines = [
        "PoC #46 Phase H v4.1 end-to-end validation (target-exits-fail gate)",
        f"input: {INPUT.name}",
        f"coastline: {COASTLINE.name}  max_snap_m={MAX_SNAP_M:g}",
        f"thresholds: alpha>={ALPHA_TARGET}  min_angle>={MIN_ANGLE_TARGET}°",
        (
            f"caps: outer_rounds={MAX_OUTER_ROUNDS}  "
            f"topology_per_round={MAX_TOPOLOGY_PER_ROUND}  "
            f"lookahead_per_round={MAX_LOOKAHEAD_PER_ROUND}"
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
    lines.append(
        f"  {'projector build_s':<24}  {'-':>14}  "
        f"{proj_build_wall:>14.1f}  -"
    )
    lines.append("")
    lines.append("  operators applied (Pass A + Pass B):")
    for op_name, n in info.get("operators_applied", {}).items():
        lines.append(f"    {op_name:<24}  applied={n:>6,}")
    lines.append("  lookahead pairs applied (Pass C):")
    pairs = info.get("lookahead_pairs_applied", {})
    if pairs:
        for pair_label, n in pairs.items():
            lines.append(f"    {pair_label:<24}  applied={n:>6,}")
    else:
        lines.append("    (none)")
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
        f"  Pass B abandoned (no op helped): {info['n_abandoned']:,}"
    )
    lines.append(
        f"  Pass C abandoned (no pair)     : "
        f"{info['n_lookahead_abandoned']:,}"
    )

    # A/B/C table vs v3 (PoC #43) and v4 (PoC #45) when available.
    def _stat_from(payload_, k):
        if payload_ is None:
            return None
        out = payload_.get("output") or {}
        if k in {"n_nodes", "n_elements", "n_fail", "wall_seconds"}:
            return out.get(k)
        return (out.get("metrics") or {}).get(k)

    v4_1_values: dict[str, object] = {
        "n_nodes": mesh_out.n_nodes,
        "n_elements": mesh_out.n_elements,
        "n_fail": n_fail_after,
        "wall_seconds": wall,
    }

    if v3_payload is not None or v4_payload is not None:
        lines.append("")
        lines.append("A/B/C — v3 (PoC #43) vs v4 (PoC #45) vs v4.1 (this PoC):")
        for k in (
            "n_elements", "alpha_mean", "alpha_p05", "min_angle_p05_deg",
            "frac_lt_20deg", "n_fail", "wall_seconds",
        ):
            v3 = _stat_from(v3_payload, k)
            v4 = _stat_from(v4_payload, k)
            v4_1 = v4_1_values.get(k, metrics_after.get(k))
            lines.append(
                f"  {k:<22}  v3={v3!s:>14}  v4={v4!s:>14}  "
                f"v4.1={v4_1!s:>14}"
            )

    SUMMARY_TXT.write_text("\n".join(lines) + "\n", encoding="utf-8")
    print()
    print("\n".join(lines))
    return 0


if __name__ == "__main__":
    raise SystemExit(main())
