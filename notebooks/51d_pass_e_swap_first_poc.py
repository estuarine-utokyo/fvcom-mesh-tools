"""PoC #51d: Pass E with edge_swap-first strategy.

PoC #51c (Pass E with full split gates) achieved C4 102 -> 86 in
108 s with 14 accepts. The remaining 86 fails are not addressable
by midpoint edge_split alone — either because every candidate
would trigger an indirect regression, or because the resulting
sub-triangles would fail C1/C2.

PoC #51d adds an ``edge_swap`` operator to Pass E. For each C4
fail edge, the round driver now tries ``edge_swap`` of the fail
edge *first* and only falls back to ``edge_split`` if the swap
cannot accept. ``edge_swap`` adds no vertices, has net-zero
valence change globally (the endpoints lose 1, the opposite
vertices gain 1), and topologically just re-triangulates the
``L ∪ S`` quad with the other diagonal. The Pass E swap
operator gates:

  1. The new diagonal CD's area_change must be at or below
     ``area_ratio_target``. If swap merely relocates the
     violation, reject — we want strict ``-1`` on the fail count.
  2. No edge in the local 5-edge block (CD + the four edges
     ``a-C``, ``a-D``, ``b-C``, ``b-D``) may newly transition
     from non-fail to fail.
  3. No element in the new block may fail C1 or C2.
  4. ``C5`` valence cap: both ``C`` and ``D`` must have
     ``valence_before < max_valence``.

Outputs:
   outputs/51d_pass_e_swap_first_optimized.14
   outputs/51d_pass_e_swap_first_summary.txt
   outputs/51d_pass_e_swap_first_summary.json
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
from fvcom_mesh_tools.diagnostics import node_valence
from fvcom_mesh_tools.io import read_fort14, write_fort14
from fvcom_mesh_tools.mesh_clean_phase_h import (
    DEFAULT_MAX_VALENCE,
    _per_edge_area_change,
    build_coastline_projector,
    phase_h_optimize,
)
from fvcom_mesh_tools.quality import compute_metrics

REPO = Path(__file__).resolve().parent.parent
INPUT = REPO / "outputs" / "50_phase_h_fvcom_full_gate_optimized.14"
COASTLINE = (
    REPO / "data" / "coastline" / "tokyo_bay"
    / "MLIT_C23" / "C23-06_TOKYOBAY.shp"
)
OUT_DIR = REPO / "outputs"
OUTPUT = OUT_DIR / "51d_pass_e_swap_first_optimized.14"
SUMMARY_TXT = OUT_DIR / "51d_pass_e_swap_first_summary.txt"
SUMMARY_JSON = OUT_DIR / "51d_pass_e_swap_first_summary.json"

POC51C_SUMMARY = OUT_DIR / "51c_pass_e_indirect_filter_summary.json"

ALPHA_TARGET = 0.95
MIN_ANGLE_TARGET = 30.0
MAX_ANGLE_TARGET = 130.0
AREA_RATIO_TARGET = 0.5
MAX_VALENCE = DEFAULT_MAX_VALENCE
MAX_OUTER_ROUNDS = 10
MAX_SNAP_M = 500.0


def _max_interior_angle(mesh) -> np.ndarray:
    p0 = mesh.nodes[mesh.elements[:, 0]]
    p1 = mesh.nodes[mesh.elements[:, 1]]
    p2 = mesh.nodes[mesh.elements[:, 2]]
    e0 = np.linalg.norm(p1 - p0, axis=1)
    e1 = np.linalg.norm(p2 - p1, axis=1)
    e2 = np.linalg.norm(p0 - p2, axis=1)

    def _ang(opp, ea, eb):
        denom = 2.0 * ea * eb
        safe = np.where(denom > 0, denom, 1.0)
        cos = np.where(
            denom > 0, (ea * ea + eb * eb - opp * opp) / safe, 0.0,
        )
        return np.arccos(np.clip(cos, -1.0, 1.0))

    return np.degrees(
        np.maximum(
            np.maximum(_ang(e1, e2, e0), _ang(e2, e0, e1)),
            _ang(e0, e1, e2),
        ),
    )


def _c4_fail_count(mesh, target: float = AREA_RATIO_TARGET) -> int:
    _u, _p, ac = _per_edge_area_change(mesh.nodes, mesh.elements)
    return int((ac > target).sum())


def _c5_fail_count(mesh, cap: int = MAX_VALENCE) -> int:
    val = node_valence(mesh.elements, mesh.n_nodes)
    return int((val > cap).sum())


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
        f"[51d] input: NP={mesh_in.n_nodes:,}  NE={mesh_in.n_elements:,}",
        flush=True,
    )
    m_before = min_interior_angle(mesh_in)
    M_before = _max_interior_angle(mesh_in)
    n_c1_before = int((m_before < MIN_ANGLE_TARGET).sum())
    n_c2_before = int((M_before > MAX_ANGLE_TARGET).sum())
    n_c4_before = _c4_fail_count(mesh_in)
    n_c5_before = _c5_fail_count(mesh_in)
    print(
        f"[51d] FVCOM residuals before Pass E: "
        f"C1={n_c1_before:,}  C2={n_c2_before:,}  "
        f"C4={n_c4_before:,}  C5={n_c5_before:,}",
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
    print(f"[51d] projector built in {proj_build_wall:.1f} s", flush=True)

    print(
        f"[51d] running Pass E only with swap-first strategy "
        f"(area_ratio_target={AREA_RATIO_TARGET}, "
        f"max_valence={MAX_VALENCE}) ...",
        flush=True,
    )
    t0 = time.perf_counter()
    mesh_out, info = phase_h_optimize(
        mesh_in,
        alpha_target=ALPHA_TARGET,
        min_angle_target=MIN_ANGLE_TARGET,
        max_angle_target=MAX_ANGLE_TARGET,
        max_smooth_sweeps=0,
        max_topology_per_round=0,
        max_outer_rounds=MAX_OUTER_ROUNDS,
        coastline_projector=projector,
        lookahead_enabled=False,
        patch_recdt_enabled=False,
        pass_e_enabled=True,
        pass_e_area_ratio_target=AREA_RATIO_TARGET,
        pass_e_max_valence=MAX_VALENCE,
        max_pass_e_splits_per_round=10_000,
    )
    wall = time.perf_counter() - t0
    print(
        f"[51d] Phase H (Pass E only) wall: {wall:.1f} s  "
        f"pass_e_accepts={info.get('pass_e_accepts', 0):,}  "
        f"pass_e_rejected={info.get('pass_e_rejected', 0):,}",
        flush=True,
    )

    m_after = min_interior_angle(mesh_out)
    M_after = _max_interior_angle(mesh_out)
    n_c1_after = int((m_after < MIN_ANGLE_TARGET).sum())
    n_c2_after = int((M_after > MAX_ANGLE_TARGET).sum())
    n_c4_after = _c4_fail_count(mesh_out)
    n_c5_after = _c5_fail_count(mesh_out)
    write_fort14(mesh_out, OUTPUT)
    print(
        f"[51d] FVCOM residuals after Pass E: "
        f"C1={n_c1_after:,}  C2={n_c2_after:,}  "
        f"C4={n_c4_after:,}  C5={n_c5_after:,}",
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

    poc51c_payload = _load_json(POC51C_SUMMARY)

    payload = {
        "config": {
            "alpha_target": ALPHA_TARGET,
            "min_angle_target": MIN_ANGLE_TARGET,
            "max_angle_target": MAX_ANGLE_TARGET,
            "area_ratio_target": AREA_RATIO_TARGET,
            "max_valence": MAX_VALENCE,
            "max_outer_rounds": MAX_OUTER_ROUNDS,
            "max_snap_m": MAX_SNAP_M,
            "coastline": str(COASTLINE.resolve()),
            "pass_e_enabled": True,
            "operator_change_vs_poc51c": (
                "swap-first: try _apply_pass_e_swap before "
                "_apply_pass_e_split on each C4 fail edge"
            ),
        },
        "input": {
            "path": str(INPUT.resolve()),
            "n_nodes": int(mesh_in.n_nodes),
            "n_elements": int(mesh_in.n_elements),
            "n_c1_violations": n_c1_before,
            "n_c2_violations": n_c2_before,
            "n_c4_violations": n_c4_before,
            "n_c5_violations": n_c5_before,
            "metrics": metrics_before,
        },
        "output": {
            "path": str(OUTPUT.resolve()),
            "n_nodes": int(mesh_out.n_nodes),
            "n_elements": int(mesh_out.n_elements),
            "n_c1_violations": n_c1_after,
            "n_c2_violations": n_c2_after,
            "n_c4_violations": n_c4_after,
            "n_c5_violations": n_c5_after,
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
        "poc51c_reference": poc51c_payload,
    }
    SUMMARY_JSON.write_text(
        json.dumps(payload, indent=2, default=str), encoding="utf-8",
    )

    lines = [
        "PoC #51d — Pass E with edge_swap-first strategy",
        f"input: {INPUT.name}",
        f"coastline: {COASTLINE.name}  max_snap_m={MAX_SNAP_M:g}",
        (
            f"thresholds: alpha>={ALPHA_TARGET}  "
            f"min_angle>={MIN_ANGLE_TARGET}°  "
            f"max_angle<={MAX_ANGLE_TARGET}°  "
            f"area_change<={AREA_RATIO_TARGET}  "
            f"valence<={MAX_VALENCE}"
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
        f"  {'C4: area_change>0.5':<24}  {n_c4_before!s:>14}  "
        f"{n_c4_after!s:>14}  {n_c4_after - n_c4_before:+}"
    )
    lines.append(
        f"  {'C5: valence>8':<24}  {n_c5_before!s:>14}  "
        f"{n_c5_after!s:>14}  {n_c5_after - n_c5_before:+}"
    )
    lines.append(
        f"  {'wall_seconds':<24}  {'-':>14}  {wall:>14.1f}  -"
    )
    lines.append("")
    pe_acc = int(info.get("pass_e_accepts", 0))
    pe_rej = int(info.get("pass_e_rejected", 0))
    lines.append(
        f"  Pass E accepts                 : {pe_acc:,}"
    )
    lines.append(
        f"  Pass E rejected                : {pe_rej:,}"
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
