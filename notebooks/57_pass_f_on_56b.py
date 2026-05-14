"""PoC #57: Pass F (C4-aware smoothing) on PoC #56b output.

PoC #56c showed the residual 70 violations on
``56b_phase_h_min5_optimized.14`` are:

  * C1 = 2: stuck thin triangles with 1 boundary edge, unfixable by
    single-operator under count-gate.
  * C4 = 68: 61 small isolated clusters (88 % pairs, none larger
    than 4 elements), worst area_ratio = 0.685, most in [0.5, 0.6].

Pass F (new) addresses the C4 residual via Laplacian / boundary-
tangent smoothing of nodes that touch any current C4 fail edge,
with a count-comparison gate: C1/C2 must not increase, C4 must
strictly decrease in the local block. C5 is invariant under
smoothing so it is not gated.

This PoC runs:

  Stage 1: Pass F ONLY (operator_order=()) on PoC #56b input.
           Isolates the marginal effect of Pass F on its own.
  Stage 2: A+B+E+F full Phase H (with Pass F newly enabled).
           Measures the integrated effect.

Outputs:
    outputs/57_phase_h_pass_f_only.14
    outputs/57_phase_h_abef.14
    outputs/57_summary.{txt,json}
"""
from __future__ import annotations

import json
import time
from pathlib import Path

import numpy as np

from fvcom_mesh_tools.algorithms.quality import min_interior_angle
from fvcom_mesh_tools.diagnostics import node_valence
from fvcom_mesh_tools.io import read_fort14, write_fort14
from fvcom_mesh_tools.mesh_clean_phase_h import (
    _per_edge_area_change,
    build_coastline_projector,
    phase_h_optimize,
)

REPO = Path(__file__).resolve().parent.parent
INPUT = REPO / "outputs" / "56b_phase_h_min5_optimized.14"
COASTLINE = (
    REPO / "data" / "coastline" / "tokyo_bay"
    / "MLIT_C23" / "C23-06_TOKYOBAY.shp"
)
OUT_DIR = REPO / "outputs"
OUT_F_ONLY = OUT_DIR / "57_phase_h_pass_f_only.14"
OUT_ABEF = OUT_DIR / "57_phase_h_abef.14"
SUMMARY_TXT = OUT_DIR / "57_summary.txt"
SUMMARY_JSON = OUT_DIR / "57_summary.json"

ALPHA_TARGET = 0.95
MIN_ANGLE_TARGET = 30.0
MAX_ANGLE_TARGET = 130.0
AREA_RATIO_TARGET = 0.5
MAX_VALENCE = 8


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


def _metrics(mesh):
    m = min_interior_angle(mesh)
    M = _max_interior_angle(mesh)
    _u, _p, ac = _per_edge_area_change(mesh.nodes, mesh.elements)
    val = node_valence(mesh.elements, mesh.n_nodes)
    return {
        "NP": int(mesh.n_nodes),
        "NE": int(mesh.n_elements),
        "C1": int((m < MIN_ANGLE_TARGET).sum()),
        "C2": int((M > MAX_ANGLE_TARGET).sum()),
        "C4": int((ac > AREA_RATIO_TARGET).sum()),
        "C5": int((val > MAX_VALENCE).sum()),
    }


def _print_metrics(label: str, m: dict) -> None:
    print(
        f"[57] {label}: NP={m['NP']:,} NE={m['NE']:,} "
        f"C1={m['C1']} C2={m['C2']} C4={m['C4']} C5={m['C5']} "
        f"(total {m['C1'] + m['C2'] + m['C4'] + m['C5']})",
        flush=True,
    )


def main() -> int:
    if not INPUT.exists():
        raise SystemExit(f"input missing: {INPUT}")
    OUT_DIR.mkdir(parents=True, exist_ok=True)

    mesh = read_fort14(INPUT)
    before = _metrics(mesh)
    _print_metrics("input (56b)", before)

    projector = build_coastline_projector(
        [COASTLINE],
        max_snap_distance_m=500.0,
        mean_latitude_deg=float(mesh.nodes[:, 1].mean()),
    )

    # ====================================================================
    # Stage 1: Pass F only (operator_order=()) — isolate the F contribution
    # ====================================================================
    print()
    print("[57] Stage 1: Pass F only on PoC #56b input")
    t0 = time.perf_counter()
    mesh_f, info_f = phase_h_optimize(
        mesh,
        alpha_target=ALPHA_TARGET,
        min_angle_target=MIN_ANGLE_TARGET,
        max_angle_target=MAX_ANGLE_TARGET,
        max_outer_rounds=5,
        operator_order=(),
        coastline_projector=projector,
        pass_f_enabled=True,
        pass_f_area_ratio_target=AREA_RATIO_TARGET,
        max_pass_f_sweeps_per_round=100,
    )
    wall_f = time.perf_counter() - t0
    after_f = _metrics(mesh_f)
    _print_metrics(f"Stage 1 (Pass F only, wall {wall_f:.1f} s)", after_f)
    print(
        f"     Pass F accepts={info_f['pass_f_accepts']}, "
        f"sweeps={info_f['pass_f_sweeps']}, "
        f"outer_rounds={info_f['n_outer_rounds']}",
        flush=True,
    )
    write_fort14(mesh_f, OUT_F_ONLY)

    # ====================================================================
    # Stage 2: A+B+E+F full Phase H on PoC #56b input
    # ====================================================================
    print()
    print("[57] Stage 2: A+B+E+F full Phase H on PoC #56b input")
    t0 = time.perf_counter()
    mesh_abef, info_abef = phase_h_optimize(
        mesh,
        alpha_target=ALPHA_TARGET,
        min_angle_target=MIN_ANGLE_TARGET,
        max_angle_target=MAX_ANGLE_TARGET,
        max_smooth_sweeps=200,
        max_topology_per_round=10_000,
        max_outer_rounds=10,
        coastline_projector=projector,
        pass_e_enabled=True,
        pass_e_area_ratio_target=AREA_RATIO_TARGET,
        pass_e_max_valence=MAX_VALENCE,
        max_pass_e_splits_per_round=10_000,
        pass_f_enabled=True,
        pass_f_area_ratio_target=AREA_RATIO_TARGET,
        max_pass_f_sweeps_per_round=100,
    )
    wall_abef = time.perf_counter() - t0
    after_abef = _metrics(mesh_abef)
    _print_metrics(f"Stage 2 (A+B+E+F, wall {wall_abef:.1f} s)", after_abef)
    print(
        f"     Pass A accepts={info_abef['operators_applied'].get('smooth_node', 0)} "
        f"Pass B accepts={sum(v for k, v in info_abef['operators_applied'].items() if k != 'smooth_node')} "
        f"Pass E accepts={info_abef['pass_e_accepts']} "
        f"(swap={info_abef['pass_e_swap_accepts']}, "
        f"split={info_abef['pass_e_split_accepts']}) "
        f"Pass F accepts={info_abef['pass_f_accepts']}, "
        f"sweeps={info_abef['pass_f_sweeps']}, "
        f"outer_rounds={info_abef['n_outer_rounds']}",
        flush=True,
    )
    write_fort14(mesh_abef, OUT_ABEF)

    # ====================================================================
    # Summary
    # ====================================================================
    payload = {
        "input": str(INPUT.resolve()),
        "before": before,
        "stage1_pass_f_only": {
            "output": str(OUT_F_ONLY.resolve()),
            "after": after_f,
            "wall_seconds": wall_f,
            "info": {
                k: (dict(v) if hasattr(v, "items") else v)
                for k, v in info_f.items()
                if not k.startswith("_")
            },
        },
        "stage2_abef": {
            "output": str(OUT_ABEF.resolve()),
            "after": after_abef,
            "wall_seconds": wall_abef,
            "info": {
                k: (dict(v) if hasattr(v, "items") else v)
                for k, v in info_abef.items()
                if not k.startswith("_")
            },
        },
    }
    SUMMARY_JSON.write_text(
        json.dumps(payload, indent=2, default=str), encoding="utf-8",
    )

    total_before = before["C1"] + before["C2"] + before["C4"] + before["C5"]
    total_f = after_f["C1"] + after_f["C2"] + after_f["C4"] + after_f["C5"]
    total_abef = (
        after_abef["C1"] + after_abef["C2"]
        + after_abef["C4"] + after_abef["C5"]
    )
    lines = [
        "PoC #57 — Pass F (C4-aware smoothing) on PoC #56b output",
        f"input : {INPUT.name}",
        "",
        f"  {'stage':<28} | {'NP':>6} | {'NE':>6} | "
        f"{'C1':>4} | {'C2':>4} | {'C4':>4} | {'C5':>3} | total",
        "  " + "-" * 82,
        f"  {'PoC #56b (input)':<28} | "
        f"{before['NP']:>6,} | {before['NE']:>6,} | "
        f"{before['C1']:>4} | {before['C2']:>4} | "
        f"{before['C4']:>4} | {before['C5']:>3} | {total_before}",
        f"  {'PoC #57 Stage 1 (F only)':<28} | "
        f"{after_f['NP']:>6,} | {after_f['NE']:>6,} | "
        f"{after_f['C1']:>4} | {after_f['C2']:>4} | "
        f"{after_f['C4']:>4} | {after_f['C5']:>3} | {total_f}",
        f"  {'PoC #57 Stage 2 (A+B+E+F)':<28} | "
        f"{after_abef['NP']:>6,} | {after_abef['NE']:>6,} | "
        f"{after_abef['C1']:>4} | {after_abef['C2']:>4} | "
        f"{after_abef['C4']:>4} | {after_abef['C5']:>3} | {total_abef}",
        "",
        f"  Stage 1 wall: {wall_f:.1f} s",
        f"  Stage 2 wall: {wall_abef:.1f} s",
        "",
        f"  delta vs 56b input:",
        f"    Stage 1: {total_f - total_before:+d} total "
        f"(C4 {after_f['C4'] - before['C4']:+d})",
        f"    Stage 2: {total_abef - total_before:+d} total "
        f"(C4 {after_abef['C4'] - before['C4']:+d})",
    ]
    SUMMARY_TXT.write_text("\n".join(lines) + "\n", encoding="utf-8")
    print()
    print("\n".join(lines))
    print()
    print(f"wrote {OUT_F_ONLY}")
    print(f"wrote {OUT_ABEF}")
    print(f"wrote {SUMMARY_TXT}")
    print(f"wrote {SUMMARY_JSON}")
    return 0


if __name__ == "__main__":
    raise SystemExit(main())
