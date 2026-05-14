"""PoC #54c: g=0.10 raw + Phase A-G with thin_chain_mode='none' +
   Phase H A+B+E.

PoC #54b's step-by-step diagnostic identified Phase C
(``repair_thin_chains`` in ``widen`` mode) as the operator that
re-creates C1/C2/C4 fails on the tighter g=0.10 raw mesh:

  step                  NP       NE      C1   C2   C4
  ---------------------------------------------------------
  raw g=0.10        50,922   89,112   240   19  395
  + A keep_components  same    same   240   19  395  (no-op)
  + B trim_dead_ends 49,476   87,681   184   12  366  (-29 C4 ✓)
  + C repair_thin    49,697   88,123   656  144  427  (+61 C4 ✗)
  + D-F-G                              653  144  426  (recovers some)

Phase B alone reduces C4 by 29. Phase C then *adds* 61 C4 fails by
splitting / widening thin channels into many sliver triangles. The
g=0.15 baseline didn't show this issue because g=0.15 produces
fewer narrow channels to begin with.

PoC #54c bypasses Phase C by calling ``clean_mesh()`` directly with
``thin_chain_mode='none'``, then runs Phase H A+B+E (PoC #52b
configuration) on the result.

Hypothesis: starting Phase H from a ~366 C4 input (vs the previous
~353 baseline) but without the +61 thin-chain induced fails should
let Phase E drive C4 well below the previous 80 ceiling.

Wall budget: Phase A-G ~ 30 s; Phase H ~ 60-90 min on the larger
~83 k-element mesh. Elapse 2 h.

Outputs:
    outputs/54c_phase_g_no_thin_chain.14
    outputs/54c_phase_h_optimized.14
    outputs/54c_summary.txt
    outputs/54c_summary.json
"""
from __future__ import annotations

import json
import time
from pathlib import Path

import numpy as np

from fvcom_mesh_tools.algorithms.quality import min_interior_angle
from fvcom_mesh_tools.diagnostics import node_valence
from fvcom_mesh_tools.io import read_fort14, write_fort14
from fvcom_mesh_tools.mesh_clean import (
    DEFAULT_SKEWED_MAX_ANGLE_DEG,
    DEFAULT_SKEWED_MIN_ANGLE_DEG,
    DEFAULT_SMOOTH_LAPLACIAN_ITERS,
    DEFAULT_SMOOTH_LAPLACIAN_TOL,
    clean_mesh,
)
from fvcom_mesh_tools.mesh_clean_phase_h import (
    DEFAULT_MAX_VALENCE,
    _per_edge_area_change,
    build_coastline_projector,
    phase_h_optimize,
)

REPO = Path(__file__).resolve().parent.parent
RAW_INPUT = REPO / "outputs" / "54a_tokyo_bay_oceanmesh_g010.14"
COASTLINE = (
    REPO / "data" / "coastline" / "tokyo_bay"
    / "MLIT_C23" / "C23-06_TOKYOBAY.shp"
)
OUT_DIR = REPO / "outputs"
PHASE_G_OUT = OUT_DIR / "54c_phase_g_no_thin_chain.14"
PHASE_H_OUT = OUT_DIR / "54c_phase_h_optimized.14"
SUMMARY_TXT = OUT_DIR / "54c_summary.txt"
SUMMARY_JSON = OUT_DIR / "54c_summary.json"

# Pipeline rung-equivalent kwargs (A+B+D+F+G; thin_chain OFF).
BBOX = (139.46, 34.99, 140.10, 35.74)
BBOX_TOL_M = 150
LAND_IBTYPE = 20
OPEN_MERGE_COAST_GAP = 50

# Phase H I1+I2+E target gates (PoC #52b configuration).
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


def _metrics(mesh) -> dict:
    m = min_interior_angle(mesh)
    M = _max_interior_angle(mesh)
    _u, _p, ac = _per_edge_area_change(mesh.nodes, mesh.elements)
    val = node_valence(mesh.elements, mesh.n_nodes)
    return {
        "NP": int(mesh.n_nodes),
        "NE": int(mesh.n_elements),
        "C1": int((m < 30.0).sum()),
        "C2": int((M > 130.0).sum()),
        "C4": int((ac > 0.5).sum()),
        "C5": int((val > 8).sum()),
        "max_valence": int(val.max()),
    }


def main() -> int:
    if not RAW_INPUT.exists():
        raise SystemExit(f"input missing: {RAW_INPUT}")
    if not COASTLINE.exists():
        raise SystemExit(f"coastline missing: {COASTLINE}")
    OUT_DIR.mkdir(parents=True, exist_ok=True)

    mesh_raw = read_fort14(RAW_INPUT)
    raw_metrics = _metrics(mesh_raw)
    print(
        f"[54c] raw input: NP={raw_metrics['NP']:,}  NE={raw_metrics['NE']:,}  "
        f"C1={raw_metrics['C1']}  C2={raw_metrics['C2']}  "
        f"C4={raw_metrics['C4']}  C5={raw_metrics['C5']}",
        flush=True,
    )

    print(
        "[54c] Stage 1: Phase A-G (thin_chain_mode='none') via clean_mesh",
        flush=True,
    )
    t0 = time.perf_counter()
    phase_g_mesh, phase_g_info = clean_mesh(
        mesh_raw,
        bbox=BBOX,
        bbox_tol_m=BBOX_TOL_M,
        land_ibtype=LAND_IBTYPE,
        open_merge_coast_gap=OPEN_MERGE_COAST_GAP,
        remove_disjoint=True,
        min_component_elements=0,
        require_open_boundary=False,
        trim_dead_ends_iters=10,
        thin_chain_mode="none",  # <-- the key change vs PoC #54a
        repair_overconnected_iters=20,
        max_nbr_elem=8,
        overconn_min_angle_floor_deg=0.0,
        repair_skewed=True,
        repair_skewed_min_angle_deg=DEFAULT_SKEWED_MIN_ANGLE_DEG,
        repair_skewed_max_angle_deg=DEFAULT_SKEWED_MAX_ANGLE_DEG,
        smooth_laplacian=True,
        smooth_laplacian_iters=DEFAULT_SMOOTH_LAPLACIAN_ITERS,
        smooth_laplacian_tol=DEFAULT_SMOOTH_LAPLACIAN_TOL,
        smooth_repair_flipped=True,
        smooth_max_repair_passes=10,
    )
    phase_g_wall = time.perf_counter() - t0
    write_fort14(phase_g_mesh, PHASE_G_OUT)
    phase_g_metrics = _metrics(phase_g_mesh)
    print(
        f"[54c] Phase G done in {phase_g_wall:.1f} s  "
        f"NP={phase_g_metrics['NP']:,}  NE={phase_g_metrics['NE']:,}  "
        f"C1={phase_g_metrics['C1']}  C2={phase_g_metrics['C2']}  "
        f"C4={phase_g_metrics['C4']}  C5={phase_g_metrics['C5']}",
        flush=True,
    )

    print(
        f"[54c] Stage 2: build coastline projector ({COASTLINE.name})",
        flush=True,
    )
    t0 = time.perf_counter()
    projector = build_coastline_projector(
        [COASTLINE],
        max_snap_distance_m=MAX_SNAP_M,
        mean_latitude_deg=float(phase_g_mesh.nodes[:, 1].mean()),
    )
    proj_build_wall = time.perf_counter() - t0
    if projector is None:
        raise SystemExit("coastline projector built no polylines")
    print(f"[54c] projector built in {proj_build_wall:.1f} s", flush=True)

    print(
        "[54c] Stage 3: Phase H A+B+E "
        f"(alpha>={ALPHA_TARGET}, min_ang>={MIN_ANGLE_TARGET}°, "
        f"max_ang<={MAX_ANGLE_TARGET}°, area_change<={AREA_RATIO_TARGET}, "
        f"max_valence<={MAX_VALENCE})",
        flush=True,
    )
    t0 = time.perf_counter()
    phase_h_mesh, h_info = phase_h_optimize(
        phase_g_mesh,
        alpha_target=ALPHA_TARGET,
        min_angle_target=MIN_ANGLE_TARGET,
        max_angle_target=MAX_ANGLE_TARGET,
        max_smooth_sweeps=200,
        max_topology_per_round=10_000,
        max_outer_rounds=MAX_OUTER_ROUNDS,
        coastline_projector=projector,
        lookahead_enabled=False,
        patch_recdt_enabled=False,
        pass_e_enabled=True,
        pass_e_area_ratio_target=AREA_RATIO_TARGET,
        pass_e_max_valence=MAX_VALENCE,
        max_pass_e_splits_per_round=10_000,
    )
    phase_h_wall = time.perf_counter() - t0
    write_fort14(phase_h_mesh, PHASE_H_OUT)
    phase_h_metrics = _metrics(phase_h_mesh)
    print(
        f"[54c] Phase H done in {phase_h_wall:.1f} s  "
        f"rounds={h_info['n_outer_rounds']:,}  "
        f"smooths={h_info['n_smooth_sweeps']:,}  "
        f"iters={h_info['n_iters']:,}  "
        f"pass_e_acc={h_info.get('pass_e_accepts', 0):,} "
        f"(swap={h_info.get('pass_e_swap_accepts', 0):,}, "
        f"split={h_info.get('pass_e_split_accepts', 0):,})  "
        f"pass_e_rej={h_info.get('pass_e_rejected', 0):,}",
        flush=True,
    )
    print(
        f"[54c] final residuals: "
        f"NP={phase_h_metrics['NP']:,}  NE={phase_h_metrics['NE']:,}  "
        f"C1={phase_h_metrics['C1']}  C2={phase_h_metrics['C2']}  "
        f"C4={phase_h_metrics['C4']}  C5={phase_h_metrics['C5']}",
        flush=True,
    )

    payload = {
        "config": {
            "gradation": 0.10,
            "thin_chain_mode": "none",
            "alpha_target": ALPHA_TARGET,
            "min_angle_target": MIN_ANGLE_TARGET,
            "max_angle_target": MAX_ANGLE_TARGET,
            "area_ratio_target": AREA_RATIO_TARGET,
            "max_valence": MAX_VALENCE,
            "max_outer_rounds": MAX_OUTER_ROUNDS,
            "passes_enabled": ["A", "B", "D", "F", "G", "Phase_H_A+B+E"],
        },
        "raw_input": {
            "path": str(RAW_INPUT.resolve()),
            "metrics": raw_metrics,
        },
        "phase_g_output": {
            "path": str(PHASE_G_OUT.resolve()),
            "wall_seconds": phase_g_wall,
            "metrics": phase_g_metrics,
        },
        "phase_h_output": {
            "path": str(PHASE_H_OUT.resolve()),
            "wall_seconds": phase_h_wall,
            "metrics": phase_h_metrics,
            "phase_h_info_summary": {
                "n_outer_rounds": int(h_info["n_outer_rounds"]),
                "n_iters": int(h_info["n_iters"]),
                "n_smooth_sweeps": int(h_info["n_smooth_sweeps"]),
                "pass_e_accepts": int(h_info.get("pass_e_accepts", 0)),
                "pass_e_swap_accepts": int(
                    h_info.get("pass_e_swap_accepts", 0),
                ),
                "pass_e_split_accepts": int(
                    h_info.get("pass_e_split_accepts", 0),
                ),
                "pass_e_rejected": int(h_info.get("pass_e_rejected", 0)),
            },
        },
    }
    SUMMARY_JSON.write_text(
        json.dumps(payload, indent=2, default=str), encoding="utf-8",
    )

    lines = [
        "PoC #54c — g=0.10 + thin_chain=none + Phase H A+B+E",
        f"raw input: {RAW_INPUT.name}",
        f"coastline: {COASTLINE.name}",
        "",
        f"  {'stage':<18} | {'NP':>6} | {'NE':>6} | "
        f"{'C1':>4} | {'C2':>4} | {'C4':>4} | {'C5':>3}",
        "  " + "-" * 76,
        f"  {'raw g=0.10':<18} | "
        f"{raw_metrics['NP']:>6,} | {raw_metrics['NE']:>6,} | "
        f"{raw_metrics['C1']:>4} | {raw_metrics['C2']:>4} | "
        f"{raw_metrics['C4']:>4} | {raw_metrics['C5']:>3}",
        f"  {'Phase G (C=off)':<18} | "
        f"{phase_g_metrics['NP']:>6,} | {phase_g_metrics['NE']:>6,} | "
        f"{phase_g_metrics['C1']:>4} | {phase_g_metrics['C2']:>4} | "
        f"{phase_g_metrics['C4']:>4} | {phase_g_metrics['C5']:>3}",
        f"  {'Phase H A+B+E':<18} | "
        f"{phase_h_metrics['NP']:>6,} | {phase_h_metrics['NE']:>6,} | "
        f"{phase_h_metrics['C1']:>4} | {phase_h_metrics['C2']:>4} | "
        f"{phase_h_metrics['C4']:>4} | {phase_h_metrics['C5']:>3}",
        "",
        f"  Phase G wall  : {phase_g_wall:.1f} s",
        f"  Phase H wall  : {phase_h_wall:.1f} s",
        f"  Phase H rounds: {h_info['n_outer_rounds']}",
        f"  Pass E acc    : {h_info.get('pass_e_accepts', 0)} "
        f"(swap={h_info.get('pass_e_swap_accepts', 0)}, "
        f"split={h_info.get('pass_e_split_accepts', 0)})",
        f"  Pass E rej    : {h_info.get('pass_e_rejected', 0)}",
        "",
        "  vs PoC #52b baseline (g=0.15 + thin_chain=widen):",
        "    final: C1=9  C2=1  C4=80  C5=0",
    ]
    SUMMARY_TXT.write_text("\n".join(lines) + "\n", encoding="utf-8")
    print()
    print("\n".join(lines))
    print()
    print(f"wrote {PHASE_G_OUT}")
    print(f"wrote {PHASE_H_OUT}")
    print(f"wrote {SUMMARY_TXT}")
    print(f"wrote {SUMMARY_JSON}")
    return 0


if __name__ == "__main__":
    raise SystemExit(main())
