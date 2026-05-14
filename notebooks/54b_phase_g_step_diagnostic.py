"""PoC #54b: Phase A-G rung-by-rung C4 (and C1/C2/C5) tracking.

PoC #54a discovered that tightening ``--om-gradation`` 0.15 -> 0.10
reduces raw C4 (502 -> 395) but the Phase G output's C4 goes
UP (+9) compared to the g=0.15 baseline. The reduction at raw is
*cancelled inside Phase A-G*. This PoC isolates which of the six
phases (A keep_components, B trim_dead_ends, C repair_thin_chains,
D repair_overconnected_nodes, F repair_skewed_elements, G
smooth_mesh_laplacian) is responsible for re-creating C4 fails on
the tighter-grading input.

Method: load the raw g=0.10 mesh from PoC #54a, apply each phase in
the order the pipeline applies them, and measure ``C1 / C2 / C4 /
C5 / NP / NE`` after every step. Per-step delta isolates the
culprit.

This is read-only diagnostic so it runs on the login node in under
a minute.

Outputs:
   outputs/54b_phase_g_step_diagnostic.txt
   outputs/54b_phase_g_step_diagnostic.json
"""
from __future__ import annotations

import json
from pathlib import Path

import numpy as np

from fvcom_mesh_tools.algorithms.quality import min_interior_angle
from fvcom_mesh_tools.diagnostics import DEFAULT_MAX_NBR_ELEM, node_valence
from fvcom_mesh_tools.io import read_fort14
from fvcom_mesh_tools.mesh_clean import (
    DEFAULT_SKEWED_MAX_ANGLE_DEG,
    DEFAULT_SKEWED_MIN_ANGLE_DEG,
    DEFAULT_SMOOTH_LAPLACIAN_ITERS,
    DEFAULT_SMOOTH_LAPLACIAN_TOL,
    _deg_per_metre,
    keep_components,
    rebuild_boundaries,
    repair_overconnected_nodes,
    repair_skewed_elements,
    repair_thin_chains,
    smooth_mesh_laplacian,
    trim_dead_ends,
)
from fvcom_mesh_tools.mesh_clean_phase_h import _per_edge_area_change

REPO = Path(__file__).resolve().parent.parent
INPUT = REPO / "outputs" / "54a_tokyo_bay_oceanmesh_g010.14"
OUT_DIR = REPO / "outputs"
SUMMARY_TXT = OUT_DIR / "54b_phase_g_step_diagnostic.txt"
SUMMARY_JSON = OUT_DIR / "54b_phase_g_step_diagnostic.json"

# Pipeline rung settings used by PoC #54a (defaults of fmesh-mesh-pipeline).
BBOX = (139.46, 34.99, 140.10, 35.74)
BBOX_TOL_M = 150
LAND_IBTYPE = 20
OPEN_MERGE_COAST_GAP = 50

TRIM_DEAD_ENDS_ITERS = 10  # clean_mesh default
THIN_CHAIN_MODE = "widen"  # default for repair_thin_chains usage in pipeline
MIN_THIN_CHAIN = 4  # commonly used default
REPAIR_OVERCONNECTED_ITERS = 20  # fmesh-mesh-pipeline default
MAX_NBR_ELEM = DEFAULT_MAX_NBR_ELEM
OVERCONN_MIN_ANGLE_FLOOR = 0.0
SKEW_MIN_ANGLE = DEFAULT_SKEWED_MIN_ANGLE_DEG
SKEW_MAX_ANGLE = DEFAULT_SKEWED_MAX_ANGLE_DEG
SMOOTH_ITERS = DEFAULT_SMOOTH_LAPLACIAN_ITERS
SMOOTH_TOL = DEFAULT_SMOOTH_LAPLACIAN_TOL


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


def _tol_deg(mesh) -> float:
    if mesh.n_nodes == 0:
        return 0.0
    return BBOX_TOL_M * _deg_per_metre(float(mesh.nodes[:, 1].mean()))


def main() -> int:
    if not INPUT.exists():
        raise SystemExit(f"input missing: {INPUT}")

    mesh = read_fort14(INPUT)
    history = []
    history.append({"step": "0_raw_g010", **_metrics(mesh)})
    print(f"[54b] step 0  raw g=0.10            "
          f"NP={history[-1]['NP']:>6,}  NE={history[-1]['NE']:>6,}  "
          f"C1={history[-1]['C1']:>4}  C2={history[-1]['C2']:>4}  "
          f"C4={history[-1]['C4']:>4}  C5={history[-1]['C5']:>3}",
          flush=True)

    # Step A: keep_components
    cur, _ = keep_components(
        mesh, min_elements=0, require_open_boundary=False,
    )
    cur = rebuild_boundaries(
        cur, bbox=BBOX, tol_deg=_tol_deg(cur),
        land_ibtype=LAND_IBTYPE,
        open_merge_coast_gap=OPEN_MERGE_COAST_GAP,
    )
    history.append({"step": "1_A_keep_components", **_metrics(cur)})
    print(f"[54b] step 1  A keep_components     "
          f"NP={history[-1]['NP']:>6,}  NE={history[-1]['NE']:>6,}  "
          f"C1={history[-1]['C1']:>4}  C2={history[-1]['C2']:>4}  "
          f"C4={history[-1]['C4']:>4}  C5={history[-1]['C5']:>3}",
          flush=True)

    # Step B: trim_dead_ends
    cur, _ = trim_dead_ends(
        cur,
        max_iters=TRIM_DEAD_ENDS_ITERS,
        bbox=BBOX, tol_deg=_tol_deg(cur),
        land_ibtype=LAND_IBTYPE,
        open_merge_coast_gap=OPEN_MERGE_COAST_GAP,
    )
    history.append({"step": "2_B_trim_dead_ends", **_metrics(cur)})
    print(f"[54b] step 2  B trim_dead_ends      "
          f"NP={history[-1]['NP']:>6,}  NE={history[-1]['NE']:>6,}  "
          f"C1={history[-1]['C1']:>4}  C2={history[-1]['C2']:>4}  "
          f"C4={history[-1]['C4']:>4}  C5={history[-1]['C5']:>3}",
          flush=True)

    # Step C: repair_thin_chains
    cur, _ = repair_thin_chains(
        cur,
        mode=THIN_CHAIN_MODE,
        min_chain_length=MIN_THIN_CHAIN,
        bbox=BBOX, tol_deg=_tol_deg(cur),
        land_ibtype=LAND_IBTYPE,
        open_merge_coast_gap=OPEN_MERGE_COAST_GAP,
    )
    history.append({"step": "3_C_repair_thin_chains", **_metrics(cur)})
    print(f"[54b] step 3  C repair_thin_chains  "
          f"NP={history[-1]['NP']:>6,}  NE={history[-1]['NE']:>6,}  "
          f"C1={history[-1]['C1']:>4}  C2={history[-1]['C2']:>4}  "
          f"C4={history[-1]['C4']:>4}  C5={history[-1]['C5']:>3}",
          flush=True)

    # Step D: repair_overconnected_nodes
    cur, _ = repair_overconnected_nodes(
        cur,
        max_nbr_elem=MAX_NBR_ELEM,
        max_iters=REPAIR_OVERCONNECTED_ITERS,
        min_angle_floor_deg=OVERCONN_MIN_ANGLE_FLOOR,
    )
    history.append(
        {"step": "4_D_repair_overconnected_nodes", **_metrics(cur)},
    )
    print(f"[54b] step 4  D repair_overconn     "
          f"NP={history[-1]['NP']:>6,}  NE={history[-1]['NE']:>6,}  "
          f"C1={history[-1]['C1']:>4}  C2={history[-1]['C2']:>4}  "
          f"C4={history[-1]['C4']:>4}  C5={history[-1]['C5']:>3}",
          flush=True)

    # Step F: repair_skewed_elements
    cur, _ = repair_skewed_elements(
        cur,
        min_angle_deg=SKEW_MIN_ANGLE,
        max_angle_deg=SKEW_MAX_ANGLE,
        bbox=BBOX, tol_deg=_tol_deg(cur),
        land_ibtype=LAND_IBTYPE,
        open_merge_coast_gap=OPEN_MERGE_COAST_GAP,
    )
    history.append({"step": "5_F_repair_skewed", **_metrics(cur)})
    print(f"[54b] step 5  F repair_skewed       "
          f"NP={history[-1]['NP']:>6,}  NE={history[-1]['NE']:>6,}  "
          f"C1={history[-1]['C1']:>4}  C2={history[-1]['C2']:>4}  "
          f"C4={history[-1]['C4']:>4}  C5={history[-1]['C5']:>3}",
          flush=True)

    # Step G: smooth_mesh_laplacian
    cur, _ = smooth_mesh_laplacian(
        cur,
        max_iter=SMOOTH_ITERS,
        tol=SMOOTH_TOL,
        repair_flipped=True,
        max_repair_passes=10,
    )
    history.append({"step": "6_G_smooth_laplacian", **_metrics(cur)})
    print(f"[54b] step 6  G smooth_laplacian    "
          f"NP={history[-1]['NP']:>6,}  NE={history[-1]['NE']:>6,}  "
          f"C1={history[-1]['C1']:>4}  C2={history[-1]['C2']:>4}  "
          f"C4={history[-1]['C4']:>4}  C5={history[-1]['C5']:>3}",
          flush=True)

    # Per-step deltas
    deltas = []
    for i in range(1, len(history)):
        d = {
            "step": history[i]["step"],
            "dNP": history[i]["NP"] - history[i - 1]["NP"],
            "dNE": history[i]["NE"] - history[i - 1]["NE"],
            "dC1": history[i]["C1"] - history[i - 1]["C1"],
            "dC2": history[i]["C2"] - history[i - 1]["C2"],
            "dC4": history[i]["C4"] - history[i - 1]["C4"],
            "dC5": history[i]["C5"] - history[i - 1]["C5"],
        }
        deltas.append(d)

    payload = {
        "input": str(INPUT.resolve()),
        "history": history,
        "deltas": deltas,
    }
    SUMMARY_JSON.write_text(
        json.dumps(payload, indent=2, default=str), encoding="utf-8",
    )

    lines = [
        "PoC #54b — Phase A-G rung-by-rung C4 tracking",
        f"input: {INPUT.name} (raw oceanmesh g=0.10)",
        "",
        f"  {'step':<32} | {'NP':>6} | {'NE':>6} | "
        f"{'C1':>4} | {'C2':>4} | {'C4':>4} | {'C5':>3}",
        "  " + "-" * 88,
    ]
    for h in history:
        lines.append(
            f"  {h['step']:<32} | "
            f"{h['NP']:>6,} | {h['NE']:>6,} | "
            f"{h['C1']:>4} | {h['C2']:>4} | {h['C4']:>4} | "
            f"{h['C5']:>3}"
        )
    lines.append("")
    lines.append("  Per-step deltas (vs previous step):")
    lines.append(
        f"  {'step':<32} | {'dNP':>6} | {'dNE':>6} | "
        f"{'dC1':>4} | {'dC2':>4} | {'dC4':>4} | {'dC5':>3}"
    )
    lines.append("  " + "-" * 88)
    for d in deltas:
        lines.append(
            f"  {d['step']:<32} | "
            f"{d['dNP']:>+6} | {d['dNE']:>+6} | "
            f"{d['dC1']:>+4} | {d['dC2']:>+4} | {d['dC4']:>+4} | "
            f"{d['dC5']:>+3}"
        )
    SUMMARY_TXT.write_text("\n".join(lines) + "\n", encoding="utf-8")
    print()
    print("\n".join(lines))
    print()
    print(f"wrote {SUMMARY_TXT}")
    print(f"wrote {SUMMARY_JSON}")
    return 0


if __name__ == "__main__":
    raise SystemExit(main())
