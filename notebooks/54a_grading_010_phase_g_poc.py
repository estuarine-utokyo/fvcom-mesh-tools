"""PoC #54a: tighten oceanmesh gradation 0.15 -> 0.10 and re-run
   the raw oceanmesh generation + Phase A-G pipeline.

PoC #53 diagnosed that the 90 residual FVCOM violations after
PoC #52b (full Phase H A+B+E) are dominated by:

  - 80 / 80 C4 fails: 74 % barely failing (ac 0.50-0.55, ratio
    ~2.0-2.2:1), 61 % endpoint-on-boundary, 73 isolated small
    clusters (size 1-8)
  - 9 / 10 C1+C2 fails: 100 % boundary-touching

Both groups point at the **build-time mesh-size grading**, not at
the Phase H operator inventory. The raw oceanmesh output already
carries 502 C4 fails (out of ~80k internal edges, 0.6 %) and 281
C1 fails. Phase A-G + Phase H can only reduce these, not eliminate
them; the auto-pipeline ceiling we keep hitting is the build-time
ceiling.

PoC #54a tests the build-time hypothesis by tightening the
``--om-gradation`` argument to ``oceanmesh.enforce_mesh_gradation``
from 0.15 (current default) to 0.10. A tighter gradation forces
the per-element size to change less per edge, which directly
attacks the 2:1 area ratio mode that dominates the C4 fail set.

This first PoC is **diagnostic-only**: it runs only the raw mesh
generation and the existing Phase A-G pipeline (the same chain
PoC #19 -> PoC #33 already validated). The FVCOM compliance
metrics at raw / Phase G stages are reported and compared against
the g=0.15 baseline. If C4 drops materially, PoC #54b will run the
full Phase H I1+I2+E on the new mesh.

Outputs:
   outputs/54a_tokyo_bay_oceanmesh_g010.14
   outputs/54a_pipeline_passing_g010.14
   outputs/54a_grading_010_summary.txt
   outputs/54a_grading_010_summary.json
"""
from __future__ import annotations

import json
import time
from pathlib import Path

import numpy as np

from fvcom_mesh_tools.algorithms.quality import min_interior_angle
from fvcom_mesh_tools.cli.buildmesh import main as buildmesh_main
from fvcom_mesh_tools.cli.meshpipeline import main as pipeline_main
from fvcom_mesh_tools.diagnostics import node_valence
from fvcom_mesh_tools.io import read_fort14
from fvcom_mesh_tools.mesh_clean_phase_h import _per_edge_area_change

REPO = Path(__file__).resolve().parent.parent
DEM = REPO / "data" / "bathymetry" / "tokyo_bay" / "dem_00_01_change.nc"
COASTLINE = (
    REPO / "data" / "coastline" / "tokyo_bay"
    / "MLIT_C23" / "C23-06_TOKYOBAY.shp"
)
RIVERS = REPO / "data" / "rivers" / "tokyo_bay" / "tokyo_bay_rivers.csv"
OUT_DIR = REPO / "outputs"

RAW_OUT = OUT_DIR / "54a_tokyo_bay_oceanmesh_g010.14"
PIPELINE_OUT = OUT_DIR / "54a_pipeline_passing_g010.14"
PIPELINE_SUMMARY = OUT_DIR / "54a_pipeline_summary_g010.json"
SUMMARY_TXT = OUT_DIR / "54a_grading_010_summary.txt"
SUMMARY_JSON = OUT_DIR / "54a_grading_010_summary.json"

# Baselines (g=0.15) for comparison — paths to existing outputs.
BASELINE_RAW = OUT_DIR / "19_tokyo_bay_oceanmesh.14"
BASELINE_PIPELINE = OUT_DIR / "33_pipeline_passing.14"

GRADATION = 0.10  # 0.15 was the baseline; tighter = smoother size transitions
HMIN_M = 200.0
HMAX_M = 5000.0
SLOPE_PARAMETER = 20.0  # match PoC #19
MAX_ITER = 50


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


def _measure_compliance(mesh):
    m = min_interior_angle(mesh)
    M = _max_interior_angle(mesh)
    _u, _p, ac = _per_edge_area_change(mesh.nodes, mesh.elements)
    val = node_valence(mesh.elements, mesh.n_nodes)
    return {
        "n_nodes": int(mesh.n_nodes),
        "n_elements": int(mesh.n_elements),
        "C1_min_ang_lt_30": int((m < 30.0).sum()),
        "C2_max_ang_gt_130": int((M > 130.0).sum()),
        "C4_area_change_gt_0_5": int((ac > 0.5).sum()),
        "C5_valence_gt_8": int((val > 8).sum()),
        "max_valence": int(val.max()),
    }


def main() -> int:
    for p in (DEM, COASTLINE, RIVERS):
        if not p.exists():
            raise SystemExit(f"required input missing: {p}")
    OUT_DIR.mkdir(parents=True, exist_ok=True)

    print(
        f"[54a] Stage 1: raw oceanmesh generation with "
        f"--om-gradation {GRADATION:g} (was 0.15)",
        flush=True,
    )
    t0 = time.perf_counter()
    rc = buildmesh_main([
        str(DEM), str(RAW_OUT),
        "--engine", "oceanmesh",
        "--hmin", str(HMIN_M), "--hmax", str(HMAX_M), "--zmax", "0.0",
        "--interp-method", "linear",
        "--coastline", str(COASTLINE),
        "--om-slope-parameter", str(SLOPE_PARAMETER),
        "--om-gradation", str(GRADATION),
        "--om-max-iter", str(MAX_ITER),
        "--om-seed", "0",
        "--open-merge-coast-gap", "50",
        "--river-inflow-points", str(RIVERS),
        "--river-segment-nodes", "5",
        "--river-ibtype", "21",
        "--river-snap-tol-m", "5000",
        "--quality-pass", "0",
        "--refine-min-angle", "0",
        "--land-ibtype", "20",
        "--perpfix-iters", "1",
    ])
    raw_wall = time.perf_counter() - t0
    print(
        f"[54a] Stage 1 wall: {raw_wall:.1f} s  (exit code {rc})",
        flush=True,
    )
    if rc != 0:
        raise SystemExit(f"buildmesh failed with exit {rc}")

    raw_mesh = read_fort14(RAW_OUT)
    raw_metrics = _measure_compliance(raw_mesh)
    print(
        f"[54a] raw compliance: "
        f"NP={raw_metrics['n_nodes']:,}  NE={raw_metrics['n_elements']:,}  "
        f"C1={raw_metrics['C1_min_ang_lt_30']}  "
        f"C2={raw_metrics['C2_max_ang_gt_130']}  "
        f"C4={raw_metrics['C4_area_change_gt_0_5']}  "
        f"C5={raw_metrics['C5_valence_gt_8']}",
        flush=True,
    )

    print(
        "\n[54a] Stage 2: Phase A-G pipeline "
        "(fmesh-mesh-pipeline) on the new raw mesh",
        flush=True,
    )
    t0 = time.perf_counter()
    rc = pipeline_main([
        str(RAW_OUT), str(PIPELINE_OUT),
        "--bbox", "139.46", "34.99", "140.10", "35.74",
        "--bbox-tol-m", "150",
        "--open-merge-coast-gap", "50",
        "--min-alpha", "0.95",
        "--max-frac-lt-20deg", "0.005",
        "--max-valence", "8",
        "--max-flipped", "0",
        "--max-disjoint-elems", "0",
        "--summary", str(PIPELINE_SUMMARY),
    ])
    pipeline_wall = time.perf_counter() - t0
    print(
        f"[54a] Stage 2 wall: {pipeline_wall:.1f} s  (exit code {rc})",
        flush=True,
    )
    if rc != 0:
        raise SystemExit(f"pipeline failed with exit {rc}")

    phaseg_mesh = read_fort14(PIPELINE_OUT)
    phaseg_metrics = _measure_compliance(phaseg_mesh)
    print(
        f"[54a] Phase G compliance: "
        f"NP={phaseg_metrics['n_nodes']:,}  NE={phaseg_metrics['n_elements']:,}  "
        f"C1={phaseg_metrics['C1_min_ang_lt_30']}  "
        f"C2={phaseg_metrics['C2_max_ang_gt_130']}  "
        f"C4={phaseg_metrics['C4_area_change_gt_0_5']}  "
        f"C5={phaseg_metrics['C5_valence_gt_8']}",
        flush=True,
    )

    # Baseline comparison (g=0.15)
    baseline_raw = _measure_compliance(read_fort14(BASELINE_RAW))
    baseline_phaseg = _measure_compliance(read_fort14(BASELINE_PIPELINE))

    payload = {
        "config": {
            "gradation": GRADATION,
            "baseline_gradation": 0.15,
            "hmin_m": HMIN_M, "hmax_m": HMAX_M,
            "slope_parameter": SLOPE_PARAMETER,
            "max_iter": MAX_ITER,
            "coastline": str(COASTLINE.resolve()),
            "dem": str(DEM.resolve()),
        },
        "stage_1_raw": {
            "path": str(RAW_OUT.resolve()),
            "wall_seconds": raw_wall,
            "metrics": raw_metrics,
        },
        "stage_2_phase_g": {
            "path": str(PIPELINE_OUT.resolve()),
            "wall_seconds": pipeline_wall,
            "metrics": phaseg_metrics,
        },
        "baseline_g015": {
            "raw_path": str(BASELINE_RAW.resolve()),
            "raw_metrics": baseline_raw,
            "phase_g_path": str(BASELINE_PIPELINE.resolve()),
            "phase_g_metrics": baseline_phaseg,
        },
    }
    SUMMARY_JSON.write_text(
        json.dumps(payload, indent=2, default=str), encoding="utf-8",
    )

    lines = [
        f"PoC #54a — oceanmesh gradation 0.15 -> {GRADATION:g}",
        f"DEM:       {DEM.name}",
        f"coastline: {COASTLINE.name}",
        "",
        "Comparison table (g=0.15 baseline vs g=0.10):",
        "",
        f"  {'stage':<18} | "
        f"{'NP':>8} | {'NE':>8} | "
        f"{'C1':>5} | {'C2':>5} | {'C4':>5} | {'C5':>3}",
        "  " + "-" * 78,
        f"  {'raw  g=0.15':<18} | "
        f"{baseline_raw['n_nodes']:>8,} | {baseline_raw['n_elements']:>8,} | "
        f"{baseline_raw['C1_min_ang_lt_30']:>5} | "
        f"{baseline_raw['C2_max_ang_gt_130']:>5} | "
        f"{baseline_raw['C4_area_change_gt_0_5']:>5} | "
        f"{baseline_raw['C5_valence_gt_8']:>3}",
        f"  {'raw  g=' + f'{GRADATION:.2f}':<18} | "
        f"{raw_metrics['n_nodes']:>8,} | {raw_metrics['n_elements']:>8,} | "
        f"{raw_metrics['C1_min_ang_lt_30']:>5} | "
        f"{raw_metrics['C2_max_ang_gt_130']:>5} | "
        f"{raw_metrics['C4_area_change_gt_0_5']:>5} | "
        f"{raw_metrics['C5_valence_gt_8']:>3}",
        f"  {'Phase G g=0.15':<18} | "
        f"{baseline_phaseg['n_nodes']:>8,} | "
        f"{baseline_phaseg['n_elements']:>8,} | "
        f"{baseline_phaseg['C1_min_ang_lt_30']:>5} | "
        f"{baseline_phaseg['C2_max_ang_gt_130']:>5} | "
        f"{baseline_phaseg['C4_area_change_gt_0_5']:>5} | "
        f"{baseline_phaseg['C5_valence_gt_8']:>3}",
        f"  {'Phase G g=' + f'{GRADATION:.2f}':<18} | "
        f"{phaseg_metrics['n_nodes']:>8,} | "
        f"{phaseg_metrics['n_elements']:>8,} | "
        f"{phaseg_metrics['C1_min_ang_lt_30']:>5} | "
        f"{phaseg_metrics['C2_max_ang_gt_130']:>5} | "
        f"{phaseg_metrics['C4_area_change_gt_0_5']:>5} | "
        f"{phaseg_metrics['C5_valence_gt_8']:>3}",
        "",
        f"raw walltime:     {raw_wall:.1f} s",
        f"pipeline walltime: {pipeline_wall:.1f} s",
        "",
    ]
    # Deltas
    raw_c4_delta = (
        raw_metrics["C4_area_change_gt_0_5"]
        - baseline_raw["C4_area_change_gt_0_5"]
    )
    phaseg_c4_delta = (
        phaseg_metrics["C4_area_change_gt_0_5"]
        - baseline_phaseg["C4_area_change_gt_0_5"]
    )
    lines.append(f"  C4 delta at raw:     {raw_c4_delta:+}")
    lines.append(f"  C4 delta at Phase G: {phaseg_c4_delta:+}")

    SUMMARY_TXT.write_text("\n".join(lines) + "\n", encoding="utf-8")
    print()
    print("\n".join(lines))
    print()
    print(f"wrote {RAW_OUT}")
    print(f"wrote {PIPELINE_OUT}")
    print(f"wrote {SUMMARY_TXT}")
    print(f"wrote {SUMMARY_JSON}")
    return 0


if __name__ == "__main__":
    raise SystemExit(main())
