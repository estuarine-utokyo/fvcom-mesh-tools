"""PoC #55: coastline simplification before mesh generation.

PoC #54d structural analysis of the 77 residual violations after
PoC #54c (g=0.10 + thin_chain=none + Phase H A+B+E):

  - C1 8 fails: 100 % isolated size-1 boundary triangles
  - C4 69 fails: 92 % cluster boundary-touching, 75 % boundary
    endpoint, ratio extremes already cleaned (worst 3.17:1)

The remaining failures are no longer a mesh-quality problem; they
are a *coastline-polyline* problem. The MLIT_C23 shapefile has
65,590 vertices across 1,325 LineStrings — many sharp inflection
points sub-100 m that force skewed boundary triangles regardless
of the mesh-size strategy.

PoC #55 preprocesses the coastline with ``shapely.simplify
(tolerance)`` using the Douglas-Peucker algorithm. Vertices that
sit within ``tolerance`` of the simplified line are removed,
preserving the overall shape while smoothing out micro-detail.
With ``hmin = 200 m`` already in oceanmesh, features below ~100 m
are not represented in the mesh anyway, so simplifying at this
scale removes geometry the mesh cannot honour.

Tolerance: 1e-3° ≈ 111 m at lat 35°. Conservative — keeps every
feature that the mesh can actually represent.

Pipeline (mirrors PoC #54c with the simplified coastline):
  1. Simplify MLIT_C23 -> outputs/55_coastline_simplified.shp
  2. buildmesh (oceanmesh, g=0.10, simplified coastline)
  3. clean_mesh (thin_chain_mode='none', A+B+D+F+G)
  4. Phase H A+B+E
  5. Compare residuals against PoC #54c

Outputs:
   outputs/55_coastline_simplified.shp   (and .shx / .dbf / .prj)
   outputs/55_tokyo_bay_oceanmesh.14
   outputs/55_phase_g.14
   outputs/55_phase_h_optimized.14
   outputs/55_summary.{txt,json}
"""
from __future__ import annotations

import json
import time
from pathlib import Path

import geopandas as gpd
import numpy as np

from fvcom_mesh_tools.algorithms.quality import min_interior_angle
from fvcom_mesh_tools.cli.buildmesh import main as buildmesh_main
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
DEM = REPO / "data" / "bathymetry" / "tokyo_bay" / "dem_00_01_change.nc"
COASTLINE_RAW = (
    REPO / "data" / "coastline" / "tokyo_bay"
    / "MLIT_C23" / "C23-06_TOKYOBAY.shp"
)
RIVERS = REPO / "data" / "rivers" / "tokyo_bay" / "tokyo_bay_rivers.csv"
OUT_DIR = REPO / "outputs"
COASTLINE_SIMPLIFIED = OUT_DIR / "55_coastline_simplified.shp"
RAW_OUT = OUT_DIR / "55_tokyo_bay_oceanmesh.14"
PHASE_G_OUT = OUT_DIR / "55_phase_g.14"
PHASE_H_OUT = OUT_DIR / "55_phase_h_optimized.14"
SUMMARY_TXT = OUT_DIR / "55_summary.txt"
SUMMARY_JSON = OUT_DIR / "55_summary.json"

# Coastline simplification tolerance. PoC #55-original tried 1e-3°
# (~111 m) but that triggered a shapely TopologyException inside
# oceanmesh's ``_classify_shoreline`` (the simplified coastline
# self-intersected at lat 35.406 / lon 139.694). The fix: drop
# aggression to 1e-4° (~11 m), still well below the mesh-
# representable scale (hmin = 200 m) so sharp sub-11 m corners
# are still smoothed, but multi-decametre islands / inlets keep
# enough vertices to stay topologically valid.
SIMPLIFY_TOLERANCE_DEG = 1.0e-4

# oceanmesh args (match PoC #54c configuration except the coastline path).
HMIN_M = 200.0
HMAX_M = 5000.0
SLOPE_PARAMETER = 20.0
GRADATION = 0.10
MAX_ITER = 50

# Phase A-G config (PoC #54c style)
BBOX = (139.46, 34.99, 140.10, 35.74)
BBOX_TOL_M = 150
LAND_IBTYPE = 20
OPEN_MERGE_COAST_GAP = 50

# Phase H A+B+E config (PoC #52b / #54c style)
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


def _vertex_count(gdf: gpd.GeoDataFrame) -> int:
    total = 0
    for geom in gdf.geometry:
        if geom.geom_type == "LineString":
            total += len(geom.coords)
        else:
            for part in geom.geoms:
                total += len(part.coords)
    return total


def main() -> int:
    for p in (DEM, COASTLINE_RAW, RIVERS):
        if not p.exists():
            raise SystemExit(f"required input missing: {p}")
    OUT_DIR.mkdir(parents=True, exist_ok=True)

    # Stage 0: coastline simplification
    print(
        f"[55] Stage 0: simplify coastline "
        f"(tolerance={SIMPLIFY_TOLERANCE_DEG:.0e} deg "
        f"≈ {SIMPLIFY_TOLERANCE_DEG * 111_000:.0f} m at lat 35)",
        flush=True,
    )
    t0 = time.perf_counter()
    gdf_raw = gpd.read_file(COASTLINE_RAW)
    if gdf_raw.crs is None:
        # The MLIT shapefile lacks CRS metadata; the data is WGS84.
        gdf_raw = gdf_raw.set_crs(epsg=4326)
    n_lines_raw = len(gdf_raw)
    n_vertices_raw = _vertex_count(gdf_raw)
    gdf_simp = gdf_raw.copy()
    gdf_simp["geometry"] = gdf_simp.geometry.simplify(
        SIMPLIFY_TOLERANCE_DEG, preserve_topology=True,
    )
    # Drop any geometries that became degenerate (single point, etc.)
    valid = gdf_simp.geometry.is_valid & ~gdf_simp.geometry.is_empty
    gdf_simp = gdf_simp.loc[valid].copy()
    n_lines_simp = len(gdf_simp)
    n_vertices_simp = _vertex_count(gdf_simp)
    gdf_simp.to_file(COASTLINE_SIMPLIFIED)
    simp_wall = time.perf_counter() - t0
    print(
        f"[55] simplification done in {simp_wall:.1f} s — "
        f"{n_lines_raw:,} -> {n_lines_simp:,} LineStrings, "
        f"{n_vertices_raw:,} -> {n_vertices_simp:,} vertices "
        f"({n_vertices_simp / max(n_vertices_raw, 1):.1%} kept)",
        flush=True,
    )

    # Stage 1: raw oceanmesh generation
    print(
        f"[55] Stage 1: raw oceanmesh (g={GRADATION:g}, simplified coastline)",
        flush=True,
    )
    t0 = time.perf_counter()
    rc = buildmesh_main([
        str(DEM), str(RAW_OUT),
        "--engine", "oceanmesh",
        "--hmin", str(HMIN_M), "--hmax", str(HMAX_M), "--zmax", "0.0",
        "--interp-method", "linear",
        "--coastline", str(COASTLINE_SIMPLIFIED),
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
    if rc != 0:
        raise SystemExit(f"buildmesh failed with exit {rc}")
    raw_mesh = read_fort14(RAW_OUT)
    raw_metrics = _metrics(raw_mesh)
    print(
        f"[55] raw built in {raw_wall:.1f} s — "
        f"NP={raw_metrics['NP']:,} NE={raw_metrics['NE']:,} "
        f"C1={raw_metrics['C1']} C2={raw_metrics['C2']} "
        f"C4={raw_metrics['C4']} C5={raw_metrics['C5']}",
        flush=True,
    )

    # Stage 2: Phase A-G (thin_chain off)
    print(
        "[55] Stage 2: clean_mesh (thin_chain_mode='none')",
        flush=True,
    )
    t0 = time.perf_counter()
    phase_g_mesh, _ = clean_mesh(
        raw_mesh,
        bbox=BBOX, bbox_tol_m=BBOX_TOL_M,
        land_ibtype=LAND_IBTYPE,
        open_merge_coast_gap=OPEN_MERGE_COAST_GAP,
        remove_disjoint=True,
        min_component_elements=0,
        require_open_boundary=False,
        trim_dead_ends_iters=10,
        thin_chain_mode="none",
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
        f"[55] Phase G done in {phase_g_wall:.1f} s — "
        f"NP={phase_g_metrics['NP']:,} NE={phase_g_metrics['NE']:,} "
        f"C1={phase_g_metrics['C1']} C2={phase_g_metrics['C2']} "
        f"C4={phase_g_metrics['C4']} C5={phase_g_metrics['C5']}",
        flush=True,
    )

    # Stage 3: Phase H A+B+E
    print(
        f"[55] Stage 3: build coastline projector ({COASTLINE_SIMPLIFIED.name})",
        flush=True,
    )
    t0 = time.perf_counter()
    projector = build_coastline_projector(
        [COASTLINE_SIMPLIFIED],
        max_snap_distance_m=MAX_SNAP_M,
        mean_latitude_deg=float(phase_g_mesh.nodes[:, 1].mean()),
    )
    proj_build_wall = time.perf_counter() - t0
    if projector is None:
        raise SystemExit("coastline projector built no polylines")
    print(f"[55] projector built in {proj_build_wall:.1f} s", flush=True)

    print("[55] Stage 4: Phase H A+B+E", flush=True)
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
        f"[55] Phase H done in {phase_h_wall:.1f} s "
        f"rounds={h_info['n_outer_rounds']} "
        f"pass_e_acc={h_info.get('pass_e_accepts', 0)} "
        f"(swap={h_info.get('pass_e_swap_accepts', 0)}, "
        f"split={h_info.get('pass_e_split_accepts', 0)}) "
        f"pass_e_rej={h_info.get('pass_e_rejected', 0)}",
        flush=True,
    )
    print(
        f"[55] final residuals: "
        f"NP={phase_h_metrics['NP']:,} NE={phase_h_metrics['NE']:,} "
        f"C1={phase_h_metrics['C1']} C2={phase_h_metrics['C2']} "
        f"C4={phase_h_metrics['C4']} C5={phase_h_metrics['C5']}",
        flush=True,
    )

    payload = {
        "config": {
            "simplify_tolerance_deg": SIMPLIFY_TOLERANCE_DEG,
            "gradation": GRADATION,
            "thin_chain_mode": "none",
            "alpha_target": ALPHA_TARGET,
            "min_angle_target": MIN_ANGLE_TARGET,
            "max_angle_target": MAX_ANGLE_TARGET,
            "area_ratio_target": AREA_RATIO_TARGET,
            "max_valence": MAX_VALENCE,
        },
        "coastline_simplification": {
            "raw_path": str(COASTLINE_RAW.resolve()),
            "simplified_path": str(COASTLINE_SIMPLIFIED.resolve()),
            "n_lines_raw": n_lines_raw,
            "n_lines_simplified": n_lines_simp,
            "n_vertices_raw": n_vertices_raw,
            "n_vertices_simplified": n_vertices_simp,
            "wall_seconds": simp_wall,
        },
        "raw_mesh": {
            "path": str(RAW_OUT.resolve()),
            "wall_seconds": raw_wall,
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
        "poc54c_baseline": {
            "C1": 8, "C2": 0, "C4": 69, "C5": 0,
            "NP": 49666, "NE": 87022,
            "note": "PoC #54c — same pipeline with raw coastline",
        },
    }
    SUMMARY_JSON.write_text(
        json.dumps(payload, indent=2, default=str), encoding="utf-8",
    )

    lines = [
        f"PoC #55 — coastline simplification "
        f"(tol={SIMPLIFY_TOLERANCE_DEG:.0e}°) + g=0.10 + thin_chain=none "
        f"+ Phase H A+B+E",
        f"raw coastline: {COASTLINE_RAW.name} "
        f"({n_lines_raw:,} lines, {n_vertices_raw:,} vertices)",
        f"simplified:    {COASTLINE_SIMPLIFIED.name} "
        f"({n_lines_simp:,} lines, {n_vertices_simp:,} vertices, "
        f"{n_vertices_simp / max(n_vertices_raw, 1):.1%} kept)",
        "",
        f"  {'stage':<22} | {'NP':>6} | {'NE':>6} | "
        f"{'C1':>4} | {'C2':>4} | {'C4':>4} | {'C5':>3}",
        "  " + "-" * 80,
        f"  {'raw (g=0.10, simp)':<22} | "
        f"{raw_metrics['NP']:>6,} | {raw_metrics['NE']:>6,} | "
        f"{raw_metrics['C1']:>4} | {raw_metrics['C2']:>4} | "
        f"{raw_metrics['C4']:>4} | {raw_metrics['C5']:>3}",
        f"  {'Phase G (C=off)':<22} | "
        f"{phase_g_metrics['NP']:>6,} | {phase_g_metrics['NE']:>6,} | "
        f"{phase_g_metrics['C1']:>4} | {phase_g_metrics['C2']:>4} | "
        f"{phase_g_metrics['C4']:>4} | {phase_g_metrics['C5']:>3}",
        f"  {'Phase H A+B+E':<22} | "
        f"{phase_h_metrics['NP']:>6,} | {phase_h_metrics['NE']:>6,} | "
        f"{phase_h_metrics['C1']:>4} | {phase_h_metrics['C2']:>4} | "
        f"{phase_h_metrics['C4']:>4} | {phase_h_metrics['C5']:>3}",
        "",
        f"  raw walltime     : {raw_wall:.1f} s",
        f"  Phase G walltime : {phase_g_wall:.1f} s",
        f"  Phase H walltime : {phase_h_wall:.1f} s",
        f"  Phase H rounds   : {h_info['n_outer_rounds']}",
        f"  Pass E acc       : {h_info.get('pass_e_accepts', 0)} "
        f"(swap={h_info.get('pass_e_swap_accepts', 0)}, "
        f"split={h_info.get('pass_e_split_accepts', 0)})",
        f"  Pass E rej       : {h_info.get('pass_e_rejected', 0)}",
        "",
        "  PoC #54c baseline (raw coastline):",
        "    final: C1=8  C2=0  C4=69  C5=0   total=77",
        "  PoC #55 (simplified coastline):",
        f"    final: C1={phase_h_metrics['C1']}  "
        f"C2={phase_h_metrics['C2']}  "
        f"C4={phase_h_metrics['C4']}  "
        f"C5={phase_h_metrics['C5']}   total="
        f"{phase_h_metrics['C1'] + phase_h_metrics['C2'] + phase_h_metrics['C4'] + phase_h_metrics['C5']}",
    ]
    SUMMARY_TXT.write_text("\n".join(lines) + "\n", encoding="utf-8")
    print()
    print("\n".join(lines))
    print()
    print(f"wrote {COASTLINE_SIMPLIFIED}")
    print(f"wrote {RAW_OUT}")
    print(f"wrote {PHASE_G_OUT}")
    print(f"wrote {PHASE_H_OUT}")
    print(f"wrote {SUMMARY_TXT}")
    print(f"wrote {SUMMARY_JSON}")
    return 0


if __name__ == "__main__":
    raise SystemExit(main())
