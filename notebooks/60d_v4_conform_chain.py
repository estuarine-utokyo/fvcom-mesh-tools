"""PoC #60d — v4 (300 m, OSM) conform chain: the two user-flagged
defects fixed at coarse resolution.

1. Fragmented sub-grid channels (Obitsu at Banzu, etc.): **exclusion
   criterion** — elements with channel w/h < 2 are removed (a 300 m
   mesh cannot carry a < 600 m channel; rivers stay represented as
   mouth-inflow boundary segments per TB-FVCOM practice). Detection =
   the QA channel metric; deletion + component cleanup follow.
2. Boundary conformity: land-boundary nodes are SNAPPED exactly onto
   the OSM true-land polygon boundaries (cap 0.6 x local edge — smooth
   coast becomes exact, sub-grid structures stay simplified); the
   west/south artificial cuts are snapped onto fitted straight lines
   so the open boundary is a clean line, then the standard repair
   ladder (junction deletes -> projector finish -> local perp fixer)
   re-converges the QA gate. The finish-stage projector uses the SAME
   OSM polylines, so repairs keep snapped nodes ON the coastline.

Output: outputs/60d_v4_gate.14 (+ QA json, cycle log,
fvcom_inputs_v4/tokyo_bay_v4_* with OBC sponge).
"""

from __future__ import annotations

import json
import time
from pathlib import Path

import numpy as np
from pyproj import Transformer

from fvcom_mesh_tools.algorithms import (
    align_open_boundary_local,
    snap_boundary_to_polylines,
    snap_nodes_to_segment,
)
from fvcom_mesh_tools.algorithms.boundary import (
    boundary_edges_from_tris,
    chain_edges_to_loops,
    outer_loop,
)
from fvcom_mesh_tools.diagnostics import under_resolved_channels_flag
from fvcom_mesh_tools.io import Fort14Mesh, read_fort14, write_fort14
from fvcom_mesh_tools.io.fvcom_native import export_fvcom_case
from fvcom_mesh_tools.mesh_clean import (
    compact_nodes,
    keep_components,
    rebuild_boundaries,
    remove_elements,
)
from fvcom_mesh_tools.mesh_clean_phase_h import (
    build_coastline_projector,
    phase_h_finish,
)
from fvcom_mesh_tools.qa import (
    _edge_topology,
    format_report,
    fvcom_boundary_element_flags,
    run_qa,
)

REPO = Path(__file__).resolve().parents[1]
SRC = REPO / "outputs" / "60c_v4_raw_300m_osm.14"
OSM_SHP = REPO / "outputs" / "osm_shoreline" / "osm_true_land_tokyo_bay.shp"
OSM_UTM_DIR = REPO / "outputs" / "osm_shoreline_utm_nocrs"
OUT = REPO / "outputs" / "60d_v4_gate.14"
QA_JSON = REPO / "outputs" / "60d_v4_gate_qa.json"
LOG_JSON = REPO / "outputs" / "60d_cycle_log.json"
EXPORT_DIR = REPO / "outputs" / "fvcom_inputs_v4"
CASENAME = "tokyo_bay_v4"

BBOX = (139.565, 35.10, 140.172, 35.86)
TOL_DEG = 450.0 / 111_195.0
MIN_DEPTH_M = 2.0
MIN_WH = 2.0                 # channel exclusion criterion (kickoff §7.4)
BAND_DEG = 0.012             # artificial-cut band for a 300 m boundary
LAND_IBTYPE = 20
MAX_CYCLES = 30
QUALITY_IDS = {"c1_min_angle", "c2_max_angle", "c4_area_change", "c5_valence"}
JUNCTION_IDS = {"r4_mixed_boundary", "isbce2_authentic", "obc_interior_neighbor",
                "obc_chain_adjacency", "obc_on_boundary"}
FATAL_OR_QUALITY = {
    "node_index_valid", "ccw_all_elements", "no_isolated_elements",
    "r4_mixed_boundary", "manifold_boundary", "no_duplicate_nodes",
    "no_orphan_nodes", "no_tiny_area", "isbce2_authentic",
    "obc_on_boundary", "obc_chain_adjacency", "obc_interior_neighbor",
    "obc_ordering", "c1_min_angle", "c2_max_angle", "c4_area_change",
    "c5_valence", "single_component", "obc_reachable", "min_depth_clip",
}
_DEG_PER_M_AT_EQ = 1.0 / 111_194.92664455873

_TO_LL = Transformer.from_crs("EPSG:32654", "EPSG:4326", always_xy=True)
_TO_M = Transformer.from_crs("EPSG:4326", "EPSG:32654", always_xy=True)


def _rebuild_ll(mesh):
    return rebuild_boundaries(mesh, bbox=BBOX, tol_deg=TOL_DEG,
                              land_ibtype=LAND_IBTYPE, open_merge_coast_gap=50)


def _osm_lines_utm():
    import geopandas as gpd

    gdf = gpd.read_file(OSM_SHP).to_crs(32654)
    lines = []
    for g in gdf.geometry:
        if g is None or g.is_empty:
            continue
        b = g.boundary
        if b.geom_type == "LineString":
            lines.append(b)
        else:
            lines.extend(list(b.geoms))
    return lines


def _osm_projector():
    import geopandas as gpd

    OSM_UTM_DIR.mkdir(parents=True, exist_ok=True)
    shp = OSM_UTM_DIR / "osm_true_land_utm54.shp"
    if not shp.exists():
        gdf = gpd.read_file(OSM_SHP).to_crs(32654)
        gdf = gdf.set_crs(None, allow_override=True)
        gdf.to_file(shp)
    return build_coastline_projector(
        [shp],
        max_snap_distance_m=500.0 / _DEG_PER_M_AT_EQ,
        mean_latitude_deg=0.0,
    )


def _arc_bands(mesh):
    """Outer-ring artificial arc plus the per-band (west / south)
    node lists for straight-line snapping."""
    lon, lat = _TO_LL.transform(mesh.nodes[:, 0], mesh.nodes[:, 1])
    loops = chain_edges_to_loops(boundary_edges_from_tris(mesh.elements))
    outer = outer_loop(loops, mesh.nodes)
    ring = outer[:-1]
    rlon, rlat = lon[ring], lat[ring]
    south_b = rlat <= rlat.min() + BAND_DEG
    west_b = rlon <= rlon.min() + BAND_DEG
    mask = south_b | west_b
    idx = np.where(mask)[0]
    runs = []
    s = p = int(idx[0])
    for q in idx[1:]:
        q = int(q)
        if q == p + 1:
            p = q
        else:
            runs.append((s, p))
            s = p = q
    runs.append((s, p))
    a, b = max(runs, key=lambda r: r[1] - r[0])
    arc_pos = np.arange(a, b + 1)
    arc_nodes = ring[a:b + 1]
    west_nodes = arc_nodes[west_b[arc_pos]]
    south_nodes = arc_nodes[south_b[arc_pos] & ~west_b[arc_pos]]
    ring = np.roll(ring, -a)
    islands = [lp[:-1].copy() for lp in loops if lp is not outer]
    return ring, islands, b - a, west_nodes, south_nodes


def _apply_obc(mesh, ring, islands, open_end, trim=1):
    lo, hi = trim, open_end - trim
    if hi - lo < 2:
        raise SystemExit("[60d] open arc trimmed away")
    open_seg = ring[lo:hi + 1].copy()
    land_seg = np.concatenate([ring[hi:], ring[:lo + 1]])
    land = [(LAND_IBTYPE, land_seg)] + [(LAND_IBTYPE, i.copy()) for i in islands]
    return Fort14Mesh(
        title=mesh.title, nodes=mesh.nodes, depths=mesh.depths,
        elements=mesh.elements, open_boundaries=[open_seg],
        land_boundaries=land,
    )


def _pinch_elements(mesh):
    topo = _edge_topology(mesh.elements, mesh.n_nodes)
    buv = topo.uv[topo.counts == 1]
    cnt = np.zeros(mesh.n_nodes, dtype=np.int64)
    if buv.size:
        np.add.at(cnt, buv.ravel(), 1)
    pinch = np.where(cnt > 2)[0]
    bad = np.zeros(mesh.n_elements, dtype=bool)
    if pinch.size:
        bad |= np.isin(mesh.elements, pinch).any(axis=1)
    for u, v in topo.uv[topo.counts > 2]:
        bad |= ((mesh.elements == u).any(axis=1)
                & (mesh.elements == v).any(axis=1))
    return bad


def _fit_chord(nodes_xy):
    """End-to-end chord of an ordered node run (robust for straight
    artificial cuts)."""
    return tuple(nodes_xy[0]), tuple(nodes_xy[-1])


def _failed(report):
    return {c.check_id for c in report.checks
            if c.gate and not c.skipped and not c.passed}


def main() -> int:
    t0 = time.perf_counter()
    mesh = read_fort14(SRC)
    print(f"[60d] input: NP={mesh.n_nodes:,} NE={mesh.n_elements:,}", flush=True)

    # -- structural + channel exclusion (lon/lat) ---------------------------
    mesh, kinfo = keep_components(mesh)
    mesh = _rebuild_ll(mesh)
    flag, _info = under_resolved_channels_flag(mesh, min_w_h=MIN_WH)
    n_ch = int(flag.sum())
    if n_ch:
        mesh = remove_elements(mesh, ~flag)
        mesh, _ = keep_components(mesh)
        mesh = _rebuild_ll(mesh)
    print(f"[60d] channel exclusion (w/h < {MIN_WH:g}): removed {n_ch} "
          "elements (sub-grid channels; rivers stay as mouth inflows)",
          flush=True)
    for _rnd in range(30):
        flags = fvcom_boundary_element_flags(mesh)
        bad = flags["r4_mask"] | flags["fake_open_mask"]
        if not bad.any():
            break
        mesh = remove_elements(mesh, ~bad)
        mesh, _ = keep_components(mesh)
        mesh = _rebuild_ll(mesh)
    mesh.depths[:] = np.maximum(mesh.depths, MIN_DEPTH_M)

    # -- project + snap ------------------------------------------------------
    x, y = _TO_M.transform(mesh.nodes[:, 0], mesh.nodes[:, 1])
    mesh.nodes = np.column_stack([x, y])
    lines = _osm_lines_utm()
    ring, islands, open_end, west_n, south_n = _arc_bands(mesh)
    arc_all = set(int(v) for v in ring[:open_end + 1])
    mesh, sinfo = snap_boundary_to_polylines(
        mesh, lines, exclude_nodes=list(arc_all),
    )
    print(f"[60d] OSM snap: {sinfo}", flush=True)
    for nodes_run, label in ((west_n, "west"), (south_n, "south")):
        if nodes_run.size >= 2:
            p0, p1 = _fit_chord(mesh.nodes[nodes_run])
            mesh, li = snap_nodes_to_segment(
                mesh, [int(v) for v in nodes_run], p0, p1, max_move=600.0,
            )
            print(f"[60d] OBC {label} line snap: {li}", flush=True)

    mesh = _apply_obc(mesh, ring, islands, open_end)
    projector = _osm_projector()

    # -- repair ladder --------------------------------------------------------
    log = []
    fingerprints = []
    report = run_qa(mesh, name=OUT.name, path=OUT)
    for cycle in range(MAX_CYCLES):
        failed = _failed(report)
        counts = {c.check_id: int(c.n_violations) for c in report.checks
                  if c.gate and not c.skipped and c.n_violations}
        print(f"[60d] cycle {cycle}: failed = {sorted(failed)} {counts}",
              flush=True)
        log.append({"cycle": cycle, "failed": sorted(failed), "counts": counts})
        if not failed:
            break
        fp = json.dumps(sorted(counts.items()))
        fingerprints.append(fp)
        if fingerprints.count(fp) >= 4:
            print("[60d] stagnation — stopping", flush=True)
            break
        if "manifold_boundary" in failed:
            bad = _pinch_elements(mesh)
            mesh = remove_elements(mesh, ~bad)
            mesh, _ = keep_components(mesh)
            ring, islands, open_end, _w, _s = _arc_bands(mesh)
            mesh = _apply_obc(mesh, ring, islands, open_end)
            print(f"[60d] cycle {cycle}: deleted {int(bad.sum())} pinch elems",
                  flush=True)
        elif failed & JUNCTION_IDS:
            flags = fvcom_boundary_element_flags(mesh)
            bad = flags["r4_mask"] | flags["fake_open_mask"]
            if bad.any():
                mesh = remove_elements(mesh, ~bad)
                mesh, _ = keep_components(mesh)
                ring, islands, open_end, _w, _s = _arc_bands(mesh)
                mesh = _apply_obc(mesh, ring, islands, open_end)
                print(f"[60d] cycle {cycle}: deleted {int(bad.sum())} R4/fake",
                      flush=True)
            else:
                mesh = _apply_obc(mesh, ring, islands, open_end, trim=2)
                print(f"[60d] cycle {cycle}: junction residual — deeper trim",
                      flush=True)
        elif failed & QUALITY_IDS:
            mesh, finfo = phase_h_finish(
                mesh, seed=300 + cycle, coastline_projector=projector,
            )
            mesh, cinfo = compact_nodes(mesh)
            print(f"[60d] cycle {cycle}: finish {finfo.get('before')} -> "
                  f"{finfo.get('after')} (compacted "
                  f"{cinfo['n_orphans_removed']})", flush=True)
            if cinfo["n_orphans_removed"]:
                ring, islands, open_end, _w, _s = _arc_bands(mesh)
                mesh = _apply_obc(mesh, ring, islands, open_end)
        elif "obc_perpendicularity" in failed:
            mesh, pinfo = align_open_boundary_local(mesh, seed=9300 + cycle)
            print(f"[60d] cycle {cycle}: perp fixer accepted="
                  f"{pinfo['accepted_total']} remaining="
                  f"{len(pinfo['remaining'])}", flush=True)
            if pinfo["accepted_total"] == 0:
                report = run_qa(mesh, name=OUT.name, path=OUT)
                break
        else:
            print(f"[60d] cycle {cycle}: unhandled — stopping", flush=True)
            break
        report = run_qa(mesh, name=OUT.name, path=OUT)

    mesh.title = "PoC 60d UTM54N Tokyo Bay v4 (300 m, OSM-conformal)"
    write_fort14(mesh, OUT)
    QA_JSON.write_text(json.dumps(report.to_dict(), indent=2, ensure_ascii=False),
                       encoding="utf-8")
    LOG_JSON.write_text(json.dumps(log, indent=2, ensure_ascii=False),
                        encoding="utf-8")
    print(format_report(report, lang="ja"), flush=True)

    failed = _failed(report)
    if not (failed & FATAL_OR_QUALITY):
        _lon, lat = _TO_LL.transform(mesh.nodes[:, 0], mesh.nodes[:, 1])
        obc = mesh.open_boundaries[0]
        sponge = [(int(v), 3000.0, 0.001) for v in obc]
        written = export_fvcom_case(
            mesh, EXPORT_DIR, CASENAME, obc_type=1, cor=lat, sponge=sponge,
        )
        for k, p in written.items():
            print(f"[60d] export {k}: {p}", flush=True)
    else:
        print(f"[60d] EXPORT SKIPPED: {sorted(failed & FATAL_OR_QUALITY)}",
              flush=True)
    print(f"[60d] wall: {time.perf_counter() - t0:.1f} s", flush=True)
    return 0 if not failed else 1


if __name__ == "__main__":
    raise SystemExit(main())
