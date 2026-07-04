"""PoC #74 — session v6: adds the targeted per-site edit stage (E).

Stage E (each step runs ONCE, no loops):
  E1  structural deletions on the current state (pinch / R4 / fake
      open) — deterministic one-shot cleanup;
  E2  insert_node_on_line on wide boundary edges next to the
      snap-deferred coordinates (the spacing slack that pure node
      motion lacks);
  E3  second chain-snap pass over the modified geometry.

v1 (#69) snapped all reachable nodes first and ground the damage
afterwards — the one bulk optimize round then removed 29% of the
nodes and still left 99 sites. v2 applies the SMS discipline at
the snap itself: each node move is accepted only if its 1-ring
does not gain C1/C2/C4 violations; rejected nodes are deferred to
the per-site list instead of poisoning the mesh.

User directive 2026-07-04: iterate-until-gates-pass is the wrong
model — SMS never converges automatically either; a human inspects
each flagged element and applies a situation-specific local edit.
This script is that workflow with the agent as the editor:

  A. bounded one-shot preprocessing (each stage runs ONCE, with a
     checkpoint): channel exclusion (w/h < 2), R4/fake deletion
     (deterministic, hard-capped), 2 m clip, UTM, snap-to-raw-OSM
     (frac 1.2), weld + tiny-element drop, arc OBC + line snap;
  B. ONE bulk grinding round: phase_h_optimize(max_outer_rounds=1,
     Pass A/B + F/G, OSM projector) — a pre-step, not a loop;
  C. ONE targeted pass: the toolkit local perp fixer (per-violation,
     the 59e model);
  D. ONE final QA; residual offenders are CLUSTERED into sites and
     REPORTED (JSON + per-site zoom figures = the SMS colouring
     equivalent) for per-site AI/human judgment — never looped.
"""

from __future__ import annotations

import json
import time
from pathlib import Path

import numpy as np
import shapely
from pyproj import Transformer
from shapely.strtree import STRtree

from fvcom_mesh_tools.algorithms import (
    align_open_boundary_local,
    snap_boundary_chains,
    snap_nodes_to_segment,
)
from fvcom_mesh_tools.algorithms.boundary import (
    boundary_edges_from_tris,
    chain_edges_to_loops,
    outer_loop,
)
from fvcom_mesh_tools.algorithms.site_edits import insert_node_on_line
from fvcom_mesh_tools.diagnostics import under_resolved_channels_flag
from fvcom_mesh_tools.io import Fort14Mesh, read_fort14, write_fort14
from fvcom_mesh_tools.io.fvcom_native import export_fvcom_case
from fvcom_mesh_tools.mesh_clean import (
    compact_nodes,
    keep_components,
    rebuild_boundaries,
    remove_elements,
    weld_close_nodes,
)
from fvcom_mesh_tools.mesh_clean_phase_h import (
    build_coastline_projector,
    phase_h_optimize,
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
CKPT_A = REPO / "outputs" / "74_stageA_prepped.14"
CKPT_B = REPO / "outputs" / "74_stageB_ground.14"
OUT = REPO / "outputs" / "74_v4_session.14"
QA_JSON = REPO / "outputs" / "74_v4_session_qa.json"
SITES_JSON = REPO / "outputs" / "74_residual_sites.json"
FIG_DIR = REPO / "outputs" / "figures" / "74_sites"
EXPORT_DIR = REPO / "outputs" / "fvcom_inputs_v4"
CASENAME = "tokyo_bay_v4"

BBOX = (139.565, 35.10, 140.172, 35.86)
TOL_DEG = 450.0 / 111_195.0
BAND_DEG = 0.012
LAND_IBTYPE = 20
MIN_DEPTH_M = 2.0
MIN_WH = 2.0
SITE_RADIUS_M = 900.0
MAX_SITE_FIGS = 24
QUALITY_IDS = {"c1_min_angle", "c2_max_angle", "c4_area_change", "c5_valence"}
FATAL = {
    "node_index_valid", "ccw_all_elements", "no_isolated_elements",
    "r4_mixed_boundary", "manifold_boundary", "no_duplicate_nodes",
    "no_orphan_nodes", "no_tiny_area", "isbce2_authentic",
    "obc_on_boundary", "obc_chain_adjacency", "obc_interior_neighbor",
    "obc_ordering", "single_component", "obc_reachable", "min_depth_clip",
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
        lines.extend(list(b.geoms) if hasattr(b, "geoms") else [b])
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
        [shp], max_snap_distance_m=500.0 / _DEG_PER_M_AT_EQ,
        mean_latitude_deg=0.0,
    )


def _arc_bands(mesh):
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


def _signed_areas_m(mesh):
    p0 = mesh.nodes[mesh.elements[:, 0]]
    p1 = mesh.nodes[mesh.elements[:, 1]]
    p2 = mesh.nodes[mesh.elements[:, 2]]
    return 0.5 * ((p1[:, 0] - p0[:, 0]) * (p2[:, 1] - p0[:, 1])
                  - (p1[:, 1] - p0[:, 1]) * (p2[:, 0] - p0[:, 0]))


def main() -> int:
    import logging

    logging.basicConfig(level=logging.INFO,
                        format="%(asctime)s %(message)s",
                        datefmt="%H:%M:%S")
    t0 = time.perf_counter()

    # ---- A. bounded one-shot preprocessing --------------------------------
    mesh = read_fort14(SRC)
    print(f"[74] A0 input: NP={mesh.n_nodes:,} NE={mesh.n_elements:,}",
          flush=True)
    mesh, _ = keep_components(mesh)
    mesh = _rebuild_ll(mesh)
    flag, _ = under_resolved_channels_flag(mesh, min_w_h=MIN_WH)
    if flag.any():
        mesh = remove_elements(mesh, ~flag)
        mesh, _ = keep_components(mesh)
        mesh = _rebuild_ll(mesh)
    print(f"[74] A1 channel exclusion: -{int(flag.sum())} elems", flush=True)
    for rnd in range(12):
        f = fvcom_boundary_element_flags(mesh)
        bad = f["r4_mask"] | f["fake_open_mask"]
        if not bad.any():
            break
        mesh = remove_elements(mesh, ~bad)
        mesh, _ = keep_components(mesh)
        mesh = _rebuild_ll(mesh)
    print(f"[74] A2 R4/fake deletion done (rounds={rnd + 1})", flush=True)
    mesh.depths[:] = np.maximum(mesh.depths, MIN_DEPTH_M)
    x, y = _TO_M.transform(mesh.nodes[:, 0], mesh.nodes[:, 1])
    mesh.nodes = np.column_stack([x, y])

    lines = _osm_lines_utm()
    ring, islands, open_end, west_n, south_n = _arc_bands(mesh)
    arc_all = set(int(v) for v in ring[:open_end + 1])
    mesh, sinfo = snap_boundary_chains(
        mesh, lines, max_snap_frac=1.2, exclude_nodes=list(arc_all),
    )
    # Record deferred nodes as COORDINATES: stage B compacts node ids.
    deferred_snap_xy = [
        (float(mesh.nodes[v, 0]), float(mesh.nodes[v, 1]))
        for c in sinfo.pop("deferred_chains") for v in c
    ]
    print(f"[74] A3 snap: {sinfo}", flush=True)
    for run_nodes, label in ((west_n, "west"), (south_n, "south")):
        if run_nodes.size >= 2:
            p0 = tuple(mesh.nodes[run_nodes[0]])
            p1 = tuple(mesh.nodes[run_nodes[-1]])
            mesh, li = snap_nodes_to_segment(
                mesh, [int(v) for v in run_nodes], p0, p1, max_move=600.0,
            )
            print(f"[74] A4 OBC {label} line: {li}", flush=True)
    mesh, winfo = weld_close_nodes(mesh, tol=2.0)
    print(f"[74] A5 weld: {winfo}", flush=True)
    sa = _signed_areas_m(mesh)
    tiny = np.abs(sa) < 100.0
    if tiny.any():
        mesh = remove_elements(mesh, ~tiny)
        mesh, _ = keep_components(mesh)
        print(f"[74] A6 tiny drop: -{int(tiny.sum())}", flush=True)
    ring, islands, open_end, _w, _s = _arc_bands(mesh)
    mesh = _apply_obc(mesh, ring, islands, open_end)
    write_fort14(mesh, CKPT_A)
    print(f"[74] A done ({time.perf_counter() - t0:.0f} s) -> {CKPT_A.name}",
          flush=True)

    # ---- B. ONE bulk grinding round ---------------------------------------
    projector = _osm_projector()
    t1 = time.perf_counter()
    mesh, _oinfo = phase_h_optimize(
        mesh,
        min_angle_target=30.0,
        max_angle_target=130.0,
        pass_f_enabled=True,
        pass_g_enabled=True,
        pass_g_min_angle_target=30.0,
        max_outer_rounds=1,
        coastline_projector=projector,
        time_budget_s=600.0,
    )
    if _oinfo.get("budget_exhausted"):
        print("[74] B: 10-min budget hit — continuing to targeted "
              "editing with the current state", flush=True)
    mesh, cinfo = compact_nodes(mesh)
    if cinfo["n_orphans_removed"]:
        ring, islands, open_end, _w, _s = _arc_bands(mesh)
        mesh = _apply_obc(mesh, ring, islands, open_end)
    write_fort14(mesh, CKPT_B)
    print(f"[74] B optimize x1 done ({time.perf_counter() - t1:.0f} s) -> "
          f"{CKPT_B.name}", flush=True)

    # ---- C. ONE targeted perp pass ----------------------------------------
    mesh, pinfo = align_open_boundary_local(mesh, seed=7400, max_outer=1)
    print(f"[74] C perp pass: accepted={pinfo['accepted_total']} "
          f"remaining={len(pinfo['remaining'])}", flush=True)

    # ---- E. targeted per-site edits (each step once) ----------------------
    # E1: structural deletions on the current state.
    for label, mask_fn in (
        ("pinch", _pinch_elements),
        ("r4fake", lambda m: (lambda f: f["r4_mask"] | f["fake_open_mask"])(
            fvcom_boundary_element_flags(m))),
    ):
        bad = mask_fn(mesh)
        if bad.any():
            mesh = remove_elements(mesh, ~bad)
            mesh, _ = keep_components(mesh)
            print(f"[74] E1 {label}: -{int(bad.sum())} elems", flush=True)
    ring, islands, open_end, _w, _s = _arc_bands(mesh)
    mesh = _apply_obc(mesh, ring, islands, open_end)
    arc_all = set(int(v) for v in ring[:open_end + 1])

    # E2: on-line inserts on wide boundary edges near deferred coords.
    from scipy.spatial import cKDTree as _KD

    from fvcom_mesh_tools.algorithms.boundary_snap import (
        _boundary_topology as _bt,
    )

    n_ins = 0
    if deferred_snap_xy:
        def_tree = _KD(np.asarray(deferred_snap_xy))
        buv, _n2e = _bt(mesh.elements, mesh.n_nodes)
        emid = 0.5 * (mesh.nodes[buv[:, 0]] + mesh.nodes[buv[:, 1]])
        elen = np.linalg.norm(
            mesh.nodes[buv[:, 0]] - mesh.nodes[buv[:, 1]], axis=1,
        )
        med = float(np.median(elen))
        d_def, _ = def_tree.query(emid, workers=-1)
        cand = np.where((d_def < 300.0) & (elen > 1.3 * med))[0]
        line_tree = STRtree(lines)
        garr = np.array(lines, dtype=object)
        order = np.argsort(-elen[cand])
        for k in cand[order]:
            u, v = int(buv[k, 0]), int(buv[k, 1])
            if u in arc_all or v in arc_all:
                continue
            mid = 0.5 * (mesh.nodes[u] + mesh.nodes[v])
            gi = line_tree.nearest(shapely.points(mid[0], mid[1]))
            res = insert_node_on_line(mesh, u, v, garr[gi])
            if res is not None:
                mesh = res[0]
                n_ins += 1
        print(f"[74] E2 inserts: {n_ins} new on-line nodes "
              f"({cand.size} candidate edges)", flush=True)
        if n_ins:
            ring, islands, open_end, _w, _s = _arc_bands(mesh)
            mesh = _apply_obc(mesh, ring, islands, open_end)
            arc_all = set(int(v) for v in ring[:open_end + 1])

    # E3: second chain-snap pass over the modified geometry.
    mesh, sinfo2 = snap_boundary_chains(
        mesh, lines, max_snap_frac=1.2, exclude_nodes=list(arc_all),
    )
    deferred_snap_xy = [
        (float(mesh.nodes[v, 0]), float(mesh.nodes[v, 1]))
        for c in sinfo2.pop("deferred_chains") for v in c
    ]
    print(f"[74] E3 snap pass 2: {sinfo2}", flush=True)

    # ---- D. final QA + residual-site report -------------------------------
    report = run_qa(mesh, name=OUT.name, path=OUT, max_offenders=100000)
    mesh.title = "PoC 74 v4 session mesh (session + per-site edits)"
    write_fort14(mesh, OUT)
    QA_JSON.write_text(json.dumps(report.to_dict(), indent=2,
                                  ensure_ascii=False), encoding="utf-8")
    print(format_report(
        report, lang="ja",
    ).split("違反箇所")[0], flush=True)

    # Cluster residual offenders into sites.
    pts, recs = [], []
    for x, y in deferred_snap_xy:
        pts.append((x, y))
        recs.append({"check": "snap_deferred", "x": x, "y": y})
    for c in report.checks:
        if not c.gate or c.skipped or c.passed:
            continue
        for o in c.offenders:
            pts.append((o["x"], o["y"]))
            recs.append({"check": c.check_id, **o})
    sites = []
    if pts:
        from scipy.cluster.hierarchy import fcluster, linkage
        from scipy.spatial.distance import pdist

        P = np.asarray(pts)
        if len(P) > 1:
            lab = fcluster(linkage(pdist(P), method="single"),
                           t=SITE_RADIUS_M, criterion="distance")
        else:
            lab = np.array([1])
        for k in sorted(set(lab)):
            sel = [r for r, m in zip(recs, lab == k) if m]
            cx = float(np.mean([r["x"] for r in sel]))
            cy = float(np.mean([r["y"] for r in sel]))
            sites.append({
                "site": int(k), "x": cx, "y": cy, "n_offenders": len(sel),
                "checks": sorted({r["check"] for r in sel}),
                "offenders": sel,
            })
        sites.sort(key=lambda s: -s["n_offenders"])
    SITES_JSON.write_text(json.dumps(sites, indent=2, ensure_ascii=False),
                          encoding="utf-8")
    print(f"[74] residual sites: {len(sites)} "
          f"(offenders total {len(recs)})", flush=True)

    # Zoom figure per site (SMS-colouring equivalent).
    import os

    os.environ.setdefault("MPLBACKEND", "Agg")
    from fvcom_mesh_tools.plotting import plot_mesh_overview

    FIG_DIR.mkdir(parents=True, exist_ok=True)
    for s in sites[:MAX_SITE_FIGS]:
        r = 700.0
        zb = (s["x"] - r, s["y"] - r, s["x"] + r, s["y"] + r)
        try:
            plot_mesh_overview(
                mesh, FIG_DIR / f"site{s['site']:03d}.png",
                crs="EPSG:32654", cell_m=None, coast=None,
                zoom=zb, dpi=200,
                title=f"site {s['site']}: {','.join(s['checks'])} "
                      f"x{s['n_offenders']}",
            )
        except Exception as e:  # keep reporting even if one figure fails
            print(f"[74] site fig {s['site']} failed: {e}", flush=True)
    print(f"[74] site figures -> {FIG_DIR}", flush=True)

    failed = {c.check_id for c in report.checks
              if c.gate and not c.skipped and not c.passed}
    if not (failed & FATAL):
        _lon, lat = _TO_LL.transform(mesh.nodes[:, 0], mesh.nodes[:, 1])
        obc = mesh.open_boundaries[0]
        sponge = [(int(v), 3000.0, 0.001) for v in obc]
        written = export_fvcom_case(mesh, EXPORT_DIR, CASENAME, obc_type=1,
                                    cor=lat, sponge=sponge)
        for k, pth in written.items():
            print(f"[74] export {k}: {pth}", flush=True)
    print(f"[74] wall: {time.perf_counter() - t0:.0f} s", flush=True)
    return 0


if __name__ == "__main__":
    raise SystemExit(main())
