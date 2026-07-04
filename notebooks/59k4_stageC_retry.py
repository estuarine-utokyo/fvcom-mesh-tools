"""PoC #59k-4 — stage C retry with a corrected convergence loop.

#59k-3's loop had three defects its log exposed: the junction branch
pre-empted every cycle (quality C1=38/C2=2/C4=20 and perp=167 were
never handled), manifold pinch nodes (3, created by the deletion
rounds) had NO handler at all, and end-trims cannot fix mid-arc
junction elements (counts sat constant while deletions converged
18 -> 6 -> 3). This retry starts from the stage-B checkpoint and
runs a strict priority ladder, one action per cycle:

  1. manifold  — delete every element incident to a pinch node or an
                 over-shared edge, keep largest component, re-derive
                 the arc;
  2. junction  — delete R4 / fake-ISBCE=2 elements directly (one
                 initial end-trim only);
  3. quality   — ``phase_h_finish`` + ``compact_nodes``;
  4. perp      — toolkit local fixer.

A stagnation detector (same failed-set fingerprint 4x) breaks out
instead of burning cycles.
"""

from __future__ import annotations

import json
import time
from pathlib import Path

import numpy as np
from pyproj import Transformer

from fvcom_mesh_tools.algorithms import align_open_boundary_local
from fvcom_mesh_tools.algorithms.boundary import (
    boundary_edges_from_tris,
    chain_edges_to_loops,
    outer_loop,
)
from fvcom_mesh_tools.io import Fort14Mesh, read_fort14, write_fort14
from fvcom_mesh_tools.io.fvcom_native import export_fvcom_case
from fvcom_mesh_tools.mesh_clean import compact_nodes, keep_components, remove_elements
from fvcom_mesh_tools.mesh_clean_phase_h import phase_h_finish
from fvcom_mesh_tools.qa import (
    _edge_topology,
    format_report,
    fvcom_boundary_element_flags,
    run_qa,
)

REPO = Path(__file__).resolve().parents[1]
SRC = REPO / "outputs" / "59k_stageB_optimized.14"
OUT = REPO / "outputs" / "59k_v3_gate_passed.14"
QA_JSON = REPO / "outputs" / "59k_v3_gate_passed_qa.json"
LOG_JSON = REPO / "outputs" / "59k4_cycle_log.json"
EXPORT_DIR = REPO / "outputs" / "fvcom_inputs_v3"
CASENAME = "tokyo_bay_v3"

BAND_DEG = 0.008
LAND_IBTYPE = 20
MAX_CYCLES = 40
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

_TO_LL = Transformer.from_crs("EPSG:32654", "EPSG:4326", always_xy=True)


def _arc(mesh):
    lon, lat = _TO_LL.transform(mesh.nodes[:, 0], mesh.nodes[:, 1])
    loops = chain_edges_to_loops(boundary_edges_from_tris(mesh.elements))
    outer = outer_loop(loops, mesh.nodes)
    ring = outer[:-1]
    rlon, rlat = lon[ring], lat[ring]
    mask = (rlat <= rlat.min() + BAND_DEG) | (rlon <= rlon.min() + BAND_DEG)
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
    ring = np.roll(ring, -a)
    islands = [lp[:-1].copy() for lp in loops if lp is not outer]
    return ring, islands, b - a


def _apply_obc(mesh, ring, islands, open_end, trim=1):
    """trim=1 by default: drop one node at each arc end so junction
    corner triangles get an interior third node."""
    lo, hi = trim, open_end - trim
    if hi - lo < 2:
        raise SystemExit("[59k4] open arc trimmed away")
    open_seg = ring[lo : hi + 1].copy()
    land_seg = np.concatenate([ring[hi:], ring[: lo + 1]])
    land = [(LAND_IBTYPE, land_seg)] + [(LAND_IBTYPE, i.copy()) for i in islands]
    return Fort14Mesh(
        title=mesh.title, nodes=mesh.nodes, depths=mesh.depths,
        elements=mesh.elements, open_boundaries=[open_seg],
        land_boundaries=land,
    )


def _pinch_elements(mesh) -> np.ndarray:
    """Mask of elements incident to a >2-boundary-edge node or an
    edge shared by >2 elements."""
    topo = _edge_topology(mesh.elements, mesh.n_nodes)
    buv = topo.uv[topo.counts == 1]
    cnt = np.zeros(mesh.n_nodes, dtype=np.int64)
    if buv.size:
        np.add.at(cnt, buv.ravel(), 1)
    pinch_nodes = np.where(cnt > 2)[0]
    over_uv = topo.uv[topo.counts > 2]
    bad = np.zeros(mesh.n_elements, dtype=bool)
    if pinch_nodes.size:
        bad |= np.isin(mesh.elements, pinch_nodes).any(axis=1)
    for u, v in over_uv:
        has_u = (mesh.elements == u).any(axis=1)
        has_v = (mesh.elements == v).any(axis=1)
        bad |= has_u & has_v
    return bad


def _failed(report):
    return {c.check_id for c in report.checks
            if c.gate and not c.skipped and not c.passed}


def _reobc(mesh):
    ring, islands, open_end = _arc(mesh)
    return _apply_obc(mesh, ring, islands, open_end)


def main() -> int:
    t0 = time.perf_counter()
    mesh = read_fort14(SRC)
    print(f"[59k4] stage C input: NP={mesh.n_nodes:,} NE={mesh.n_elements:,}",
          flush=True)
    mesh = _reobc(mesh)
    seg = mesh.open_boundaries[0]
    print(f"[59k4] arc OBC: {len(seg)} nodes", flush=True)

    log = []
    fingerprints: list[str] = []
    report = run_qa(mesh, name=OUT.name, path=OUT)
    for cycle in range(MAX_CYCLES):
        failed = _failed(report)
        counts = {c.check_id: int(c.n_violations) for c in report.checks
                  if c.gate and not c.skipped and c.n_violations}
        print(f"[59k4] cycle {cycle}: failed = {sorted(failed)} {counts}",
              flush=True)
        log.append({"cycle": cycle, "failed": sorted(failed), "counts": counts})
        if not failed:
            break
        fp = json.dumps(sorted(counts.items()))
        fingerprints.append(fp)
        if fingerprints.count(fp) >= 4:
            print("[59k4] stagnation detected — stopping", flush=True)
            break

        if "manifold_boundary" in failed:
            bad = _pinch_elements(mesh)
            mesh = remove_elements(mesh, ~bad)
            mesh, _ = keep_components(mesh)
            mesh = _reobc(mesh)
            print(f"[59k4] cycle {cycle}: deleted {int(bad.sum())} pinch "
                  "elements", flush=True)
        elif failed & JUNCTION_IDS:
            flags = fvcom_boundary_element_flags(mesh)
            bad = flags["r4_mask"] | flags["fake_open_mask"]
            if bad.any():
                mesh = remove_elements(mesh, ~bad)
                mesh, _ = keep_components(mesh)
                mesh = _reobc(mesh)
                print(f"[59k4] cycle {cycle}: deleted {int(bad.sum())} R4/fake",
                      flush=True)
            else:
                # obc_interior_neighbor etc. without R4 elements: the
                # necked node's ring is fully open — trim deeper.
                mesh = _reobc(mesh)
                print(f"[59k4] cycle {cycle}: junction residual w/o R4 — "
                      "re-applied arc", flush=True)
        elif failed & QUALITY_IDS:
            mesh, finfo = phase_h_finish(mesh, seed=100 + cycle)
            mesh, cinfo = compact_nodes(mesh)
            print(f"[59k4] cycle {cycle}: finish {finfo.get('before')} -> "
                  f"{finfo.get('after')}; compacted "
                  f"{cinfo['n_orphans_removed']}", flush=True)
            if cinfo["n_orphans_removed"]:
                mesh = _reobc(mesh)
        elif "obc_perpendicularity" in failed:
            mesh, pinfo = align_open_boundary_local(mesh, seed=9100 + cycle)
            print(f"[59k4] cycle {cycle}: perp fixer accepted="
                  f"{pinfo['accepted_total']} remaining="
                  f"{len(pinfo['remaining'])}", flush=True)
            if pinfo["accepted_total"] == 0:
                report = run_qa(mesh, name=OUT.name, path=OUT)
                log.append({"cycle": f"{cycle}-final",
                            "failed": sorted(_failed(report))})
                break
        else:
            print(f"[59k4] cycle {cycle}: unhandled — stopping", flush=True)
            break
        report = run_qa(mesh, name=OUT.name, path=OUT)

    mesh.title = "PoC 59k UTM54N Tokyo Bay v3 (100 m coastline-fidelity)"
    write_fort14(mesh, OUT)
    QA_JSON.write_text(json.dumps(report.to_dict(), indent=2, ensure_ascii=False),
                       encoding="utf-8")
    LOG_JSON.write_text(json.dumps(log, indent=2, ensure_ascii=False),
                        encoding="utf-8")
    print(format_report(report, lang="ja"), flush=True)

    failed = _failed(report)
    if not (failed & FATAL_OR_QUALITY):
        lon, lat = _TO_LL.transform(mesh.nodes[:, 0], mesh.nodes[:, 1])
        obc = mesh.open_boundaries[0]
        sponge = [(int(v), 3000.0, 0.001) for v in obc]
        written = export_fvcom_case(
            mesh, EXPORT_DIR, CASENAME, obc_type=1, cor=lat, sponge=sponge,
        )
        for k, p in written.items():
            print(f"[59k4] export {k}: {p}", flush=True)
    else:
        print(f"[59k4] EXPORT SKIPPED: {sorted(failed & FATAL_OR_QUALITY)}",
              flush=True)
    print(f"[59k4] wall: {time.perf_counter() - t0:.1f} s", flush=True)
    return 0 if not failed else 1


if __name__ == "__main__":
    raise SystemExit(main())
