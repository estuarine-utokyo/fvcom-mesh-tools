"""PoC #81 — per-site manual edits on the 16 residual sites of #80.

The endgame of the AI manual-editing workflow: each residual site is
treated individually with the operator an SMS editor would pick,
accept-or-rollback per site, all decisions logged.

* quality sites (C1/C2/C4): patch re-CDT — collect the offending
  elements plus their node-sharing neighbours, keep the rim fixed,
  re-triangulate the interior (first WITHOUT the old interior nodes —
  sliver apexes simply vanish — retrying WITH them), accept only if
  the site's local violation count drops;
* perpendicularity sites: slide the OBC node ALONG its boundary
  chord (whitelist: on-line movement only), 17-sample deterministic
  line search gated on the local quality patch.

One pass, no loops; whatever remains is reported with fresh figures.
"""

from __future__ import annotations

import json
import time
from pathlib import Path

import numpy as np
from pyproj import Transformer

from fvcom_mesh_tools.algorithms.perp_local import _tri_quality
from fvcom_mesh_tools.io import Fort14Mesh, read_fort14, write_fort14
from fvcom_mesh_tools.io.fvcom_native import export_fvcom_case
from fvcom_mesh_tools.mesh_clean import (
    _patch_rim_polygon,
    _retriangulate_patch,
    compact_nodes,
)
from fvcom_mesh_tools.qa import format_report, run_qa

REPO = Path(__file__).resolve().parents[1]
SRC = REPO / "outputs" / "80_v4_session.14"
SITES = REPO / "outputs" / "80_residual_sites.json"
OUT = REPO / "outputs" / "81_v4_final.14"
QA_JSON = REPO / "outputs" / "81_v4_final_qa.json"
LOG_JSON = REPO / "outputs" / "81_edit_log.json"
EXPORT_DIR = REPO / "outputs" / "fvcom_inputs_v4"
CASENAME = "tokyo_bay_v4"

MIN_ANGLE = 30.0
MAX_ANGLE = 130.0
MAX_AC = 0.5
BAND_DEG = 0.012
LAND_IBTYPE = 20
_TO_LL = Transformer.from_crs("EPSG:32654", "EPSG:4326", always_xy=True)

FATAL = {
    "node_index_valid", "ccw_all_elements", "no_isolated_elements",
    "r4_mixed_boundary", "manifold_boundary", "no_duplicate_nodes",
    "no_orphan_nodes", "no_tiny_area", "isbce2_authentic",
    "obc_on_boundary", "obc_chain_adjacency", "obc_interior_neighbor",
    "obc_ordering", "single_component", "obc_reachable", "min_depth_clip",
}


def _region_violations(nodes, elements, region_e) -> int:
    """C1/C2/flip violations plus C4 pairs inside ``region_e``."""
    tri = elements[region_e]
    mn, mx, twice = _tri_quality(nodes[tri])
    n_bad = int((mn < MIN_ANGLE).sum() + (mx > MAX_ANGLE).sum()
                + (twice <= 0).sum())
    areas = 0.5 * np.abs(twice)
    pos = {int(e): k for k, e in enumerate(region_e)}
    edges: dict[tuple[int, int], list[int]] = {}
    for k, e in enumerate(region_e):
        a, b, c = (int(x) for x in elements[e])
        for u, v in ((a, b), (b, c), (c, a)):
            edges.setdefault((min(u, v), max(u, v)), []).append(int(e))
    for pair in edges.values():
        if len(pair) == 2:
            a1, a2 = areas[pos[pair[0]]], areas[pos[pair[1]]]
            hi = max(a1, a2)
            if hi > 0 and (hi - min(a1, a2)) / hi > MAX_AC:
                n_bad += 1
    return n_bad


def _neighbours_by_node(elements, n_nodes, elem_ids):
    sel = np.zeros(n_nodes, dtype=bool)
    sel[np.unique(elements[elem_ids].ravel())] = True
    return np.where(sel[elements].any(axis=1))[0]


def _try_patch_recdt(mesh, patch_e, protected: set[int]):
    """Re-CDT ``patch_e`` with the rim fixed. Returns (elements2,
    dropped_nodes, tag) or None."""
    elements = mesh.elements
    rim = _patch_rim_polygon(elements, patch_e)
    if rim is None:
        return None
    rim = np.asarray([int(v) for v in rim])
    patch_nodes = np.unique(elements[patch_e].ravel())
    interior = [int(v) for v in patch_nodes if v not in set(rim)]
    if any(v in protected for v in interior):
        keep_interior = interior
        variants = [("keep", keep_interior)]
    else:
        variants = [("drop", []), ("keep", interior)]
    for tag, spine_ids in variants:
        spine_xy = (mesh.nodes[spine_ids]
                    if spine_ids else np.empty((0, 2)))
        tris, reason = _retriangulate_patch(
            mesh.nodes[rim], spine_xy, len(rim),
        )
        if tris is None:
            continue
        id_map = np.concatenate([rim, np.asarray(spine_ids, dtype=int)]) \
            if spine_ids else rim
        new_tris = id_map[tris]
        elements2 = np.vstack([
            np.delete(elements, patch_e, axis=0), new_tris,
        ])
        dropped = [] if spine_ids else interior
        return elements2, dropped, tag
    return None


def main() -> int:
    t0 = time.perf_counter()
    mesh = read_fort14(SRC)
    sites = json.loads(SITES.read_text())
    print(f"[81] input: NP={mesh.n_nodes:,} NE={mesh.n_elements:,} "
          f"sites={len(sites)}", flush=True)

    protected = set(int(v) for s in mesh.open_boundaries for v in s)
    edit_log = []

    # ---- quality sites: patch re-CDT --------------------------------------
    for site in sorted(sites, key=lambda s: s["site"]):
        checks = {r["check"] for r in site["offenders"]}
        if checks == {"obc_perpendicularity"}:
            continue
        elem_ids = sorted({int(r["element"]) for r in site["offenders"]
                           if "element" in r})
        if not elem_ids:
            edit_log.append({"site": site["site"], "action": "skip",
                             "reason": "no element offenders"})
            continue
        patch_e = _neighbours_by_node(
            mesh.elements, mesh.n_nodes, np.asarray(elem_ids),
        )
        region_e = _neighbours_by_node(
            mesh.elements, mesh.n_nodes, patch_e,
        )
        before = _region_violations(mesh.nodes, mesh.elements, region_e)
        res = _try_patch_recdt(mesh, patch_e, protected)
        if res is None:
            edit_log.append({"site": site["site"], "action": "recdt",
                             "accepted": False, "reason": "no valid patch"})
            continue
        elements2, dropped, tag = res
        # Region ids shifted by the deletion — recompute on the new
        # array via the same offender nodes.
        sel_nodes = np.unique(mesh.elements[patch_e].ravel())
        mask = np.zeros(mesh.n_nodes, dtype=bool)
        mask[sel_nodes] = True
        region2 = np.where(mask[elements2].any(axis=1))[0]
        after = _region_violations(mesh.nodes, elements2, region2)
        if after < before:
            mesh = Fort14Mesh(
                title=mesh.title, nodes=mesh.nodes, depths=mesh.depths,
                elements=elements2,
                open_boundaries=[s.copy() for s in mesh.open_boundaries],
                land_boundaries=[(ib, s.copy())
                                 for ib, s in mesh.land_boundaries],
            )
            edit_log.append({"site": site["site"], "action": "recdt",
                             "variant": tag, "accepted": True,
                             "violations": [before, after],
                             "n_dropped_nodes": len(dropped)})
        else:
            edit_log.append({"site": site["site"], "action": "recdt",
                             "variant": tag, "accepted": False,
                             "violations": [before, after]})
        print(f"[81] site {site['site']}: recdt[{tag}] "
              f"{before}->{after} "
              f"{'ACCEPT' if after < before else 'rollback'}", flush=True)

    # ---- perp sites: slide the OBC node along its chord -------------------
    obc = mesh.open_boundaries[0]
    pos_in_obc = {int(v): i for i, v in enumerate(obc)}
    report0 = run_qa(mesh, name="tmp", path=SRC, max_offenders=1000)
    perp_nodes = []
    for c in report0.checks:
        if c.check_id == "obc_perpendicularity" and not c.passed:
            perp_nodes = [int(o["node"]) for o in c.offenders]
    for v in perp_nodes:
        i = pos_in_obc.get(v)
        if i is None or i == 0 or i == len(obc) - 1:
            edit_log.append({"site": f"perp:{v}", "action": "slide",
                             "accepted": False, "reason": "endpoint"})
            continue
        p_prev, p_next = mesh.nodes[obc[i - 1]], mesh.nodes[obc[i + 1]]
        ring_e = np.where((mesh.elements == v).any(axis=1))[0]
        nbr_e = _neighbours_by_node(mesh.elements, mesh.n_nodes, ring_e)
        before = _region_violations(mesh.nodes, mesh.elements, nbr_e)
        best = None
        old = mesh.nodes[v].copy()
        for t in np.linspace(0.2, 0.8, 17):
            cand = (1 - t) * p_prev + t * p_next
            mesh.nodes[v] = cand
            tri = mesh.elements[ring_e]
            _mn, _mx, twice = _tri_quality(mesh.nodes[tri])
            if (twice <= 0).any():
                continue
            dev = _perp_dev(mesh, obc, i)
            reg = _region_violations(mesh.nodes, mesh.elements, nbr_e)
            if reg <= before and (best is None or dev < best[0]):
                best = (dev, cand.copy())
        mesh.nodes[v] = old
        dev0 = _perp_dev(mesh, obc, i)
        if best is not None and best[0] < dev0:
            mesh.nodes[v] = best[1]
            edit_log.append({"site": f"perp:{v}", "action": "slide",
                             "accepted": True,
                             "deviation_deg": [dev0, best[0]]})
            print(f"[81] perp node {v}: {dev0:.1f} -> {best[0]:.1f} deg",
                  flush=True)
        else:
            edit_log.append({"site": f"perp:{v}", "action": "slide",
                             "accepted": False,
                             "deviation_deg": [dev0, None]})
            print(f"[81] perp node {v}: no on-line improvement "
                  f"({dev0:.1f} deg)", flush=True)

    # ---- finalize ----------------------------------------------------------
    mesh, cinfo = compact_nodes(mesh)
    print(f"[81] compact: {cinfo}", flush=True)
    report = run_qa(mesh, name=OUT.name, path=OUT, max_offenders=10000)
    mesh.title = "PoC 81 v4 final (per-site edited)"
    write_fort14(mesh, OUT)
    QA_JSON.write_text(json.dumps(report.to_dict(), indent=2,
                                  ensure_ascii=False), encoding="utf-8")
    LOG_JSON.write_text(json.dumps(edit_log, indent=2,
                                   ensure_ascii=False), encoding="utf-8")
    print(format_report(report, lang="ja").split("違反箇所")[0], flush=True)

    failed = {c.check_id for c in report.checks
              if c.gate and not c.skipped and not c.passed}
    if not (failed & FATAL):
        _lon, lat = _TO_LL.transform(mesh.nodes[:, 0], mesh.nodes[:, 1])
        obc2 = mesh.open_boundaries[0]
        sponge = [(int(v), 3000.0, 0.001) for v in obc2]
        written = export_fvcom_case(mesh, EXPORT_DIR, CASENAME, obc_type=1,
                                    cor=lat, sponge=sponge)
        for k, pth in written.items():
            print(f"[81] export {k}: {pth}", flush=True)
    print(f"[81] wall: {time.perf_counter() - t0:.0f} s", flush=True)
    return 0


def _perp_dev(mesh, obc, i) -> float:
    """Best-edge perpendicularity deviation (deg) at OBC node i."""
    v = int(obc[i])
    tangent = mesh.nodes[obc[min(i + 1, len(obc) - 1)]] \
        - mesh.nodes[obc[max(i - 1, 0)]]
    tangent = tangent / (np.linalg.norm(tangent) or 1.0)
    ring = np.where((mesh.elements == v).any(axis=1))[0]
    obc_set = set(int(x) for x in obc)
    best = 90.0
    for e in ring:
        for w in mesh.elements[e]:
            w = int(w)
            if w == v or w in obc_set:
                continue
            d = mesh.nodes[w] - mesh.nodes[v]
            n = np.linalg.norm(d)
            if n == 0:
                continue
            dev = abs(np.degrees(np.arcsin(
                abs(float(np.dot(d / n, tangent)))
            )))
            best = min(best, dev)
    return best


if __name__ == "__main__":
    raise SystemExit(main())
