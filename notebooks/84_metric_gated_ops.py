"""PoC #84 — metric-gated operators, two deterministic passes.

PoC #83 cut C1 18->6 and cleared C2, but every C4 split was
rejected: one split improves 0.84 -> ~0.68 — real progress the
violation COUNT cannot see. Both operators now accept
metric-improving edits (gate_metric), and the driver makes TWO
deterministic worst-first passes (a bounded second pass, not a
convergence loop). The stubborn perpendicularity node gets a
facing-edge split to create a perpendicular neighbour.

Profile of the 39 residuals on 81_v4_final (see session log):

* 14/18 C1 offenders are boundary "needle ears" (2 boundary edges, a
  10-300 m short edge vs ~2 long ones) left by snapping/extrusion —
  remedy: :func:`collapse_edge` of the shortest edge (survivor
  projected onto the coastline, whitelist-conform);
* most C4 offenders pair a needle with its neighbour (die with it);
* the interior C1/C4 cluster (Futtsu long chords, 324-1431 m edges)
  needs grading — remedy: :func:`split_edge_pair` at the offending
  shared edge / the element's longest interior edge.

Single worst-first pass, per-edit accept-or-rollback (built into the
operators), node-id-based re-resolution (element ids shift after
collapses), full edit log. No loops.
"""

from __future__ import annotations

import json
import time
from pathlib import Path

import numpy as np
from pyproj import Transformer

from fvcom_mesh_tools.algorithms.boundary_snap import load_polylines
from fvcom_mesh_tools.algorithms.perp_local import _tri_quality
from fvcom_mesh_tools.algorithms.site_edits import (
    collapse_edge,
    split_edge_pair,
)
from fvcom_mesh_tools.io import read_fort14, write_fort14
from fvcom_mesh_tools.io.fvcom_native import export_fvcom_case
from fvcom_mesh_tools.mesh_clean import compact_nodes
from fvcom_mesh_tools.qa import format_report, run_qa

REPO = Path(__file__).resolve().parents[1]
SRC = REPO / "outputs" / "83_v4_final.14"
OSM_SHP = REPO / "outputs" / "osm_shoreline" / "osm_true_land_tokyo_bay.shp"
OUT = REPO / "outputs" / "84_v4_final.14"
QA_JSON = REPO / "outputs" / "84_v4_final_qa.json"
LOG_JSON = REPO / "outputs" / "84_edit_log.json"
EXPORT_DIR = REPO / "outputs" / "fvcom_inputs_v4"
CASENAME = "tokyo_bay_v4"

_TO_LL = Transformer.from_crs("EPSG:32654", "EPSG:4326", always_xy=True)
FATAL = {
    "node_index_valid", "ccw_all_elements", "no_isolated_elements",
    "r4_mixed_boundary", "manifold_boundary", "no_duplicate_nodes",
    "no_orphan_nodes", "no_tiny_area", "isbce2_authentic",
    "obc_on_boundary", "obc_chain_adjacency", "obc_interior_neighbor",
    "obc_ordering", "single_component", "obc_reachable", "min_depth_clip",
}


def _boundary_edge_codes(mesh):
    els = mesh.elements
    raw = np.vstack([els[:, [0, 1]], els[:, [1, 2]], els[:, [2, 0]]])
    raw.sort(axis=1)
    codes = raw[:, 0].astype(np.int64) * mesh.n_nodes + raw[:, 1]
    uniq, counts = np.unique(codes, return_counts=True)
    return set(uniq[counts == 1].tolist())


def _find_element(mesh, triple: frozenset):
    for k, tri in enumerate(mesh.elements):
        if frozenset(int(x) for x in tri) == triple:
            return k
    return None




def _one_pass(mesh, lines, protected, edit_log, pass_no):
    report0 = run_qa(mesh, name="pass", path=SRC, max_offenders=10000)
    c1_off, c4_off = [], []
    for c in report0.checks:
        if c.check_id == "c1_min_angle" and not c.passed:
            c1_off = sorted(c.offenders, key=lambda o: o["min_angle_deg"])
        if c.check_id == "c4_area_change" and not c.passed:
            c4_off = sorted(c.offenders, key=lambda o: -o["area_change"])
    print(f"[84] pass {pass_no}: C1={len(c1_off)} C4={len(c4_off)}",
          flush=True)

    triples = [frozenset(int(x) for x in mesh.elements[int(o["id"])])
               for o in c1_off]
    for o, triple in zip(c1_off, triples):
        e = _find_element(mesh, triple)
        if e is None:
            continue
        tri = [int(x) for x in mesh.elements[e]]
        P = mesh.nodes[tri]
        mn, _mx, _tw = _tri_quality(P[None, :, :])
        if mn[0] >= 30.0:
            continue
        bcodes = _boundary_edge_codes(mesh)
        pairs = [(tri[0], tri[1]), (tri[1], tri[2]), (tri[2], tri[0])]
        lens = [float(np.linalg.norm(mesh.nodes[a] - mesh.nodes[b]))
                for a, b in pairs]
        n_bedges = sum(
            1 for a, b in pairs
            if min(a, b) * mesh.n_nodes + max(a, b) in bcodes
        )
        order = np.argsort(lens)
        applied = "unresolved"
        if n_bedges >= 1:
            for k in order:
                a, b = pairs[int(k)]
                res = collapse_edge(mesh, a, b, lines=lines,
                                    protected=protected, gate_metric=True)
                if res is not None:
                    mesh, info = res
                    applied = "collapse"
                    edit_log.append({"pass": pass_no, "target": tri,
                                     "op": "collapse", **info})
                    break
        if applied == "unresolved":
            for k in order[::-1]:
                a, b = pairs[int(k)]
                res = split_edge_pair(mesh, a, b, gate_metric=True)
                if res is not None:
                    mesh, info = res
                    applied = "split"
                    edit_log.append({"pass": pass_no, "target": tri,
                                     "op": "split", **info})
                    break
        if applied == "unresolved":
            edit_log.append({"pass": pass_no, "target": tri,
                             "op": "unresolved",
                             "min_angle_deg": float(mn[0])})
        print(f"[84] p{pass_no} C1 {tri} ({o['min_angle_deg']:.1f}): "
              f"{applied}", flush=True)

    for o in c4_off:
        u, v = (int(x) for x in o["id"])
        carriers = np.where(
            (mesh.elements == u).any(axis=1)
            & (mesh.elements == v).any(axis=1)
        )[0]
        if carriers.size != 2:
            continue
        tri2 = mesh.elements[carriers]
        _mn, _mx, twice = _tri_quality(mesh.nodes[tri2])
        a1, a2 = 0.5 * np.abs(twice)
        hi = max(a1, a2)
        if hi <= 0 or (hi - min(a1, a2)) / hi <= 0.5:
            continue
        applied = "unresolved"
        res = split_edge_pair(mesh, u, v, gate_metric=True)
        if res is not None:
            mesh, info = res
            applied = "split"
            edit_log.append({"pass": pass_no, "target": [u, v],
                             "op": "split", **info})
        else:
            res = collapse_edge(mesh, u, v, lines=lines,
                                protected=protected, gate_metric=True)
            if res is not None:
                mesh, info = res
                applied = "collapse"
                edit_log.append({"pass": pass_no, "target": [u, v],
                                 "op": "collapse", **info})
        if applied == "unresolved":
            edit_log.append({"pass": pass_no, "target": [u, v],
                             "op": "unresolved",
                             "area_change": float(o["area_change"])})
        print(f"[84] p{pass_no} C4 ({u},{v}) ({o['area_change']:.2f}): "
              f"{applied}", flush=True)
    return mesh


def _fix_perp_node(mesh, edit_log):
    """Create a perpendicular neighbour for the stubborn OBC node by
    splitting the interior edge between its two flanking non-OBC
    neighbours."""
    report = run_qa(mesh, name="perp", path=SRC, max_offenders=1000)
    nodes_bad = []
    for c in report.checks:
        if c.check_id == "obc_perpendicularity" and not c.passed:
            nodes_bad = [int(o["id"]) for o in c.offenders
                         if o.get("kind") == "node"]
    obc_set = set(int(v) for s in mesh.open_boundaries for v in s)
    for v in nodes_bad:
        ring = np.where((mesh.elements == v).any(axis=1))[0]
        nbrs = sorted({int(x) for e in ring for x in mesh.elements[e]
                       if int(x) != v and int(x) not in obc_set})
        applied = "unresolved"
        for a in nbrs:
            for b in nbrs:
                if b <= a:
                    continue
                shared = np.where(
                    (mesh.elements == a).any(axis=1)
                    & (mesh.elements == b).any(axis=1)
                )[0]
                if shared.size != 2:
                    continue
                if not any(v in mesh.elements[e] for e in shared):
                    continue
                res = split_edge_pair(mesh, a, b, gate_metric=True)
                if res is not None:
                    mesh, info = res
                    applied = "facing-split"
                    edit_log.append({"target": f"perp:{v}",
                                     "op": "facing-split", **info})
                    break
            if applied != "unresolved":
                break
        edit_log.append({"target": f"perp:{v}", "op": applied})
        print(f"[84] perp node {v}: {applied}", flush=True)
    return mesh


def main() -> int:
    t0 = time.perf_counter()
    mesh = read_fort14(SRC)
    lines = load_polylines(OSM_SHP, to_crs=32654)
    protected = [int(v) for s in mesh.open_boundaries for v in s]
    print(f"[84] input: NP={mesh.n_nodes:,} NE={mesh.n_elements:,}",
          flush=True)

    edit_log = []
    for pass_no in (1, 2):
        mesh = _one_pass(mesh, lines, protected, edit_log, pass_no)

    mesh = _fix_perp_node(mesh, edit_log)

    # ---- finalize -----------------------------------------------------------
    mesh, cinfo = compact_nodes(mesh)
    print(f"[84] compact: {cinfo}", flush=True)
    report = run_qa(mesh, name=OUT.name, path=OUT, max_offenders=10000)
    mesh.title = "PoC 83 v4 final (dedicated operators)"
    write_fort14(mesh, OUT)
    QA_JSON.write_text(json.dumps(report.to_dict(), indent=2,
                                  ensure_ascii=False), encoding="utf-8")
    LOG_JSON.write_text(json.dumps(edit_log, indent=2,
                                   ensure_ascii=False), encoding="utf-8")
    print(format_report(report, lang="ja").split("違反箇所")[0], flush=True)

    ops = {}
    for r in edit_log:
        ops[r["op"]] = ops.get(r["op"], 0) + 1
    print(f"[84] ops: {ops}", flush=True)

    failed = {c.check_id for c in report.checks
              if c.gate and not c.skipped and not c.passed}
    if not (failed & FATAL):
        _lon, lat = _TO_LL.transform(mesh.nodes[:, 0], mesh.nodes[:, 1])
        obc2 = mesh.open_boundaries[0]
        sponge = [(int(v), 3000.0, 0.001) for v in obc2]
        written = export_fvcom_case(mesh, EXPORT_DIR, CASENAME, obc_type=1,
                                    cor=lat, sponge=sponge)
        for k, pth in written.items():
            print(f"[84] export {k}: {pth}", flush=True)
    print(f"[84] wall: {time.perf_counter() - t0:.0f} s", flush=True)
    return 0


if __name__ == "__main__":
    raise SystemExit(main())
