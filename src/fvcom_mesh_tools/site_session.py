"""Per-site operator pass (package version of PoC #83/#86).

One or two deterministic worst-first passes over the C1/C4 QA
offenders, applying the operators that the PoC series validated:

* C1 elements: ``equalize_pair`` on each edge first (interior-apex
  rebalancing), then ``collapse_edge`` (boundary needles) shortest
  edge first, then ``split_edge_pair`` for interior elements;
* C4 pairs: ``equalize_pair`` first, ``split_edge_pair`` /
  ``collapse_edge`` as fallbacks.

Every decision is returned in an edit log. No convergence loops.
"""

from __future__ import annotations

from typing import Any

import numpy as np

from fvcom_mesh_tools.algorithms.perp_local import _tri_quality
from fvcom_mesh_tools.algorithms.site_edits import (
    collapse_edge,
    equalize_pair,
    split_edge_pair,
)
from fvcom_mesh_tools.io import Fort14Mesh
from fvcom_mesh_tools.qa import run_qa

__all__ = ["apply_site_operators", "collapse_short_edges"]


def _boundary_codes(mesh: Fort14Mesh) -> set[int]:
    els = mesh.elements
    raw = np.vstack([els[:, [0, 1]], els[:, [1, 2]], els[:, [2, 0]]])
    raw.sort(axis=1)
    codes = raw[:, 0].astype(np.int64) * mesh.n_nodes + raw[:, 1]
    uniq, counts = np.unique(codes, return_counts=True)
    return set(uniq[counts == 1].tolist())


def _find_element(mesh: Fort14Mesh, triple: frozenset):
    for k, tri in enumerate(mesh.elements):
        if frozenset(int(x) for x in tri) == triple:
            return k
    return None


def collapse_short_edges(
    mesh: Fort14Mesh,
    lines,
    *,
    min_len: float,
    log=print,
) -> tuple[Fort14Mesh, int]:
    """Single ascending sweep collapsing edges shorter than
    ``min_len`` (gated collapse_edge; OBC nodes protected). Sub-scale
    edges are a SOLVER hazard, not just a QA item: the v5 M2 test
    went NaN in T/S with a 24 m minimum edge at hmin=300."""
    protected = [int(v) for s in mesh.open_boundaries for v in s]
    els = mesh.elements
    raw = np.vstack([els[:, [0, 1]], els[:, [1, 2]], els[:, [2, 0]]])
    raw.sort(axis=1)
    uv = np.unique(raw, axis=0)
    lens = np.linalg.norm(mesh.nodes[uv[:, 0]] - mesh.nodes[uv[:, 1]],
                          axis=1)
    order = np.argsort(lens)
    n_done = 0
    for k in order:
        if lens[k] >= min_len:
            break
        u, v = int(uv[k, 0]), int(uv[k, 1])
        if u >= mesh.n_nodes or v >= mesh.n_nodes:
            continue
        d = float(np.linalg.norm(mesh.nodes[u] - mesh.nodes[v]))
        if d >= min_len:
            continue
        res = collapse_edge(mesh, u, v, lines=lines,
                            protected=protected, gate_metric=True)
        if res is not None:
            mesh, _ = res
            n_done += 1
    log(f"[siteops] short-edge collapse (<{min_len:g} m): {n_done}")
    return mesh, n_done


def apply_site_operators(
    mesh: Fort14Mesh,
    lines,
    *,
    passes: int = 2,
    log=print,
) -> tuple[Fort14Mesh, list[dict[str, Any]]]:
    """Run ``passes`` worst-first operator passes. ``lines`` are the
    engineered-shoreline polylines in mesh coordinates (survivor
    positions of boundary collapses are projected onto them)."""
    protected = [int(v) for s in mesh.open_boundaries for v in s]
    edit_log: list[dict[str, Any]] = []

    for pass_no in range(1, passes + 1):
        report = run_qa(mesh, name=f"siteops-p{pass_no}", path=None,
                        max_offenders=10000)
        c1_off, c4_off = [], []
        for c in report.checks:
            if c.check_id == "c1_min_angle" and not c.passed:
                c1_off = sorted(c.offenders,
                                key=lambda o: o["min_angle_deg"])
            if c.check_id == "c4_area_change" and not c.passed:
                c4_off = sorted(c.offenders,
                                key=lambda o: -o["area_change"])
        log(f"[siteops] pass {pass_no}: C1={len(c1_off)} "
            f"C4={len(c4_off)}")
        if not c1_off and not c4_off:
            break

        triples = [
            frozenset(int(x) for x in mesh.elements[int(o["id"])])
            for o in c1_off
        ]
        for o, triple in zip(c1_off, triples):
            e = _find_element(mesh, triple)
            if e is None:
                continue
            tri = [int(x) for x in mesh.elements[e]]
            mn, _mx, _tw = _tri_quality(mesh.nodes[tri][None, :, :])
            if mn[0] >= 30.0:
                continue
            pairs = [(tri[0], tri[1]), (tri[1], tri[2]),
                     (tri[2], tri[0])]
            lens = [float(np.linalg.norm(mesh.nodes[x] - mesh.nodes[y]))
                    for x, y in pairs]
            bcodes = _boundary_codes(mesh)
            n_bedges = sum(
                1 for x, y in pairs
                if min(x, y) * mesh.n_nodes + max(x, y) in bcodes
            )
            applied = "unresolved"
            for x, y in pairs:
                res = equalize_pair(mesh, x, y)
                if res is not None:
                    mesh, inf = res
                    applied = "equalize"
                    edit_log.append({"pass": pass_no, "target": tri,
                                     "op": "equalize", **inf})
                    break
            if applied == "unresolved" and n_bedges >= 1:
                for k in np.argsort(lens):
                    x, y = pairs[int(k)]
                    res = collapse_edge(mesh, x, y, lines=lines,
                                        protected=protected,
                                        gate_metric=True)
                    if res is not None:
                        mesh, inf = res
                        applied = "collapse"
                        edit_log.append({"pass": pass_no,
                                         "target": tri,
                                         "op": "collapse", **inf})
                        break
            if applied == "unresolved" and n_bedges == 0:
                for k in np.argsort(lens)[::-1]:
                    x, y = pairs[int(k)]
                    res = split_edge_pair(mesh, x, y, gate_metric=True)
                    if res is not None:
                        mesh, inf = res
                        applied = "split"
                        edit_log.append({"pass": pass_no,
                                         "target": tri,
                                         "op": "split", **inf})
                        break
            if applied == "unresolved":
                edit_log.append({"pass": pass_no, "target": tri,
                                 "op": "unresolved",
                                 "min_angle_deg": float(mn[0])})

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
            res = equalize_pair(mesh, u, v)
            if res is not None:
                mesh, inf = res
                applied = "equalize"
                edit_log.append({"pass": pass_no, "target": [u, v],
                                 "op": "equalize", **inf})
            if applied == "unresolved":
                res = split_edge_pair(mesh, u, v, gate_metric=True)
                if res is not None:
                    mesh, inf = res
                    applied = "split"
                    edit_log.append({"pass": pass_no, "target": [u, v],
                                     "op": "split", **inf})
            if applied == "unresolved":
                res = collapse_edge(mesh, u, v, lines=lines,
                                    protected=protected,
                                    gate_metric=True)
                if res is not None:
                    mesh, inf = res
                    applied = "collapse"
                    edit_log.append({"pass": pass_no, "target": [u, v],
                                     "op": "collapse", **inf})
                else:
                    edit_log.append({"pass": pass_no, "target": [u, v],
                                     "op": "unresolved",
                                     "area_change":
                                         float(o["area_change"])})
    from fvcom_mesh_tools.mesh_clean import keep_components

    mesh, kinfo = keep_components(mesh)
    if kinfo.get("n_removed_elements"):
        edit_log.append({"op": "keep_components", **kinfo})
    return mesh, edit_log
