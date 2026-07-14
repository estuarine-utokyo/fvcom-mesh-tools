"""Human-judgment MESH edits (owner 2026-07-15: "a human sees at
a glance that adding one cell on top of the red pair completes
OW05 -- can we do manual editing based on such judgment?").

The mesh-level analogue of ``recipes/edits`` (polygon edits): a
version-controlled ledger of small, per-site operations applied
AFTER the automatic finishing chain. Every operation is
coordinate-addressed (fort.14 ids change with every realization),
fully quality-gated, and FAILS LOUDLY -- no silent fallbacks.

Operations
----------
flip
    ``{"op": "flip", "edge_lonlat": [[lon,lat],[lon,lat]]}``
    Flip the interior edge whose endpoints are nearest the two
    coordinates. Gates: both cells exist, the new edge does not
    already exist, post-flip angles >= 30 deg (or strictly
    improving) and <= 130 deg, seam area ratio <= 0.5, receiving
    valences stay <= 8. Typical use: valence relief before a
    split (OW05: node 2215 sat at valence 8, blocking the choke
    transaction).

widen_split
    ``{"op": "widen_split", "edge_lonlat": [[...],[...]]}``
    Run the standard widen-then-split choke transaction
    (:func:`widen_choke_sections`) restricted to that single
    edge -- all of its gates (wall-thickness rays, angles, exact
    C4, CFL sub-triangles, valence) apply unchanged.
"""

from __future__ import annotations

from typing import Any

import numpy as np

from fvcom_mesh_tools.io import Fort14Mesh

__all__ = ["apply_mesh_edits"]

_SNAP_M = 120.0     # max node-locator distance; farther = FAIL


def _nearest_node(nodes, p, exclude=()):
    d = np.hypot(nodes[:, 0] - p[0], nodes[:, 1] - p[1])
    for j in exclude:
        d[j] = np.inf
    j = int(np.argmin(d))
    return j, float(d[j])


def _flip(mesh: Fort14Mesh, a: int, b: int) -> dict[str, Any]:
    from collections import defaultdict

    nodes, els = mesh.nodes, mesh.elements
    ee = np.vstack([els[:, [0, 1]], els[:, [1, 2]],
                    els[:, [2, 0]]])
    ee.sort(axis=1)
    edge_cells = defaultdict(list)
    for irow, (x, y) in enumerate(ee):
        edge_cells[(int(x), int(y))].append(irow % len(els))
    node_els = defaultdict(list)
    for j, t3 in enumerate(els):
        for v in t3:
            node_els[int(v)].append(j)
    key = (min(a, b), max(a, b))
    cc = edge_cells[key]
    if len(cc) != 2:
        return {"ok": False,
                "reason": f"edge ({a+1},{b+1}) is not an "
                          f"interior edge (cells={len(cc)})"}
    P = [v for v in (int(x) for x in els[cc[0]])
         if v not in (a, b)][0]
    Q = [v for v in (int(x) for x in els[cc[1]])
         if v not in (a, b)][0]
    if edge_cells.get((min(P, Q), max(P, Q))):
        return {"ok": False,
                "reason": f"new edge ({P+1},{Q+1}) exists"}
    if len(node_els[P]) >= 8 or len(node_els[Q]) >= 8:
        return {"ok": False, "reason": "receiving valence >= 8"}

    def _tri(pa, pb, pc):
        q = nodes[[pa, pb, pc]]
        ar = ((q[1, 0] - q[0, 0]) * (q[2, 1] - q[0, 1])
              - (q[2, 0] - q[0, 0]) * (q[1, 1] - q[0, 1]))
        if ar < 0:
            pa, pb, pc = pa, pc, pb
            q = nodes[[pa, pb, pc]]
            ar = -ar
        ang = []
        for k in range(3):
            u = q[(k + 1) % 3] - q[k]
            v3 = q[(k + 2) % 3] - q[k]
            c3 = np.dot(u, v3) / (np.linalg.norm(u)
                                  * np.linalg.norm(v3) + 1e-300)
            ang.append(float(np.degrees(
                np.arccos(np.clip(c3, -1, 1)))))
        return [pa, pb, pc], 0.5 * ar, ang

    def _cur_ang(ci):
        q = nodes[els[ci]]
        ang = []
        for k in range(3):
            u = q[(k + 1) % 3] - q[k]
            v3 = q[(k + 2) % 3] - q[k]
            c3 = np.dot(u, v3) / (np.linalg.norm(u)
                                  * np.linalg.norm(v3) + 1e-300)
            ang.append(float(np.degrees(
                np.arccos(np.clip(c3, -1, 1)))))
        return ang

    cur = _cur_ang(cc[0]) + _cur_ang(cc[1])
    t1, A1, g1 = _tri(P, Q, a)
    t2, A2, g2 = _tri(Q, P, b)
    if A1 <= 0 or A2 <= 0:
        return {"ok": False, "reason": "degenerate flip"}
    ang = g1 + g2
    if not ((min(ang) >= 30.0 or min(ang) > min(cur))
            and max(ang) <= 130.0):
        return {"ok": False,
                "reason": f"angles {min(ang):.1f}/{max(ang):.1f}"}
    if abs(A1 - A2) / max(A1, A2) > 0.5:
        return {"ok": False,
                "reason": f"seam C4 "
                          f"{abs(A1-A2)/max(A1,A2):.2f}"}
    els[cc[0]] = t1
    els[cc[1]] = t2
    return {"ok": True, "cells": [int(cc[0]) + 1,
                                  int(cc[1]) + 1],
            "new_edge": [P + 1, Q + 1],
            "min_angle": round(min(ang), 1)}


def apply_mesh_edits(
    mesh: Fort14Mesh,
    edits: list[dict[str, Any]],
    land_union,
    to_mesh_xy,
) -> tuple[Fort14Mesh, dict[str, Any]]:
    """Apply the ledger. ``to_mesh_xy(lon, lat) -> (x, y)``
    converts edit coordinates into the mesh CRS. Every failure is
    reported in ``info["results"]`` with its reason -- loudly, and
    the run log must show it."""
    from fvcom_mesh_tools.algorithms.obc_finish import (
        widen_choke_sections,
    )

    info: dict[str, Any] = {"applied": 0, "failed": 0,
                            "results": [], "widen_ops": []}
    for ed in edits:
        op = ed.get("op")
        pts = [to_mesh_xy(*p) for p in ed["edge_lonlat"]]
        a, da = _nearest_node(mesh.nodes, pts[0])
        b, db = _nearest_node(mesh.nodes, pts[1], exclude=(a,))
        res: dict[str, Any] = {
            "id": ed.get("id", "?"), "op": op,
            "nodes_f14": [a + 1, b + 1],
            "locator_dist_m": [round(da, 1), round(db, 1)]}
        if max(da, db) > _SNAP_M:
            res.update(ok=False,
                       reason=f"node locator missed "
                              f"({da:.0f}/{db:.0f} m > "
                              f"{_SNAP_M:.0f})")
            info["failed"] += 1
            info["results"].append(res)
            continue
        if op == "flip":
            out = _flip(mesh, a, b)
        elif op == "widen_split":
            mesh, out2 = widen_choke_sections(
                mesh, land_union,
                only_edges=[(min(a, b), max(a, b))])
            out = {"ok": out2["widened"] > 0}
            if out["ok"]:
                out.update(out2["ops"][0])
                info["widen_ops"].extend(out2["ops"])
            else:
                out["reason"] = ("widen transaction refused by "
                                 "its gates at this edge")
        else:
            out = {"ok": False, "reason": f"unknown op {op!r}"}
        res.update(out)
        info["applied" if out.get("ok") else "failed"] += 1
        info["results"].append(res)
    return mesh, info
