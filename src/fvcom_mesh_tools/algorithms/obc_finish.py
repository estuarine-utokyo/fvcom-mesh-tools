"""Open-boundary finishing operators (promoted from the Tokyo-Bay
sample-reproduction runner, 2026-07-11).

All operators are topology-local, quality-gated, deterministic and
never move a node that lies on the open boundary; positions of the
user-supplied OBC line are sacrosanct throughout.

``finish_obc_mesh`` chains them with the 59e-proven movers
(``align_open_boundary_local`` + ``phase_h_finish`` with the OBC
frozen) into the standard finishing sequence for meshes built with
a constrained OBC line.
"""

from __future__ import annotations

from typing import Any

import numpy as np

from fvcom_mesh_tools.algorithms.perpendicularity import boundary_tangents
from fvcom_mesh_tools.io import Fort14Mesh

__all__ = [
    "prune_one_wide_protected",
    "flip_for_obc_perp",
    "fix_r4",
    "split_r4_end_cells",
    "flip_c4_edges",
    "collapse_short_boundary_edges",
    "widen_choke_sections",
    "finish_obc_mesh",
]


def _tri_angles(p0, p1, p2):
    out = []
    for x, y, z in ((p0, p1, p2), (p1, p2, p0), (p2, p0, p1)):
        u, v = y - x, z - x
        c = np.dot(u, v) / (np.linalg.norm(u) * np.linalg.norm(v))
        out.append(np.degrees(np.arccos(np.clip(c, -1, 1))))
    return out


def _area(p0, p1, p2):
    return 0.5 * ((p1[0] - p0[0]) * (p2[1] - p0[1])
                  - (p1[1] - p0[1]) * (p2[0] - p0[0]))


def _boundary_mask(nodes: np.ndarray, els: np.ndarray) -> np.ndarray:
    ee = np.vstack([els[:, [0, 1]], els[:, [1, 2]], els[:, [2, 0]]])
    ee.sort(axis=1)
    uq, ct = np.unique(ee, axis=0, return_counts=True)
    bnd = np.zeros(len(nodes), bool)
    bnd[uq[ct == 1].ravel()] = True
    return bnd


def prune_one_wide_protected(
    points: np.ndarray,
    cells: np.ndarray,
    protected_pts: np.ndarray,
) -> tuple[np.ndarray, np.ndarray]:
    """Iteratively delete faces with a single face-neighbour
    (1-element-wide dead-end strips, e.g. upstream river reaches),
    except faces containing a protected (constrained) node — the
    unprotected version eats OBC band end cells and orphans the
    constrained line."""
    from collections import defaultdict

    from scipy.spatial import cKDTree

    _, pidx = cKDTree(points).query(np.asarray(protected_pts, float))
    prot = np.zeros(len(points), bool)
    prot[pidx] = True
    t = np.asarray(cells)
    while True:
        ee = np.vstack([t[:, [0, 1]], t[:, [1, 2]], t[:, [2, 0]]])
        ee.sort(axis=1)
        ef = defaultdict(list)
        for k, (a, b) in enumerate(map(tuple, ee)):
            ef[(a, b)].append(k % len(t))
        nnb = np.zeros(len(t), int)
        for fs in ef.values():
            if len(fs) == 2:
                nnb[fs[0]] += 1
                nnb[fs[1]] += 1
        kill = (nnb <= 1) & ~prot[t].any(axis=1)
        if not kill.any():
            break
        t = t[~kill]
    used = np.unique(t)
    remap = np.full(len(points), -1, dtype=np.int64)
    remap[used] = np.arange(used.size)
    return points[used], remap[t]


def flip_for_obc_perp(
    mesh: Fort14Mesh,
    dev_max: float = 20.0,
    min_angle: float = 30.0,
    max_angle: float = 130.0,
    max_area_change: float = 0.5,
) -> dict[str, Any]:
    """IN-PLACE diagonal flips at OBC nodes whose best incident
    edge is > ``dev_max`` from perpendicular: flip the opposite
    edge of a 1-ring element so its far vertex connects to the OBC
    node (a wedge split 2->3 — the sample's own corner pattern; a
    pure node move cannot reach it). Quality-gated; node positions
    never change."""
    nodes, els = mesh.nodes, mesh.elements
    seg = np.asarray(mesh.open_boundaries[0], int)
    tang = boundary_tangents(nodes[seg])
    fixed, unfixed = [], []
    for k, v in enumerate(seg):
        v = int(v)
        t_v = tang[k]

        def _dev_to(w):
            e = nodes[w] - nodes[v]
            e = e / np.linalg.norm(e)
            return abs(90.0 - np.degrees(
                np.arccos(min(1.0, abs(float(np.dot(e, t_v)))))))

        ring = np.where((els == v).any(axis=1))[0]
        devs = [_dev_to(int(w)) for ei in ring for w in els[ei]
                if int(w) != v]
        if not devs or min(devs) <= dev_max:
            continue
        best = None
        for ei in ring:
            a, b = [int(x) for x in els[ei] if int(x) != v]
            nb = [int(ej) for ej in np.where(
                ((els == a).any(axis=1)) & ((els == b).any(axis=1))
            )[0] if int(ej) != int(ei)]
            if not nb:
                continue
            ej = nb[0]
            m = int([x for x in els[ej] if int(x) not in (a, b)][0])
            if m == v:
                continue
            dev = _dev_to(m)
            if dev > dev_max - 1.0:
                continue
            t1, t2 = [v, a, m], [v, m, b]
            if _area(*nodes[t1]) < 0:
                t1 = [v, m, a]
            if _area(*nodes[t2]) < 0:
                t2 = [v, b, m]
            A1, A2 = _area(*nodes[t1]), _area(*nodes[t2])
            if A1 <= 0 or A2 <= 0:
                continue
            ang = _tri_angles(*nodes[t1]) + _tri_angles(*nodes[t2])
            if min(ang) < min_angle or max(ang) > max_angle:
                continue
            if abs(A1 - A2) / max(A1, A2) > max_area_change:
                continue
            if best is None or dev < best[0]:
                best = (dev, int(ei), int(ej), t1, t2)
        if best is None:
            unfixed.append(v)
            continue
        _, ei, ej, t1, t2 = best
        els[ei] = t1
        els[ej] = t2
        fixed.append((v, round(best[0], 1)))
    return {"fixed": fixed, "unfixed": unfixed}


def fix_r4(
    mesh: Fort14Mesh,
    min_angle: float = 30.0,
    max_angle: float = 130.0,
    regression_deg: float = 2.5,
) -> dict[str, Any]:
    """IN-PLACE repair of R4 violations (an element carrying an OBC
    edge whose third node is ALSO a boundary node — FVCOM-fatal
    ISONB sum 5): flip an internal edge so a neighbour's interior
    node replaces the boundary third-node. R4 outranks C1, so a
    flip is accepted with up to ``regression_deg`` of min-angle
    regression (polish afterwards with a moves-only pass)."""
    nodes, els = mesh.nodes, mesh.elements
    bnd = _boundary_mask(nodes, els)
    ob = set(int(v) for v in np.asarray(mesh.open_boundaries[0]))
    fixed, unfixed = [], []
    for ei in range(len(els)):
        tri = [int(x) for x in els[ei]]
        if not all(bnd[v] for v in tri):
            continue
        obe = [(a, b) for a, b in ((tri[0], tri[1]),
                                   (tri[1], tri[2]),
                                   (tri[2], tri[0]))
               if a in ob and b in ob]
        if not obe:
            continue
        o1, o2 = obe[0]
        w = [v for v in tri if v not in (o1, o2)][0]
        done = False
        for oo in (o1, o2):
            other = o2 if oo == o1 else o1
            nb = [int(ej) for ej in np.where(
                ((els == oo).any(axis=1)) & ((els == w).any(axis=1))
            )[0] if ej != ei]
            if not nb:
                continue
            ej = nb[0]
            m = int([x for x in els[ej] if int(x) not in (oo, w)][0])
            if bnd[m]:
                continue
            t1, t2 = [other, oo, m], [other, m, w]
            if _area(*nodes[t1]) < 0:
                t1 = [other, m, oo]
            if _area(*nodes[t2]) < 0:
                t2 = [other, w, m]
            A1, A2 = _area(*nodes[t1]), _area(*nodes[t2])
            if A1 <= 0 or A2 <= 0:
                continue
            ang = _tri_angles(*nodes[t1]) + _tri_angles(*nodes[t2])
            cur = _tri_angles(*nodes[tri]) + _tri_angles(
                *nodes[[int(x) for x in els[ej]]])
            if (min(ang) < min(min_angle, min(cur) - regression_deg)
                    or max(ang) > max(max_angle, max(cur) + 5.0)):
                continue
            els[ei] = t1
            els[ej] = t2
            fixed.append(ei)
            done = True
            break
        if not done:
            unfixed.append(ei)
    return {"fixed": fixed, "unfixed": unfixed}


def flip_c4_edges(
    mesh: Fort14Mesh,
    max_area_change: float = 0.5,
    min_angle: float = 30.0,
    max_angle: float = 130.0,
) -> dict[str, Any]:
    """IN-PLACE: flip internal edges whose adjacent-element area
    ratio exceeds the C4 bound — the swapped diagonal mixes the
    small and large triangle, rebalancing areas. Quality-gated."""
    nodes, els = mesh.nodes, mesh.elements
    fixed, unfixed = [], []
    ee = np.vstack([els[:, [0, 1]], els[:, [1, 2]], els[:, [2, 0]]])
    ee.sort(axis=1)
    uq, inv, ct = np.unique(ee, axis=0, return_inverse=True,
                            return_counts=True)
    for k in np.where(ct == 2)[0]:
        eids = np.where(inv == k)[0] % len(els)
        ei, ej = int(eids[0]), int(eids[1])
        A1 = abs(_area(*nodes[[int(x) for x in els[ei]]]))
        A2 = abs(_area(*nodes[[int(x) for x in els[ej]]]))
        if abs(A1 - A2) / max(A1, A2) <= max_area_change:
            continue
        a, b = [int(v) for v in uq[k]]
        m1 = int([x for x in els[ei] if int(x) not in (a, b)][0])
        m2 = int([x for x in els[ej] if int(x) not in (a, b)][0])
        t1, t2 = [m1, a, m2], [m1, m2, b]
        if _area(*nodes[t1]) < 0:
            t1 = [m1, m2, a]
        if _area(*nodes[t2]) < 0:
            t2 = [m1, b, m2]
        B1, B2 = _area(*nodes[t1]), _area(*nodes[t2])
        if B1 <= 0 or B2 <= 0:
            unfixed.append((a, b))
            continue
        ang = _tri_angles(*nodes[t1]) + _tri_angles(*nodes[t2])
        if (min(ang) < min_angle or max(ang) > max_angle
                or abs(B1 - B2) / max(B1, B2) > max_area_change):
            unfixed.append((a, b))
            continue
        els[ei] = t1
        els[ej] = t2
        fixed.append((a, b))
    return {"fixed": fixed, "unfixed": unfixed}




def split_r4_end_cells(
    mesh: Fort14Mesh, elem_ids: list[int],
) -> tuple[Fort14Mesh, dict[str, Any]]:
    """R4 cells whose only flip candidate fails the angle gate
    (single-internal-edge end cells): split that internal edge at
    the fraction (30-75%) that maximises min angle subject to the
    C4 bound over the 4 new sub-triangles (internal seams AND
    external neighbours). Deterministic; run it after all
    node-moving stages so nothing can undo the insertion."""
    import dataclasses

    nodes2, els2, dep2 = mesh.nodes, mesh.elements, mesh.depths
    ee2 = np.vstack([els2[:, [0, 1]], els2[:, [1, 2]],
                     els2[:, [2, 0]]])
    ee2.sort(axis=1)
    uq2, ct2 = np.unique(ee2, axis=0, return_counts=True)
    ob2 = set(int(v) for v in np.asarray(mesh.open_boundaries[0]))
    done, failed = [], []
    for ei in list(elem_ids):
        tri = [int(x) for x in els2[ei]]
        obe = [(a3, b3) for a3, b3 in ((tri[0], tri[1]),
                                       (tri[1], tri[2]),
                                       (tri[2], tri[0]))
               if a3 in ob2 and b3 in ob2]
        if not obe:
            failed.append(ei)
            continue
        o1, o2 = obe[0]
        w = [v for v in tri if v not in (o1, o2)][0]
        cand = None
        for oo in (o1, o2):
            nb = [int(ej) for ej in np.where(
                ((els2 == oo).any(axis=1))
                & ((els2 == w).any(axis=1)))[0] if ej != ei]
            if nb:
                cand = (oo, o2 if oo == o1 else o1, nb[0])
        if cand is None:
            failed.append(ei)
            continue
        oo, other, ej = cand
        mfar = int([x for x in els2[ej]
                    if int(x) not in (oo, w)][0])

        def _ext_area(a4, b4):
            nb4 = [int(ek) for ek in np.where(
                ((els2 == a4).any(axis=1))
                & ((els2 == b4).any(axis=1)))[0]
                if ek not in (ei, ej)]
            if not nb4:
                return None
            t4 = [int(x) for x in els2[nb4[0]]]
            return abs(_area(*nodes2[t4]))

        ext = {(other, oo): _ext_area(other, oo),
               (other, w): _ext_area(other, w),
               (oo, mfar): _ext_area(oo, mfar),
               (mfar, w): _ext_area(mfar, w)}
        best = None
        for fr in np.linspace(0.30, 0.75, 19):
            sN = (1 - fr) * nodes2[oo] + fr * nodes2[w]
            tris = [[other, oo, -1], [other, -1, w],
                    [oo, mfar, -1] if _area(
                        nodes2[oo], nodes2[mfar], sN) > 0
                    else [mfar, oo, -1],
                    [mfar, w, -1] if _area(
                        nodes2[mfar], nodes2[w], sN) > 0
                    else [w, mfar, -1]]
            angs = []
            areas4 = []
            ok = True
            for t3 in tris:
                P3 = [nodes2[v] if v >= 0 else sN for v in t3]
                A3 = _area(*P3)
                if A3 <= 0:
                    ok = False
                    break
                areas4.append(A3)
                angs += _tri_angles(*P3)
            if not ok:
                continue
            c4v = []
            for (ia, ib) in ((0, 1), (0, 2), (1, 3), (2, 3)):
                c4v.append(abs(areas4[ia] - areas4[ib])
                           / max(areas4[ia], areas4[ib]))
            for k4, key in enumerate(((other, oo), (other, w),
                                      (oo, mfar), (mfar, w))):
                Ae = ext[key]
                if Ae is not None:
                    c4v.append(abs(areas4[k4] - Ae)
                               / max(areas4[k4], Ae))
            feas = (min(angs) >= 30.0 and max(angs) <= 130.0
                    and max(c4v) <= 0.5)
            score = (1 if feas else 0, min(angs) - max(c4v))
            if best is None or score > best[0]:
                best = (score, fr, sN, tris)
        if best is None:
            failed.append(ei)
            continue
        score, fr, sN, tris = best
        si = len(nodes2)
        nodes2 = np.vstack([nodes2, sN[None, :]])
        dep2 = np.append(dep2, 0.5 * (dep2[oo] + dep2[w]))
        tt = [[si if v < 0 else v for v in t3] for t3 in tris]
        els2 = els2.copy()
        els2[ei] = tt[0]
        els2[ej] = tt[2]
        els2 = np.vstack([els2, [tt[1]], [tt[3]]])
        done.append((int(ei), round(float(fr), 2)))
        mesh = dataclasses.replace(mesh, nodes=nodes2,
                                   depths=dep2, elements=els2)
        nodes2, els2, dep2 = mesh.nodes, mesh.elements, mesh.depths
    return mesh, {"split": done, "failed": failed}


def split_choke_edges(
    mesh: Fort14Mesh,
) -> tuple[Fort14Mesh, dict[str, Any]]:
    """Split bank-to-bank CHOKE edges: an interior edge whose two
    endpoints both lie on the LAND boundary without being
    near-neighbours along it throttles a channel to a single edge
    (owner 2026-07-12, Haneda I8-a2: "the mesh must be increased
    to keep the waterway"). Insert the midpoint and split both
    adjacent cells, so the cross-section carries two cells.
    Deterministic; run before the final polish."""
    import dataclasses
    from collections import defaultdict

    nodes, els, dep = mesh.nodes, mesh.elements, mesh.depths
    n_nodes = len(nodes)
    ee = np.vstack([els[:, [0, 1]], els[:, [1, 2]],
                    els[:, [2, 0]]])
    ee.sort(axis=1)
    uq, ct = np.unique(ee, axis=0, return_counts=True)
    bnode = np.zeros(n_nodes, bool)
    badj = defaultdict(set)
    for a, b in uq[ct == 1]:
        bnode[a] = bnode[b] = True
        badj[int(a)].add(int(b))
        badj[int(b)].add(int(a))
    obc = set()
    for ob in mesh.open_boundaries:
        obc.update(int(v) for v in np.asarray(ob))

    def _hops(a, b, max_hops=3, skip_direct=False):
        """Boundary-graph hop distance a->b; ``skip_direct``
        ignores the immediate a-b edge itself (needed when the
        candidate IS a boundary edge, e.g. the terminal wedge of
        a carved corridor)."""
        seen = {a}
        first = badj[a] - ({b} if skip_direct else set())
        front = set(first) - seen
        if b in front:
            return 1
        seen |= front
        for hop in range(2, max_hops + 1):
            front = set().union(*(badj[v] for v in front)) - seen
            if b in front:
                return hop
            seen |= front
        return None

    edge_cells = defaultdict(list)
    for irow, (a, b) in enumerate(ee):
        edge_cells[(int(a), int(b))].append(irow % len(els))
    new_nodes, new_dep = [nodes], [dep]
    replaced: dict[int, list[list[int]]] = {}
    used_cells: set[int] = set()
    split_edges = []
    for (a, b), nshare in zip(uq, ct):
        a, b = int(a), int(b)
        if not (bnode[a] and bnode[b]) or a in obc or b in obc:
            continue
        cells = edge_cells[(a, b)]
        if len(cells) != int(nshare) or used_cells & set(cells):
            continue
        if nshare == 1:
            # BOUNDARY edges: excluding the edge itself from the
            # hop test qualifies nearly every coastline edge (the
            # only alternative boundary path is the whole loop --
            # 346 bogus splits, run 6186463). Restrict to the
            # TERMINAL-WEDGE signature: single all-boundary cell
            # whose spanning edge is clearly its longest side.
            tri = [int(x) for x in els[cells[0]]]
            w0 = [v for v in tri if v not in (a, b)][0]
            e_ab = float(np.hypot(*(nodes[a] - nodes[b])))
            e_aw = float(np.hypot(*(nodes[a] - nodes[w0])))
            e_bw = float(np.hypot(*(nodes[b] - nodes[w0])))
            if e_ab < 1.3 * 0.5 * (e_aw + e_bw):
                continue
        else:
            if _hops(a, b) is not None:
                continue

        # QUALITY GATE (run 6185030: blind midpoint splits made 16
        # C1 violations): search the split fraction that maximises
        # the worst angle over the 4 sub-triangles; skip the split
        # if even the best is a sliver. Skipped chokes stay in the
        # one-wide ledger.
        def _min_angle(tris, pos):
            worst = 180.0
            for t3 in tris:
                q = np.array([pos[v] for v in t3])
                for k3 in range(3):
                    u = q[(k3 + 1) % 3] - q[k3]
                    v3 = q[(k3 + 2) % 3] - q[k3]
                    c3 = np.dot(u, v3) / (
                        np.linalg.norm(u) * np.linalg.norm(v3)
                        + 1e-300)
                    worst = min(worst, np.degrees(
                        np.arccos(np.clip(c3, -1, 1))))
            return worst

        def _tri_area(t3, pos):
            q = np.array([pos[v] for v in t3])
            return 0.5 * abs(
                (q[1, 0] - q[0, 0]) * (q[2, 1] - q[0, 1])
                - (q[2, 0] - q[0, 0]) * (q[1, 1] - q[0, 1]))

        def _c4_pred(frac):
            """Worst predicted area-change ratio of the 4 sub-
            triangles against their EXTERNAL neighbours and the
            internal split seams (QA C4 definition)."""
            worst = 0.0
            for ci in cells:
                tri = [int(x) for x in els[ci]]
                w3 = [v for v in tri if v not in (a, b)][0]
                q = nodes[[tri[0], tri[1], tri[2]]]
                A = 0.5 * abs(
                    (q[1, 0] - q[0, 0]) * (q[2, 1] - q[0, 1])
                    - (q[2, 0] - q[0, 0]) * (q[1, 1] - q[0, 1]))
                subs = {(a, w3): frac * A,
                        (b, w3): (1 - frac) * A}
                worst = max(worst, abs(2 * frac - 1)
                            / max(frac, 1 - frac))
                for (v1, v2), As in subs.items():
                    lo3, hi3 = sorted((v1, v2))
                    for cj in edge_cells[(lo3, hi3)]:
                        if cj == ci:
                            continue
                        tj = [int(x) for x in els[cj]]
                        qj = nodes[[tj[0], tj[1], tj[2]]]
                        Aj = 0.5 * abs(
                            (qj[1, 0] - qj[0, 0])
                            * (qj[2, 1] - qj[0, 1])
                            - (qj[2, 0] - qj[0, 0])
                            * (qj[1, 1] - qj[0, 1]))
                        worst = max(worst, abs(As - Aj)
                                    / max(As, Aj, 1e-300))
            return worst

        # CURRENT quality of the un-split pair (a split that
        # strictly improves an already-violating site is accepted
        # even when the absolute targets stay out of reach --
        # element 4797, run 6185775: the end-wedge cell violates
        # C4 at 0.71 and the absolute-only gate refused to touch
        # it)
        cur_pos = {}
        cur_tris = []
        for ci in cells:
            tri = [int(x) for x in els[ci]]
            for v in tri:
                cur_pos[v] = nodes[v]
            cur_tris.append(tri)
        cur_ma = _min_angle(cur_tris, cur_pos)
        cur_c4 = _c4_pred(1.0 - 1e-9)   # frac->1: subs ~ originals
        best = (None, -1.0, np.inf)
        for frac in (0.5, 0.4, 0.6, 0.35, 0.65):
            c4p = _c4_pred(frac)
            mpos = (1 - frac) * nodes[a] + frac * nodes[b]
            pos = {a: nodes[a], b: nodes[b], -1: mpos}
            tris4 = []
            for ci in cells:
                tri = [int(x) for x in els[ci]]
                w3 = [v for v in tri if v not in (a, b)][0]
                pos[w3] = nodes[w3]
                tris4.append([(-1 if v == b else v)
                              for v in tri])
                tris4.append([(-1 if v == a else v)
                              for v in tri])
            ma = _min_angle(tris4, pos)
            if ma > best[1]:
                best = (frac, ma, c4p)
        ok_abs = best[1] >= 30.0 and best[2] <= 0.55
        ok_impr = (max(cur_c4, 0.0) > 0.5
                   and best[2] < cur_c4 - 0.02
                   and best[1] >= min(30.0, cur_ma - 1.0))
        if best[0] is None or not (ok_abs or ok_impr):
            continue
        frac = best[0]
        mid = n_nodes + sum(len(x) for x in new_nodes[1:])
        new_nodes.append(
            ((1 - frac) * nodes[a] + frac * nodes[b])[None, :])
        new_dep.append(np.array(
            [(1 - frac) * dep[a] + frac * dep[b]]))
        for ci in cells:
            tri = [int(x) for x in els[ci]]
            # preserve orientation: replace a->mid and b->mid
            t1 = [mid if v == b else v for v in tri]
            t2 = [mid if v == a else v for v in tri]
            replaced[ci] = [t1, t2]
            used_cells.add(ci)
        split_edges.append((a, b))
    if not split_edges:
        return mesh, {"split": 0}
    out_els = []
    for ci in range(len(els)):
        if ci in replaced:
            out_els.extend(replaced[ci])
        else:
            out_els.append([int(x) for x in els[ci]])
    mesh2 = dataclasses.replace(
        mesh,
        nodes=np.vstack(new_nodes),
        depths=np.concatenate(new_dep),
        elements=np.asarray(out_els, dtype=els.dtype))
    return mesh2, {"split": len(split_edges),
                   "edges": [(int(a), int(b))
                             for a, b in split_edges]}


def widen_choke_sections(
    mesh: Fort14Mesh,
    land_union,
    *,
    target_rows_h: float = 1.85,
    max_push_frac: float = 0.6,
    margin_frac: float = 0.4,
) -> tuple[Fort14Mesh, dict[str, Any]]:
    """Owner 2026-07-14: when a one-wide section is detected,
    WIDEN the channel locally so it carries two cells. For every
    bank-to-bank choke edge (both endpoints on the land boundary,
    no boundary path within 3 hops -- the split_choke_edges
    eligibility) the two shoreline endpoints are pushed APART
    along the edge direction, into land, until the spanning edge
    reaches ``target_rows_h`` x the local ambient edge; then the
    standard quality-gated midpoint split makes the section two
    cells across.

    Safety, per push side:
    - WALL-THICKNESS GUARD: the segment from the old position to
      the pushed position PLUS a ``margin_frac``*h margin must
      stay inside ``land_union`` (mesh CRS). A thin levee or a
      protected parallel basin blocks the push -- the site then
      stays in the one-wide ledger, never punched through.
    - per-side push capped at ``max_push_frac``*h;
    - OBC nodes are never moved;
    - all-or-nothing per edge: the split must pass the same gates
      as split_choke_edges (local min angle >= 30 deg or strictly
      improving, all affected elements CCW), else positions are
      restored.

    Every applied push is returned in ``info["ops"]`` (old/new
    node coordinates and the swept quadrilateral) so the
    connectivity comparator and the over-resolution gate classify
    the moved banks as INTENDED widening."""
    import dataclasses
    from collections import defaultdict

    import shapely

    nodes = mesh.nodes.copy()
    els = mesh.elements
    dep = mesh.depths
    n_nodes = len(nodes)
    ee = np.vstack([els[:, [0, 1]], els[:, [1, 2]],
                    els[:, [2, 0]]])
    ee.sort(axis=1)
    uq, ct = np.unique(ee, axis=0, return_counts=True)
    bnode = np.zeros(n_nodes, bool)
    badj = defaultdict(set)
    for a, b in uq[ct == 1]:
        bnode[a] = bnode[b] = True
        badj[int(a)].add(int(b))
        badj[int(b)].add(int(a))
    obc = set()
    for ob in mesh.open_boundaries:
        obc.update(int(v) for v in np.asarray(ob))

    def _hops(a, b, max_hops=3):
        seen = {a}
        front = set(badj[a]) - seen
        if b in front:
            return 1
        seen |= front
        for hop in range(2, max_hops + 1):
            front = set().union(*(badj[v] for v in front)) - seen
            if b in front:
                return hop
            seen |= front
        return None

    edge_cells = defaultdict(list)
    for irow, (a, b) in enumerate(ee):
        edge_cells[(int(a), int(b))].append(irow % len(els))
    node_els = defaultdict(list)
    for j, t3 in enumerate(els):
        for v in t3:
            node_els[int(v)].append(j)

    def _worst_angle(rows, pos):
        worst = 180.0
        for j in rows:
            q = pos[els[j]]
            for k3 in range(3):
                u = q[(k3 + 1) % 3] - q[k3]
                v3 = q[(k3 + 2) % 3] - q[k3]
                c3 = np.dot(u, v3) / (
                    np.linalg.norm(u) * np.linalg.norm(v3)
                    + 1e-300)
                worst = min(worst, float(np.degrees(
                    np.arccos(np.clip(c3, -1, 1)))))
        return worst

    def _all_ccw(rows, pos):
        for j in rows:
            q = pos[els[j]]
            if ((q[1, 0] - q[0, 0]) * (q[2, 1] - q[0, 1])
                    - (q[2, 0] - q[0, 0])
                    * (q[1, 1] - q[0, 1])) <= 1e-6:
                return False
        return True

    ops = []
    new_pts: list = []
    new_dep: list = []
    repl_rows: list = []
    touched: set[int] = set()
    for a, b in uq[ct == 2]:
        a, b = int(a), int(b)
        if not (bnode[a] and bnode[b]) or a in obc or b in obc:
            continue
        if a in touched or b in touched:
            continue
        if _hops(a, b) is not None:
            continue
        # local ambient edge scale from the two adjacent cells
        cells2 = edge_cells[(min(a, b), max(a, b))]
        h_loc = float(np.median([
            np.hypot(*(nodes[t3[(k + 1) % 3]] - nodes[t3[k]]))
            for ci in cells2 for t3 in [els[ci]]
            for k in range(3)]))
        L = float(np.hypot(*(nodes[a] - nodes[b])))
        need = target_rows_h * h_loc - L
        if need <= 0.05 * h_loc:
            continue                # long enough; split refused
            #                         for other reasons -- skip
        u = (nodes[a] - nodes[b]) / max(L, 1e-9)
        margin = margin_frac * h_loc
        cap = max_push_frac * h_loc

        def _avail(v, sgn):
            """Feasible push distance along the outward ray.
            Mesh boundary nodes float 0-110 m off the land
            polygon (run 6195753), so scan the ray in 20 m
            steps. POLYGON WATER passes freely: DistMesh retreats
            from water strips thinner than ~0.5 cell, and pushing
            the bank back onto the terrain line just reclaims
            water the shoreline already grants (owner 2026-07-14,
            OW13 Urayasu: 120-240 m unmeshed fringes made the
            land-within-170 m rule refuse exactly the corners a
            human would mesh). Land positions are allowed only
            inside the FIRST land run, keeping ``margin`` of wall
            beyond the new bank -- a second land run past a water
            gap is never entered (no spit-hopping)."""
            step = 20.0
            ts = np.arange(step, cap + margin + 190.0, step)
            pts = shapely.points(
                nodes[v][0] + sgn * u[0] * ts,
                nodes[v][1] + sgn * u[1] * ts)
            inl = shapely.covers(land_union, pts)
            if not bool(inl.any()):
                return float(cap)     # open water all the way
            i0 = int(np.argmax(inl))
            i1 = i0
            while i1 + 1 < len(ts) and inl[i1 + 1]:
                i1 += 1
            # up to the land-run end minus the wall margin; the
            # leading water fringe [0, ts[i0]) is always allowed
            d_land_ok = ts[i1] - margin
            if d_land_ok < ts[i0]:
                d_land_ok = ts[i0] - step   # water fringe only
            return float(max(0.0, min(cap, d_land_ok)))

        d_a = _avail(a, +1.0)
        d_b = _avail(b, -1.0)
        # foreign-mesh guard: a push across a water fringe must
        # stop short of any other mesh node (an island coast
        # meshed across the strip would otherwise overlap)
        from scipy.spatial import cKDTree
        if not hasattr(_avail, "_tree"):
            _avail._tree = cKDTree(nodes)
        for v0, sgn0, dv in ((a, +1.0, d_a), (b, -1.0, d_b)):
            dcur = dv
            while dcur > 0:
                p_t = nodes[v0] + sgn0 * u * dcur
                dd, jj = _avail._tree.query(p_t, k=3)
                dd = [d2 for d2, j2 in zip(np.atleast_1d(dd),
                                           np.atleast_1d(jj))
                      if j2 not in (a, b)]
                if dd and min(dd) < 0.45 * h_loc:
                    dcur -= 40.0
                    continue
                break
            if v0 == a:
                d_a = max(0.0, dcur)
            else:
                d_b = max(0.0, dcur)
        if d_a + d_b < min(need, 0.4 * h_loc):
            continue             # not enough room to matter
        frac = min(1.0, need / max(d_a + d_b, 1e-9))
        # PUSH + SPLIT evaluated TOGETHER (run 6195753: a push
        # alone flattens the two adjacent cells and fails any
        # honest angle gate -- only the midpoint split restores
        # the proportions, so the pair is one transaction). The
        # push amount is searched downward (full need, then 70%,
        # then 45%): a full push can over-stretch the ring's
        # boundary cells, while a partial widening that still
        # passes the gates beats no widening at all.
        rows = sorted(set(node_els[a]) | set(node_els[b]))
        old_worst = _worst_angle(rows, nodes)
        saved = {v: nodes[v].copy() for v in (a, b)}
        w12 = []
        for ci in cells2:
            w12.append(int([x for x in els[ci]
                            if int(x) not in (a, b)][0]))
        # a midpoint split adds +1 valence to both off-edge
        # nodes (run 6195779: C5=9)
        if any(len(node_els[w0]) >= 8 for w0 in w12):
            continue
        best = None

        def _tri_area_now(t3, pm=None):
            q = np.array([pm if v == -1 else nodes[v]
                          for v in t3])
            return 0.5 * abs(
                (q[1, 0] - q[0, 0]) * (q[2, 1] - q[0, 1])
                - (q[2, 0] - q[0, 0]) * (q[1, 1] - q[0, 1]))

        def _ring_c4_ok():
            for j in rows:
                if j in cells2:
                    continue
                Aj = _tri_area_now(els[j])
                t3 = els[j]
                for e2 in ((t3[0], t3[1]), (t3[1], t3[2]),
                           (t3[2], t3[0])):
                    lo3, hi3 = sorted((int(e2[0]), int(e2[1])))
                    for k2 in edge_cells[(lo3, hi3)]:
                        if k2 == j or k2 in cells2:
                            continue
                        Ak = _tri_area_now(els[k2])
                        if (abs(Aj - Ak)
                                / max(Aj, Ak, 1e-9)) > 0.5:
                            return False
            return True

        # GATE-FIRST dense search (runs 6197628 replay: the old
        # best-angle-then-gate selection discarded candidates
        # that pass EVERY gate at a slightly smaller angle --
        # OW23 had 31.8 deg with only sub-C4 over, while another
        # f held sub-C4 and the angle both). Largest passing
        # pscale wins (most widening), then the best angle.
        for pscale in np.arange(1.0, 0.29, -0.05):
            for v, sgn2, dv in ((a, +1.0, d_a), (b, -1.0, d_b)):
                nodes[v] = saved[v] + sgn2 * u * (
                    dv * frac * pscale)
            if not _all_ccw([j for j in rows
                             if j not in cells2], nodes):
                continue
            if not _ring_c4_ok():
                continue
            wring = _worst_angle(
                [j for j in rows if j not in cells2], nodes)
            # split-point search: along the edge (f) AND a small
            # normal offset (owner analysis 2026-07-14, OW08: all
            # on-line candidates capped at 29.6 deg -- 0.4 deg
            # under the C1 bar; an off-line midpoint redistributes
            # the sub angles). Orientation-consistency guarded.
            _nrm = np.array([-u[1], u[0]])
            for f, toff in [(float(f2), float(t2))
                            for f2 in np.arange(0.30, 0.7001,
                                                0.025)
                            for t2 in (0.0, 0.15, -0.15)]:
                if abs(2 * f - 1) / max(f, 1 - f) > 0.5:
                    continue
                pm = (nodes[a] * (1 - f) + nodes[b] * f
                      + _nrm * (toff * h_loc))
                worst_f = 180.0
                ok_f = True
                sub_c4 = 0.0
                sub_area = {}     # EXACT areas (run 6198270:
                #                   the f*A approximation is
                #                   wrong off-line -> C4 0.586)
                for ci, w0 in zip(cells2, w12):
                    qp = np.array([nodes[a], nodes[b],
                                   nodes[w0]])
                    ar_par = ((qp[1, 0] - qp[0, 0])
                              * (qp[2, 1] - qp[0, 1])
                              - (qp[2, 0] - qp[0, 0])
                              * (qp[1, 1] - qp[0, 1]))
                    for t3 in ([a, -1, w0], [-1, b, w0]):
                        q = np.array([pm if v == -1 else nodes[v]
                                      for v in t3])
                        ar = ((q[1, 0] - q[0, 0])
                              * (q[2, 1] - q[0, 1])
                              - (q[2, 0] - q[0, 0])
                              * (q[1, 1] - q[0, 1]))
                        if abs(ar) < 1e-6 or ar * ar_par < 0:
                            ok_f = False
                            break
                        As = 0.5 * abs(ar)
                        v0 = t3[0] if t3[0] != -1 else t3[1]
                        sub_area[(int(v0), int(w0))] = As
                        for k3 in range(3):
                            u3 = q[(k3 + 1) % 3] - q[k3]
                            v3 = q[(k3 + 2) % 3] - q[k3]
                            c3 = np.dot(u3, v3) / (
                                np.linalg.norm(u3)
                                * np.linalg.norm(v3) + 1e-300)
                            worst_f = min(worst_f, float(
                                np.degrees(np.arccos(
                                    np.clip(c3, -1, 1)))))
                        lo3, hi3 = sorted((int(v0), int(w0)))
                        for cj in edge_cells[(lo3, hi3)]:
                            if cj in cells2:
                                continue
                            Ajx = _tri_area_now(els[cj])
                            sub_c4 = max(
                                sub_c4, abs(As - Ajx)
                                / max(As, Ajx, 1e-9))
                    if not ok_f:
                        break
                if not ok_f:
                    continue
                # internal seam pairs, exact areas: (a,w)-(b,w)
                # across (m,w) per side, and (a,w1)-(a,w2) /
                # (b,w1)-(b,w2) across (a,m)/(m,b)
                for k1, k2 in (
                        ((a, w12[0]), (b, w12[0])),
                        ((a, w12[1]), (b, w12[1])),
                        ((a, w12[0]), (a, w12[1])),
                        ((b, w12[0]), (b, w12[1]))):
                    A1 = sub_area.get((int(k1[0]), int(k1[1])))
                    A2 = sub_area.get((int(k2[0]), int(k2[1])))
                    if A1 is None or A2 is None:
                        continue
                    sub_c4 = max(sub_c4, abs(A1 - A2)
                                 / max(A1, A2, 1e-9))
                if sub_c4 > 0.5:
                    continue
                wtot = min(worst_f, wring)
                if not (wtot >= 30.0 or wtot > old_worst):
                    continue
                if best is None or wtot > best[1]:
                    best = (f, wtot, toff)
            if best is not None:
                break
        accept = best is not None
        if not accept:
            for v, p in saved.items():
                nodes[v] = p
            continue
        f = best[0]
        pm = (nodes[a] * (1 - f) + nodes[b] * f
              + np.array([-u[1], u[0]]) * (best[2] * h_loc))
        mid = n_nodes + len(new_pts)
        new_pts.append(pm)
        new_dep.append(0.5 * float(dep[a] + dep[b]))
        for ci, w0 in zip(cells2, w12):
            repl_rows.append((ci, [[a, mid, w0], [mid, b, w0]]))
        # freeze the 2-RING (runs 6198321/6198347: a later
        # transaction 600 m away moved a node it shared with an
        # already-gated pair of this one -- C4 0.616 post hoc;
        # freezing only the quad was not enough). Every node of
        # every cell touching a/b/w1/w2 becomes immovable for
        # later transactions in this pass.
        for v2 in (a, b, w12[0], w12[1]):
            for j2 in node_els[v2]:
                touched.update(int(x) for x in els[j2])
        ops.append({
            "nodes": [int(a), int(b)],
            "old": [saved[a].tolist(), saved[b].tolist()],
            "new": [nodes[a].tolist(), nodes[b].tolist()],
            "h_loc": round(h_loc, 1),
            "split_frac": round(f, 2),
            "local_worst_deg": round(best[1], 2),
        })
    info: dict[str, Any] = {"widened": len(ops), "ops": ops,
                            "split": len(repl_rows)}
    if not ops:
        return mesh, info
    out_els = []
    repl = {ci: subs for ci, subs in repl_rows}
    for j in range(len(els)):
        if j in repl:
            for t3 in repl[j]:
                # enforce CCW
                q = np.array([nodes[v] if v < n_nodes
                              else new_pts[v - n_nodes]
                              for v in t3])
                ar = ((q[1, 0] - q[0, 0]) * (q[2, 1] - q[0, 1])
                      - (q[2, 0] - q[0, 0])
                      * (q[1, 1] - q[0, 1]))
                out_els.append(t3 if ar > 0
                               else [t3[0], t3[2], t3[1]])
        else:
            out_els.append([int(x) for x in els[j]])
    mesh2 = dataclasses.replace(
        mesh,
        nodes=np.vstack([nodes] + [np.asarray(p2)[None, :]
                                   for p2 in new_pts]),
        depths=np.concatenate([dep, np.asarray(new_dep)]),
        elements=np.asarray(out_els, dtype=els.dtype))
    return mesh2, info


def collapse_short_boundary_edges(
    mesh: Fort14Mesh,
    ratio: float = 0.5,
) -> tuple[Fort14Mesh, dict[str, Any]]:
    """Collapse a boundary edge much shorter than its flanking
    boundary edges (run 6191386, G8-c5: a 137 m carve-bank step
    between 360/398 m neighbours forces a 20 deg sliver on the
    triangle spanning it -- no node move can fix it because both
    endpoints are shoreline nodes). The lower-valence endpoint is
    merged into the other; quality-gated: every affected element
    must stay CCW, the local worst angle must improve (or reach
    30 deg), and the survivor valence stays <= 8. OBC nodes are
    never victims."""
    import dataclasses
    from collections import defaultdict

    nodes = mesh.nodes.copy()
    els = mesh.elements.copy()
    dep = mesh.depths.copy()
    ee = np.vstack([els[:, [0, 1]], els[:, [1, 2]],
                    els[:, [2, 0]]])
    ee.sort(axis=1)
    uq, ct = np.unique(ee, axis=0, return_counts=True)
    badj = defaultdict(set)
    for a, b in uq[ct == 1]:
        badj[int(a)].add(int(b))
        badj[int(b)].add(int(a))
    obc = set()
    for ob in mesh.open_boundaries:
        obc.update(int(v) for v in np.asarray(ob))
    node_els = defaultdict(list)
    for j, t3 in enumerate(els):
        for v in t3:
            node_els[int(v)].append(j)

    def _elen(a, b):
        return float(np.hypot(*(nodes[a] - nodes[b])))

    def _worst_angle(rows, pos):
        worst = 180.0
        for j in rows:
            q = pos[els[j]]
            if abs((q[1, 0] - q[0, 0]) * (q[2, 1] - q[0, 1])
                   - (q[2, 0] - q[0, 0])
                   * (q[1, 1] - q[0, 1])) < 1e-9:
                continue
            for k3 in range(3):
                u = q[(k3 + 1) % 3] - q[k3]
                v3 = q[(k3 + 2) % 3] - q[k3]
                c3 = np.dot(u, v3) / (
                    np.linalg.norm(u) * np.linalg.norm(v3)
                    + 1e-300)
                worst = min(worst, float(np.degrees(
                    np.arccos(np.clip(c3, -1, 1)))))
        return worst

    cand = []
    for a, b in uq[ct == 1]:
        a, b = int(a), int(b)
        ln = _elen(a, b)
        flank = [_elen(a, x) for x in badj[a] if x != b] \
            + [_elen(b, x) for x in badj[b] if x != a]
        if flank and ln < ratio * min(flank):
            cand.append((ln, a, b))
    cand.sort()
    touched: set[int] = set()
    collapsed = []
    dead_els: set[int] = set()
    for ln, a, b in cand:
        if a in touched or b in touched:
            continue
        pick = [(len(node_els[a]), a, b), (len(node_els[b]),
                                           b, a)]
        pick.sort()
        victim, survivor = None, None
        for _, v0, s0 in pick:
            if v0 not in obc:
                victim, survivor = v0, s0
                break
        if victim is None:
            continue
        rows = sorted(set(node_els[victim])
                      | set(node_els[survivor]) - dead_els)
        rows = [j for j in rows if j not in dead_els]
        old_worst = _worst_angle(rows, nodes)
        trial = els.copy()
        deg = []
        for j in rows:
            t3 = trial[j]
            t3[t3 == victim] = survivor
            if len(set(int(x) for x in t3)) < 3:
                deg.append(j)
        keep_rows = [j for j in rows if j not in deg]
        # CCW check on kept affected elements
        ok = True
        for j in keep_rows:
            q = nodes[trial[j]]
            ar = ((q[1, 0] - q[0, 0]) * (q[2, 1] - q[0, 1])
                  - (q[2, 0] - q[0, 0]) * (q[1, 1] - q[0, 1]))
            if ar <= 1e-6:
                ok = False
                break
        if not ok or not deg:
            continue
        # quality gate: local worst angle must improve or pass
        saved = els[rows].copy()
        els[rows] = trial[rows]
        new_worst = _worst_angle(keep_rows, nodes)
        if not (new_worst > old_worst or new_worst >= 30.0):
            els[rows] = saved
            continue
        # valence gate on the survivor
        n_srv = sum(1 for j in keep_rows
                    if survivor in els[j]) \
            + sum(1 for j in node_els[survivor]
                  if j not in rows and j not in dead_els)
        if n_srv > 8:
            els[rows] = saved
            continue
        dead_els.update(deg)
        for j in keep_rows:
            if j not in node_els[survivor]:
                node_els[survivor].append(j)
        touched.update((victim, survivor))
        collapsed.append((victim, survivor, round(ln, 1)))
    if not collapsed:
        return mesh, {"collapsed": 0, "edges": []}
    keep = np.array([j for j in range(len(els))
                     if j not in dead_els])
    els = els[keep]
    used = np.zeros(len(nodes), bool)
    used[els.ravel()] = True
    for ob in mesh.open_boundaries:
        used[np.asarray(ob, int)] = True
    remap = -np.ones(len(nodes), int)
    remap[used] = np.arange(int(used.sum()))
    els = remap[els]
    obs = [remap[np.asarray(ob, int)]
           for ob in mesh.open_boundaries]
    lbs = []
    for ibt, ids in mesh.land_boundaries:
        ids2 = remap[np.asarray(ids, int)]
        lbs.append((ibt, ids2[ids2 >= 0]))
    mesh2 = dataclasses.replace(
        mesh, nodes=nodes[used], depths=dep[used],
        elements=els, open_boundaries=obs,
        land_boundaries=lbs)
    return mesh2, {"collapsed": len(collapsed),
                   "edges": collapsed}


def split_c4_edges(
    mesh: Fort14Mesh,
) -> tuple[Fort14Mesh, dict[str, Any]]:
    """Residual C4 fixer: for a neighbour pair whose area change
    exceeds 0.5 after flips/polish, split the LARGER cell along
    its longest non-shared edge (halving it fixes the ratio;
    splitting the shared edge would preserve it). Quality-gated
    like the choke splitter; one deterministic pass."""
    import dataclasses
    from collections import defaultdict

    nodes, els, dep = mesh.nodes, mesh.elements, mesh.depths
    ee = np.vstack([els[:, [0, 1]], els[:, [1, 2]],
                    els[:, [2, 0]]])
    ee.sort(axis=1)
    edge_cells = defaultdict(list)
    for irow, (a, b) in enumerate(ee):
        edge_cells[(int(a), int(b))].append(irow % len(els))

    def _area(ci):
        q = nodes[els[ci]]
        return 0.5 * abs(
            (q[1, 0] - q[0, 0]) * (q[2, 1] - q[0, 1])
            - (q[2, 0] - q[0, 0]) * (q[1, 1] - q[0, 1]))

    def _min_angle_tris(tris, pos):
        worst = 180.0
        for t3 in tris:
            q = np.array([pos[v] for v in t3])
            for k3 in range(3):
                u = q[(k3 + 1) % 3] - q[k3]
                v3 = q[(k3 + 2) % 3] - q[k3]
                c3 = np.dot(u, v3) / (
                    np.linalg.norm(u) * np.linalg.norm(v3)
                    + 1e-300)
                worst = min(worst, np.degrees(
                    np.arccos(np.clip(c3, -1, 1))))
        return worst

    obc = set()
    for ob in mesh.open_boundaries:
        obc.update(int(v) for v in np.asarray(ob))
    new_nodes, new_dep = [nodes], [dep]
    replaced: dict[int, list[list[int]]] = {}
    used: set[int] = set()
    n_nodes = len(nodes)
    n_split = 0
    for (a, b), cells in list(edge_cells.items()):
        if len(cells) != 2:
            continue
        A0, A1 = _area(cells[0]), _area(cells[1])
        big, small = ((cells[0], cells[1]) if A0 >= A1
                      else (cells[1], cells[0]))
        Ab, As = max(A0, A1), min(A0, A1)
        if abs(Ab - As) / max(Ab, 1e-300) <= 0.5:
            continue
        if used & {big, small}:
            continue
        tri = [int(x) for x in els[big]]
        # longest edge of BIG that is not the shared edge
        cand = []
        for k in range(3):
            u2, v2 = tri[k], tri[(k + 1) % 3]
            if {u2, v2} == {a, b}:
                continue
            cand.append(((u2, v2),
                         float(np.hypot(*(nodes[u2]
                                          - nodes[v2])))))
        (u2, v2), _len = max(cand, key=lambda t: t[1])
        if u2 in obc or v2 in obc:
            continue
        w3 = [v for v in tri if v not in (u2, v2)][0]
        mid = 0.5 * (nodes[u2] + nodes[v2])
        pos = {u2: nodes[u2], v2: nodes[v2], w3: nodes[w3],
               -1: mid}
        t1 = [-1 if v == v2 else v for v in tri]
        t2 = [-1 if v == u2 else v for v in tri]
        if _min_angle_tris([t1, t2], pos) < 28.0:
            continue
        # the split must also help the OTHER side of (u2, v2)
        others = [c for c in edge_cells[tuple(sorted((u2, v2)))]
                  if c != big]
        ok = True
        for c2 in others:
            Ao = _area(c2)
            if abs(Ab / 2 - Ao) / max(Ab / 2, Ao, 1e-300) > 0.55:
                ok = False
        if not ok:
            continue
        midx = n_nodes + sum(len(x) for x in new_nodes[1:])
        new_nodes.append(mid[None, :])
        new_dep.append(np.array([0.5 * (dep[u2] + dep[v2])]))
        t1 = [midx if v == v2 else v for v in tri]
        t2 = [midx if v == u2 else v for v in tri]
        replaced[big] = [t1, t2]
        used.add(big)
        # subdivide the neighbour across (u2,v2) too, so the new
        # node stays conforming
        for c2 in others:
            tj = [int(x) for x in els[c2]]
            replaced[c2] = [[midx if v == v2 else v for v in tj],
                            [midx if v == u2 else v for v in tj]]
            used.add(c2)
        n_split += 1
    if not n_split:
        return mesh, {"split": 0}
    out_els = []
    for ci in range(len(els)):
        if ci in replaced:
            out_els.extend(replaced[ci])
        else:
            out_els.append([int(x) for x in els[ci]])
    mesh2 = dataclasses.replace(
        mesh, nodes=np.vstack(new_nodes),
        depths=np.concatenate(new_dep),
        elements=np.asarray(out_els, dtype=els.dtype))
    return mesh2, {"split": n_split}


def finish_obc_mesh(
    mesh: Fort14Mesh,
    *,
    seed: int = 42,
    verify_tol_m: float = 1e-3,
    land_union=None,
) -> tuple[Fort14Mesh, dict[str, Any]]:
    """Standard finishing chain for a mesh built with a constrained
    OBC line: perp-local moves -> phase_h (OBC frozen) -> compact ->
    perp flips -> R4 flips -> moves-only polish -> C4 flips.
    STOPS (raises) if any OBC node moved more than
    ``verify_tol_m``."""
    from fvcom_mesh_tools.algorithms.perp_local import (
        align_open_boundary_local,
    )
    from fvcom_mesh_tools.mesh_clean import compact_nodes
    from fvcom_mesh_tools.mesh_clean_phase_h import (
        _stochastic_local_fix_round,
        phase_h_finish,
    )

    info: dict[str, Any] = {}
    arc0 = mesh.nodes[np.asarray(mesh.open_boundaries[0], int)].copy()
    mesh, info["perp_local"] = align_open_boundary_local(mesh)
    mesh, hinfo = phase_h_finish(mesh, seed=seed,
                                 freeze_open_boundary=True)
    info["phase_h"] = {k: v for k, v in hinfo.items()
                       if not isinstance(v, (list, dict))}
    mesh, info["compact"] = compact_nodes(mesh)
    info["perp_flips"] = flip_for_obc_perp(mesh)
    info["r4_flips"] = fix_r4(mesh)
    mesh, hinfo2 = phase_h_finish(mesh, seed=seed,
                                  freeze_open_boundary=True)
    info["phase_h_2"] = {k: v for k, v in hinfo2.items()
                         if not isinstance(v, (list, dict))}
    mesh, info["compact_2"] = compact_nodes(mesh)
    info["r4_recheck"] = fix_r4(mesh)
    if info["r4_recheck"]["unfixed"]:
        mesh, info["r4_split"] = split_r4_end_cells(
            mesh, info["r4_recheck"]["unfixed"])
        info["polish"] = _stochastic_local_fix_round(
            mesh, np.random.default_rng(seed * 101),
            min_angle_target=30.0, max_angle_target=130.0,
            area_ratio_target=0.5, max_valence=8,
            max_tries_per_fail=500, perturbation_sigma=0.3,
            max_outer_passes=5, coastline_projector=None,
            freeze_open_boundary=True)
    info["c4_flips"] = flip_c4_edges(mesh)
    # CHOKE-EDGE split (owner 2026-07-12): a widened channel can
    # still be throttled to one bank-to-bank edge; insert the
    # midpoint so the section carries two cells, then polish.
    mesh, info["choke_split"] = split_choke_edges(mesh)
    mesh, info["c4_split"] = split_c4_edges(mesh)
    # SHORT BOUNDARY EDGES (run 6191386): carve-bank steps left
    # sub-half-cell boundary edges whose spanning triangles are
    # unfixable slivers -- collapse them before the final polish
    # smooths the merged region.
    mesh, info["bnd_collapse"] = collapse_short_boundary_edges(
        mesh)
    # FINAL single-pass stochastic polish (seeded, OBC frozen):
    # widened-corridor geometry can leave a 1-2 element C1/C4 tail
    # that the earlier passes miss (element 3894, run 6184675).
    # One pass, not a convergence loop.
    info["polish_final"] = _stochastic_local_fix_round(
        mesh, np.random.default_rng(seed * 211),
        min_angle_target=30.0, max_angle_target=130.0,
        area_ratio_target=0.5, max_valence=8,
        max_tries_per_fail=500, perturbation_sigma=0.3,
        max_outer_passes=5, coastline_projector=None,
        freeze_open_boundary=True)
    # NOTE (run 6190672): re-running split_choke_edges HERE, after
    # polish, was tried and REVERTED -- 17 late splits fixed no C1
    # (a midpoint split cannot open an acute angle AT an endpoint
    # of the split edge) and added 15 C4 violations because no
    # smoothing follows. Late chokes stay in the one-wide ledger.
    # WIDEN-then-split (owner 2026-07-14) runs HERE instead: the
    # ledger chokes only exist in the post-polish configuration
    # (mid-chain placement saw none, run 6195749), and unlike the
    # bare late split the bank push gives the section real room,
    # so the gated split makes two well-proportioned cells.
    if land_union is not None:
        mesh, info["choke_widen"] = widen_choke_sections(
            mesh, land_union)
        # second pass: the 2-ring freeze defers any transaction
        # that neighbours an accepted one (C4 post-hoc safety,
        # run 6198347); a fresh pass on the updated topology
        # picks those up with clean gates (run 6198387: G8-e4)
        if info["choke_widen"]["widened"]:
            mesh, info["choke_widen_2"] = widen_choke_sections(
                mesh, land_union)
    # FINAL perp pass: phase_h moves can re-tilt an OBC node's
    # best edge after the first alignment (node 1319, run 6184643:
    # perp_local had fixed it, phase_h re-broke it, flips could
    # not help). perp_local moves only INTERIOR neighbours, so the
    # constrained line stays fixed.
    mesh, info["perp_local_final"] = align_open_boundary_local(mesh)
    arc1 = mesh.nodes[np.asarray(mesh.open_boundaries[0], int)]
    mv = float(np.hypot(*(arc1 - arc0).T).max())
    info["obc_displacement_m"] = mv
    if mv > verify_tol_m:
        raise RuntimeError(
            f"finishing moved constrained OBC nodes (max {mv:.3f} "
            "m); the input line must stay fixed. Inspect the chain "
            "stages before trusting this mesh.")
    return mesh, info
