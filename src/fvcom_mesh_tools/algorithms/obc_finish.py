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

    def _hops(a, b, max_hops=3):
        seen = {a}
        front = {a}
        for hop in range(1, max_hops + 1):
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
    for a, b in uq[ct == 2]:
        a, b = int(a), int(b)
        if not (bnode[a] and bnode[b]) or a in obc or b in obc:
            continue
        if _hops(a, b) is not None:
            continue
        cells = edge_cells[(a, b)]
        if len(cells) != 2 or used_cells & set(cells):
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

        best = (None, -1.0)
        for frac in (0.5, 0.4, 0.6, 0.35, 0.65):
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
                best = (frac, ma)
        if best[1] < 28.0:
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


def finish_obc_mesh(
    mesh: Fort14Mesh,
    *,
    seed: int = 42,
    verify_tol_m: float = 1e-3,
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
