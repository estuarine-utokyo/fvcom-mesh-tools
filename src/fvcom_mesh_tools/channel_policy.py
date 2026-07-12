"""Narrow-channel policy (owner 2026-07-11): what to do with
channels resolved at one element width (``w/h < 1``).

Decision rule
-------------
Remove the flagged throat elements from the element dual graph and
look at what each throat cluster connects:

* both sides reach the MAIN water body (through-channel, e.g. the
  Keihin canal) -> **widen** (centroid fan; keeps the hydraulic
  connection at two cells across);
* the far side is a basin with at least ``min_basin_elements``
  elements (a real port, e.g. Funabashi class) -> **widen**;
* the far side is smaller than ``min_basin_elements`` or nothing
  (dead-end inlet, tiny harbour) -> **delete** the throat AND the
  small basin (the goto2023 sample meshes none of these).

Elements touching an open-boundary node are never flagged (the
width metric is documented-unreliable at open/land junctions).

Pure numpy/scipy; no oceanmesh import (license policy).
"""

from __future__ import annotations

from typing import Any

import numpy as np

from fvcom_mesh_tools.diagnostics import under_resolved_channels_flag
from fvcom_mesh_tools.io import Fort14Mesh
from fvcom_mesh_tools.mesh_clean import (
    compact_nodes,
    remove_elements,
    widen_thin_elements_at_centroid,
)

__all__ = ["resolve_narrow_channels"]


def _edge_key(a: np.ndarray, b: np.ndarray, n: int) -> np.ndarray:
    lo = np.minimum(a, b).astype(np.int64)
    hi = np.maximum(a, b).astype(np.int64)
    return lo * n + hi


def _dual_adjacency(els: np.ndarray, n_nodes: int):
    """(K,2) pairs of element ids sharing an edge."""
    ne = len(els)
    keys = np.concatenate([
        _edge_key(els[:, 0], els[:, 1], n_nodes),
        _edge_key(els[:, 1], els[:, 2], n_nodes),
        _edge_key(els[:, 2], els[:, 0], n_nodes),
    ])
    eids = np.tile(np.arange(ne), 3)
    order = np.argsort(keys, kind="stable")
    keys, eids = keys[order], eids[order]
    same = keys[1:] == keys[:-1]
    return np.column_stack([eids[:-1][same], eids[1:][same]])


def _components(pairs: np.ndarray, n: int, active: np.ndarray):
    """Connected components over ``active`` items (inactive = -1)."""
    import scipy.sparse as sp
    from scipy.sparse.csgraph import connected_components

    keep = active[pairs[:, 0]] & active[pairs[:, 1]]
    pr = pairs[keep]
    g = sp.coo_matrix(
        (np.ones(len(pr)), (pr[:, 0], pr[:, 1])), shape=(n, n))
    _, lab = connected_components(g + g.T, directed=False)
    lab = lab.copy()
    lab[~active] = -1
    # renumber active labels densely
    uniq = np.unique(lab[active])
    remap = {int(u): k for k, u in enumerate(uniq)}
    out = np.full(n, -1, dtype=np.int64)
    for i in np.where(active)[0]:
        out[i] = remap[int(lab[i])]
    return out, len(uniq)


def _boundary_loops(els: np.ndarray, n_nodes: int) -> list[np.ndarray]:
    from collections import defaultdict

    keys = np.concatenate([
        _edge_key(els[:, 0], els[:, 1], n_nodes),
        _edge_key(els[:, 1], els[:, 2], n_nodes),
        _edge_key(els[:, 2], els[:, 0], n_nodes),
    ])
    uq, ct = np.unique(keys, return_counts=True)
    bkeys = set(uq[ct == 1].tolist())
    nxt = defaultdict(list)
    unused = set()
    for a, b in ((0, 1), (1, 2), (2, 0)):
        for x, y in zip(els[:, a], els[:, b]):
            k = _edge_key(np.asarray(x), np.asarray(y), n_nodes)
            if int(k) in bkeys:
                nxt[int(x)].append(int(y))
                nxt[int(y)].append(int(x))
                unused.add((min(int(x), int(y)), max(int(x), int(y))))
    loops = []
    while unused:
        a, b = next(iter(unused))
        unused.discard((a, b))
        loop = [a, b]
        while True:
            cur, prev = loop[-1], loop[-2]
            cands = [w for w in nxt[cur] if w != prev
                     and (min(cur, w), max(cur, w)) in unused]
            if not cands:
                break
            w = cands[0]
            unused.discard((min(cur, w), max(cur, w)))
            if w == loop[0]:
                break
            loop.append(w)
        loops.append(np.asarray(loop, dtype=np.int64))
    return loops


def _fix_pinch_nodes(els: np.ndarray) -> tuple[np.ndarray, int]:
    """Boundary-manifold repair: a manifold boundary node has
    exactly 2 boundary edges. Cap deletions can leave two fans
    touching at a single node (pinch); delete every fan except the
    largest until the boundary is manifold again."""
    from collections import defaultdict

    n_del = 0
    while True:
        ee = np.vstack([els[:, [0, 1]], els[:, [1, 2]],
                        els[:, [2, 0]]])
        ee.sort(axis=1)
        uq, ct = np.unique(ee, axis=0, return_counts=True)
        bdeg = defaultdict(int)
        for a, b in uq[ct == 1]:
            bdeg[int(a)] += 1
            bdeg[int(b)] += 1
        pinch = [v for v, d in bdeg.items() if d > 2]
        if not pinch:
            return els, n_del
        v = pinch[0]
        inc = np.where((els == v).any(axis=1))[0]
        # split the incident elements into fans edge-connected
        # through v
        parent = {int(e): int(e) for e in inc}

        def find(x):
            while parent[x] != x:
                parent[x] = parent[parent[x]]
                x = parent[x]
            return x

        for i2 in range(len(inc)):
            for j2 in range(i2 + 1, len(inc)):
                sh = set(els[inc[i2]]) & set(els[inc[j2]])
                if len(sh) == 2 and v in sh:
                    parent[find(int(inc[i2]))] = find(int(inc[j2]))
        fans = defaultdict(list)
        for e in inc:
            fans[find(int(e))].append(int(e))
        keep_fan = max(fans.values(), key=len)
        kill = [e for fan in fans.values() if fan is not keep_fan
                for e in fan]
        els = np.delete(els, kill, axis=0)
        n_del += len(kill)


def resolve_narrow_channels(
    mesh: Fort14Mesh,
    *,
    min_basin_elements: int = 6,
    min_w_h: float = 1.0,
    coords: str = "metric",
    analyze_only: bool = False,
    apply_widen: bool = True,
    small_cluster_delete: int = 0,
) -> tuple[Fort14Mesh, dict[str, Any]]:
    """Apply the owner's narrow-channel policy. Returns
    ``(new_mesh, info)``; raises if the open boundary cannot be
    preserved through a deletion. With ``analyze_only`` the mesh is
    returned unchanged and ``info["clusters"]`` carries the decision
    per cluster plus centroids/widths -- callers use this to lay
    REFINEMENT corridors (width/2) for the widen clusters and
    REGENERATE, because a true 2-cell cross-section in a narrow
    canal needs local h ~ width/2 (centroid fans violate C1 there
    and get undone by quality finishing). For the same reason pass
    ``apply_widen=False`` in a FINISHING context (after the two-pass
    refinement): widen clusters are then only reported.
    ``small_cluster_delete``: flagged clusters with at most this
    many members are DELETED regardless of classification -- after
    the geometry-stage policy every real channel is >= 2 cells
    wide, so tiny leftovers are junction corner caps (the bridging
    triangle DistMesh drops where a widened passage opens into
    wide water); nibbling them off is the sample's own look."""
    flag, chinfo = under_resolved_channels_flag(
        mesh, min_w_h=min_w_h, coords=coords)
    ob_nodes = set(
        int(v) for s in mesh.open_boundaries
        for v in np.asarray(s).ravel())
    touches_ob = np.isin(mesh.elements,
                         np.asarray(sorted(ob_nodes))).any(axis=1)
    flag = flag & ~touches_ob
    info: dict[str, Any] = {"n_flagged": int(flag.sum())}
    if not flag.any():
        info.update(n_widened=0, n_deleted_elements=0, clusters=[])
        return mesh, info

    els = mesh.elements
    n_el = len(els)
    pairs = _dual_adjacency(els, mesh.n_nodes)
    # water-body components with throats removed
    lab, _ = _components(pairs, n_el, ~flag)
    sizes = np.bincount(lab[lab >= 0])
    # main = the component holding the open boundary
    ob_elems = np.isin(els, np.asarray(sorted(ob_nodes))).any(axis=1)
    main_labels = lab[ob_elems & (lab >= 0)]
    if main_labels.size == 0:
        raise RuntimeError(
            "no unflagged element touches the open boundary -- "
            "cannot identify the main water body")
    main = int(np.bincount(main_labels).argmax())

    # throat clusters
    clab, n_cl = _components(pairs, n_el, flag)
    delete = np.zeros(n_el, dtype=bool)
    widen = np.zeros(n_el, dtype=bool)
    clusters = []
    for c in range(n_cl):
        members = np.where(clab == c)[0]
        nbr = set()
        for pa, pb in pairs:
            if clab[pa] == c and lab[pb] >= 0:
                nbr.add(int(lab[pb]))
            elif clab[pb] == c and lab[pa] >= 0:
                nbr.add(int(lab[pa]))
        nonmain = [l for l in nbr if l != main]
        big = [l for l in nonmain if sizes[l] >= min_basin_elements]
        small = [l for l in nonmain if sizes[l] < min_basin_elements]
        if (0 < small_cluster_delete
                and len(members) <= small_cluster_delete):
            # junction corner caps: after the geometry-stage
            # policy every real channel is >= 2 cells wide, so
            # tiny leftovers are bridging triangles at passage
            # mouths -- nibble them off (sample look)
            action = "delete"
            delete[members] = True
            for l in small:
                delete[lab == l] = True
        elif (main in nbr and not nonmain) or big:
            action = "widen"
            widen[members] = True
        else:
            action = "delete"
            delete[members] = True
            for l in small:
                delete[lab == l] = True
        wvals = chinfo["channel_width_m"][members]
        wvals = wvals[np.isfinite(wvals)]
        clusters.append({
            "n_members": int(len(members)),
            "action": action,
            "neighbor_sizes": sorted(
                int(sizes[l]) for l in nonmain),
            "centroids": mesh.nodes[els[members]].mean(axis=1),
            "width_m": float(np.median(wvals)) if wvals.size
            else float("nan"),
        })
    info["clusters"] = clusters
    if analyze_only:
        info.update(n_widened=0, n_deleted_elements=0)
        return mesh, info

    if delete.any():
        arc_xy = {i: mesh.nodes[np.asarray(s, int)].copy()
                  for i, s in enumerate(mesh.open_boundaries)}
        land_ibs = [ib for ib, ids in mesh.land_boundaries]
        island_ib = None
        land_ib = None
        for ib, ids in mesh.land_boundaries:
            ids = np.asarray(ids)
            if len(ids) > 2 and ids[0] == ids[-1]:
                island_ib = ib
            else:
                land_ib = ib
        land_ib = land_ib if land_ib is not None else (
            land_ibs[0] if land_ibs else 20)
        island_ib = island_ib if island_ib is not None else land_ib

        keep = ~delete
        widen = widen[keep]
        mesh = remove_elements(mesh, keep)
        mesh, _ = compact_nodes(mesh)
        # cap deletions can pinch the boundary (two fans touching
        # at one node -- FVCOM-fatal); repair before rebuilding
        els_fixed, n_pinch_del = _fix_pinch_nodes(mesh.elements)
        if n_pinch_del:
            import dataclasses
            keep2 = np.ones(len(mesh.elements), dtype=bool)
            # _fix_pinch_nodes returns the surviving element array;
            # rebuild the widen mask by matching rows
            mesh = dataclasses.replace(mesh, elements=els_fixed)
            mesh, _ = compact_nodes(mesh)
            info["n_pinch_elements_deleted"] = int(n_pinch_del)
            widen = np.zeros(len(mesh.elements), dtype=bool)
        # drop any fragments disconnected from the main body
        pairs2 = _dual_adjacency(mesh.elements, mesh.n_nodes)
        lab2, n2 = _components(
            pairs2, len(mesh.elements),
            np.ones(len(mesh.elements), dtype=bool))
        if n2 > 1:
            from scipy.spatial import cKDTree
            _, ob_idx = cKDTree(mesh.nodes).query(arc_xy[0])
            ob_el = np.isin(mesh.elements, ob_idx).any(axis=1)
            mainc = int(np.bincount(lab2[ob_el]).argmax())
            frag = lab2 != mainc
            info["n_disconnected_dropped"] = int(frag.sum())
            widen = widen[~frag]
            mesh = remove_elements(mesh, ~frag)
            mesh, _ = compact_nodes(mesh)
        # rebuild boundaries: open strings re-resolved by COORDS
        # (deletion never touches them -- verified), land/islands
        # re-derived from the surviving topology
        from scipy.spatial import cKDTree
        opens = []
        for i in sorted(arc_xy):
            d, idx = cKDTree(mesh.nodes).query(arc_xy[i])
            if (d > 1e-6).any():
                raise RuntimeError(
                    "narrow-channel deletion moved/removed open-"
                    f"boundary nodes (max offset {d.max():.3g}); "
                    "this must never happen -- inspect the flags")
            opens.append(idx.astype(np.int64))
        loops = _boundary_loops(mesh.elements, mesh.n_nodes)

        def _loop_area(lp):
            x, y = mesh.nodes[lp, 0], mesh.nodes[lp, 1]
            return 0.5 * float(np.dot(x, np.roll(y, -1))
                               - np.dot(y, np.roll(x, -1)))

        outer = max(loops, key=lambda lp: abs(_loop_area(lp)))
        ob0 = opens[0]
        start = int(np.where(outer == ob0[0])[0][0])
        ring = np.roll(outer, -start)
        if ring[1] not in set(ob0.tolist()):
            ring = np.roll(ring[::-1], 1)
        stop = int(np.where(ring == ob0[-1])[0][0])
        land = np.append(ring[stop:], ring[0])
        lands = [(land_ib, land)]
        for lp in loops:
            if lp is not outer:
                lands.append((island_ib,
                              np.append(lp, lp[0])))
        import dataclasses
        mesh = dataclasses.replace(
            mesh, open_boundaries=opens, land_boundaries=lands)
    info["n_deleted_elements"] = int(delete.sum())

    if widen.any() and apply_widen:
        mesh, winfo = widen_thin_elements_at_centroid(mesh, widen)
        info["n_widened"] = int(widen.sum())
        info["widen_info"] = {k: v for k, v in winfo.items()
                              if isinstance(v, (int, float))}
    else:
        info["n_widened"] = 0
        info["n_widen_skipped"] = int(widen.sum())
    return mesh, info
