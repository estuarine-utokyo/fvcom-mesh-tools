"""Targeted per-site edit operators (the AI manual-editing toolbox).

These are the SMS-editor moves the session workflow applies to one
site at a time, honouring the boundary-edit whitelist: a node on the
coastline / open-boundary line may only move ALONG that line or be
deleted; new boundary nodes are created ON the line.

Currently:

* :func:`insert_node_on_line` — split a boundary edge, placing the
  new node exactly on the reference polyline. This creates the
  spacing slack that pure node motion lacks (PoC #72/#73: ~half of
  the nearshore nodes cannot conform without it).
* :func:`extrude_boundary_strip` — append a row of well-shaped
  triangles between a retreated boundary chain and the coastline,
  new nodes exactly ON the line. Reverses the post-generation
  cleanup retreat (PoC #77/#78: the boundary sits ~0.7 h inside the
  SDF zero line at every resolution, so no amount of node motion
  can conform a true-coarse mesh).
"""

from __future__ import annotations

from typing import Any

import numpy as np

from fvcom_mesh_tools.io import Fort14Mesh


def _signed_area(p0, p1, p2) -> float:
    return 0.5 * float(
        (p1[0] - p0[0]) * (p2[1] - p0[1])
        - (p1[1] - p0[1]) * (p2[0] - p0[0])
    )


def insert_node_on_line(
    mesh: Fort14Mesh,
    u: int,
    v: int,
    line,
    *,
    min_angle: float = 30.0,
    max_angle: float = 130.0,
) -> tuple[Fort14Mesh, dict[str, Any]] | None:
    """Split boundary edge ``(u, v)`` by a new node ON ``line``.

    The new node is the projection of the edge midpoint onto the
    shapely ``line``. The single triangle carrying the boundary edge
    is split in two; the land-boundary segment listing ``u, v``
    consecutively gains the new node between them. Returns ``None``
    (mesh untouched) when the edge is not a boundary edge, the split
    would flip/degenerate a triangle, or the patch would end up with
    more C1/C2 violations than before.
    """
    import shapely

    from fvcom_mesh_tools.algorithms.perp_local import _tri_quality

    u, v = int(u), int(v)
    elements = mesh.elements
    # The boundary edge belongs to exactly one triangle.
    has_u = (elements == u).any(axis=1)
    has_v = (elements == v).any(axis=1)
    carriers = np.where(has_u & has_v)[0]
    if carriers.size != 1:
        return None
    eid = int(carriers[0])
    tri = [int(x) for x in elements[eid]]

    mid = 0.5 * (mesh.nodes[u] + mesh.nodes[v])
    q = line.interpolate(line.project(shapely.Point(mid[0], mid[1])))
    p_new = np.array([q.x, q.y])

    # Keep each half's orientation identical to the original triangle:
    # replace one endpoint by the new node in the original vertex
    # order, so CCW is preserved by construction.
    def _replace(seq, a, b):
        return [b if x == a else x for x in seq]

    n_new = mesh.n_nodes
    tri_a = _replace(tri, v, n_new)
    tri_b = _replace(tri, u, n_new)
    nodes2 = np.vstack([mesh.nodes, p_new[None, :]])
    for t3 in (tri_a, tri_b):
        pts = [nodes2[t3[0]], nodes2[t3[1]], nodes2[t3[2]]]
        if _signed_area(*pts) <= 0:
            return None

    # Patch quality must not get worse (C1/C2 on the carrier before
    # vs the two halves after).
    before_tri = np.asarray([tri])
    mn_b, mx_b, _tw = _tri_quality(mesh.nodes[before_tri])
    before_bad = int((mn_b < min_angle).sum() + (mx_b > max_angle).sum())
    after_tri = np.asarray([tri_a, tri_b])
    mn_a, mx_a, _tw = _tri_quality(nodes2[after_tri])
    after_bad = int((mn_a < min_angle).sum() + (mx_a > max_angle).sum())
    if after_bad > before_bad:
        return None

    elements2 = np.vstack([elements, np.asarray([tri_b])])
    elements2[eid] = tri_a

    # Insert the new node between u and v in the land boundary that
    # lists them consecutively (either direction).
    land2 = []
    inserted = False
    for ib, seg in mesh.land_boundaries:
        seg = np.asarray(seg)
        if not inserted:
            for k in range(len(seg) - 1):
                a, b = int(seg[k]), int(seg[k + 1])
                if {a, b} == {u, v}:
                    seg = np.concatenate(
                        [seg[:k + 1], [n_new], seg[k + 1:]]
                    )
                    inserted = True
                    break
        land2.append((ib, seg))

    depths2 = np.concatenate([
        mesh.depths, [0.5 * (mesh.depths[u] + mesh.depths[v])],
    ])
    out = Fort14Mesh(
        title=mesh.title,
        nodes=nodes2,
        depths=depths2,
        elements=elements2,
        open_boundaries=[s.copy() for s in mesh.open_boundaries],
        land_boundaries=land2,
    )
    return out, {
        "new_node": int(n_new),
        "carrier_element": eid,
        "boundary_updated": inserted,
    }


def extrude_boundary_strip(
    mesh: Fort14Mesh,
    polylines,
    *,
    d_lo_frac: float = 0.25,
    d_hi_frac: float = 1.6,
    min_run: int = 3,
    min_angle: float = 30.0,
    max_angle: float = 130.0,
    exclude_nodes=None,
) -> tuple[Fort14Mesh, dict[str, Any]]:
    """Fill the cleanup-retreat gap: for boundary chains lying between
    ``d_lo_frac`` and ``d_hi_frac`` local-edge-lengths from a coastline
    polyline, add one row of triangles whose outer nodes sit exactly
    ON that line.

    Per run of consecutive eligible boundary nodes (same nearest line,
    monotonic footpoints, length >= ``min_run``), every quad
    ``(v_i, v_{i+1}, N_{i+1}, N_i)`` is split along the better
    diagonal; a quad whose best split violates the C1/C2 gates ends
    the run segment (the strip is trimmed, never forced). New nodes
    inherit the parent node's depth. Boundary lists are NOT updated —
    rebuild them afterwards (the session re-derives arcs anyway).
    """
    import shapely
    from shapely.strtree import STRtree

    from fvcom_mesh_tools.algorithms.perp_local import _tri_quality

    nodes = mesh.nodes
    elements = mesh.elements
    n_nodes = mesh.n_nodes

    raw = np.vstack([
        elements[:, [0, 1]], elements[:, [1, 2]], elements[:, [2, 0]],
    ])
    raw.sort(axis=1)
    codes = raw[:, 0].astype(np.int64) * n_nodes + raw[:, 1]
    uniq, counts = np.unique(codes, return_counts=True)
    b = uniq[counts == 1]
    boundary_uv = np.column_stack([b // n_nodes, b % n_nodes])

    excl = set(int(v) for s in mesh.open_boundaries for v in s)
    if exclude_nodes is not None:
        excl.update(int(v) for v in exclude_nodes)

    elen = np.linalg.norm(
        nodes[boundary_uv[:, 0]] - nodes[boundary_uv[:, 1]], axis=1,
    )
    lsum = np.zeros(n_nodes)
    lcnt = np.zeros(n_nodes)
    np.add.at(lsum, boundary_uv[:, 0], elen)
    np.add.at(lsum, boundary_uv[:, 1], elen)
    np.add.at(lcnt, boundary_uv[:, 0], 1)
    np.add.at(lcnt, boundary_uv[:, 1], 1)
    local_h = np.divide(lsum, np.maximum(lcnt, 1))

    bnodes = np.unique(boundary_uv.ravel())
    tree = STRtree(list(polylines))
    geoms = np.array(list(polylines), dtype=object)
    pts = shapely.points(nodes[bnodes, 0], nodes[bnodes, 1])
    nearest_idx = tree.nearest(pts)
    d = shapely.distance(pts, geoms[nearest_idx])

    eligible: dict[int, int] = {}
    for v, gi, dist in zip(bnodes, nearest_idx, d):
        v = int(v)
        if v in excl:
            continue
        if d_lo_frac * local_h[v] <= dist <= d_hi_frac * local_h[v]:
            eligible[v] = int(gi)

    # Consecutive runs along the boundary with a common line.
    nbr: dict[int, list[int]] = {}
    for u, v in boundary_uv:
        u, v = int(u), int(v)
        if u in eligible and v in eligible and eligible[u] == eligible[v]:
            nbr.setdefault(u, []).append(v)
            nbr.setdefault(v, []).append(u)
    runs: list[list[int]] = []
    seen: set[int] = set()
    starts = sorted(v for v in nbr if len(nbr[v]) == 1)
    for v0 in starts:
        if v0 in seen:
            continue
        run = [v0]
        seen.add(v0)
        prev, cur = None, v0
        while True:
            nxt = [w for w in nbr[cur] if w != prev and w not in seen]
            if not nxt:
                break
            prev, cur = cur, nxt[0]
            run.append(cur)
            seen.add(cur)
        runs.append(run)

    new_nodes: list[np.ndarray] = []
    new_depths: list[float] = []
    new_tris: list[list[int]] = []
    n_strips = n_trimmed = 0

    def _tri_ccw(a_id, b_id, c_id, coords):
        pa, pb, pc = coords[a_id], coords[b_id], coords[c_id]
        if _signed_area(pa, pb, pc) > 0:
            return [a_id, b_id, c_id]
        return [a_id, c_id, b_id]

    def _tri_ok(tri, coords):
        arr = np.asarray([[coords[i] for i in tri]])
        mn, mx, twice = _tri_quality(arr)
        return bool(
            twice[0] > 0 and mn[0] >= min_angle and mx[0] <= max_angle
        )

    for run in runs:
        if len(run) < 2:
            continue
        line = geoms[eligible[run[0]]]
        s_vals = [float(line.project(shapely.Point(*nodes[v])))
                  for v in run]
        foot = [line.interpolate(s) for s in s_vals]
        coords = {v: nodes[v] for v in run}

        # Per-quad evaluation; a defective pair CLOSES the current
        # sub-segment instead of killing the whole run (the same
        # narrowing lesson as the bisecting snap, PoC #71 -> #72).
        seg_quads: list[tuple[int, list[list[int]] | None]] = []
        direction = 0.0
        for k in range(len(run) - 1):
            ds = s_vals[k + 1] - s_vals[k]
            good = abs(ds) >= 0.2 * local_h[run[k]]
            if good and direction != 0.0 and ds * direction < 0:
                good = False
            pick = None
            if good:
                if direction == 0.0:
                    direction = np.sign(ds)
                vi, vj = run[k], run[k + 1]
                ni, nj = -(k + 1), -(k + 2)  # placeholder ids
                coords[ni] = np.array([foot[k].x, foot[k].y])
                coords[nj] = np.array([foot[k + 1].x, foot[k + 1].y])
                cand_a = [_tri_ccw(vi, vj, nj, coords),
                          _tri_ccw(vi, nj, ni, coords)]
                cand_b = [_tri_ccw(vi, vj, ni, coords),
                          _tri_ccw(vj, nj, ni, coords)]
                for cand in (cand_a, cand_b):
                    if all(_tri_ok(c, coords) for c in cand):
                        pick = cand
                        break
            if pick is None:
                direction = 0.0
                n_trimmed += 1
            seg_quads.append((k, pick))

        # Emit contiguous stretches of >= (min_run - 1) good quads.
        stretch: list[tuple[int, list[list[int]]]] = []

        def _flush():
            nonlocal n_strips
            if len(stretch) < max(1, min_run - 1):
                stretch.clear()
                return
            ks = [k for k, _q in stretch]
            node_ks = list(range(ks[0], ks[-1] + 2))
            base = n_nodes + len(new_nodes)
            id_of = {kk: base + i for i, kk in enumerate(node_ks)}
            for kk in node_ks:
                new_nodes.append(
                    np.array([foot[kk].x, foot[kk].y])
                )
                new_depths.append(float(mesh.depths[run[kk]]))
            for k, quad in stretch:
                for tri in quad:
                    remap = [
                        id_of[k] if i == -(k + 1)
                        else id_of[k + 1] if i == -(k + 2)
                        else i
                        for i in tri
                    ]
                    new_tris.append(remap)
            n_strips += 1
            stretch.clear()

        for k, pick in seg_quads:
            if pick is None:
                _flush()
            else:
                stretch.append((k, pick))
        _flush()

    if not new_nodes:
        return mesh, {"n_strips": 0, "n_new_nodes": 0,
                      "n_new_elements": 0, "n_trimmed": n_trimmed}

    out = Fort14Mesh(
        title=mesh.title,
        nodes=np.vstack([nodes, np.asarray(new_nodes)]),
        depths=np.concatenate([mesh.depths, np.asarray(new_depths)]),
        elements=np.vstack([elements, np.asarray(new_tris)]),
        open_boundaries=[s.copy() for s in mesh.open_boundaries],
        land_boundaries=[(ib, s.copy())
                         for ib, s in mesh.land_boundaries],
    )
    return out, {
        "n_strips": n_strips,
        "n_new_nodes": len(new_nodes),
        "n_new_elements": len(new_tris),
        "n_trimmed": n_trimmed,
    }


def collapse_edge(
    mesh: Fort14Mesh,
    u: int,
    v: int,
    *,
    lines=None,
    min_angle: float = 30.0,
    max_angle: float = 130.0,
    protected=None,
) -> tuple[Fort14Mesh, dict[str, Any]] | None:
    """Collapse edge ``(u, v)`` into ``u`` (the sliver-needle remedy).

    Position rules honour the boundary whitelist: if both endpoints
    are mesh-boundary nodes and ``lines`` is given, the survivor
    moves to the midpoint PROJECTED onto the nearest line (it stays
    on the coastline); if exactly one endpoint is a boundary node,
    that node survives in place; two interior nodes meet at the
    midpoint. Rejected (returns ``None``) when an endpoint is
    ``protected`` (OBC), when a flip would occur, or when the 1-ring
    violation count (C1/C2/flips) does not decrease.
    """
    import shapely

    from fvcom_mesh_tools.algorithms.perp_local import _tri_quality

    u, v = int(u), int(v)
    protected = set(protected or ())
    if u in protected or v in protected:
        return None
    elements = mesh.elements
    n_nodes = mesh.n_nodes

    raw = np.vstack([
        elements[:, [0, 1]], elements[:, [1, 2]], elements[:, [2, 0]],
    ])
    raw.sort(axis=1)
    codes = raw[:, 0].astype(np.int64) * n_nodes + raw[:, 1]
    uniq, counts = np.unique(codes, return_counts=True)
    bnd_nodes = set(
        np.column_stack([
            uniq[counts == 1] // n_nodes, uniq[counts == 1] % n_nodes,
        ]).ravel().tolist()
    )

    u_b, v_b = u in bnd_nodes, v in bnd_nodes
    if u_b and not v_b:
        keep, drop, pos = u, v, mesh.nodes[u].copy()
    elif v_b and not u_b:
        keep, drop, pos = v, u, mesh.nodes[v].copy()
    else:
        keep, drop = u, v
        pos = 0.5 * (mesh.nodes[u] + mesh.nodes[v])
        if u_b and v_b and lines is not None:
            from shapely.strtree import STRtree

            tree = STRtree(list(lines))
            pt = shapely.Point(pos[0], pos[1])
            line = list(lines)[int(tree.nearest(pt))]
            q = line.interpolate(line.project(pt))
            pos = np.array([q.x, q.y])

    ring = np.where((elements == keep).any(axis=1)
                    | (elements == drop).any(axis=1))[0]

    def _fails(nodes_arr, elems_arr, region):
        tri = elems_arr[region]
        mn, mx, twice = _tri_quality(nodes_arr[tri])
        return int((mn < min_angle).sum() + (mx > max_angle).sum()
                   + (twice <= 0).sum())

    before = _fails(mesh.nodes, elements, ring)

    nodes2 = mesh.nodes.copy()
    nodes2[keep] = pos
    elements2 = np.where(elements == drop, keep, elements)
    degen = (
        (elements2[:, 0] == elements2[:, 1])
        | (elements2[:, 1] == elements2[:, 2])
        | (elements2[:, 2] == elements2[:, 0])
    )
    elements2 = elements2[~degen]
    ring2 = np.where((elements2 == keep).any(axis=1))[0]
    if ring2.size == 0:
        return None
    after = _fails(nodes2, elements2, ring2)
    if after >= before:
        return None

    def _remap(seg):
        s = np.where(np.asarray(seg) == drop, keep, np.asarray(seg))
        keep_m = np.ones(s.size, dtype=bool)
        keep_m[1:] = s[1:] != s[:-1]
        return s[keep_m]

    out = Fort14Mesh(
        title=mesh.title, nodes=nodes2, depths=mesh.depths,
        elements=elements2,
        open_boundaries=[_remap(s) for s in mesh.open_boundaries],
        land_boundaries=[(ib, _remap(s))
                         for ib, s in mesh.land_boundaries],
    )
    return out, {"kept": keep, "dropped": drop,
                 "violations": [before, after],
                 "n_elements_removed": int(degen.sum())}


def split_edge_pair(
    mesh: Fort14Mesh,
    u: int,
    v: int,
    *,
    min_angle: float = 30.0,
    max_angle: float = 130.0,
    max_area_change: float = 0.5,
    relax_iters: int = 3,
) -> tuple[Fort14Mesh, dict[str, Any]] | None:
    """Split interior edge ``(u, v)`` at its midpoint (2 -> 4
    elements) and relax the new node — the long-chord grading remedy
    for C4 jumps. Gated on the C1/C2/C4/flip count of the affected
    ring; returns ``None`` if the edge is not interior or the gate
    fails.
    """
    from fvcom_mesh_tools.algorithms.perp_local import _tri_quality

    u, v = int(u), int(v)
    elements = mesh.elements
    has_u = (elements == u).any(axis=1)
    has_v = (elements == v).any(axis=1)
    carriers = np.where(has_u & has_v)[0]
    if carriers.size != 2:
        return None

    def _replace(seq, a, b):
        return [b if x == a else x for x in seq]

    w = mesh.n_nodes
    nodes2 = np.vstack([
        mesh.nodes, (0.5 * (mesh.nodes[u] + mesh.nodes[v]))[None, :],
    ])
    new_rows = []
    for eid in carriers:
        tri = [int(x) for x in elements[eid]]
        new_rows.append((eid, _replace(tri, v, w), _replace(tri, u, w)))
    elements2 = np.vstack([elements, np.zeros((2, 3), dtype=elements.dtype)])
    for k2, (eid, tri_a, tri_b) in enumerate(new_rows):
        elements2[eid] = tri_a
        elements2[mesh.n_elements + k2] = tri_b

    region = np.where((elements2 == u).any(axis=1)
                      | (elements2 == v).any(axis=1)
                      | (elements2 == w).any(axis=1))[0]

    def _region_bad(nodes_arr):
        tri = elements2[region]
        mn, mx, twice = _tri_quality(nodes_arr[tri])
        n_bad = int((mn < min_angle).sum() + (mx > max_angle).sum()
                    + (twice <= 0).sum())
        areas = 0.5 * np.abs(twice)
        pos_of = {int(e): k3 for k3, e in enumerate(region)}
        edges: dict[tuple[int, int], list[int]] = {}
        for k3, e in enumerate(region):
            a, b, c = (int(x) for x in elements2[e])
            for p_, q_ in ((a, b), (b, c), (c, a)):
                edges.setdefault((min(p_, q_), max(p_, q_)), []).append(int(e))
        for pair in edges.values():
            if len(pair) == 2:
                a1, a2 = areas[pos_of[pair[0]]], areas[pos_of[pair[1]]]
                hi = max(a1, a2)
                if hi > 0 and (hi - min(a1, a2)) / hi > max_area_change:
                    n_bad += 1
        return n_bad

    # Pre-split baseline on the same region topology is not defined;
    # gate against the PRE-EDIT ring of (u, v) instead.
    ring0 = np.where(has_u | has_v)[0]
    tri0 = elements[ring0]
    mn0, mx0, twice0 = _tri_quality(mesh.nodes[tri0])
    before = int((mn0 < min_angle).sum() + (mx0 > max_angle).sum()
                 + (twice0 <= 0).sum())
    areas0 = 0.5 * np.abs(twice0)
    pos0 = {int(e): k3 for k3, e in enumerate(ring0)}
    edges0: dict[tuple[int, int], list[int]] = {}
    for k3, e in enumerate(ring0):
        a, b, c = (int(x) for x in elements[e])
        for p_, q_ in ((a, b), (b, c), (c, a)):
            edges0.setdefault((min(p_, q_), max(p_, q_)), []).append(int(e))
    for pair in edges0.values():
        if len(pair) == 2:
            a1, a2 = areas0[pos0[pair[0]]], areas0[pos0[pair[1]]]
            hi = max(a1, a2)
            if hi > 0 and (hi - min(a1, a2)) / hi > max_area_change:
                before += 1

    # Gated centroid relax of the new node.
    nbr = set()
    for e in region:
        tri = elements2[e]
        if w in tri:
            nbr.update(int(x) for x in tri if x != w)
    best_bad = _region_bad(nodes2)
    for _ in range(relax_iters):
        old_pos = nodes2[w].copy()
        nodes2[w] = nodes2[sorted(nbr)].mean(axis=0)
        bad = _region_bad(nodes2)
        if bad > best_bad:
            nodes2[w] = old_pos
            break
        best_bad = bad

    if best_bad >= before:
        return None
    depths2 = np.concatenate([
        mesh.depths, [0.5 * (mesh.depths[u] + mesh.depths[v])],
    ])
    out = Fort14Mesh(
        title=mesh.title, nodes=nodes2, depths=depths2,
        elements=elements2,
        open_boundaries=[s.copy() for s in mesh.open_boundaries],
        land_boundaries=[(ib, s.copy())
                         for ib, s in mesh.land_boundaries],
    )
    return out, {"new_node": int(w), "violations": [before, best_bad]}
