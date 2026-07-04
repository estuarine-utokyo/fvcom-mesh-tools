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
        if len(run) < min_run:
            continue
        line = geoms[eligible[run[0]]]
        s_vals = [float(line.project(shapely.Point(*nodes[v])))
                  for v in run]
        diffs = np.diff(s_vals)
        if not ((diffs > 0).all() or (diffs < 0).all()):
            n_trimmed += 1
            continue
        if (np.abs(diffs) < 0.2 * np.array(
                [local_h[v] for v in run[:-1]])).any():
            n_trimmed += 1
            continue
        foot = [line.interpolate(s) for s in s_vals]
        base = n_nodes + len(new_nodes)
        ids = list(range(base, base + len(run)))
        coords = {v: nodes[v] for v in run}
        for k, f in enumerate(foot):
            coords[ids[k]] = np.array([f.x, f.y])
        strip: list[list[int]] = []
        ok = True
        for k in range(len(run) - 1):
            vi, vj = run[k], run[k + 1]
            ni, nj = ids[k], ids[k + 1]
            # Two candidate diagonal splits of quad (vi, vj, nj, ni).
            cand_a = [_tri_ccw(vi, vj, nj, coords),
                      _tri_ccw(vi, nj, ni, coords)]
            cand_b = [_tri_ccw(vi, vj, ni, coords),
                      _tri_ccw(vj, nj, ni, coords)]
            pick = None
            for cand in (cand_a, cand_b):
                if all(_tri_ok(c, coords) for c in cand):
                    pick = cand
                    break
            if pick is None:
                ok = False
                break
            strip.extend(pick)
        if not ok:
            n_trimmed += 1
            continue
        for k, f in enumerate(foot):
            new_nodes.append(np.array([f.x, f.y]))
            new_depths.append(float(mesh.depths[run[k]]))
        new_tris.extend(strip)
        n_strips += 1

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
