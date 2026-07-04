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
