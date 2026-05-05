"""Delaunay-style edge swapping for unstructured triangle meshes.

For every interior edge ``(i, j)`` shared by two triangles ``T1, T2``,
the swap replaces the diagonal: ``T1, T2 -> T1', T2'`` where ``T1'`` and
``T2'`` share the new diagonal connecting the two "opposite" vertices.
The swap is accepted if it strictly improves the worst minimum-interior
angle of the two triangles. This is the classical sliver-removal step
that pure Laplacian smoothing (PoC #8) cannot perform - slivers are
topology, not metric, problems.

Within one pass:

* All swap candidates are evaluated in a single vectorised pass.
* Candidates that produce a degenerate or flipped triangle are skipped.
* Accepted candidates are applied in order of largest improvement,
  greedily skipping any whose two triangles have already been
  consumed by an earlier (better) swap. This avoids cascade
  conflicts within one pass.

Multiple passes converge quickly: in typical use only 3-6 passes are
required before no further improvement is possible.
"""

from __future__ import annotations

from dataclasses import replace

import numpy as np

from fvcom_mesh_tools.io import Fort14Mesh

# ---------------------------------------------------------------------------
# Per-triangle vectorised helpers (independent of the package-level metrics
# in algorithms.quality, so this module stays self-contained and we can
# evaluate candidate triangles without round-tripping through Fort14Mesh).
# ---------------------------------------------------------------------------


def _signed_area(nodes: np.ndarray, tris: np.ndarray) -> np.ndarray:
    p0 = nodes[tris[:, 0]]
    p1 = nodes[tris[:, 1]]
    p2 = nodes[tris[:, 2]]
    return 0.5 * (
        (p1[:, 0] - p0[:, 0]) * (p2[:, 1] - p0[:, 1])
        - (p1[:, 1] - p0[:, 1]) * (p2[:, 0] - p0[:, 0])
    )


def _min_angle(nodes: np.ndarray, tris: np.ndarray) -> np.ndarray:
    p0 = nodes[tris[:, 0]]
    p1 = nodes[tris[:, 1]]
    p2 = nodes[tris[:, 2]]
    l01 = np.linalg.norm(p1 - p0, axis=1)
    l12 = np.linalg.norm(p2 - p1, axis=1)
    l20 = np.linalg.norm(p0 - p2, axis=1)
    a, b, c = l12, l20, l01

    def _ang(opp: np.ndarray, e1: np.ndarray, e2: np.ndarray) -> np.ndarray:
        denom = np.where(e1 * e2 == 0, 1.0, 2.0 * e1 * e2)
        cos = (e1 ** 2 + e2 ** 2 - opp ** 2) / denom
        return np.arccos(np.clip(cos, -1.0, 1.0))

    return np.degrees(np.minimum(np.minimum(_ang(a, b, c), _ang(b, c, a)), _ang(c, a, b)))


# ---------------------------------------------------------------------------
# Interior-edge / triangle-pair extraction
# ---------------------------------------------------------------------------


def _interior_edge_pairs(elements: np.ndarray) -> tuple[np.ndarray, np.ndarray]:
    """For every interior edge, return its two endpoint nodes and its two
    incident triangle indices.

    Returns
    -------
    edges:
        ``(M, 2)`` int array of (sorted) endpoint node indices.
    tri_pairs:
        ``(M, 2)`` int array of triangle indices into ``elements``.
    """
    ne = elements.shape[0]
    flat = np.vstack(
        [elements[:, [0, 1]], elements[:, [1, 2]], elements[:, [2, 0]]]
    )
    flat_sorted = np.sort(flat, axis=1)
    uniq, inverse, counts = np.unique(
        flat_sorted, axis=0, return_inverse=True, return_counts=True,
    )
    is_interior = counts == 2

    # For each unique edge, group its flat-row indices together. After
    # sorting by `inverse`, all entries with the same edge id are
    # adjacent; the cumsum of `counts` indexes their start.
    order = np.argsort(inverse, kind="stable")
    starts = np.zeros(uniq.shape[0] + 1, dtype=np.int64)
    np.cumsum(counts, out=starts[1:])

    # ``np.vstack([elements[:, [0,1]], elements[:, [1,2]], elements[:, [2,0]]])``
    # lays out 3 NE-blocks: rows ``[0, NE)`` come from edge01 of every
    # triangle, rows ``[NE, 2 NE)`` from edge12, etc. So the triangle id
    # of flat row ``k`` is ``k % NE``, not ``k // 3``.
    interior_starts = starts[:-1][is_interior]
    flat_tri_id = order % ne
    tri_pairs = np.stack(
        [flat_tri_id[interior_starts], flat_tri_id[interior_starts + 1]], axis=1
    )
    return uniq[is_interior], tri_pairs


# ---------------------------------------------------------------------------
# Edge-swap pass and full driver
# ---------------------------------------------------------------------------


def _swap_pass(elements: np.ndarray, nodes: np.ndarray) -> int:
    """Apply one pass of non-conflicting min-angle edge swaps in place.

    Returns the number of swaps applied.
    """
    edges, tri_pairs = _interior_edge_pairs(elements)
    if edges.size == 0:
        return 0

    t1_idx = tri_pairs[:, 0]
    t2_idx = tri_pairs[:, 1]
    T1 = elements[t1_idx]
    T2 = elements[t2_idx]
    i, j = edges[:, 0], edges[:, 1]

    # Third vertex of each triangle (the "opposite" node to the edge).
    k_nodes = T1.sum(axis=1) - i - j
    m_nodes = T2.sum(axis=1) - i - j

    # Candidate triangles after the diagonal flip.
    cand1 = np.stack([i, m_nodes, k_nodes], axis=1)
    cand2 = np.stack([j, k_nodes, m_nodes], axis=1)

    # Normalise orientation to CCW (positive signed area). If a candidate
    # has zero or negative area in its preferred orientation, the quad
    # is non-convex and the swap is invalid; mark it so.
    for cand in (cand1, cand2):
        sa = _signed_area(nodes, cand)
        flip = sa < 0
        if flip.any():
            cand[flip] = cand[flip][:, [0, 2, 1]]

    sa1 = _signed_area(nodes, cand1)
    sa2 = _signed_area(nodes, cand2)
    valid = (sa1 > 0) & (sa2 > 0)

    # Quality before / after.
    q_before = np.minimum(_min_angle(nodes, T1), _min_angle(nodes, T2))
    q_after = np.minimum(_min_angle(nodes, cand1), _min_angle(nodes, cand2))
    improves = valid & (q_after > q_before + 1e-9)

    if not improves.any():
        return 0

    # Sort improving swaps by delta-quality (descending) and apply
    # non-conflicting ones greedily. A "conflict" is two swaps that
    # share a triangle.
    delta = np.where(improves, q_after - q_before, -1.0)
    order = np.argsort(-delta, kind="stable")

    used = np.zeros(elements.shape[0], dtype=bool)
    n_applied = 0
    for k in order:
        if not improves[k]:
            break
        t1 = int(t1_idx[k])
        t2 = int(t2_idx[k])
        if used[t1] or used[t2]:
            continue
        elements[t1] = cand1[k]
        elements[t2] = cand2[k]
        used[t1] = True
        used[t2] = True
        n_applied += 1

    return n_applied


def swap_edges_for_quality(
    mesh: Fort14Mesh,
    *,
    max_iters: int = 20,
) -> tuple[Fort14Mesh, dict]:
    """Iterate edge swaps until no further improvement is possible.

    Boundary edges (those incident on only one triangle) cannot be
    swapped, so coast and open arcs are preserved by construction.

    Parameters
    ----------
    mesh:
        Input mesh.
    max_iters:
        Hard upper bound on the number of passes (default 20). Real
        meshes converge in 3-6 passes; the cap is just a safety net.

    Returns
    -------
    (swapped_mesh, info)
        ``info`` keys: ``"max_iters"``, ``"swaps_per_iter"``,
        ``"total_swaps"``.
    """
    elements = mesh.elements.copy()
    swaps_per_iter: list[int] = []
    for _ in range(max_iters):
        n = _swap_pass(elements, mesh.nodes)
        swaps_per_iter.append(n)
        if n == 0:
            break

    out = replace(mesh, elements=elements)
    return out, {
        "max_iters": max_iters,
        "swaps_per_iter": swaps_per_iter,
        "total_swaps": int(sum(swaps_per_iter)),
    }
