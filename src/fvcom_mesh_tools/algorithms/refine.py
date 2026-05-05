"""Local refinement of bad-quality triangles via longest-edge bisection.

PoC #14 left ``frac<20deg ~ 2.84 %`` after coastline-aware sizing,
quality-pass smoothing and edge swapping. The remaining slivers are
local features that no in-place rearrangement of the existing nodes
can fix.

This module performs a *topological* refinement: for every triangle
whose minimum interior angle is below a user-supplied threshold, the
*longest* edge of that triangle is bisected by inserting a midpoint.
Each affected triangle is replaced by 2-4 sub-triangles depending on
how many of its edges happen to be bisected in the same pass.

Centroid insertion (an obvious alternative) was tried first and
empirically *increased* sliver count on Tokyo Bay - it splits the
sliver into one wedge that retains the sliver's long base.

Boundary edges are never bisected (would split a boundary segment in
fort.14). A bad triangle whose longest edge is on the boundary is
left untouched in this pass; the edge-swap pass that runs after may
flip a different edge to fix it.
"""

from __future__ import annotations

from dataclasses import replace

import numpy as np

from fvcom_mesh_tools.algorithms.edge_swap import swap_edges_for_quality
from fvcom_mesh_tools.algorithms.quality import min_interior_angle
from fvcom_mesh_tools.io import Fort14Mesh


def _split_subtriangles(tri: np.ndarray, m: list[int | None]) -> list[list[int]]:
    """Return CCW sub-triangles for a single parent triangle whose edges
    are bisected at ``m[0..2]`` (None when not bisected).

    Edge indexing matches ``elements`` columns:
        edge 0 = (tri[0], tri[1])
        edge 1 = (tri[1], tri[2])
        edge 2 = (tri[2], tri[0])
    """
    a, b, c = int(tri[0]), int(tri[1]), int(tri[2])
    m0, m1, m2 = m
    n = (m0 is not None) + (m1 is not None) + (m2 is not None)

    if n == 0:
        return [[a, b, c]]

    if n == 3:
        # Triforce split: 4 sub-triangles, all CCW.
        return [
            [a, m0, m2],
            [m0, b, m1],
            [m0, m1, m2],
            [m2, m1, c],
        ]

    if n == 1:
        if m0 is not None:
            # Bisect (a,b) at m0; c is the apex.
            return [[a, m0, c], [m0, b, c]]
        if m1 is not None:
            return [[b, m1, a], [m1, c, a]]
        return [[c, m2, b], [m2, a, b]]

    # n == 2: 3 sub-triangles. The unsplit edge defines the "tip"
    # vertex pair; the apex is the vertex opposite to the tip.
    if m0 is None:
        # Edges (b,c) and (c,a) split at m1 and m2.
        return [
            [c, m2, m1],   # apex with the two new midpoints
            [a, b, m1],    # base side, includes unsplit edge (a,b)
            [a, m1, m2],
        ]
    if m1 is None:
        # Edges (a,b) and (c,a) split at m0 and m2.
        return [
            [a, m0, m2],
            [b, c, m0],
            [c, m2, m0],
        ]
    # m2 is None: edges (a,b) and (b,c) split at m0 and m1.
    return [
        [b, m1, m0],
        [c, a, m1],
        [a, m0, m1],
    ]


def _bisect_longest_edges(
    mesh: Fort14Mesh, bad_mask: np.ndarray
) -> tuple[Fort14Mesh, int]:
    """One pass of longest-edge bisection on every triangle flagged in
    ``bad_mask``. Boundary edges are not split; bad triangles whose
    longest edge is on the boundary are left in place this pass.

    Returns ``(new_mesh, n_edges_split)``.
    """
    nodes = mesh.nodes
    depths = mesh.depths
    elements = mesh.elements
    ne = elements.shape[0]

    # Edge -> count adjacency: any edge with count == 1 is on the
    # boundary and must not be split.
    flat = np.vstack(
        [elements[:, [0, 1]], elements[:, [1, 2]], elements[:, [2, 0]]]
    )
    flat_sorted = np.sort(flat, axis=1)
    uniq, inverse, counts = np.unique(
        flat_sorted, axis=0, return_inverse=True, return_counts=True,
    )
    is_interior_edge = counts == 2

    # Per-triangle edge lengths.
    p0 = nodes[elements[:, 0]]
    p1 = nodes[elements[:, 1]]
    p2 = nodes[elements[:, 2]]
    lens = np.column_stack([
        np.linalg.norm(p1 - p0, axis=1),
        np.linalg.norm(p2 - p1, axis=1),
        np.linalg.norm(p0 - p2, axis=1),
    ])
    longest_pos = np.argmax(lens, axis=1)  # (NE,) in {0, 1, 2}

    # Global edge id for each triangle's longest edge. The vstack
    # layout means flat row ``k`` belongs to triangle ``k % NE`` and
    # edge ``k // NE``; equivalently global edge id of triangle t's
    # edge p is ``inverse[p * NE + t]``.
    bad_idx = np.where(bad_mask)[0]
    bad_longest_global = inverse[longest_pos[bad_idx] * ne + bad_idx]
    interior_ok = is_interior_edge[bad_longest_global]
    edges_to_split = np.unique(bad_longest_global[interior_ok])
    if edges_to_split.size == 0:
        return mesh, 0

    # Insert midpoints + interpolate depths.
    edge_pairs = uniq[edges_to_split]
    midpoints = 0.5 * (nodes[edge_pairs[:, 0]] + nodes[edge_pairs[:, 1]])
    new_depths_seg = 0.5 * (
        depths[edge_pairs[:, 0]] + depths[edge_pairs[:, 1]]
    )
    new_node_ids = np.arange(
        nodes.shape[0], nodes.shape[0] + edges_to_split.size, dtype=np.int64
    )

    midpoint_for: dict[int, int] = {
        int(eg): int(nid) for eg, nid in zip(edges_to_split, new_node_ids)
    }

    # Per-triangle: which of its three edges (by position) was split.
    # inverse_pos[k, p] = global edge id of triangle k's edge p.
    inverse_per_pos = inverse.reshape(3, ne).T  # (NE, 3)

    new_elements: list[list[int]] = []
    for t in range(ne):
        eg = inverse_per_pos[t]
        m = [midpoint_for.get(int(eg[p])) for p in range(3)]
        new_elements.extend(_split_subtriangles(elements[t], m))

    new_nodes = np.vstack([nodes, midpoints])
    new_depths = np.concatenate([depths, new_depths_seg])
    new_elements_np = np.asarray(new_elements, dtype=np.int64)

    out = replace(
        mesh, nodes=new_nodes, depths=new_depths, elements=new_elements_np,
    )
    return out, int(edges_to_split.size)


def refine_bad_triangles(
    mesh: Fort14Mesh,
    *,
    min_angle_threshold: float = 20.0,
    max_passes: int = 5,
    swap_after_each_pass: bool = True,
) -> tuple[Fort14Mesh, dict]:
    """Iterative longest-edge bisection of triangles below the
    ``min_angle_threshold`` (degrees). Each pass:

    1. Identify bad triangles.
    2. Bisect their longest interior edge (skip if boundary).
    3. Run :func:`swap_edges_for_quality` to clean up resulting
       slivers (skippable via ``swap_after_each_pass=False``).

    Stops early when no triangle is below threshold or no edge is
    eligible for bisection.

    Boundary segments (open + land) are unchanged: only interior
    nodes are added (mid-edge of strictly interior edges).
    """
    if min_angle_threshold <= 0:
        raise ValueError("min_angle_threshold must be > 0 deg.")
    if max_passes < 1:
        raise ValueError("max_passes must be >= 1.")

    best = mesh
    best_bad = int((min_interior_angle(mesh) < min_angle_threshold).sum())
    out = mesh
    history: list[dict] = []
    nodes_inserted = 0
    swaps_total = 0
    passes_used = 0
    stop_reason = "max_passes"

    for p in range(max_passes):
        ma = min_interior_angle(out)
        bad = ma < min_angle_threshold
        n_bad_before = int(bad.sum())
        if n_bad_before == 0:
            stop_reason = "no_bad_triangles"
            break

        candidate, n_edges_split = _bisect_longest_edges(out, bad)
        if n_edges_split == 0:
            stop_reason = "only_boundary_edges_eligible"
            break

        swap_count = 0
        if swap_after_each_pass:
            candidate, swap_info = swap_edges_for_quality(candidate, max_iters=10)
            swap_count = int(swap_info["total_swaps"])

        ma_after = min_interior_angle(candidate)
        n_bad_after = int((ma_after < min_angle_threshold).sum())
        history.append({
            "pass": p + 1,
            "n_bad_before": n_bad_before,
            "n_bad_after": n_bad_after,
            "edges_split": n_edges_split,
            "swaps": swap_count,
            "frac_below_after": float((ma_after < min_angle_threshold).mean()),
        })

        # Early stop: if this pass would *increase* bad-triangle count,
        # back out and return the best mesh seen so far. Longest-edge
        # bisection without proper Rivara propagation can oscillate
        # once the easy slivers are gone.
        if n_bad_after >= best_bad:
            stop_reason = "regression_rolled_back"
            break

        out = candidate
        nodes_inserted += n_edges_split
        swaps_total += swap_count
        best = candidate
        best_bad = n_bad_after
        passes_used = p + 1

    info = {
        "min_angle_threshold": min_angle_threshold,
        "passes": passes_used,
        "stop_reason": stop_reason,
        "per_pass": history,
        "total_nodes_inserted": nodes_inserted,
        "total_swaps": swaps_total,
    }
    return best, info
