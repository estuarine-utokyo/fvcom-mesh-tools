"""Mesh smoothing for interior nodes.

Currently exposes one routine, :func:`laplacian_smooth`, which performs
vectorised Jacobi-style Laplacian smoothing on every interior node
while pinning every boundary node (open + land). Triangle flips are
detected after each iteration and the offending nodes are reverted, so
the returned mesh is always strictly valid.
"""

from __future__ import annotations

from dataclasses import replace

import numpy as np

from fvcom_mesh_tools.algorithms.perpendicularity import (
    fixed_node_mask,
    signed_areas,
    unique_edges,
)
from fvcom_mesh_tools.io import Fort14Mesh


def _neighbor_means(nodes: np.ndarray, edges: np.ndarray) -> tuple[np.ndarray, np.ndarray]:
    n = nodes.shape[0]
    sums = np.zeros((n, 2), dtype=np.float64)
    deg = np.zeros(n, dtype=np.int64)
    np.add.at(sums, edges[:, 0], nodes[edges[:, 1]])
    np.add.at(sums, edges[:, 1], nodes[edges[:, 0]])
    np.add.at(deg, edges[:, 0], 1)
    np.add.at(deg, edges[:, 1], 1)
    safe = np.maximum(deg, 1)
    return sums / safe[:, None], deg


def laplacian_smooth(
    mesh: Fort14Mesh,
    *,
    n_iters: int = 10,
    alpha: float = 0.5,
    prevent_flips: bool = True,
) -> tuple[Fort14Mesh, dict]:
    """Damped Laplacian smoothing of every interior node.

    Each iteration moves every movable node toward the centroid of its
    one-ring neighbours:

        new_pos[i] = (1 - alpha) * pos[i] + alpha * mean(neighbours[i])

    Boundary nodes (open + land) are kept fixed. With
    ``prevent_flips=True`` (default), any node whose new position would
    cause an incident triangle to flip is reverted within that
    iteration. The returned mesh is therefore always strictly valid.

    Parameters
    ----------
    mesh:
        Input mesh.
    n_iters:
        Number of smoothing iterations (default 10).
    alpha:
        Per-iteration damping factor in (0, 1] (default 0.5). Higher
        values converge faster but are more likely to require a revert
        on slivers.
    prevent_flips:
        Revert any node move that would flip an incident triangle.

    Returns
    -------
    (smoothed_mesh, info)
        ``info`` keys: ``"n_iters"``, ``"alpha"``, ``"moved"`` (final
        boolean mask of nodes whose position changed), ``"reverts"``
        (per-iteration revert counts), ``"degenerate_remaining"``
        (count of triangles still <= 0 area).
    """
    if not 0.0 < alpha <= 1.0:
        raise ValueError("alpha must lie in (0, 1].")
    if n_iters < 1:
        raise ValueError("n_iters must be >= 1.")

    nodes = mesh.nodes.copy()
    elements = mesh.elements
    fixed = fixed_node_mask(mesh)
    movable = ~fixed
    edges = unique_edges(elements)

    reverts: list[int] = []
    for _ in range(n_iters):
        means, deg = _neighbor_means(nodes, edges)
        target = nodes.copy()
        valid = movable & (deg > 0)
        target[valid] = (1.0 - alpha) * nodes[valid] + alpha * means[valid]
        candidate = target.copy()

        if prevent_flips:
            # Compute signed areas with the candidate update applied.
            p0 = candidate[elements[:, 0]]
            p1 = candidate[elements[:, 1]]
            p2 = candidate[elements[:, 2]]
            sa = 0.5 * (
                (p1[:, 0] - p0[:, 0]) * (p2[:, 1] - p0[:, 1])
                - (p1[:, 1] - p0[:, 1]) * (p2[:, 0] - p0[:, 0])
            )
            bad = sa <= 0
            n_revert = 0
            if bool(bad.any()):
                # Revert every movable node that participates in a bad
                # triangle. One iteration of "revert and rerun the
                # signed-area test" is normally enough; loop a few
                # times to handle cascades.
                for _ in range(3):
                    bad_nodes = np.zeros(nodes.shape[0], dtype=bool)
                    bad_nodes[np.unique(elements[bad].ravel())] = True
                    revert_mask = bad_nodes & movable & np.any(
                        candidate != nodes, axis=1
                    )
                    n_revert += int(revert_mask.sum())
                    if not revert_mask.any():
                        break
                    candidate[revert_mask] = nodes[revert_mask]
                    p0 = candidate[elements[:, 0]]
                    p1 = candidate[elements[:, 1]]
                    p2 = candidate[elements[:, 2]]
                    sa = 0.5 * (
                        (p1[:, 0] - p0[:, 0]) * (p2[:, 1] - p0[:, 1])
                        - (p1[:, 1] - p0[:, 1]) * (p2[:, 0] - p0[:, 0])
                    )
                    bad = sa <= 0
                    if not bad.any():
                        break
            reverts.append(n_revert)
        else:
            reverts.append(0)

        nodes = candidate

    out = replace(mesh, nodes=nodes)
    final_sa = signed_areas(out)
    info = {
        "n_iters": n_iters,
        "alpha": alpha,
        "moved": np.any(nodes != mesh.nodes, axis=1),
        "reverts": reverts,
        "degenerate_remaining": int((final_sa <= 0).sum()),
    }
    return out, info
