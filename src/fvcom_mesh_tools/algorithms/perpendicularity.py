"""Make first-ring interior edges perpendicular to the open-boundary tangent.

The user-facing routine is :func:`align_open_boundary_first_ring`. The pre-
existing helpers :func:`unique_edges`, :func:`signed_areas`, and
:func:`open_bdy_perpendicularity` are also exposed for use in metrics and
tests.
"""

from __future__ import annotations

from dataclasses import dataclass

import numpy as np

from fvcom_mesh_tools.io import Fort14Mesh

# ---------------------------------------------------------------------------
# Geometry helpers
# ---------------------------------------------------------------------------

def unique_edges(elements: np.ndarray) -> np.ndarray:
    """Return ``(M, 2)`` int array of unique mesh edges (sorted node indices)."""
    e = np.vstack([
        elements[:, [0, 1]], elements[:, [1, 2]], elements[:, [2, 0]],
    ])
    e.sort(axis=1)
    return np.unique(e, axis=0)


def boundary_tangents(bdy_xy: np.ndarray) -> np.ndarray:
    """Per-node unit tangent along a boundary polyline.

    Central difference for interior nodes, one-sided at the two ends.
    """
    tangents = np.empty_like(bdy_xy)
    tangents[1:-1] = bdy_xy[2:] - bdy_xy[:-2]
    tangents[0] = bdy_xy[1] - bdy_xy[0]
    tangents[-1] = bdy_xy[-1] - bdy_xy[-2]
    norms = np.linalg.norm(tangents, axis=1, keepdims=True)
    return tangents / np.where(norms == 0, 1.0, norms)


def signed_areas(mesh: Fort14Mesh) -> np.ndarray:
    """Per-element signed area (negative => triangle is flipped)."""
    p0 = mesh.nodes[mesh.elements[:, 0]]
    p1 = mesh.nodes[mesh.elements[:, 1]]
    p2 = mesh.nodes[mesh.elements[:, 2]]
    return 0.5 * (
        (p1[:, 0] - p0[:, 0]) * (p2[:, 1] - p0[:, 1])
        - (p1[:, 1] - p0[:, 1]) * (p2[:, 0] - p0[:, 0])
    )


def fixed_node_mask(mesh: Fort14Mesh) -> np.ndarray:
    """Boolean mask of nodes whose position must be preserved.

    By policy: every open and land boundary node is fixed so that boundary
    geometry (open arc, coastline) is invariant under any algorithm exposed
    in this module.
    """
    fixed = np.zeros(mesh.n_nodes, dtype=bool)
    for ids in mesh.open_boundaries:
        fixed[ids] = True
    for _ibtype, ids in mesh.land_boundaries:
        fixed[ids] = True
    return fixed


# ---------------------------------------------------------------------------
# Metrics
# ---------------------------------------------------------------------------

def open_bdy_perpendicularity(
    mesh: Fort14Mesh,
    segment_index: int = 0,
) -> np.ndarray:
    """Per-edge deviation from 90 deg between every interior edge incident on
    the selected open-boundary segment and the local boundary tangent.

    Returns an empty array if the mesh has no open boundary segments.
    """
    if not mesh.open_boundaries:
        return np.array([])

    bdy = np.asarray(mesh.open_boundaries[segment_index], dtype=np.int64)
    bdy_xy = mesh.nodes[bdy]
    tangents = boundary_tangents(bdy_xy)

    inv_map = np.full(mesh.n_nodes, -1, dtype=np.int64)
    inv_map[bdy] = np.arange(len(bdy))

    edges = unique_edges(mesh.elements)
    a_in = inv_map[edges[:, 0]] >= 0
    b_in = inv_map[edges[:, 1]] >= 0
    incident = a_in ^ b_in
    inc = edges[incident]
    inc_a_in = a_in[incident]
    bdy_node = np.where(inc_a_in, inc[:, 0], inc[:, 1])
    int_node = np.where(inc_a_in, inc[:, 1], inc[:, 0])

    edge_vec = mesh.nodes[int_node] - mesh.nodes[bdy_node]
    edge_norms = np.linalg.norm(edge_vec, axis=1, keepdims=True)
    edge_vec = edge_vec / np.where(edge_norms == 0, 1.0, edge_norms)
    edge_tangent = tangents[inv_map[bdy_node]]
    cos_angle = np.clip((edge_vec * edge_tangent).sum(axis=1), -1.0, 1.0)
    angle = np.degrees(np.arccos(np.abs(cos_angle)))
    return 90.0 - angle


# ---------------------------------------------------------------------------
# Algorithm
# ---------------------------------------------------------------------------

@dataclass
class _IncidenceCache:
    bdy: np.ndarray
    tangents: np.ndarray
    perp: np.ndarray
    inv_map: np.ndarray
    bdy_node: np.ndarray
    int_node: np.ndarray
    edge_len_orig: np.ndarray
    sign_orig: np.ndarray


def _build_incidence(mesh: Fort14Mesh, segment_index: int) -> _IncidenceCache:
    bdy = np.asarray(mesh.open_boundaries[segment_index], dtype=np.int64)
    bdy_xy = mesh.nodes[bdy]
    tangents = boundary_tangents(bdy_xy)
    perp = np.column_stack([-tangents[:, 1], tangents[:, 0]])

    inv_map = np.full(mesh.n_nodes, -1, dtype=np.int64)
    inv_map[bdy] = np.arange(len(bdy))

    edges = unique_edges(mesh.elements)
    a_in = inv_map[edges[:, 0]] >= 0
    b_in = inv_map[edges[:, 1]] >= 0
    incident = a_in ^ b_in
    inc = edges[incident]
    inc_a_in = a_in[incident]
    bdy_node = np.where(inc_a_in, inc[:, 0], inc[:, 1])
    int_node = np.where(inc_a_in, inc[:, 1], inc[:, 0])

    edge_vec_orig = mesh.nodes[int_node] - mesh.nodes[bdy_node]
    edge_len_orig = np.linalg.norm(edge_vec_orig, axis=1)
    perp_at_bdy = perp[inv_map[bdy_node]]
    sign_orig = np.sign((edge_vec_orig * perp_at_bdy).sum(axis=1))
    sign_orig = np.where(sign_orig == 0, 1.0, sign_orig)

    return _IncidenceCache(
        bdy=bdy,
        tangents=tangents,
        perp=perp,
        inv_map=inv_map,
        bdy_node=bdy_node,
        int_node=int_node,
        edge_len_orig=edge_len_orig,
        sign_orig=sign_orig,
    )


def _apply_perp_step(
    nodes: np.ndarray,
    cache: _IncidenceCache,
    fixed: np.ndarray,
    alpha: float,
) -> np.ndarray:
    """One damped step toward per-edge perpendicular targets.

    Targets are computed using the *original* edge length and the *original*
    sign (which side of the boundary the interior node sits on); the boundary
    tangents are also fixed because boundary nodes never move. Iterating with
    ``alpha < 1`` therefore relaxes interior nodes monotonically toward the
    same fixed point that ``alpha == 1`` reaches in a single step.
    """
    bdy_pos = nodes[cache.bdy_node]
    perp_at_bdy = cache.perp[cache.inv_map[cache.bdy_node]]
    target_pos = (
        bdy_pos
        + cache.edge_len_orig[:, None]
        * (perp_at_bdy * cache.sign_orig[:, None])
    )

    accum = np.zeros_like(nodes)
    cnt = np.zeros(len(nodes), dtype=np.int64)
    np.add.at(accum, cache.int_node, target_pos)
    np.add.at(cnt, cache.int_node, 1)

    new_nodes = nodes.copy()
    moved_idx = np.where(cnt > 0)[0]
    if moved_idx.size:
        avg_target = accum[moved_idx] / cnt[moved_idx, None]
        new_nodes[moved_idx] = (
            (1.0 - alpha) * nodes[moved_idx] + alpha * avg_target
        )

    new_nodes[fixed] = nodes[fixed]
    return new_nodes


def _laplacian_smooth_second_ring(
    nodes: np.ndarray,
    elements: np.ndarray,
    first_ring: np.ndarray,
    fixed: np.ndarray,
    alpha: float,
) -> np.ndarray:
    """One Laplacian smoothing pass over the second-ring nodes.

    "Second ring" = nodes that share an edge with a first-ring node but are
    themselves neither boundary nor first-ring. This pass absorbs the
    disturbance that the perpendicular projection introduces beyond the
    immediate ring.
    """
    edges = unique_edges(elements)
    in_first = np.zeros(len(nodes), dtype=bool)
    in_first[first_ring] = True

    a_in = in_first[edges[:, 0]]
    b_in = in_first[edges[:, 1]]
    second_ring_edges = edges[a_in ^ b_in]
    second_ring_node = np.where(
        a_in[a_in ^ b_in],
        second_ring_edges[:, 1],
        second_ring_edges[:, 0],
    )
    candidate = np.unique(second_ring_node)
    candidate = candidate[~fixed[candidate] & ~in_first[candidate]]

    # For each candidate, its Laplacian neighbours are all nodes connected by
    # any edge in the full mesh.
    n0, n1 = edges[:, 0], edges[:, 1]
    new_nodes = nodes.copy()
    for c in candidate:
        nbr = np.concatenate([n1[n0 == c], n0[n1 == c]])
        if nbr.size == 0:
            continue
        avg = nodes[nbr].mean(axis=0)
        new_nodes[c] = (1.0 - alpha) * nodes[c] + alpha * avg
    return new_nodes


def align_open_boundary_first_ring(
    mesh: Fort14Mesh,
    *,
    alpha: float = 1.0,
    n_iters: int = 1,
    smooth_iters: int = 0,
    smooth_alpha: float = 0.3,
    segment_index: int = 0,
) -> tuple[Fort14Mesh, dict]:
    """Move the first-ring interior neighbours of every open-boundary node
    so that each incident interior edge is as perpendicular to the local
    boundary tangent as possible while preserving the original edge length.

    Boundary nodes (open and land) are never moved. When an interior node is
    shared between several boundary parents, the per-edge targets are
    averaged; the resulting edge is therefore not perfectly perpendicular to
    any single tangent but represents a best length-preserving compromise.

    Parameters
    ----------
    mesh:
        Input mesh. Topology and boundary lists are copied unchanged.
    alpha:
        Per-iteration damping factor. ``alpha=1`` jumps directly to the
        target; ``alpha<1`` relaxes monotonically toward the same fixed
        point.
    n_iters:
        Number of perpendicular-projection iterations.
    smooth_iters:
        If ``>0``, run that many Laplacian smoothing passes on second-ring
        interior nodes after the projection passes. Helps absorb the
        disturbance into the rest of the mesh; does not change the
        first-ring perpendicularity result.
    smooth_alpha:
        Damping factor for the Laplacian smoothing pass.
    segment_index:
        Which open-boundary segment to operate on (default: the first).

    Returns
    -------
    (Fort14Mesh, dict)
        A new mesh with updated node coordinates, plus an info dict
        containing parameter echo, count of moved interior nodes, and the
        per-parent-count breakdown of incident edges.
    """
    if not mesh.open_boundaries:
        return mesh, {"moved": 0, "note": "no open boundaries"}

    cache = _build_incidence(mesh, segment_index)
    fixed = fixed_node_mask(mesh)

    nodes = mesh.nodes.copy()
    for _ in range(n_iters):
        nodes = _apply_perp_step(nodes, cache, fixed, alpha)

    # Per-interior-node parent counts (1 = single-parent, >=2 = multi-parent).
    cnt = np.zeros(mesh.n_nodes, dtype=np.int64)
    np.add.at(cnt, cache.int_node, 1)
    first_ring = np.where(cnt > 0)[0]
    movable_first_ring = first_ring[~fixed[first_ring]]

    def _by_count(idx: np.ndarray) -> dict[int, int]:
        counts = cnt[idx]
        return {int(k): int((counts == k).sum()) for k in np.unique(counts)}

    by_parent_count = _by_count(first_ring)
    movable_by_parent_count = _by_count(movable_first_ring)

    if smooth_iters > 0:
        for _ in range(smooth_iters):
            nodes = _laplacian_smooth_second_ring(
                nodes, mesh.elements, first_ring, fixed, smooth_alpha,
            )

    out = Fort14Mesh(
        title=mesh.title,
        nodes=nodes,
        depths=mesh.depths.copy(),
        elements=mesh.elements.copy(),
        open_boundaries=[a.copy() for a in mesh.open_boundaries],
        land_boundaries=[(ib, a.copy()) for (ib, a) in mesh.land_boundaries],
    )
    info = {
        "alpha": alpha,
        "n_iters": n_iters,
        "smooth_iters": smooth_iters,
        "smooth_alpha": smooth_alpha,
        "segment_index": segment_index,
        "moved": int(movable_first_ring.size),
        "first_ring_by_parent_count": by_parent_count,
        "movable_first_ring_by_parent_count": movable_by_parent_count,
    }
    return out, info
