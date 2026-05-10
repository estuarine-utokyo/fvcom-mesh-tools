"""Phase H — per-element greedy quality optimiser.

Phase H is the planned automation of the SMS manual mesh-edit
workflow: visit each element that fails a strict per-element gate
(``alpha >= alpha_target`` ∧ ``min_angle >= min_angle_target``) and
try a sequence of local-edit operators until one improves the
surrounding 1-ring without making any neighbour worse.

The operator inventory:

* :func:`_apply_smooth_node` — Gauss-Seidel-style move of an
  interior vertex to its 1-ring centroid; cheap, applied in batch
  in :func:`_batch_smooth_sweep` before any topology operator runs.
* :func:`_apply_edge_swap` — Lawson swap on an internal edge with
  the alpha-driven acceptance criterion.
* :func:`_apply_edge_split_interior` — insert a midpoint on an
  interior edge; the two incident triangles become four sub-
  triangles. NP +1, NE +2.
* :func:`_apply_vertex_remove` — remove an interior vertex, gather
  its 1-ring, and re-triangulate the resulting polygon via Delaunay
  pruned by the rim (single-element variant of the Stage 2 medial-
  axis re-mesh).

Boundary handling — v1 is conservative: the operators refuse to
move or insert a node that lies on (or would lie on) an open / land
boundary. v2 will add a coastline-projecting boundary edge_split
and a boundary-tangent smooth.

:func:`phase_h_optimize` runs in two passes for performance: first
a sequence of in-place batch smooth sweeps (Gauss-Seidel; each
sweep visits every interior vertex once and accepts the centroid
move iff the per-1-ring penalty strictly decreases without flipping
a triangle), then a per-element greedy loop that pops fail elements
by descending penalty and tries the topology-changing operators.
The smooth pass is asymptotically O(NP × ring) per sweep and does
not rebuild the topology aux dicts; the topology pass costs an aux
rebuild per accept but accepts are few (the residual after the
smooth pass is small).
"""
from __future__ import annotations

import heapq
from collections import defaultdict
from typing import Any

import numpy as np

from fvcom_mesh_tools.algorithms.quality import (
    alpha_quality,
    min_interior_angle,
)
from fvcom_mesh_tools.io import Fort14Mesh

# Re-use the Stage 2 retriangulation helpers.
from fvcom_mesh_tools.mesh_clean import (
    _patch_rim_polygon,
    _retriangulate_patch,
)

EARTH_R_M: float = 6_371_000.0


# ---------------------------------------------------------------------------
# Penalty / quality helpers
# ---------------------------------------------------------------------------


def _per_element_quality(
    nodes: np.ndarray, elements: np.ndarray,
) -> tuple[np.ndarray, np.ndarray]:
    """Per-element ``(alpha, min_angle_deg)`` arrays."""
    if elements.size == 0:
        return np.empty(0), np.empty(0)
    mesh = Fort14Mesh(
        title="phaseh", nodes=nodes,
        depths=np.zeros(nodes.shape[0]),
        elements=elements,
        open_boundaries=[], land_boundaries=[],
    )
    return alpha_quality(mesh), min_interior_angle(mesh)


def _inline_quality(
    p0: np.ndarray, p1: np.ndarray, p2: np.ndarray,
) -> tuple[np.ndarray, np.ndarray, np.ndarray]:
    """Vectorised per-triangle ``(alpha, min_angle_deg, twice_signed)``.

    Avoids the ``Fort14Mesh`` constructor overhead by operating on raw
    coordinate arrays. Used by the hot batch-smooth sweep.
    """
    twice_signed = (
        (p1[:, 0] - p0[:, 0]) * (p2[:, 1] - p0[:, 1])
        - (p1[:, 1] - p0[:, 1]) * (p2[:, 0] - p0[:, 0])
    )
    area = 0.5 * np.abs(twice_signed)
    e0 = np.linalg.norm(p1 - p0, axis=1)  # opp v2
    e1 = np.linalg.norm(p2 - p1, axis=1)  # opp v0
    e2 = np.linalg.norm(p0 - p2, axis=1)  # opp v1
    sum_sq = e0 * e0 + e1 * e1 + e2 * e2
    alpha = np.where(
        sum_sq > 0, 4.0 * np.sqrt(3.0) * area / np.where(sum_sq > 0, sum_sq, 1.0),
        0.0,
    )
    # Min interior angle via law of cosines on each vertex.
    def _ang(opp, ea, eb):
        denom = 2.0 * ea * eb
        safe_denom = np.where(denom > 0, denom, 1.0)
        cos = np.where(
            denom > 0, (ea * ea + eb * eb - opp * opp) / safe_denom, 0.0,
        )
        return np.arccos(np.clip(cos, -1.0, 1.0))
    angle_v0 = _ang(e1, e2, e0)
    angle_v1 = _ang(e2, e0, e1)
    angle_v2 = _ang(e0, e1, e2)
    min_ang_rad = np.minimum(np.minimum(angle_v0, angle_v1), angle_v2)
    return alpha, np.degrees(min_ang_rad), twice_signed


def _penalty(
    alpha: np.ndarray, min_ang: np.ndarray,
    *, alpha_target: float, min_angle_target: float,
) -> np.ndarray:
    """Element penalty: zero iff both gates met. Squared deficits with
    the angle term scaled by 1/100 to keep both contributions in the
    same range."""
    a_pen = np.maximum(0.0, alpha_target - alpha) ** 2
    g_pen = np.maximum(0.0, min_angle_target - min_ang) ** 2 / 100.0
    return a_pen + g_pen


def _is_fail(
    alpha: np.ndarray, min_ang: np.ndarray,
    *, alpha_target: float, min_angle_target: float,
) -> np.ndarray:
    return (alpha < alpha_target) | (min_ang < min_angle_target)


# ---------------------------------------------------------------------------
# Topology / boundary helpers
# ---------------------------------------------------------------------------


def _node_to_elements(elements: np.ndarray, n_nodes: int
                      ) -> dict[int, np.ndarray]:
    """Map node-id → array of incident element ids."""
    # ``elements.ravel()`` lays out vertices in row-major order:
    # position ``k*3 + i`` is vertex ``i`` of element ``k``. The
    # element-id at each position is therefore ``k = pos // 3``,
    # i.e. ``np.repeat(np.arange(NE), 3)`` (not ``np.tile``).
    rows = elements.ravel()
    cols = np.repeat(np.arange(elements.shape[0]), 3)
    order = np.argsort(rows, kind="stable")
    rows = rows[order]
    cols = cols[order]
    boundaries = np.searchsorted(rows, np.arange(n_nodes + 1))
    out: dict[int, np.ndarray] = {}
    for n in range(n_nodes):
        s, e = boundaries[n], boundaries[n + 1]
        if s < e:
            out[n] = cols[s:e].copy()
    return out


def _edge_use_counts(elements: np.ndarray) -> dict[tuple[int, int], list[int]]:
    """For every undirected edge, return the list of incident element
    ids. Boundary edges have len 1; interior edges have len 2.
    """
    out: dict[tuple[int, int], list[int]] = defaultdict(list)
    for k, tri in enumerate(elements):
        for i in range(3):
            a = int(tri[i])
            b = int(tri[(i + 1) % 3])
            key = (min(a, b), max(a, b))
            out[key].append(k)
    return out


def _boundary_node_mask(mesh: Fort14Mesh) -> np.ndarray:
    mask = np.zeros(mesh.n_nodes, dtype=bool)
    for seg in mesh.open_boundaries:
        mask[np.asarray(seg, dtype=np.int64)] = True
    for _ib, seg in mesh.land_boundaries:
        mask[np.asarray(seg, dtype=np.int64)] = True
    return mask


# ---------------------------------------------------------------------------
# Operator: smooth_node (move 1 interior vertex to its 1-ring centroid)
# ---------------------------------------------------------------------------


def _apply_smooth_node(
    mesh: Fort14Mesh, vertex_id: int, ring_elem_ids: np.ndarray,
    *, alpha_target: float, min_angle_target: float,
    boundary_node_mask: np.ndarray,
) -> tuple[Fort14Mesh, dict[str, Any]] | None:
    """Greedy smooth: move ``vertex_id`` to its 1-ring centroid if
    that strictly reduces the local penalty AND keeps every signed
    area positive. Returns the updated mesh + info, or ``None`` if
    the move is rejected.
    """
    if boundary_node_mask[vertex_id]:
        return None
    elem_block = mesh.elements[ring_elem_ids]
    neighbours = np.unique(elem_block.ravel())
    neighbours = neighbours[neighbours != vertex_id]
    if neighbours.size == 0:
        return None

    new_xy = mesh.nodes[neighbours].mean(axis=0)
    nodes_proposed = mesh.nodes.copy()
    nodes_proposed[vertex_id] = new_xy

    a_before, m_before = _per_element_quality(mesh.nodes, elem_block)
    a_after, m_after = _per_element_quality(nodes_proposed, elem_block)
    p_before = _penalty(
        a_before, m_before,
        alpha_target=alpha_target, min_angle_target=min_angle_target,
    ).sum()
    p_after = _penalty(
        a_after, m_after,
        alpha_target=alpha_target, min_angle_target=min_angle_target,
    ).sum()
    if p_after + 1e-12 >= p_before:
        return None

    # Signed-area check.
    p0 = nodes_proposed[elem_block[:, 0]]
    p1 = nodes_proposed[elem_block[:, 1]]
    p2 = nodes_proposed[elem_block[:, 2]]
    cross = (
        (p1[:, 0] - p0[:, 0]) * (p2[:, 1] - p0[:, 1])
        - (p1[:, 1] - p0[:, 1]) * (p2[:, 0] - p0[:, 0])
    )
    if (cross <= 0).any():
        return None

    new_mesh = Fort14Mesh(
        title=mesh.title,
        nodes=nodes_proposed,
        depths=mesh.depths.copy(),
        elements=mesh.elements.copy(),
        open_boundaries=[np.asarray(s).copy() for s in mesh.open_boundaries],
        land_boundaries=[(int(ib), np.asarray(s).copy())
                         for ib, s in mesh.land_boundaries],
    )
    return new_mesh, {
        "operator": "smooth_node",
        "vertex": int(vertex_id),
        "penalty_before": float(p_before),
        "penalty_after": float(p_after),
        "affected_elements": [int(x) for x in ring_elem_ids],
    }


# ---------------------------------------------------------------------------
# Operator: edge_swap (Lawson, alpha-driven on the 2-element block)
# ---------------------------------------------------------------------------


def _apply_edge_swap(
    mesh: Fort14Mesh, elem_id: int, edge_local: int,
    *, alpha_target: float, min_angle_target: float,
    edge_uses: dict[tuple[int, int], list[int]],
    boundary_edge_keys: set,
) -> tuple[Fort14Mesh, dict[str, Any]] | None:
    """Swap the shared edge between ``elem_id`` and its buddy across
    edge index ``edge_local``. Accept iff the local penalty drops and
    no signed area goes negative.
    """
    a = int(mesh.elements[elem_id, edge_local])
    b = int(mesh.elements[elem_id, (edge_local + 1) % 3])
    key = (min(a, b), max(a, b))
    if key in boundary_edge_keys:
        return None
    incident = edge_uses.get(key, [])
    if len(incident) != 2:
        return None
    buddy_id = incident[0] if incident[1] == elem_id else incident[1]
    third = int(mesh.elements[elem_id, (edge_local + 2) % 3])
    buddy_set = {int(x) for x in mesh.elements[buddy_id]}
    fourth = (buddy_set - {a, b}).pop()

    block_after = np.array(
        [[a, third, fourth], [b, fourth, third]],
        dtype=mesh.elements.dtype,
    )
    block_before = mesh.elements[[elem_id, buddy_id]]
    a_b, m_b = _per_element_quality(mesh.nodes, block_before)
    a_a, m_a = _per_element_quality(mesh.nodes, block_after)
    p_before = _penalty(
        a_b, m_b,
        alpha_target=alpha_target, min_angle_target=min_angle_target,
    ).sum()
    p_after = _penalty(
        a_a, m_a,
        alpha_target=alpha_target, min_angle_target=min_angle_target,
    ).sum()
    if p_after + 1e-12 >= p_before:
        return None

    # Signed-area check.
    p0 = mesh.nodes[block_after[:, 0]]
    p1 = mesh.nodes[block_after[:, 1]]
    p2 = mesh.nodes[block_after[:, 2]]
    cross = (
        (p1[:, 0] - p0[:, 0]) * (p2[:, 1] - p0[:, 1])
        - (p1[:, 1] - p0[:, 1]) * (p2[:, 0] - p0[:, 0])
    )
    if (cross <= 0).any():
        return None

    new_elements = mesh.elements.copy()
    new_elements[elem_id] = block_after[0]
    new_elements[buddy_id] = block_after[1]
    new_mesh = Fort14Mesh(
        title=mesh.title,
        nodes=mesh.nodes.copy(),
        depths=mesh.depths.copy(),
        elements=new_elements,
        open_boundaries=[np.asarray(s).copy() for s in mesh.open_boundaries],
        land_boundaries=[(int(ib), np.asarray(s).copy())
                         for ib, s in mesh.land_boundaries],
    )
    return new_mesh, {
        "operator": "edge_swap",
        "edge": [int(min(a, b)), int(max(a, b))],
        "elements_modified": [int(elem_id), int(buddy_id)],
        "penalty_before": float(p_before),
        "penalty_after": float(p_after),
    }


# ---------------------------------------------------------------------------
# Operator: edge_split_interior (insert midpoint on internal edge)
# ---------------------------------------------------------------------------


def _split_triangle_at_edge(
    tri: np.ndarray, edge_a: int, edge_b: int, n_new: int,
) -> list[list[int]]:
    """Split CCW triangle ``tri`` at edge ``{edge_a, edge_b}`` by
    inserting ``n_new``. Returns the two CCW sub-triangles. The
    output preserves the parent's orientation regardless of which
    cyclic position the edge appears in.
    """
    v0, v1, v2 = int(tri[0]), int(tri[1]), int(tri[2])
    edge_set = {edge_a, edge_b}
    if {v0, v1} == edge_set:
        a, b, c = v0, v1, v2
    elif {v1, v2} == edge_set:
        a, b, c = v1, v2, v0
    elif {v2, v0} == edge_set:
        a, b, c = v2, v0, v1
    else:
        raise ValueError(
            f"edge ({edge_a}, {edge_b}) not in triangle ({v0}, {v1}, {v2})"
        )
    return [[a, n_new, c], [n_new, b, c]]


def _apply_edge_split_interior(
    mesh: Fort14Mesh, elem_id: int, edge_local: int,
    *, alpha_target: float, min_angle_target: float,
    edge_uses: dict[tuple[int, int], list[int]],
    boundary_edge_keys: set,
) -> tuple[Fort14Mesh, dict[str, Any]] | None:
    """Insert a node at the midpoint of edge ``edge_local`` of element
    ``elem_id``. Replace the two incident triangles (e1, e2) with
    four sub-triangles. Accept iff the per-element penalty over the
    new four is strictly less than over the original two AND no
    signed area is non-positive. Boundary edges are rejected (caller
    is expected to use a boundary-aware variant in v2).
    """
    a = int(mesh.elements[elem_id, edge_local])
    b = int(mesh.elements[elem_id, (edge_local + 1) % 3])
    key = (min(a, b), max(a, b))
    if key in boundary_edge_keys:
        return None
    incident = edge_uses.get(key, [])
    if len(incident) != 2:
        return None
    e1, e2 = incident
    if e1 == elem_id:
        e_self, e_other = e1, e2
    else:
        e_self, e_other = e2, e1

    midpoint = 0.5 * (mesh.nodes[a] + mesh.nodes[b])
    mid_depth = 0.5 * (mesh.depths[a] + mesh.depths[b])
    n_new = mesh.n_nodes  # appended at the end

    new_self = _split_triangle_at_edge(
        mesh.elements[e_self], a, b, n_new,
    )
    new_other = _split_triangle_at_edge(
        mesh.elements[e_other], a, b, n_new,
    )

    block_before = mesh.elements[[e_self, e_other]]
    block_after = np.array(new_self + new_other, dtype=mesh.elements.dtype)
    nodes_proposed = np.vstack([mesh.nodes, midpoint[None, :]])

    a_b, m_b = _per_element_quality(mesh.nodes, block_before)
    a_a, m_a = _per_element_quality(nodes_proposed, block_after)
    p_before = _penalty(
        a_b, m_b,
        alpha_target=alpha_target, min_angle_target=min_angle_target,
    ).sum()
    p_after = _penalty(
        a_a, m_a,
        alpha_target=alpha_target, min_angle_target=min_angle_target,
    ).sum()
    if p_after + 1e-12 >= p_before:
        return None

    # Signed-area check on all 4 new triangles.
    p0 = nodes_proposed[block_after[:, 0]]
    p1 = nodes_proposed[block_after[:, 1]]
    p2 = nodes_proposed[block_after[:, 2]]
    cross = (
        (p1[:, 0] - p0[:, 0]) * (p2[:, 1] - p0[:, 1])
        - (p1[:, 1] - p0[:, 1]) * (p2[:, 0] - p0[:, 0])
    )
    if (cross <= 0).any():
        return None

    # Build new mesh: keep all elements except {e_self, e_other},
    # append the 4 new triangles. The new node is appended at the end.
    keep_mask = np.ones(mesh.n_elements, dtype=bool)
    keep_mask[[e_self, e_other]] = False
    new_elements = np.vstack([
        mesh.elements[keep_mask],
        block_after,
    ])
    new_depths = np.concatenate([mesh.depths, [mid_depth]])
    new_mesh = Fort14Mesh(
        title=mesh.title,
        nodes=nodes_proposed,
        depths=new_depths,
        elements=new_elements,
        open_boundaries=[np.asarray(s).copy() for s in mesh.open_boundaries],
        land_boundaries=[(int(ib), np.asarray(s).copy())
                         for ib, s in mesh.land_boundaries],
    )
    return new_mesh, {
        "operator": "edge_split_interior",
        "edge": [int(min(a, b)), int(max(a, b))],
        "new_node": int(n_new),
        "removed_elements": [int(e_self), int(e_other)],
        "penalty_before": float(p_before),
        "penalty_after": float(p_after),
    }


# ---------------------------------------------------------------------------
# Operator: vertex_remove (remove interior node, retriangulate 1-ring)
# ---------------------------------------------------------------------------


def _apply_vertex_remove(
    mesh: Fort14Mesh, vertex_id: int, ring_elem_ids: np.ndarray,
    *, alpha_target: float, min_angle_target: float,
    boundary_node_mask: np.ndarray,
) -> tuple[Fort14Mesh, dict[str, Any]] | None:
    """Remove ``vertex_id`` (interior only) and re-triangulate the
    enclosing 1-ring polygon via Delaunay (pruned by the rim).
    Single-element variant of the Stage 2 medial-axis re-mesh.
    """
    if boundary_node_mask[vertex_id]:
        return None
    if ring_elem_ids.size < 3:
        return None  # need at least a 3-element 1-ring

    rim_node_ids = _patch_rim_polygon(mesh.elements, ring_elem_ids)
    if rim_node_ids is None:
        return None
    if int(vertex_id) in set(int(x) for x in rim_node_ids):
        # Vertex lies on the patch rim — not a true interior of its
        # own 1-ring (can happen at pinch points / degenerate
        # connectivity). Reject.
        return None

    # Convert rim coords to a metres-equivalent local frame for stable
    # Delaunay near the equator.
    rim_xy = mesh.nodes[rim_node_ids]
    lat_centre = float(rim_xy[:, 1].mean())
    deg_per_m_lat = 1.0 / (EARTH_R_M * np.pi / 180.0)
    deg_per_m_lon = deg_per_m_lat / max(np.cos(np.deg2rad(lat_centre)), 1e-6)
    rim_xy_m = np.column_stack([
        (rim_xy[:, 0] - rim_xy[0, 0]) / deg_per_m_lon,
        (rim_xy[:, 1] - rim_xy[0, 1]) / deg_per_m_lat,
    ])

    # CCW orient the rim using signed area (matches Stage 2).
    sx, sy = rim_xy_m[:, 0], rim_xy_m[:, 1]
    if 0.5 * float(np.sum(sx * np.roll(sy, -1) - np.roll(sx, -1) * sy)) < 0:
        rim_node_ids = rim_node_ids[::-1].copy()
        rim_xy = mesh.nodes[rim_node_ids]
        rim_xy_m = np.column_stack([
            (rim_xy[:, 0] - rim_xy[0, 0]) / deg_per_m_lon,
            (rim_xy[:, 1] - rim_xy[0, 1]) / deg_per_m_lat,
        ])

    # Re-triangulate (rim only — no spine).
    triangles, reason = _retriangulate_patch(
        rim_xy_m, np.empty((0, 2), dtype=float), int(rim_node_ids.size),
    )
    if triangles is None:
        return None

    # Map local triangle indices to global node IDs.
    new_block = rim_node_ids[triangles].astype(mesh.elements.dtype)

    # Score new vs old.
    block_before = mesh.elements[ring_elem_ids]
    a_b, m_b = _per_element_quality(mesh.nodes, block_before)
    a_a, m_a = _per_element_quality(mesh.nodes, new_block)
    p_before = _penalty(
        a_b, m_b,
        alpha_target=alpha_target, min_angle_target=min_angle_target,
    ).sum()
    p_after = _penalty(
        a_a, m_a,
        alpha_target=alpha_target, min_angle_target=min_angle_target,
    ).sum()
    if p_after + 1e-12 >= p_before:
        return None

    # Build new mesh: drop the 1-ring elements, append the new
    # triangulation. The node ``vertex_id`` is no longer referenced;
    # we keep it in the node array (its ID stays valid for any
    # boundary segments that referenced it) but it becomes "orphan".
    # In v1 we rely on the fact that removed interior vertices are
    # never on boundary (we rejected boundary nodes above), so no
    # boundary segment references them.
    keep_mask = np.ones(mesh.n_elements, dtype=bool)
    keep_mask[ring_elem_ids] = False
    new_elements = np.vstack([
        mesh.elements[keep_mask],
        new_block,
    ])
    new_mesh = Fort14Mesh(
        title=mesh.title,
        nodes=mesh.nodes.copy(),
        depths=mesh.depths.copy(),
        elements=new_elements,
        open_boundaries=[np.asarray(s).copy() for s in mesh.open_boundaries],
        land_boundaries=[(int(ib), np.asarray(s).copy())
                         for ib, s in mesh.land_boundaries],
    )
    return new_mesh, {
        "operator": "vertex_remove",
        "vertex": int(vertex_id),
        "rim_size": int(rim_node_ids.size),
        "removed_elements": [int(x) for x in ring_elem_ids],
        "n_new_elements": int(new_block.shape[0]),
        "penalty_before": float(p_before),
        "penalty_after": float(p_after),
    }


# ---------------------------------------------------------------------------
# Driver
# ---------------------------------------------------------------------------


DEFAULT_ALPHA_TARGET: float = 0.95
DEFAULT_MIN_ANGLE_TARGET: float = 20.0
DEFAULT_OPERATOR_ORDER: tuple[str, ...] = (
    "smooth_node",
    "edge_swap",
    "edge_split_interior",
    "vertex_remove",
)
DEFAULT_MAX_SMOOTH_SWEEPS: int = 200


def _batch_smooth_sweep(
    mesh: Fort14Mesh,
    *, alpha_target: float, min_angle_target: float,
    boundary_node_mask: np.ndarray,
    n2e: dict[int, np.ndarray],
) -> int:
    """One Gauss-Seidel pass over the interior vertices: for each
    non-boundary node ``v`` whose 1-ring exists, propose the centroid
    of its 1-ring neighbours and accept the move iff (a) no flipped
    triangle results and (b) the per-1-ring penalty strictly drops.

    Mutates ``mesh.nodes`` in place. Returns the count of accepted
    moves. The caller iterates the sweep until no accepts. Topology
    is unchanged so ``n2e`` and ``boundary_node_mask`` stay valid
    across sweeps; the caller passes them in to avoid rebuild costs.
    """
    accepts = 0
    nodes = mesh.nodes  # mutable view
    for v in np.where(~boundary_node_mask)[0]:
        ring = n2e.get(int(v))
        if ring is None:
            continue
        elem_block = mesh.elements[ring]
        nbrs = np.unique(elem_block.ravel())
        nbrs = nbrs[nbrs != v]
        if nbrs.size == 0:
            continue
        proposed_v = nodes[nbrs].mean(axis=0)

        p0 = nodes[elem_block[:, 0]]
        p1 = nodes[elem_block[:, 1]]
        p2 = nodes[elem_block[:, 2]]
        a_b, m_b, _ts_b = _inline_quality(p0, p1, p2)

        p0p = p0.copy()
        p1p = p1.copy()
        p2p = p2.copy()
        m0 = elem_block[:, 0] == v
        m1 = elem_block[:, 1] == v
        m2 = elem_block[:, 2] == v
        p0p[m0] = proposed_v
        p1p[m1] = proposed_v
        p2p[m2] = proposed_v
        a_a, m_a, ts_a = _inline_quality(p0p, p1p, p2p)
        if (ts_a <= 0).any():
            continue

        p_b = float(_penalty(
            a_b, m_b,
            alpha_target=alpha_target, min_angle_target=min_angle_target,
        ).sum())
        p_a = float(_penalty(
            a_a, m_a,
            alpha_target=alpha_target, min_angle_target=min_angle_target,
        ).sum())
        if p_a + 1e-12 >= p_b:
            continue

        nodes[v] = proposed_v
        accepts += 1
    return accepts


def _topology_round(
    cur: Fort14Mesh,
    *, alpha_target: float, min_angle_target: float,
    operator_order: tuple[str, ...],
    max_topology_accepts: int,
) -> tuple[Fort14Mesh, dict[str, int], int]:
    """Run a single Pass-B round: pop fail elements by descending
    penalty and apply topology operators. ``smooth_node`` is *not*
    in this loop's inventory (it lives in Pass A). Each accepted op
    rebuilds the aux dicts (n2e / eu / bnd_*) on the new mesh.
    Returns ``(updated_mesh, accepts_per_op, abandoned_count)``.
    """
    accepts_per_op: dict[str, int] = defaultdict(int)
    abandoned: set = set()
    accepts_total = 0
    while accepts_total < max_topology_accepts:
        a, m = _per_element_quality(cur.nodes, cur.elements)
        fail = _is_fail(
            a, m,
            alpha_target=alpha_target, min_angle_target=min_angle_target,
        )
        if not fail.any():
            break
        pen = _penalty(
            a, m,
            alpha_target=alpha_target, min_angle_target=min_angle_target,
        )
        heap: list[tuple[float, int]] = []
        for eid in np.where(fail)[0]:
            heapq.heappush(heap, (-float(pen[eid]), int(eid)))
        bnd_node = _boundary_node_mask(cur)
        n2e = _node_to_elements(cur.elements, cur.n_nodes)
        eu = _edge_use_counts(cur.elements)
        bnd_edges = {k for k, v in eu.items() if len(v) == 1}

        progress = False
        while heap:
            _neg_pen, eid = heapq.heappop(heap)
            sig = (
                int(cur.n_elements),
                frozenset(int(x) for x in cur.elements[eid]),
            )
            if sig in abandoned:
                continue

            applied: tuple[Fort14Mesh, dict[str, Any]] | None = None
            for op_name in operator_order:
                if op_name == "smooth_node":
                    continue  # handled in Pass A
                if op_name == "edge_swap":
                    for k in range(3):
                        out = _apply_edge_swap(
                            cur, int(eid), k,
                            alpha_target=alpha_target,
                            min_angle_target=min_angle_target,
                            edge_uses=eu,
                            boundary_edge_keys=bnd_edges,
                        )
                        if out is not None:
                            applied = out
                            break
                elif op_name == "edge_split_interior":
                    for k in range(3):
                        out = _apply_edge_split_interior(
                            cur, int(eid), k,
                            alpha_target=alpha_target,
                            min_angle_target=min_angle_target,
                            edge_uses=eu,
                            boundary_edge_keys=bnd_edges,
                        )
                        if out is not None:
                            applied = out
                            break
                elif op_name == "vertex_remove":
                    for v in cur.elements[eid]:
                        ring = n2e.get(int(v))
                        if ring is None:
                            continue
                        out = _apply_vertex_remove(
                            cur, int(v), ring,
                            alpha_target=alpha_target,
                            min_angle_target=min_angle_target,
                            boundary_node_mask=bnd_node,
                        )
                        if out is not None:
                            applied = out
                            break
                else:
                    raise ValueError(f"unknown operator: {op_name!r}")
                if applied is not None:
                    break

            if applied is None:
                abandoned.add(sig)
                continue

            cur = applied[0]
            accepts_per_op[applied[1]["operator"]] += 1
            accepts_total += 1
            progress = True
            break  # rebuild aux dicts

        if not progress:
            break

    return cur, dict(accepts_per_op), len(abandoned)


def phase_h_optimize(
    mesh: Fort14Mesh,
    *,
    alpha_target: float = DEFAULT_ALPHA_TARGET,
    min_angle_target: float = DEFAULT_MIN_ANGLE_TARGET,
    max_smooth_sweeps: int = DEFAULT_MAX_SMOOTH_SWEEPS,
    max_topology_per_round: int = 10_000,
    max_outer_rounds: int = 10,
    operator_order: tuple[str, ...] = DEFAULT_OPERATOR_ORDER,
) -> tuple[Fort14Mesh, dict[str, Any]]:
    """Phase H driver: alternating Pass A (batch smooth) ↔ Pass B
    (topology operators) until both stop making progress.

    **Pass A** — batch Gauss-Seidel smooth. Each sweep visits every
    interior vertex, proposes the 1-ring centroid, and accepts the
    move iff the per-1-ring penalty strictly decreases without
    flipping a triangle. Topology is fixed across the sweeps so
    the per-mesh aux structures are built once and reused.

    **Pass B** — per-element greedy with topology operators only
    (``edge_swap``, ``edge_split_interior``, ``vertex_remove``).
    Pops fail elements by descending penalty and applies the first
    operator that strictly reduces the local penalty without
    flipping. Each accept rebuilds the aux dicts on the new mesh.

    The two passes alternate: Pass A runs to exhaustion, Pass B
    runs up to ``max_topology_per_round`` accepts, then Pass A
    cleans up the local perturbations Pass B introduced. The loop
    terminates when both passes contribute zero accepts in the
    same outer round, or when ``max_outer_rounds`` is reached.

    Returns ``(new_mesh, info)`` with operator histograms, sweep
    counts, and pre/post quality summaries.
    """
    cur = Fort14Mesh(
        title=mesh.title,
        nodes=mesh.nodes.copy(),
        depths=mesh.depths.copy(),
        elements=mesh.elements.copy(),
        open_boundaries=[np.asarray(s).copy() for s in mesh.open_boundaries],
        land_boundaries=[(int(ib), np.asarray(s).copy())
                         for ib, s in mesh.land_boundaries],
    )

    info: dict[str, Any] = {
        "alpha_target": float(alpha_target),
        "min_angle_target": float(min_angle_target),
        "max_smooth_sweeps": int(max_smooth_sweeps),
        "max_topology_per_round": int(max_topology_per_round),
        "max_outer_rounds": int(max_outer_rounds),
        "operator_order": list(operator_order),
        "operators_applied": defaultdict(int),
        "n_iters": 0,
        "n_smooth_sweeps": 0,
        "n_outer_rounds": 0,
        "n_abandoned": 0,
    }

    a0, m0 = _per_element_quality(cur.nodes, cur.elements)
    fail0 = _is_fail(
        a0, m0, alpha_target=alpha_target, min_angle_target=min_angle_target,
    )
    info["n_input_fail"] = int(fail0.sum())
    info["alpha_mean_before"] = float(a0.mean()) if a0.size else float("nan")

    do_smooth = "smooth_node" in operator_order
    topology_ops = tuple(
        op for op in operator_order if op != "smooth_node"
    )

    for outer_round in range(max_outer_rounds):
        info["n_outer_rounds"] = outer_round + 1
        round_accepts = 0

        # Pass A: batch Gauss-Seidel smooth.
        if do_smooth:
            bnd_node = _boundary_node_mask(cur)
            n2e = _node_to_elements(cur.elements, cur.n_nodes)
            for _sweep in range(max_smooth_sweeps):
                n_acc = _batch_smooth_sweep(
                    cur,
                    alpha_target=alpha_target,
                    min_angle_target=min_angle_target,
                    boundary_node_mask=bnd_node,
                    n2e=n2e,
                )
                info["operators_applied"]["smooth_node"] += int(n_acc)
                info["n_iters"] += int(n_acc)
                info["n_smooth_sweeps"] += 1
                round_accepts += int(n_acc)
                if n_acc == 0:
                    break

        # Pass B: topology operators.
        if topology_ops:
            cur, topo_acc, n_aband = _topology_round(
                cur,
                alpha_target=alpha_target,
                min_angle_target=min_angle_target,
                operator_order=operator_order,
                max_topology_accepts=max_topology_per_round,
            )
            for op_name, n in topo_acc.items():
                info["operators_applied"][op_name] += int(n)
                info["n_iters"] += int(n)
                round_accepts += int(n)
            info["n_abandoned"] = n_aband

        if round_accepts == 0:
            break

    a_after, m_after = _per_element_quality(cur.nodes, cur.elements)
    fail_after = _is_fail(
        a_after, m_after,
        alpha_target=alpha_target, min_angle_target=min_angle_target,
    )
    info["n_output_fail"] = int(fail_after.sum())
    info["alpha_mean_after"] = (
        float(a_after.mean()) if a_after.size else float("nan")
    )
    info["n_nodes"] = int(cur.n_nodes)
    info["n_elements"] = int(cur.n_elements)
    info["operators_applied"] = dict(info["operators_applied"])
    return cur, info


__all__ = [
    "DEFAULT_ALPHA_TARGET",
    "DEFAULT_MIN_ANGLE_TARGET",
    "DEFAULT_OPERATOR_ORDER",
    "phase_h_optimize",
]
