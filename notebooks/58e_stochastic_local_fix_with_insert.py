"""PoC #58e: stochastic local fixer + insert operator.

Extension of PoC #58d (move-only): adds the random-barycentric
vertex insertion operator to escape the move-only fixed point
that left 5 stuck residuals (1 C1 + 4 C4) on the PoC #57 Stage 1
input.

  rng = numpy.random.default_rng(seed=42)
  for each fail element / fail edge in priority order:
      for try_i in range(max_tries):
          op = rng.choice(["move_vertex", "insert_vertex"],
                          p=[0.7, 0.3])
          if op == "move_vertex":
              propose Gaussian perturbation of a random vertex of
              the fail element; validate on the moved vertex's 1-ring
              (= every element whose geometry changes)
          else:  # insert_vertex
              propose a new node at a random barycentric position
              inside the fail element; split the element into 3
              sub-triangles; validate on the 3 sub-triangles plus
              their external buddies' C4-affected edges
          if all C1 / C2 / C4 / C5 + signed-area gates pass on the
                  patch AND the original fail element no longer
                  fails:
              accept; break
          else:
              revert and try again

Operator definitions:

  move_vertex:
    1. pick a vertex uniformly from the fail element's 3 vertices
    2. interior nodes: 2-D Gaussian perturbation, σ = 0.3 × mean
       local edge length
    3. boundary nodes: 1-D Gaussian along the tangent line clamped
       to [0.05, 0.95] of the segment, projected onto the coastline

  insert_vertex:
    1. sample barycentric (u, v, w) on the open 2-simplex with
       min(u, v, w) >= 0.1 so sub-triangles are not degenerate
    2. append the new node to ``mesh.nodes`` and depth (depth =
       barycentric blend of the 3 corner depths)
    3. overwrite the original element with sub1 = (v0, v1, new) and
       append sub2 = (v1, v2, new), sub3 = (v2, v0, new)
    4. validate the 3 new sub-triangles + check C4 on each of their
       6 incident internal edges (3 new spokes + 3 original perimeter
       edges, each against its external buddy)
    5. C5: new vertex has valence 3 (always OK), v0/v1/v2 valence
       each gain 1

On insert accept, the topology aux maps (n2e, edge_uses, valence)
are refreshed in place so the outer pass sees the new mesh
consistently.

Outputs:
    outputs/58e_stochastic_with_insert.14
    outputs/58e_summary.{txt,json}
"""
from __future__ import annotations

import json
import time
from collections import defaultdict
from pathlib import Path

import numpy as np

from fvcom_mesh_tools.algorithms.quality import min_interior_angle
from fvcom_mesh_tools.diagnostics import node_valence
from fvcom_mesh_tools.io import Fort14Mesh, read_fort14, write_fort14
from fvcom_mesh_tools.mesh_clean_phase_h import (
    BOUNDARY_TANGENT_T_MIN,
    BOUNDARY_TANGENT_T_MAX,
    _boundary_node_mask,
    _boundary_topology,
    _edge_use_counts,
    _per_edge_area_change,
    _per_element_quality,
    build_coastline_projector,
)

REPO = Path(__file__).resolve().parent.parent
INPUT = REPO / "outputs" / "57_phase_h_pass_f_only.14"
COASTLINE = (
    REPO / "data" / "coastline" / "tokyo_bay"
    / "MLIT_C23" / "C23-06_TOKYOBAY.shp"
)
OUT_DIR = REPO / "outputs"
OUTPUT = OUT_DIR / "58e_stochastic_with_insert.14"
SUMMARY_TXT = OUT_DIR / "58e_summary.txt"
SUMMARY_JSON = OUT_DIR / "58e_summary.json"

SEED = 42
MIN_ANGLE_TARGET = 30.0
MAX_ANGLE_TARGET = 130.0
AREA_RATIO_TARGET = 0.5
MAX_VALENCE = 8

MAX_TRIES_PER_FAIL = 500
PERTURBATION_SIGMA = 0.30   # fraction of mean local edge length
MAX_OUTER_PASSES = 5        # restart from worst residual N times
OPERATOR_MOVE_WEIGHT = 0.7   # 70 % move, 30 % insert
BARYCENTRIC_MIN = 0.10      # min(u, v, w) floor → avoid degenerate sub-tris


def _max_interior_angle(mesh: Fort14Mesh) -> np.ndarray:
    p0 = mesh.nodes[mesh.elements[:, 0]]
    p1 = mesh.nodes[mesh.elements[:, 1]]
    p2 = mesh.nodes[mesh.elements[:, 2]]
    e0 = np.linalg.norm(p1 - p0, axis=1)
    e1 = np.linalg.norm(p2 - p1, axis=1)
    e2 = np.linalg.norm(p0 - p2, axis=1)

    def _ang(opp, ea, eb):
        denom = 2.0 * ea * eb
        safe = np.where(denom > 0, denom, 1.0)
        cos = np.where(
            denom > 0, (ea * ea + eb * eb - opp * opp) / safe, 0.0,
        )
        return np.arccos(np.clip(cos, -1.0, 1.0))

    return np.degrees(
        np.maximum(
            np.maximum(_ang(e1, e2, e0), _ang(e2, e0, e1)),
            _ang(e0, e1, e2),
        ),
    )


def _global_metrics(mesh: Fort14Mesh) -> dict:
    m = min_interior_angle(mesh)
    M = _max_interior_angle(mesh)
    _u, _p, ac = _per_edge_area_change(mesh.nodes, mesh.elements)
    val = node_valence(mesh.elements, mesh.n_nodes)
    return {
        "NP": int(mesh.n_nodes),
        "NE": int(mesh.n_elements),
        "C1": int((m < MIN_ANGLE_TARGET).sum()),
        "C2": int((M > MAX_ANGLE_TARGET).sum()),
        "C4": int((ac > AREA_RATIO_TARGET).sum()),
        "C5": int((val > MAX_VALENCE).sum()),
    }


def _vertex_ring_patch(
    n2e: dict[int, np.ndarray], v: int,
) -> np.ndarray:
    """Return the element IDs in the 1-ring of vertex ``v`` — i.e.,
    every element containing ``v``. This is the complete set of
    elements whose geometry changes when ``v`` moves."""
    ring = n2e.get(int(v))
    if ring is None:
        return np.empty(0, dtype=np.int64)
    return np.asarray(ring, dtype=np.int64)


def _affected_internal_edge_pairs(
    mesh: Fort14Mesh, patch_eids: np.ndarray,
    edge_uses: dict[tuple[int, int], list[int]],
) -> np.ndarray:
    """Return ``(M, 2)`` array of (e_i, e_j) element pairs sharing an
    internal edge such that at least one of e_i, e_j is in
    ``patch_eids``."""
    patch_set = set(int(e) for e in patch_eids)
    pairs: list[tuple[int, int]] = []
    seen: set[tuple[int, int]] = set()
    for e in patch_eids:
        tri = mesh.elements[int(e)]
        for k in range(3):
            a = int(tri[k])
            b = int(tri[(k + 1) % 3])
            key = (min(a, b), max(a, b))
            if key in seen:
                continue
            buds = edge_uses.get(key, [])
            if len(buds) != 2:
                continue
            seen.add(key)
            if int(buds[0]) in patch_set or int(buds[1]) in patch_set:
                pairs.append((int(buds[0]), int(buds[1])))
    if not pairs:
        return np.empty((0, 2), dtype=np.int64)
    return np.asarray(pairs, dtype=np.int64)


def _local_quality_ok(
    mesh: Fort14Mesh, patch_eids: np.ndarray,
    affected_pairs: np.ndarray,
    valence: np.ndarray,
) -> tuple[bool, dict]:
    """Check C1/C2/C4/C5 + signed area on the local patch. Returns
    ``(ok, diagnostics)``."""
    if patch_eids.size == 0:
        return True, {"reason": "empty_patch"}
    block = mesh.elements[patch_eids]
    p0 = mesh.nodes[block[:, 0]]
    p1 = mesh.nodes[block[:, 1]]
    p2 = mesh.nodes[block[:, 2]]
    cross = (
        (p1[:, 0] - p0[:, 0]) * (p2[:, 1] - p0[:, 1])
        - (p1[:, 1] - p0[:, 1]) * (p2[:, 0] - p0[:, 0])
    )
    if (cross <= 0).any():
        return False, {"reason": "flipped_triangle"}
    _alpha, min_ang, max_ang = _per_element_quality(mesh.nodes, block)
    c1 = int((min_ang < MIN_ANGLE_TARGET).sum())
    c2 = int((max_ang > MAX_ANGLE_TARGET).sum())
    if affected_pairs.shape[0]:
        pair_elems = mesh.elements[affected_pairs.ravel()].reshape(
            affected_pairs.shape[0], 2, 3,
        )
        pi = mesh.nodes[pair_elems[:, 0]]
        pj = mesh.nodes[pair_elems[:, 1]]
        ai = np.abs(0.5 * (
            (pi[:, 1, 0] - pi[:, 0, 0]) * (pi[:, 2, 1] - pi[:, 0, 1])
            - (pi[:, 1, 1] - pi[:, 0, 1]) * (pi[:, 2, 0] - pi[:, 0, 0])
        ))
        aj = np.abs(0.5 * (
            (pj[:, 1, 0] - pj[:, 0, 0]) * (pj[:, 2, 1] - pj[:, 0, 1])
            - (pj[:, 1, 1] - pj[:, 0, 1]) * (pj[:, 2, 0] - pj[:, 0, 0])
        ))
        larger = np.maximum(ai, aj)
        smaller = np.minimum(ai, aj)
        ratio = (larger - smaller) / np.maximum(larger, 1e-30)
        c4 = int((ratio > AREA_RATIO_TARGET).sum())
    else:
        c4 = 0
    # C5 is topology-invariant under move; check anyway for completeness.
    patch_verts = np.unique(block.ravel())
    c5 = int((valence[patch_verts] > MAX_VALENCE).sum())
    ok = (c1 == 0 and c2 == 0 and c4 == 0 and c5 == 0)
    return ok, {
        "c1_local": c1, "c2_local": c2, "c4_local": c4, "c5_local": c5,
        "reason": "ok" if ok else "criterion_fail",
    }


def _propose_move(
    mesh: Fort14Mesh, eid: int, rng,
    *,
    boundary_node_mask: np.ndarray,
    boundary_prev: np.ndarray,
    boundary_next: np.ndarray,
    coastline_projector,
    n2e: dict[int, np.ndarray],
) -> tuple[int, np.ndarray, np.ndarray] | None:
    """Propose a random node move. Returns
    ``(vertex_id, new_pos, ring_eids)`` or ``None`` if no candidate
    move can be constructed. ``ring_eids`` is the 1-ring of the
    chosen vertex — the patch on which to validate the move."""
    tri = mesh.elements[int(eid)]
    v = int(rng.choice(tri))
    ring_eids = _vertex_ring_patch(n2e, v)
    if ring_eids.size == 0:
        return None

    # Local mean edge length (over the vertex's 1-ring) as scale.
    block = mesh.elements[ring_eids]
    pa = mesh.nodes[block[:, 0]]
    pb = mesh.nodes[block[:, 1]]
    pc = mesh.nodes[block[:, 2]]
    edge_lens = np.concatenate([
        np.linalg.norm(pb - pa, axis=1),
        np.linalg.norm(pc - pb, axis=1),
        np.linalg.norm(pa - pc, axis=1),
    ])
    local_h = float(edge_lens.mean())
    if local_h <= 0:
        return None

    sigma = PERTURBATION_SIGMA * local_h
    if boundary_node_mask[v]:
        pv = int(boundary_prev[v])
        nv = int(boundary_next[v])
        if pv < 0 or nv < 0:
            return None  # corner, skip
        a_pos = mesh.nodes[pv]
        b_pos = mesh.nodes[nv]
        ab = b_pos - a_pos
        ab_len = float(np.linalg.norm(ab))
        if ab_len < 1e-20:
            return None
        # Current parameter along the tangent line.
        t_now = float(((mesh.nodes[v] - a_pos) @ ab) / (ab @ ab))
        # Gaussian perturbation in the parameter; scale by σ / ab_len so
        # the spatial sigma is ~PERTURBATION_SIGMA × local_h.
        delta_t = float(rng.standard_normal()) * (sigma / ab_len)
        t_new = max(
            BOUNDARY_TANGENT_T_MIN,
            min(BOUNDARY_TANGENT_T_MAX, t_now + delta_t),
        )
        new_pos = a_pos + t_new * ab
        if coastline_projector is not None:
            snapped = coastline_projector(new_pos)
            if snapped is not None:
                new_pos = snapped
    else:
        delta = rng.standard_normal(2) * sigma
        new_pos = mesh.nodes[v] + delta
    return v, np.asarray(new_pos, dtype=mesh.nodes.dtype), ring_eids


def _propose_insert(
    mesh: Fort14Mesh, eid: int, rng,
    *,
    edge_uses: dict[tuple[int, int], list[int]],
) -> dict | None:
    """Propose a random-barycentric vertex insertion inside element
    ``eid``. Returns a dict describing the proposed operation:

        {
            "new_pos": np.ndarray (2,),
            "new_depth": float,
            "v0": v0, "v1": v1, "v2": v2,
            "sub1": (v0, v1, new), "sub2": (v1, v2, new),
            "sub3": (v2, v0, new),
            "external_buddies": dict mapping each original perimeter
                edge to its external buddy id (or -1 if boundary),
        }

    or ``None`` if no insertion can be constructed. Does NOT mutate
    the mesh. Caller is responsible for applying via
    :func:`_apply_insert` after validation.
    """
    tri = mesh.elements[int(eid)]
    v0, v1, v2 = int(tri[0]), int(tri[1]), int(tri[2])
    # Uniform barycentric sample with min-coordinate floor.
    while True:
        u = float(rng.random())
        v = float(rng.random())
        if u + v > 1.0:
            u, v = 1.0 - u, 1.0 - v
        w = 1.0 - u - v
        if min(u, v, w) >= BARYCENTRIC_MIN:
            break
    new_pos = (
        u * mesh.nodes[v0] + v * mesh.nodes[v1] + w * mesh.nodes[v2]
    )
    new_depth = (
        u * mesh.depths[v0] + v * mesh.depths[v1] + w * mesh.depths[v2]
    )

    # Identify external buddies across the 3 original edges of eid.
    external_buddies: dict[tuple[int, int], int] = {}
    for a, b in ((v0, v1), (v1, v2), (v2, v0)):
        key = (min(a, b), max(a, b))
        buds = edge_uses.get(key, [])
        if len(buds) == 2:
            be = buds[0] if buds[1] == int(eid) else buds[1]
            external_buddies[(a, b)] = int(be)
        else:
            external_buddies[(a, b)] = -1  # boundary edge
    return {
        "new_pos": np.asarray(new_pos, dtype=mesh.nodes.dtype),
        "new_depth": float(new_depth),
        "v0": v0, "v1": v1, "v2": v2,
        "external_buddies": external_buddies,
    }


def _insert_local_quality_ok(
    mesh: Fort14Mesh, *,
    proposal: dict,
    new_idx: int,
    valence: np.ndarray,
) -> bool:
    """Validate a proposed insert *without mutating the mesh*.

    The 3 sub-triangles are evaluated against C1 / C2 and signed
    area; the 3 original perimeter edges of the eid being split are
    evaluated against C4 between the sub-triangle on the inside and
    the unchanged external buddy on the outside; and the 3 new
    spokes (from new vertex to each original corner) are evaluated
    against C4 between the two sub-triangles they separate. C5:
    new vertex has valence 3, each of v0/v1/v2 gains 1 valence.
    """
    v0, v1, v2 = proposal["v0"], proposal["v1"], proposal["v2"]
    new_pos = proposal["new_pos"]

    # Sub-triangle vertex coordinates (proposed).
    p_v0 = mesh.nodes[v0]
    p_v1 = mesh.nodes[v1]
    p_v2 = mesh.nodes[v2]
    p_new = new_pos
    sub_corners = [
        (p_v0, p_v1, p_new),  # sub1
        (p_v1, p_v2, p_new),  # sub2
        (p_v2, p_v0, p_new),  # sub3
    ]

    # Signed area + C1/C2 on each sub-triangle.
    sub_areas: list[float] = []
    for pa, pb, pc in sub_corners:
        cross = (
            (pb[0] - pa[0]) * (pc[1] - pa[1])
            - (pb[1] - pa[1]) * (pc[0] - pa[0])
        )
        if cross <= 0:
            return False
        a_sub = float(0.5 * cross)
        sub_areas.append(a_sub)
        # Angles
        e_ab = float(np.linalg.norm(np.asarray(pb) - np.asarray(pa)))
        e_bc = float(np.linalg.norm(np.asarray(pc) - np.asarray(pb)))
        e_ca = float(np.linalg.norm(np.asarray(pa) - np.asarray(pc)))
        if min(e_ab, e_bc, e_ca) <= 0:
            return False
        # Use the law of cosines for each angle.
        def _angle(opp, ea, eb):
            denom = 2.0 * ea * eb
            if denom <= 0:
                return 180.0
            cos = (ea * ea + eb * eb - opp * opp) / denom
            cos = max(-1.0, min(1.0, cos))
            return float(np.degrees(np.arccos(cos)))
        angles = [
            _angle(e_bc, e_ab, e_ca),  # at pa
            _angle(e_ca, e_ab, e_bc),  # at pb
            _angle(e_ab, e_bc, e_ca),  # at pc
        ]
        if min(angles) < MIN_ANGLE_TARGET:
            return False
        if max(angles) > MAX_ANGLE_TARGET:
            return False

    # C4 on the 3 perimeter edges of the original element (sub vs
    # external buddy). external_buddies[(va, vb)] == -1 ⇒ boundary,
    # no C4 metric.
    sub_for_edge = {
        (v0, v1): 0, (v1, v0): 0,
        (v1, v2): 1, (v2, v1): 1,
        (v2, v0): 2, (v0, v2): 2,
    }
    for (a, b), buddy_id in proposal["external_buddies"].items():
        if buddy_id < 0:
            continue
        # Buddy area is unchanged.
        buddy_tri = mesh.elements[int(buddy_id)]
        bp = [mesh.nodes[int(buddy_tri[k])] for k in range(3)]
        cross_b = (
            (bp[1][0] - bp[0][0]) * (bp[2][1] - bp[0][1])
            - (bp[1][1] - bp[0][1]) * (bp[2][0] - bp[0][0])
        )
        a_buddy = float(abs(0.5 * cross_b))
        a_sub = sub_areas[sub_for_edge[(a, b)]]
        larger = max(a_sub, a_buddy)
        smaller = min(a_sub, a_buddy)
        ratio = (larger - smaller) / max(larger, 1e-30)
        if ratio > AREA_RATIO_TARGET:
            return False

    # C4 on the 3 spokes (new vertex to each corner). Each spoke
    # connects 2 sub-triangles.
    spoke_pairs = [
        (0, 2),  # (v0, new) → sub1 and sub3
        (0, 1),  # (v1, new) → sub1 and sub2
        (1, 2),  # (v2, new) → sub2 and sub3
    ]
    for i, j in spoke_pairs:
        larger = max(sub_areas[i], sub_areas[j])
        smaller = min(sub_areas[i], sub_areas[j])
        ratio = (larger - smaller) / max(larger, 1e-30)
        if ratio > AREA_RATIO_TARGET:
            return False

    # C5 check: each corner gains 1 valence. New vertex has valence 3.
    if int(valence[v0]) + 1 > MAX_VALENCE:
        return False
    if int(valence[v1]) + 1 > MAX_VALENCE:
        return False
    if int(valence[v2]) + 1 > MAX_VALENCE:
        return False

    return True


def _apply_insert(
    mesh: Fort14Mesh, eid: int, proposal: dict,
    *,
    n2e: dict[int, np.ndarray],
    edge_uses: dict[tuple[int, int], list[int]],
    valence: np.ndarray,
) -> tuple[Fort14Mesh, int, np.ndarray]:
    """Apply a validated insert. Mutates the topology aux maps in
    place. Returns ``(new_mesh, new_vertex_idx, updated_valence)``.

    The mesh fields are reassigned because ``np.ndarray`` cannot be
    grown in place; the caller must rebind any local mesh reference
    from the returned object.
    """
    v0, v1, v2 = proposal["v0"], proposal["v1"], proposal["v2"]
    new_pos = proposal["new_pos"]
    new_depth = proposal["new_depth"]
    new_idx = int(mesh.n_nodes)
    ne_old = int(mesh.n_elements)
    sub1 = np.array([v0, v1, new_idx], dtype=mesh.elements.dtype)
    sub2 = np.array([v1, v2, new_idx], dtype=mesh.elements.dtype)
    sub3 = np.array([v2, v0, new_idx], dtype=mesh.elements.dtype)

    new_nodes = np.vstack([mesh.nodes, new_pos[None, :]])
    new_depths = np.concatenate([mesh.depths, [new_depth]])
    new_elements = np.vstack([
        mesh.elements,
        sub2[None, :], sub3[None, :],
    ])
    new_elements[int(eid)] = sub1
    new_mesh = Fort14Mesh(
        title=mesh.title,
        nodes=new_nodes,
        depths=new_depths,
        elements=new_elements,
        open_boundaries=[np.asarray(s).copy() for s in mesh.open_boundaries],
        land_boundaries=[
            (int(ib), np.asarray(s).copy())
            for ib, s in mesh.land_boundaries
        ],
    )

    # Update n2e in place. eid stays in v0 and v1 (sub1), leaves v2.
    sub2_id = ne_old      # appended first
    sub3_id = ne_old + 1  # appended second
    n2e[v0] = np.concatenate([n2e[v0], np.array([sub3_id])])
    n2e[v1] = np.concatenate([n2e[v1], np.array([sub2_id])])
    # v2: remove eid, add sub2 and sub3
    old_v2_ring = n2e[v2]
    n2e[v2] = np.concatenate([
        old_v2_ring[old_v2_ring != int(eid)],
        np.array([sub2_id, sub3_id]),
    ])
    n2e[new_idx] = np.array(
        [int(eid), sub2_id, sub3_id], dtype=np.int64,
    )

    # Update edge_uses.
    # Original perimeter edges: (v0, v1), (v1, v2), (v2, v0).
    #   (v0, v1): sub1 contains this edge. Element id is now eid.
    #             Unchanged (eid was on this edge before; still is).
    #   (v1, v2): now in sub2 (= sub2_id), not in eid.
    #             Replace eid → sub2_id in this edge's list.
    #   (v2, v0): now in sub3 (= sub3_id), not in eid.
    #             Replace eid → sub3_id in this edge's list.
    def _replace(key_list: list[int], old: int, new: int) -> list[int]:
        return [new if x == old else x for x in key_list]
    edge_uses[(min(v1, v2), max(v1, v2))] = _replace(
        edge_uses[(min(v1, v2), max(v1, v2))], int(eid), sub2_id,
    )
    edge_uses[(min(v2, v0), max(v2, v0))] = _replace(
        edge_uses[(min(v2, v0), max(v2, v0))], int(eid), sub3_id,
    )
    # New spoke edges.
    edge_uses[(min(v0, new_idx), max(v0, new_idx))] = [
        int(eid), sub3_id,
    ]
    edge_uses[(min(v1, new_idx), max(v1, new_idx))] = [
        int(eid), sub2_id,
    ]
    edge_uses[(min(v2, new_idx), max(v2, new_idx))] = [
        sub2_id, sub3_id,
    ]

    # Update valence.
    new_valence = np.append(valence, 3)
    new_valence[v0] += 1
    new_valence[v1] += 1
    new_valence[v2] += 1
    return new_mesh, new_idx, new_valence


def _eid_still_fails(
    mesh: Fort14Mesh, eid: int,
    edge_uses: dict[tuple[int, int], list[int]],
) -> bool:
    """Element-level fail check: C1 (min_angle) + C2 (max_angle) + C4
    (area_change vs each internal-edge buddy)."""
    block = mesh.elements[[int(eid)]]
    _a, m, M = _per_element_quality(mesh.nodes, block)
    if float(m[0]) < MIN_ANGLE_TARGET:
        return True
    if float(M[0]) > MAX_ANGLE_TARGET:
        return True
    tri = mesh.elements[int(eid)]
    for k in range(3):
        a = int(tri[k])
        b = int(tri[(k + 1) % 3])
        key = (min(a, b), max(a, b))
        buds = edge_uses.get(key, [])
        if len(buds) != 2:
            continue
        be = buds[0] if buds[1] == int(eid) else buds[1]
        block2 = mesh.elements[[int(eid), int(be)]]
        p0 = mesh.nodes[block2[:, 0]]
        p1 = mesh.nodes[block2[:, 1]]
        p2 = mesh.nodes[block2[:, 2]]
        areas = np.abs(0.5 * (
            (p1[:, 0] - p0[:, 0]) * (p2[:, 1] - p0[:, 1])
            - (p1[:, 1] - p0[:, 1]) * (p2[:, 0] - p0[:, 0])
        ))
        larger = areas.max()
        smaller = areas.min()
        if (larger - smaller) / max(larger, 1e-30) > AREA_RATIO_TARGET:
            return True
    return False


def _process_fail_element(
    mesh: Fort14Mesh, eid: int, rng,
    *,
    boundary_node_mask: np.ndarray,
    boundary_prev: np.ndarray,
    boundary_next: np.ndarray,
    coastline_projector,
    valence: np.ndarray,
    n2e: dict[int, np.ndarray],
    edge_uses: dict[tuple[int, int], list[int]],
) -> tuple[Fort14Mesh, np.ndarray, bool, int, int]:
    """Try up to MAX_TRIES_PER_FAIL random operations on element
    ``eid``. Each try randomly picks "move" (probability
    ``OPERATOR_MOVE_WEIGHT``) or "insert" (the complement) and
    validates against the appropriate local patch.

    Returns ``(updated_mesh, updated_valence, fixed, n_tries_used,
    n_inserts_accepted)``. The returned mesh and valence may be new
    objects if an insert was accepted (which grows the arrays).
    """
    if not _eid_still_fails(mesh, eid, edge_uses):
        return mesh, valence, True, 0, 0

    n_inserts = 0
    for try_i in range(1, MAX_TRIES_PER_FAIL + 1):
        op = "move" if rng.random() < OPERATOR_MOVE_WEIGHT else "insert"

        if op == "move":
            proposal = _propose_move(
                mesh, eid, rng,
                boundary_node_mask=boundary_node_mask,
                boundary_prev=boundary_prev,
                boundary_next=boundary_next,
                coastline_projector=coastline_projector,
                n2e=n2e,
            )
            if proposal is None:
                continue
            v, new_pos, ring_eids = proposal
            affected_pairs = _affected_internal_edge_pairs(
                mesh, ring_eids, edge_uses,
            )
            old_pos = mesh.nodes[v].copy()
            mesh.nodes[v] = new_pos
            ok, _diag = _local_quality_ok(
                mesh, ring_eids, affected_pairs, valence,
            )
            if ok and not _eid_still_fails(mesh, eid, edge_uses):
                return mesh, valence, True, try_i, n_inserts
            mesh.nodes[v] = old_pos

        else:  # insert
            proposal = _propose_insert(
                mesh, eid, rng, edge_uses=edge_uses,
            )
            if proposal is None:
                continue
            ok = _insert_local_quality_ok(
                mesh, proposal=proposal,
                new_idx=int(mesh.n_nodes), valence=valence,
            )
            if not ok:
                continue
            # Validated — apply.
            mesh, _new_idx, valence = _apply_insert(
                mesh, eid, proposal,
                n2e=n2e, edge_uses=edge_uses, valence=valence,
            )
            n_inserts += 1
            if not _eid_still_fails(mesh, eid, edge_uses):
                return mesh, valence, True, try_i, n_inserts
            # Insert applied but target somehow still fails (e.g.
            # because C4 against a buddy regressed after re-test on
            # the actual mesh state). This is rare and we keep
            # trying with the now-updated mesh.

    return mesh, valence, False, MAX_TRIES_PER_FAIL, n_inserts


def main() -> int:
    if not INPUT.exists():
        raise SystemExit(f"input missing: {INPUT}")
    OUT_DIR.mkdir(parents=True, exist_ok=True)

    mesh = read_fort14(INPUT)
    before = _global_metrics(mesh)
    print(
        f"[58e] input : NP={before['NP']:,} NE={before['NE']:,} "
        f"C1={before['C1']} C2={before['C2']} "
        f"C4={before['C4']} C5={before['C5']}",
        flush=True,
    )

    projector = build_coastline_projector(
        [COASTLINE],
        max_snap_distance_m=500.0,
        mean_latitude_deg=float(mesh.nodes[:, 1].mean()),
    )
    rng = np.random.default_rng(SEED)

    fix_records: list[dict] = []
    t_total = time.perf_counter()

    for outer in range(1, MAX_OUTER_PASSES + 1):
        m = min_interior_angle(mesh)
        M = _max_interior_angle(mesh)
        _u, _p, ac = _per_edge_area_change(mesh.nodes, mesh.elements)
        c1_fail = m < MIN_ANGLE_TARGET
        c2_fail = M > MAX_ANGLE_TARGET
        # Map C4 fail edges -> their incident elements.
        edge_uv, elem_pair, _ac_full = _per_edge_area_change(
            mesh.nodes, mesh.elements,
        )
        c4_edge_fail = _ac_full > AREA_RATIO_TARGET
        c4_fail_elems = set()
        for pair in elem_pair[c4_edge_fail]:
            c4_fail_elems.add(int(pair[0]))
            c4_fail_elems.add(int(pair[1]))

        fail_set = (
            set(int(e) for e in np.where(c1_fail)[0])
            | set(int(e) for e in np.where(c2_fail)[0])
            | c4_fail_elems
        )
        if not fail_set:
            print(f"[58e] outer pass {outer}: all fails cleared",
                  flush=True)
            break

        # Sort fail elements deterministically (by id) for reproducibility.
        fail_list = sorted(fail_set)
        print(
            f"[58e] outer pass {outer}: {len(fail_list)} fail elements"
            f" (C1={int(c1_fail.sum())}, C2={int(c2_fail.sum())}, "
            f"C4_elems={len(c4_fail_elems)})",
            flush=True,
        )

        bnd_node = _boundary_node_mask(mesh)
        bnd_prev, bnd_next, _ = _boundary_topology(mesh)
        valence_arr = node_valence(mesh.elements, mesh.n_nodes)
        edge_uses = _edge_use_counts(mesh.elements)
        n2e_map: dict[int, np.ndarray] = defaultdict(
            lambda: np.empty(0, dtype=np.int64),
        )
        # Build node-to-element map.
        from collections import defaultdict as _defaultdict
        tmp: dict[int, list[int]] = _defaultdict(list)
        for k, tri in enumerate(mesh.elements):
            tmp[int(tri[0])].append(k)
            tmp[int(tri[1])].append(k)
            tmp[int(tri[2])].append(k)
        n2e_map = {
            int(v): np.asarray(es, dtype=np.int64)
            for v, es in tmp.items()
        }

        outer_fixed = 0
        outer_stuck = 0
        outer_tries_sum = 0
        outer_inserts = 0
        for eid in fail_list:
            mesh, valence_arr, fixed, ntries, n_ins = (
                _process_fail_element(
                    mesh, eid, rng,
                    boundary_node_mask=bnd_node,
                    boundary_prev=bnd_prev,
                    boundary_next=bnd_next,
                    coastline_projector=projector,
                    valence=valence_arr,
                    n2e=n2e_map,
                    edge_uses=edge_uses,
                )
            )
            fix_records.append({
                "outer": outer, "eid": int(eid),
                "fixed": bool(fixed), "n_tries": int(ntries),
                "n_inserts": int(n_ins),
            })
            outer_tries_sum += ntries
            outer_inserts += n_ins
            if fixed:
                outer_fixed += 1
            else:
                outer_stuck += 1

        avg_tries = (
            outer_tries_sum / max(1, outer_fixed + outer_stuck)
        )
        print(
            f"[58e] outer pass {outer}: fixed={outer_fixed} "
            f"stuck={outer_stuck} avg_tries={avg_tries:.1f} "
            f"inserts={outer_inserts}",
            flush=True,
        )
        if outer_fixed == 0:
            print(f"[58e] no progress in outer pass {outer}; "
                  "terminating", flush=True)
            break

    wall = time.perf_counter() - t_total
    after = _global_metrics(mesh)
    write_fort14(mesh, OUTPUT)
    n_inserts_total = sum(r.get("n_inserts", 0) for r in fix_records)
    print(
        f"[58e] output: NP={after['NP']:,} NE={after['NE']:,} "
        f"C1={after['C1']} C2={after['C2']} "
        f"C4={after['C4']} C5={after['C5']} "
        f"(wall {wall:.1f} s, total_inserts={n_inserts_total})",
        flush=True,
    )

    payload = {
        "input": str(INPUT.resolve()),
        "output": str(OUTPUT.resolve()),
        "seed": SEED,
        "max_tries_per_fail": MAX_TRIES_PER_FAIL,
        "perturbation_sigma": PERTURBATION_SIGMA,
        "max_outer_passes": MAX_OUTER_PASSES,
        "operator_move_weight": OPERATOR_MOVE_WEIGHT,
        "barycentric_min": BARYCENTRIC_MIN,
        "before": before,
        "after": after,
        "delta": {k: after[k] - before[k] for k in before},
        "wall_seconds": wall,
        "n_records": len(fix_records),
        "n_fixed_total": sum(1 for r in fix_records if r["fixed"]),
        "n_stuck_total": sum(1 for r in fix_records if not r["fixed"]),
        "n_inserts_total": n_inserts_total,
    }
    SUMMARY_JSON.write_text(
        json.dumps(payload, indent=2, default=str), encoding="utf-8",
    )

    total_b = before["C1"] + before["C2"] + before["C4"] + before["C5"]
    total_a = after["C1"] + after["C2"] + after["C4"] + after["C5"]
    lines = [
        f"PoC #58e — stochastic local fixer "
        f"(move {OPERATOR_MOVE_WEIGHT:.0%} / insert "
        f"{1-OPERATOR_MOVE_WEIGHT:.0%}, seed=42)",
        f"input : {INPUT.name}",
        f"output: {OUTPUT.name}",
        "",
        f"  {'stage':<28} | {'NP':>6} | {'NE':>6} | "
        f"{'C1':>4} | {'C2':>4} | {'C4':>4} | {'C5':>3} | total",
        "  " + "-" * 82,
        f"  {'PoC #57 Stage 1 (input)':<28} | "
        f"{before['NP']:>6,} | {before['NE']:>6,} | "
        f"{before['C1']:>4} | {before['C2']:>4} | "
        f"{before['C4']:>4} | {before['C5']:>3} | {total_b}",
        f"  {'PoC #58e (after fixer)':<28} | "
        f"{after['NP']:>6,} | {after['NE']:>6,} | "
        f"{after['C1']:>4} | {after['C2']:>4} | "
        f"{after['C4']:>4} | {after['C5']:>3} | {total_a}",
        "",
        f"  wall          : {wall:.2f} s",
        f"  n_fix_records : {len(fix_records)}",
        f"  n_fixed       : {payload['n_fixed_total']}",
        f"  n_stuck       : {payload['n_stuck_total']}",
        f"  n_inserts     : {n_inserts_total}",
        f"  delta total   : {total_a - total_b:+d}",
    ]
    SUMMARY_TXT.write_text("\n".join(lines) + "\n", encoding="utf-8")
    print()
    print("\n".join(lines))
    print()
    print(f"wrote {OUTPUT}")
    print(f"wrote {SUMMARY_TXT}")
    print(f"wrote {SUMMARY_JSON}")
    return 0


if __name__ == "__main__":
    raise SystemExit(main())
