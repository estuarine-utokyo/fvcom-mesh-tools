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
from pathlib import Path
from typing import Any, Callable

import numpy as np

from fvcom_mesh_tools.io import Fort14Mesh

# Re-use the Stage 2 retriangulation helpers.
from fvcom_mesh_tools.mesh_clean import (
    _patch_rim_polygon,
    _retriangulate_patch,
)

CoastlineProjector = Callable[[np.ndarray], np.ndarray | None]
"""``(xy: np.ndarray of shape (2,)) -> projected np.ndarray or None``.

Returns the input point projected onto the nearest coastline polyline
when it lies within the configured snap distance, else ``None`` so
the caller falls back to the un-projected position.
"""

EARTH_R_M: float = 6_371_000.0


def build_coastline_projector(
    coastline_paths: list[Path] | list[str] | None,
    *,
    max_snap_distance_m: float = 500.0,
    mean_latitude_deg: float | None = None,
) -> CoastlineProjector | None:
    """Construct a projector callable that snaps a point in EPSG:4326
    onto the nearest coastline polyline within ``max_snap_distance_m``.

    Loads ``coastline_paths`` (shapefile / GeoJSON / any GeoPandas-
    readable vector source) and assembles every ``LineString``,
    ``MultiLineString``, ``Polygon`` boundary, and ``MultiPolygon``
    boundary into a flat polyline list. A :class:`shapely.STRtree`
    indexes the polylines for O(log N) nearest-neighbour lookup.

    ``max_snap_distance_m`` is converted to degrees at
    ``mean_latitude_deg`` (or the mean latitude of the union of
    polyline bounding boxes if ``None``) — this is the longitudinal
    metres-per-degree, which is the shorter axis at non-equatorial
    latitudes and therefore the conservative bound on the snap
    distance the projector enforces in lon-lat space.

    Returns ``None`` when ``coastline_paths`` is empty or no
    geometries load successfully (the caller then falls back to
    midpoint behaviour).
    """
    if not coastline_paths:
        return None
    try:
        import geopandas as gpd  # noqa: I001
        from shapely import STRtree
        from shapely.geometry import Point
    except ImportError as e:
        raise RuntimeError(
            "build_coastline_projector requires geopandas + shapely; "
            f"install the [io-vector] extra. Underlying error: {e}"
        ) from e

    polylines: list = []
    for raw_path in coastline_paths:
        path = Path(raw_path)
        gdf = gpd.read_file(path)
        # Reproject to EPSG:4326 if metadata says otherwise.
        if gdf.crs is not None and gdf.crs.to_epsg() != 4326:
            gdf = gdf.to_crs(4326)
        for geom in gdf.geometry:
            if geom is None or geom.is_empty:
                continue
            gt = geom.geom_type
            if gt == "LineString":
                polylines.append(geom)
            elif gt == "MultiLineString":
                polylines.extend(list(geom.geoms))
            elif gt == "Polygon":
                polylines.append(geom.boundary)
            elif gt == "MultiPolygon":
                for sub in geom.geoms:
                    polylines.append(sub.boundary)
            # Other types (Point / MultiPoint) are ignored.

    if not polylines:
        return None

    if mean_latitude_deg is None:
        # Mean latitude across all polylines' bounding boxes.
        lat_centres = []
        for ls in polylines:
            (_minx, miny, _maxx, maxy) = ls.bounds
            lat_centres.append(0.5 * (miny + maxy))
        mean_latitude_deg = float(np.mean(lat_centres))

    deg_per_m_lon = 1.0 / (
        EARTH_R_M
        * max(np.cos(np.deg2rad(mean_latitude_deg)), 1e-6)
        * np.pi / 180.0
    )
    max_snap_deg = float(max_snap_distance_m) * deg_per_m_lon

    tree = STRtree(polylines)

    def _project(xy: np.ndarray) -> np.ndarray | None:
        pt = Point(float(xy[0]), float(xy[1]))
        idx = tree.nearest(pt)
        if idx is None:
            return None
        polyline = polylines[int(idx)]
        if polyline.distance(pt) > max_snap_deg:
            return None
        proj = polyline.interpolate(polyline.project(pt))
        return np.asarray([proj.x, proj.y], dtype=float)

    return _project


# ---------------------------------------------------------------------------
# Penalty / quality helpers
# ---------------------------------------------------------------------------


#: Default ``max_angle_target``. 180° makes the ``max_ang > target``
#: check vacuous (every triangle has max angle <= 180°), so the gate
#: is a no-op unless the caller overrides — preserves backward
#: compatibility for existing tests / PoCs while exposing the FVCOM
#: criterion ``max_angle <= 130°`` as a tunable parameter.
DEFAULT_MAX_ANGLE_TARGET: float = 180.0

#: Default per-edge area-change target for Pass E (gradation
#: refinement). The FVCOM manual lists 0.5 as the upper bound on
#: ``(max(A_i, A_j) - min(A_i, A_j)) / max(A_i, A_j)`` between
#: adjacent triangles (criterion C4). 1.0 makes the check vacuous,
#: so Pass E is a no-op until the caller lowers the target.
DEFAULT_AREA_RATIO_TARGET: float = 0.5

#: Default upper bound on per-node valence (FVCOM manual criterion
#: C5: each node may participate in at most 8 elements). Pass E
#: gates against this because ``edge_split_interior`` raises the
#: valence of the two "opposite" vertices by 1 (and
#: ``edge_split_boundary`` raises one), so unconstrained splits at
#: nodes already at valence 8 would regress C5.
DEFAULT_MAX_VALENCE: int = 8


def _per_element_quality(
    nodes: np.ndarray, elements: np.ndarray,
) -> tuple[np.ndarray, np.ndarray, np.ndarray]:
    """Per-element ``(alpha, min_angle_deg, max_angle_deg)`` arrays."""
    if elements.size == 0:
        return np.empty(0), np.empty(0), np.empty(0)
    p0 = nodes[elements[:, 0]]
    p1 = nodes[elements[:, 1]]
    p2 = nodes[elements[:, 2]]
    alpha, min_ang, max_ang, _twice = _inline_quality(p0, p1, p2)
    return alpha, min_ang, max_ang


def _inline_quality(
    p0: np.ndarray, p1: np.ndarray, p2: np.ndarray,
) -> tuple[np.ndarray, np.ndarray, np.ndarray, np.ndarray]:
    """Vectorised per-triangle
    ``(alpha, min_angle_deg, max_angle_deg, twice_signed)``.

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
    # Interior angles via law of cosines on each vertex.
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
    max_ang_rad = np.maximum(np.maximum(angle_v0, angle_v1), angle_v2)
    return (
        alpha,
        np.degrees(min_ang_rad),
        np.degrees(max_ang_rad),
        twice_signed,
    )


def _penalty(
    alpha: np.ndarray, min_ang: np.ndarray, max_ang: np.ndarray | None = None,
    *, alpha_target: float, min_angle_target: float,
    max_angle_target: float = DEFAULT_MAX_ANGLE_TARGET,
) -> np.ndarray:
    """Element penalty: zero iff every gate met. Squared deficits with
    each angle term scaled by 1/100 so all contributions stay in the
    same dynamic range. ``max_ang`` may be ``None`` for backward-
    compatible callers that did not yet supply max-angle arrays; the
    max-angle term is then omitted.
    """
    a_pen = np.maximum(0.0, alpha_target - alpha) ** 2
    g_pen = np.maximum(0.0, min_angle_target - min_ang) ** 2 / 100.0
    pen = a_pen + g_pen
    if max_ang is not None:
        G_pen = np.maximum(0.0, max_ang - max_angle_target) ** 2 / 100.0
        pen = pen + G_pen
    return pen


def _is_fail(
    alpha: np.ndarray, min_ang: np.ndarray, max_ang: np.ndarray | None = None,
    *, alpha_target: float, min_angle_target: float,
    max_angle_target: float = DEFAULT_MAX_ANGLE_TARGET,
) -> np.ndarray:
    mask = (alpha < alpha_target) | (min_ang < min_angle_target)
    if max_ang is not None:
        mask = mask | (max_ang > max_angle_target)
    return mask


def _signed_areas(nodes: np.ndarray, elements: np.ndarray) -> np.ndarray:
    """Per-element signed area (positive for CCW)."""
    if elements.size == 0:
        return np.empty(0)
    p0 = nodes[elements[:, 0]]
    p1 = nodes[elements[:, 1]]
    p2 = nodes[elements[:, 2]]
    return 0.5 * (
        (p1[:, 0] - p0[:, 0]) * (p2[:, 1] - p0[:, 1])
        - (p1[:, 1] - p0[:, 1]) * (p2[:, 0] - p0[:, 0])
    )


def _internal_edge_buddies(
    elements: np.ndarray,
) -> tuple[np.ndarray, np.ndarray]:
    """Return ``(edge_uv, elem_pair)`` for every internal edge.

    An "internal" edge is shared by exactly two elements. ``edge_uv``
    is an ``(N_int, 2)`` sorted-vertex array; ``elem_pair`` is the
    ``(N_int, 2)`` element-id pair sharing that edge. Boundary edges
    (incident on exactly one element) are skipped.
    """
    if elements.size == 0:
        empty_uv = np.empty((0, 2), dtype=np.int64)
        empty_pair = np.empty((0, 2), dtype=np.int64)
        return empty_uv, empty_pair
    NE = elements.shape[0]
    edges = np.vstack([
        elements[:, [0, 1]],
        elements[:, [1, 2]],
        elements[:, [2, 0]],
    ])
    edges_sorted = np.sort(edges, axis=1)
    elem_ids = np.tile(np.arange(NE, dtype=np.int64), 3)
    # Lex-sort by (u, v) so duplicate edges sit adjacent.
    order = np.lexsort((edges_sorted[:, 1], edges_sorted[:, 0]))
    edges_sorted = edges_sorted[order]
    elem_ids = elem_ids[order]
    # Two consecutive rows with identical (u, v) = one internal edge.
    same = (
        (edges_sorted[:-1, 0] == edges_sorted[1:, 0])
        & (edges_sorted[:-1, 1] == edges_sorted[1:, 1])
    )
    internal_idx = np.where(same)[0]
    edge_uv = edges_sorted[internal_idx]
    elem_pair = np.stack(
        [elem_ids[internal_idx], elem_ids[internal_idx + 1]], axis=1,
    )
    return edge_uv.astype(np.int64), elem_pair.astype(np.int64)


def _per_edge_area_change(
    nodes: np.ndarray, elements: np.ndarray,
) -> tuple[np.ndarray, np.ndarray, np.ndarray]:
    """Compute ``(edge_uv, elem_pair, area_change)`` for every
    internal edge.

    Formula matches PoC #48 / the FVCOM manual:
        area_change = (max(|A_i|, |A_j|) - min(|A_i|, |A_j|))
                      / max(|A_i|, |A_j|)
    where ``A_i`` and ``A_j`` are the signed areas of the two
    adjacent elements. Result is in ``[0, 1]``; 0 = same size, 1 =
    one element vanishingly small.
    """
    edge_uv, elem_pair = _internal_edge_buddies(elements)
    if edge_uv.shape[0] == 0:
        return edge_uv, elem_pair, np.empty(0)
    areas = np.abs(_signed_areas(nodes, elements))
    a_i = areas[elem_pair[:, 0]]
    a_j = areas[elem_pair[:, 1]]
    larger = np.maximum(a_i, a_j)
    smaller = np.minimum(a_i, a_j)
    area_change = (larger - smaller) / np.maximum(larger, 1e-30)
    return edge_uv, elem_pair, area_change


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


def _boundary_topology(
    mesh: Fort14Mesh,
) -> tuple[np.ndarray, np.ndarray, dict[tuple[int, int], tuple[str, int, int, int]]]:
    """Build per-node boundary tangent neighbours + per-edge segment
    map.

    Returns ``(bnd_prev, bnd_next, edge_to_segment)`` where:

    * ``bnd_prev[v]`` and ``bnd_next[v]`` are the boundary-tangent
      neighbours of node ``v`` along its segment (``-1`` if ``v`` is
      not interior to any segment — either off-boundary, at a
      segment endpoint, or repeated across multiple segments).
    * ``edge_to_segment`` maps each topological boundary edge
      ``(min(a, b), max(a, b))`` to ``(kind, seg_idx, position,
      land_ibtype)`` where ``kind`` is ``"open"`` or ``"land"``,
      ``seg_idx`` indexes into ``mesh.open_boundaries`` /
      ``mesh.land_boundaries``, ``position`` is the index ``j`` such
      that ``seg[j], seg[j+1] == a, b`` (in segment order, not the
      sorted edge key), and ``land_ibtype`` is the ibtype (``-1`` for
      open segments).
    """
    n_nodes = mesh.n_nodes
    bnd_prev = np.full(n_nodes, -1, dtype=np.int64)
    bnd_next = np.full(n_nodes, -1, dtype=np.int64)
    seen_v: dict[int, int] = {}  # node-id → count of segments it appears in
    edge_to_segment: dict[tuple[int, int], tuple[str, int, int, int]] = {}

    def _walk(seg_list, kind: str, ibtype_of_seg):
        for seg_idx, seg in enumerate(seg_list):
            ibtype = ibtype_of_seg(seg_idx)
            arr = np.asarray(seg, dtype=np.int64)
            for j in range(arr.size):
                v = int(arr[j])
                seen_v[v] = seen_v.get(v, 0) + 1
                if j > 0:
                    bnd_prev[v] = int(arr[j - 1])
                if j < arr.size - 1:
                    bnd_next[v] = int(arr[j + 1])
                    a = int(arr[j])
                    b = int(arr[j + 1])
                    key = (min(a, b), max(a, b))
                    edge_to_segment[key] = (kind, seg_idx, j, ibtype)

    _walk(mesh.open_boundaries, "open", lambda _i: -1)
    _walk([s for _ib, s in mesh.land_boundaries], "land",
          lambda i: int(mesh.land_boundaries[i][0]))

    # Nodes appearing in multiple segments are corners; clear their
    # tangent neighbours so smoothers refuse to move them.
    for v, count in seen_v.items():
        if count > 1:
            bnd_prev[v] = -1
            bnd_next[v] = -1

    return bnd_prev, bnd_next, edge_to_segment


# ---------------------------------------------------------------------------
# Operator: smooth_node (move 1 interior vertex to its 1-ring centroid)
# ---------------------------------------------------------------------------


def _apply_smooth_node(
    mesh: Fort14Mesh, vertex_id: int, ring_elem_ids: np.ndarray,
    *, alpha_target: float, min_angle_target: float,
    max_angle_target: float = DEFAULT_MAX_ANGLE_TARGET,
    boundary_node_mask: np.ndarray,
    force: bool = False,
) -> tuple[Fort14Mesh, dict[str, Any]] | None:
    """Greedy smooth: move ``vertex_id`` to its 1-ring centroid if
    that strictly reduces the local penalty AND keeps every signed
    area positive. Returns the updated mesh + info, or ``None`` if
    the move is rejected. ``force=True`` skips the penalty gate
    (validity checks remain) — used by the 2-step lookahead driver
    to evaluate op1 candidates that do not directly improve their
    local 1-ring.
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

    a_before, m_before, M_before = _per_element_quality(mesh.nodes, elem_block)
    a_after, m_after, M_after = _per_element_quality(nodes_proposed, elem_block)
    p_before = _penalty(
        a_before, m_before, M_before,
        alpha_target=alpha_target, min_angle_target=min_angle_target,
        max_angle_target=max_angle_target,
    ).sum()
    p_after = _penalty(
        a_after, m_after, M_after,
        alpha_target=alpha_target, min_angle_target=min_angle_target,
        max_angle_target=max_angle_target,
    ).sum()
    if not force and p_after + 1e-12 >= p_before:
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
        "forced": bool(force),
    }


# ---------------------------------------------------------------------------
# Operator: edge_swap (Lawson, alpha-driven on the 2-element block)
# ---------------------------------------------------------------------------


def _apply_edge_swap(
    mesh: Fort14Mesh, elem_id: int, edge_local: int,
    *, alpha_target: float, min_angle_target: float,
    max_angle_target: float = DEFAULT_MAX_ANGLE_TARGET,
    edge_uses: dict[tuple[int, int], list[int]],
    boundary_edge_keys: set,
    force: bool = False,
) -> tuple[Fort14Mesh, dict[str, Any]] | None:
    """Swap the shared edge between ``elem_id`` and its buddy across
    edge index ``edge_local``. Accept iff the local penalty drops and
    no signed area goes negative. ``force=True`` skips the penalty
    gate (validity checks remain) — see ``_apply_smooth_node``.
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
    a_b, m_b, M_b = _per_element_quality(mesh.nodes, block_before)
    a_a, m_a, M_a = _per_element_quality(mesh.nodes, block_after)
    p_before = _penalty(
        a_b, m_b, M_b,
        alpha_target=alpha_target, min_angle_target=min_angle_target,
        max_angle_target=max_angle_target,
    ).sum()
    p_after = _penalty(
        a_a, m_a, M_a,
        alpha_target=alpha_target, min_angle_target=min_angle_target,
        max_angle_target=max_angle_target,
    ).sum()
    if not force and p_after + 1e-12 >= p_before:
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
        "forced": bool(force),
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
    max_angle_target: float = DEFAULT_MAX_ANGLE_TARGET,
    edge_uses: dict[tuple[int, int], list[int]],
    boundary_edge_keys: set,
    force: bool = False,
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

    a_b, m_b, M_b = _per_element_quality(mesh.nodes, block_before)
    a_a, m_a, M_a = _per_element_quality(nodes_proposed, block_after)
    p_before = _penalty(
        a_b, m_b, M_b,
        alpha_target=alpha_target, min_angle_target=min_angle_target,
        max_angle_target=max_angle_target,
    ).sum()
    p_after = _penalty(
        a_a, m_a, M_a,
        alpha_target=alpha_target, min_angle_target=min_angle_target,
        max_angle_target=max_angle_target,
    ).sum()
    if not force and p_after + 1e-12 >= p_before:
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
        "forced": bool(force),
    }


# ---------------------------------------------------------------------------
# Operator: edge_split_boundary (insert midpoint on a boundary edge)
# ---------------------------------------------------------------------------


def _apply_edge_split_boundary(
    mesh: Fort14Mesh, elem_id: int, edge_local: int,
    *, alpha_target: float, min_angle_target: float,
    max_angle_target: float = DEFAULT_MAX_ANGLE_TARGET,
    edge_uses: dict[tuple[int, int], list[int]],
    edge_to_segment: dict[tuple[int, int], tuple[str, int, int, int]],
    coastline_projector: CoastlineProjector | None = None,
    force: bool = False,
) -> tuple[Fort14Mesh, dict[str, Any]] | None:
    """Insert a node at the midpoint of a boundary edge of ``elem_id``.
    Replace the single incident triangle with two sub-triangles and
    update the boundary segment to thread the new node between its
    two endpoints.

    When ``coastline_projector`` is supplied, the straight chord
    midpoint is projected onto the nearest coastline polyline within
    the projector's configured snap distance (v3); when projection
    is unavailable or out-of-range, the straight midpoint stands (v2
    behaviour). Accept iff per-element penalty drops and no signed
    area becomes non-positive.
    """
    a = int(mesh.elements[elem_id, edge_local])
    b = int(mesh.elements[elem_id, (edge_local + 1) % 3])
    key = (min(a, b), max(a, b))
    seg_meta = edge_to_segment.get(key)
    if seg_meta is None:
        return None  # not a topological boundary edge
    incident = edge_uses.get(key, [])
    if len(incident) != 1:
        return None  # geometric boundary disagrees — skip rather than corrupt
    if incident[0] != elem_id:
        return None  # caller passed an edge that does not include elem_id

    kind, seg_idx, position, ibtype = seg_meta
    midpoint = 0.5 * (mesh.nodes[a] + mesh.nodes[b])
    snapped_to_coastline = False
    if coastline_projector is not None:
        snapped = coastline_projector(midpoint)
        if snapped is not None:
            midpoint = snapped
            snapped_to_coastline = True
    mid_depth = 0.5 * (mesh.depths[a] + mesh.depths[b])
    n_new = mesh.n_nodes

    new_self = _split_triangle_at_edge(
        mesh.elements[elem_id], a, b, n_new,
    )

    block_before = mesh.elements[[elem_id]]
    block_after = np.array(new_self, dtype=mesh.elements.dtype)
    nodes_proposed = np.vstack([mesh.nodes, midpoint[None, :]])

    a_b, m_b, M_b = _per_element_quality(mesh.nodes, block_before)
    a_a, m_a, M_a = _per_element_quality(nodes_proposed, block_after)
    p_before = _penalty(
        a_b, m_b, M_b,
        alpha_target=alpha_target, min_angle_target=min_angle_target,
        max_angle_target=max_angle_target,
    ).sum()
    p_after = _penalty(
        a_a, m_a, M_a,
        alpha_target=alpha_target, min_angle_target=min_angle_target,
        max_angle_target=max_angle_target,
    ).sum()
    if not force and p_after + 1e-12 >= p_before:
        return None

    # Signed-area check on the new sub-triangles.
    p0 = nodes_proposed[block_after[:, 0]]
    p1 = nodes_proposed[block_after[:, 1]]
    p2 = nodes_proposed[block_after[:, 2]]
    cross = (
        (p1[:, 0] - p0[:, 0]) * (p2[:, 1] - p0[:, 1])
        - (p1[:, 1] - p0[:, 1]) * (p2[:, 0] - p0[:, 0])
    )
    if (cross <= 0).any():
        return None

    # Update the segment array: insert n_new between positions
    # ``position`` and ``position + 1``.
    new_open = [np.asarray(s).copy() for s in mesh.open_boundaries]
    new_land = [(int(ib), np.asarray(s).copy())
                for ib, s in mesh.land_boundaries]
    if kind == "open":
        seg = new_open[seg_idx]
        new_open[seg_idx] = np.concatenate(
            [seg[: position + 1], [n_new], seg[position + 1:]]
        ).astype(seg.dtype)
    else:
        ib, seg = new_land[seg_idx]
        new_seg = np.concatenate(
            [seg[: position + 1], [n_new], seg[position + 1:]]
        ).astype(seg.dtype)
        new_land[seg_idx] = (ib, new_seg)

    keep_mask = np.ones(mesh.n_elements, dtype=bool)
    keep_mask[elem_id] = False
    new_elements = np.vstack([mesh.elements[keep_mask], block_after])
    new_depths = np.concatenate([mesh.depths, [mid_depth]])
    new_mesh = Fort14Mesh(
        title=mesh.title,
        nodes=nodes_proposed,
        depths=new_depths,
        elements=new_elements,
        open_boundaries=new_open,
        land_boundaries=new_land,
    )
    return new_mesh, {
        "operator": "edge_split_boundary",
        "edge": [int(min(a, b)), int(max(a, b))],
        "new_node": int(n_new),
        "boundary_kind": kind,
        "segment_idx": int(seg_idx),
        "land_ibtype": int(ibtype),
        "removed_elements": [int(elem_id)],
        "penalty_before": float(p_before),
        "penalty_after": float(p_after),
        "snapped_to_coastline": bool(snapped_to_coastline),
        "forced": bool(force),
    }


# ---------------------------------------------------------------------------
# Operator: vertex_remove (remove interior node, retriangulate 1-ring)
# ---------------------------------------------------------------------------


def _apply_vertex_remove(
    mesh: Fort14Mesh, vertex_id: int, ring_elem_ids: np.ndarray,
    *, alpha_target: float, min_angle_target: float,
    max_angle_target: float = DEFAULT_MAX_ANGLE_TARGET,
    boundary_node_mask: np.ndarray,
    force: bool = False,
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
    a_b, m_b, M_b = _per_element_quality(mesh.nodes, block_before)
    a_a, m_a, M_a = _per_element_quality(mesh.nodes, new_block)
    p_before = _penalty(
        a_b, m_b, M_b,
        alpha_target=alpha_target, min_angle_target=min_angle_target,
        max_angle_target=max_angle_target,
    ).sum()
    p_after = _penalty(
        a_a, m_a, M_a,
        alpha_target=alpha_target, min_angle_target=min_angle_target,
        max_angle_target=max_angle_target,
    ).sum()
    if not force and p_after + 1e-12 >= p_before:
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
        "forced": bool(force),
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
    "edge_split_boundary",
    "vertex_remove",
)
DEFAULT_MAX_SMOOTH_SWEEPS: int = 200


BOUNDARY_TANGENT_T_MIN: float = 0.05
BOUNDARY_TANGENT_T_MAX: float = 0.95


def _batch_smooth_sweep(
    mesh: Fort14Mesh,
    *, alpha_target: float, min_angle_target: float,
    max_angle_target: float = DEFAULT_MAX_ANGLE_TARGET,
    boundary_node_mask: np.ndarray,
    n2e: dict[int, np.ndarray],
    boundary_prev: np.ndarray | None = None,
    boundary_next: np.ndarray | None = None,
    coastline_projector: CoastlineProjector | None = None,
) -> int:
    """One Gauss-Seidel pass over every node:

    * **Interior** (``~boundary_node_mask``): propose the 1-ring
      centroid.
    * **Boundary, segment-interior**: project the 1-ring centroid
      onto the line ``[bnd_prev[v], bnd_next[v]]`` and clamp the
      parameter to ``(BOUNDARY_TANGENT_T_MIN, BOUNDARY_TANGENT_T_MAX)``
      so the moving node cannot collapse onto either tangent
      neighbour. Requires ``boundary_prev`` / ``boundary_next``;
      otherwise boundary nodes are skipped (v1 behaviour).
    * **Corners** (no tangent neighbour, or shared by multiple
      segments): skipped.

    Each candidate is accepted iff (a) no flipped triangle results
    and (b) the per-1-ring penalty strictly drops. Mutates
    ``mesh.nodes`` in place. Topology is unchanged so the aux
    structures the caller passes in stay valid across sweeps.
    """
    accepts = 0
    nodes = mesh.nodes  # mutable view
    has_tangent = (
        boundary_prev is not None and boundary_next is not None
    )
    for v in range(mesh.n_nodes):
        ring = n2e.get(int(v))
        if ring is None:
            continue
        is_bnd = bool(boundary_node_mask[v])
        if is_bnd:
            if not has_tangent:
                continue
            pv = int(boundary_prev[v])
            nv = int(boundary_next[v])
            if pv < 0 or nv < 0:
                continue  # corner
            a_pos = nodes[pv]
            b_pos = nodes[nv]
            ab = b_pos - a_pos
            ab_len_sq = float((ab * ab).sum())
            if ab_len_sq < 1e-20:
                continue
            elem_block = mesh.elements[ring]
            nbrs = np.unique(elem_block.ravel())
            nbrs = nbrs[nbrs != v]
            if nbrs.size == 0:
                continue
            ring_centroid = nodes[nbrs].mean(axis=0)
            t = float(((ring_centroid - a_pos) @ ab) / ab_len_sq)
            t = max(BOUNDARY_TANGENT_T_MIN, min(BOUNDARY_TANGENT_T_MAX, t))
            proposed_v = a_pos + t * ab
            if coastline_projector is not None:
                snapped = coastline_projector(proposed_v)
                if snapped is not None:
                    proposed_v = snapped
        else:
            elem_block = mesh.elements[ring]
            nbrs = np.unique(elem_block.ravel())
            nbrs = nbrs[nbrs != v]
            if nbrs.size == 0:
                continue
            proposed_v = nodes[nbrs].mean(axis=0)

        p0 = nodes[elem_block[:, 0]]
        p1 = nodes[elem_block[:, 1]]
        p2 = nodes[elem_block[:, 2]]
        a_b, m_b, M_b, _ts_b = _inline_quality(p0, p1, p2)

        p0p = p0.copy()
        p1p = p1.copy()
        p2p = p2.copy()
        m0 = elem_block[:, 0] == v
        m1 = elem_block[:, 1] == v
        m2 = elem_block[:, 2] == v
        p0p[m0] = proposed_v
        p1p[m1] = proposed_v
        p2p[m2] = proposed_v
        a_a, m_a, M_a, ts_a = _inline_quality(p0p, p1p, p2p)
        if (ts_a <= 0).any():
            continue

        p_b = float(_penalty(
            a_b, m_b, M_b,
            alpha_target=alpha_target, min_angle_target=min_angle_target,
            max_angle_target=max_angle_target,
        ).sum())
        p_a = float(_penalty(
            a_a, m_a, M_a,
            alpha_target=alpha_target, min_angle_target=min_angle_target,
            max_angle_target=max_angle_target,
        ).sum())
        if p_a + 1e-12 >= p_b:
            continue

        nodes[v] = proposed_v
        accepts += 1
    return accepts


def _topology_round(
    cur: Fort14Mesh,
    *, alpha_target: float, min_angle_target: float,
    max_angle_target: float = DEFAULT_MAX_ANGLE_TARGET,
    operator_order: tuple[str, ...],
    max_topology_accepts: int,
    coastline_projector: CoastlineProjector | None = None,
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
        a, m, M = _per_element_quality(cur.nodes, cur.elements)
        fail = _is_fail(
            a, m, M,
            alpha_target=alpha_target, min_angle_target=min_angle_target,
            max_angle_target=max_angle_target,
        )
        if not fail.any():
            break
        pen = _penalty(
            a, m, M,
            alpha_target=alpha_target, min_angle_target=min_angle_target,
            max_angle_target=max_angle_target,
        )
        heap: list[tuple[float, int]] = []
        for eid in np.where(fail)[0]:
            heapq.heappush(heap, (-float(pen[eid]), int(eid)))
        bnd_node = _boundary_node_mask(cur)
        n2e = _node_to_elements(cur.elements, cur.n_nodes)
        eu = _edge_use_counts(cur.elements)
        bnd_edges = {k for k, v in eu.items() if len(v) == 1}
        _bp, _bn, e2s = _boundary_topology(cur)

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

                            max_angle_target=max_angle_target,
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

                            max_angle_target=max_angle_target,
                            edge_uses=eu,
                            boundary_edge_keys=bnd_edges,
                        )
                        if out is not None:
                            applied = out
                            break
                elif op_name == "edge_split_boundary":
                    for k in range(3):
                        out = _apply_edge_split_boundary(
                            cur, int(eid), k,
                            alpha_target=alpha_target,
                            min_angle_target=min_angle_target,

                            max_angle_target=max_angle_target,
                            edge_uses=eu,
                            edge_to_segment=e2s,
                            coastline_projector=coastline_projector,
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

                            max_angle_target=max_angle_target,
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


#: Default inventory for Pass C op1. Restricted to ``smooth_node``
#: because:
#:
#:   * Under the v4.1 ``target_exits_fail`` gate, destructive ops
#:     (``vertex_remove``, ``edge_split_*``, ``edge_swap``) erase
#:     the target's vertex set so the gate rejects them by
#:     construction. Including them only burns compute.
#:   * Under the legacy v4 ``union_penalty`` gate, PoC #44 measured
#:     ``vertex_remove + smooth_node`` to be 53 % of accepted pairs
#:     in dry-run but only 0.8 % in the iterative driver (PoC #45);
#:     the productive yield is small enough that defaulting to the
#:     narrower inventory is a better cost/benefit trade.
#:
#: Callers can opt back into the wider inventory via the
#: ``lookahead_op1_inventory`` kwarg when reproducing PoC #45.
DEFAULT_LOOKAHEAD_OP1_INVENTORY: tuple[str, ...] = ("smooth_node",)
DEFAULT_LOOKAHEAD_OP2_INVENTORY: tuple[str, ...] = ("smooth_node",)


def _affected_nodes_from_info(info: dict[str, Any]) -> set[int]:
    """Node IDs whose 1-ring may have been modified by an op. The IDs
    reference the *resulting* mesh; new-node IDs from splits are
    included. Used by the 2-step lookahead union-penalty gate.
    """
    op = info["operator"]
    if op == "smooth_node":
        return {int(info["vertex"])}
    if op == "edge_swap":
        return {int(info["edge"][0]), int(info["edge"][1])}
    if op == "edge_split_interior":
        return {
            int(info["edge"][0]),
            int(info["edge"][1]),
            int(info["new_node"]),
        }
    if op == "edge_split_boundary":
        return {
            int(info["edge"][0]),
            int(info["edge"][1]),
            int(info["new_node"]),
        }
    if op == "vertex_remove":
        return {int(info["vertex"])}
    raise ValueError(f"unknown op: {op!r}")


def _affected_elements_in_mesh(
    mesh: Fort14Mesh, info: dict[str, Any],
) -> list[int]:
    """Element IDs in the *result* mesh that overlap an op's
    modified region. Used to enumerate candidate ``e2`` for Pass C
    op2 search.
    """
    op = info["operator"]
    if op == "smooth_node":
        return [int(x) for x in info["affected_elements"]]
    if op == "edge_swap":
        return [int(x) for x in info["elements_modified"]]
    if op == "edge_split_interior":
        ne = int(mesh.n_elements)
        return list(range(ne - 4, ne))
    if op == "edge_split_boundary":
        ne = int(mesh.n_elements)
        return list(range(ne - 2, ne))
    if op == "vertex_remove":
        ne = int(mesh.n_elements)
        n_new = int(info["n_new_elements"])
        return list(range(ne - n_new, ne))
    raise ValueError(f"unknown op: {op!r}")


def _union_penalty(
    mesh: Fort14Mesh, nodes: set[int], *,
    alpha_target: float, min_angle_target: float,
    max_angle_target: float = DEFAULT_MAX_ANGLE_TARGET,
) -> float:
    """Sum penalty over every element in ``mesh`` that touches any
    node in ``nodes``. IDs outside ``mesh.n_nodes`` are silently
    skipped (a new-node ID from a future mesh has no incident
    elements in the initial mesh, contributing zero).
    """
    if not nodes:
        return 0.0
    n2e = _node_to_elements(mesh.elements, mesh.n_nodes)
    affected_eids: set[int] = set()
    for v in nodes:
        if 0 <= v < mesh.n_nodes:
            for e in n2e.get(int(v), ()):
                affected_eids.add(int(e))
    if not affected_eids:
        return 0.0
    eids = np.fromiter(affected_eids, dtype=np.int64)
    block = mesh.elements[eids]
    a, m, M = _per_element_quality(mesh.nodes, block)
    return float(_penalty(
        a, m, M,
        alpha_target=alpha_target, min_angle_target=min_angle_target,
        max_angle_target=max_angle_target,
    ).sum())


def _ctx_for_lookahead(mesh: Fort14Mesh) -> dict[str, Any]:
    """Bundle the aux dicts ``_iter_op_candidates`` and the operators
    need. One call per ``_lookahead_round`` accept (mesh changed)."""
    n2e = _node_to_elements(mesh.elements, mesh.n_nodes)
    eu = _edge_use_counts(mesh.elements)
    bnd_node = _boundary_node_mask(mesh)
    bnd_edges = {k for k, v in eu.items() if len(v) == 1}
    _bp, _bn, e2s = _boundary_topology(mesh)
    return {
        "n2e": n2e,
        "eu": eu,
        "bnd_node": bnd_node,
        "bnd_edges": bnd_edges,
        "e2s": e2s,
    }


#: Acceptance gates available to :func:`_try_lookahead_pair`. The
#: PoC #45 negative result motivates the stricter v4.1 default
#: ``"target_exits_fail"``: each accept guarantees the target
#: element exits fail status (alpha ≥ alpha_target ∧ min_angle ≥
#: min_angle_target) on the resulting mesh, or has been removed by
#: the operator chain (vertex_remove of one of its vertices).
#: ``"union_penalty"`` reproduces the v4 behaviour (accept iff
#: penalty over op1 ∪ op2 affected nodes strictly drops); kept for
#: reproducibility of PoC #45.
LOOKAHEAD_GATES: tuple[str, ...] = ("target_exits_fail", "union_penalty")
DEFAULT_LOOKAHEAD_GATE: str = "target_exits_fail"


def _target_exits_fail(
    m_after: Fort14Mesh, target_vset: frozenset[int], *,
    alpha_target: float, min_angle_target: float,
    max_angle_target: float = DEFAULT_MAX_ANGLE_TARGET,
) -> bool:
    """v4.1 acceptance gate (strict). Returns True iff the triangle
    whose vertex set equals ``target_vset`` is **present in
    ``m_after``** and passes the per-element quality gate
    (``alpha >= alpha_target ∧ min_angle >= min_angle_target``).

    The target is located by vertex set rather than by element ID
    because topology ops shift element IDs.

    **An absent vertex set returns False, not True.** An earlier
    "fixed by elimination" interpretation produced PoC #46's
    catastrophic regression: under the (smooth_node, vertex_remove)
    inventory, vertex_remove deletes E by construction, so the
    elimination branch made every valid vertex_remove auto-accept
    and ~19 600 interior vertices were stripped from the Tokyo-Bay
    mesh, dropping ``n_elements`` 47 426 → 7 050 with
    ``alpha_p05`` 0.88 → 0.14 and ``frac<20°`` 0.13 % → 39.5 %.
    The strict spec demands that ``alpha(E)`` and ``min_angle(E)``
    actually exist and pass, which is only meaningful when E
    itself survives the operator chain.
    """
    target = set(target_vset)
    if len(target) != 3:
        return False
    v0 = next(iter(target))
    mask = (m_after.elements == v0).any(axis=1)
    for eid in np.where(mask)[0]:
        if set(int(x) for x in m_after.elements[eid]) == target:
            block = m_after.elements[[int(eid)]]
            a, m, M = _per_element_quality(m_after.nodes, block)
            return bool(
                float(a[0]) >= alpha_target
                and float(m[0]) >= min_angle_target
            )
    return False


def _iter_op_candidates(
    mesh: Fort14Mesh, eid: int, op_name: str, *,
    force: bool, ctx: dict[str, Any],
    alpha_target: float, min_angle_target: float,
    max_angle_target: float = DEFAULT_MAX_ANGLE_TARGET,
    coastline_projector: CoastlineProjector | None = None,
):
    """Yield ``(new_mesh, info)`` for each local variant of
    ``op_name`` applied to ``mesh.elements[eid]``. Skips variants the
    operator rejects (validity or, when ``force=False``, the penalty
    gate).
    """
    if op_name == "smooth_node":
        for v in mesh.elements[eid]:
            ring = ctx["n2e"].get(int(v))
            if ring is None:
                continue
            out = _apply_smooth_node(
                mesh, int(v), ring,
                alpha_target=alpha_target,
                min_angle_target=min_angle_target,

                max_angle_target=max_angle_target,
                boundary_node_mask=ctx["bnd_node"],
                force=force,
            )
            if out is not None:
                yield out
    elif op_name == "edge_swap":
        for k in range(3):
            out = _apply_edge_swap(
                mesh, int(eid), k,
                alpha_target=alpha_target,
                min_angle_target=min_angle_target,

                max_angle_target=max_angle_target,
                edge_uses=ctx["eu"],
                boundary_edge_keys=ctx["bnd_edges"],
                force=force,
            )
            if out is not None:
                yield out
    elif op_name == "edge_split_interior":
        for k in range(3):
            out = _apply_edge_split_interior(
                mesh, int(eid), k,
                alpha_target=alpha_target,
                min_angle_target=min_angle_target,

                max_angle_target=max_angle_target,
                edge_uses=ctx["eu"],
                boundary_edge_keys=ctx["bnd_edges"],
                force=force,
            )
            if out is not None:
                yield out
    elif op_name == "edge_split_boundary":
        for k in range(3):
            out = _apply_edge_split_boundary(
                mesh, int(eid), k,
                alpha_target=alpha_target,
                min_angle_target=min_angle_target,

                max_angle_target=max_angle_target,
                edge_uses=ctx["eu"],
                edge_to_segment=ctx["e2s"],
                coastline_projector=coastline_projector,
                force=force,
            )
            if out is not None:
                yield out
    elif op_name == "vertex_remove":
        for v in mesh.elements[eid]:
            ring = ctx["n2e"].get(int(v))
            if ring is None:
                continue
            out = _apply_vertex_remove(
                mesh, int(v), ring,
                alpha_target=alpha_target,
                min_angle_target=min_angle_target,

                max_angle_target=max_angle_target,
                boundary_node_mask=ctx["bnd_node"],
                force=force,
            )
            if out is not None:
                yield out
    else:
        raise ValueError(f"unknown op: {op_name!r}")


def _try_lookahead_pair(
    cur: Fort14Mesh, eid: int, ctx: dict[str, Any], *,
    alpha_target: float, min_angle_target: float,
    max_angle_target: float = DEFAULT_MAX_ANGLE_TARGET,
    op1_inventory: tuple[str, ...],
    op2_inventory: tuple[str, ...],
    coastline_projector: CoastlineProjector | None,
    gate: str = DEFAULT_LOOKAHEAD_GATE,
) -> tuple[Fort14Mesh, str] | None:
    """Search for an accepting ``(op1, op2)`` chain on fail element
    ``eid``. op1 is applied with ``force=True`` (validity only);
    op2 is searched on the elements overlapping op1's affected
    region that still fail the per-element gate in ``m1``.

    Acceptance depends on ``gate``:

    * ``"target_exits_fail"`` (v4.1 default) — accept iff the
      target element ``E`` (located in ``m_after`` by vertex set)
      exits fail status, or was removed by the operator chain.
      Also tries **op1-only** first: if op1 alone makes E exit
      fail, accept without searching for op2.
    * ``"union_penalty"`` (v4 reproduction) — accept iff the
      penalty summed over elements touching
      ``op1.affected ∪ op2.affected`` nodes strictly drops between
      ``cur`` and ``m2``. No op1-only path.

    Returns ``(m_after, "op1+op2")`` (or ``"op1+none"`` for an
    op1-only accept under the ``target_exits_fail`` gate) on the
    first accept; ``None`` if every candidate is rejected.
    """
    if gate not in LOOKAHEAD_GATES:
        raise ValueError(f"unknown lookahead gate: {gate!r}")

    target_vset: frozenset[int] = frozenset(
        int(v) for v in cur.elements[eid]
    )

    for op1_name in op1_inventory:
        for m1, info1 in _iter_op_candidates(
            cur, eid, op1_name, force=True, ctx=ctx,
            alpha_target=alpha_target, min_angle_target=min_angle_target,
            coastline_projector=coastline_projector,
        ):
            affected1_nodes = _affected_nodes_from_info(info1)

            # v4.1 op1-only fast path. Only meaningful under the
            # target_exits_fail gate — under union_penalty, op1
            # alone is already covered by Pass B's strict gate.
            if gate == "target_exits_fail" and _target_exits_fail(
                m1, target_vset,
                alpha_target=alpha_target,
                min_angle_target=min_angle_target,

                max_angle_target=max_angle_target,
            ):
                return m1, f"{op1_name}+none"

            affected1_eids = _affected_elements_in_mesh(m1, info1)
            if not affected1_eids:
                continue
            block = m1.elements[affected1_eids]
            a_blk, m_blk, M_blk = _per_element_quality(m1.nodes, block)
            fail_blk = _is_fail(
                a_blk, m_blk, M_blk,
                alpha_target=alpha_target, min_angle_target=min_angle_target,
                max_angle_target=max_angle_target,
            )
            cand_e2 = [
                int(e) for e, f in zip(affected1_eids, fail_blk) if f
            ]
            if not cand_e2:
                continue
            ctx_m1 = _ctx_for_lookahead(m1)
            for e2 in cand_e2:
                for op2_name in op2_inventory:
                    for m2, info2 in _iter_op_candidates(
                        m1, e2, op2_name, force=True, ctx=ctx_m1,
                        alpha_target=alpha_target,
                        min_angle_target=min_angle_target,

                        max_angle_target=max_angle_target,
                        coastline_projector=coastline_projector,
                    ):
                        if gate == "target_exits_fail":
                            if _target_exits_fail(
                                m2, target_vset,
                                alpha_target=alpha_target,
                                min_angle_target=min_angle_target,

                                max_angle_target=max_angle_target,
                            ):
                                return m2, f"{op1_name}+{op2_name}"
                        else:  # union_penalty
                            affected2_nodes = _affected_nodes_from_info(info2)
                            union = affected1_nodes | affected2_nodes
                            pen_before = _union_penalty(
                                cur, union,
                                alpha_target=alpha_target,
                                min_angle_target=min_angle_target,

                                max_angle_target=max_angle_target,
                            )
                            pen_after = _union_penalty(
                                m2, union,
                                alpha_target=alpha_target,
                                min_angle_target=min_angle_target,

                                max_angle_target=max_angle_target,
                            )
                            if pen_after + 1e-12 < pen_before:
                                return m2, f"{op1_name}+{op2_name}"
    return None


def _lookahead_round(
    cur: Fort14Mesh, *,
    alpha_target: float, min_angle_target: float,
    max_angle_target: float = DEFAULT_MAX_ANGLE_TARGET,
    op1_inventory: tuple[str, ...],
    op2_inventory: tuple[str, ...],
    max_lookahead_accepts: int,
    coastline_projector: CoastlineProjector | None = None,
    gate: str = DEFAULT_LOOKAHEAD_GATE,
) -> tuple[Fort14Mesh, dict[str, int], int]:
    """Pass C: 2-step lookahead. Pops fail elements by descending
    penalty and applies the first accepting ``(op1, op2)`` chain
    found by :func:`_try_lookahead_pair`. Each accept rebuilds the
    aux dicts on the new mesh. ``gate`` selects the acceptance rule
    (see :func:`_try_lookahead_pair`). Returns
    ``(updated_mesh, accepts_per_pair, n_abandoned)`` where
    ``accepts_per_pair`` is keyed by ``"op1+op2"`` (or
    ``"op1+none"`` for op1-only accepts under the v4.1 gate).
    """
    accepts_per_pair: dict[str, int] = defaultdict(int)
    abandoned: set = set()
    accepts_total = 0
    while accepts_total < max_lookahead_accepts:
        a, m, M = _per_element_quality(cur.nodes, cur.elements)
        fail = _is_fail(
            a, m, M,
            alpha_target=alpha_target, min_angle_target=min_angle_target,
            max_angle_target=max_angle_target,
        )
        if not fail.any():
            break
        pen = _penalty(
            a, m, M,
            alpha_target=alpha_target, min_angle_target=min_angle_target,
            max_angle_target=max_angle_target,
        )
        heap: list[tuple[float, int]] = []
        for eid in np.where(fail)[0]:
            heapq.heappush(heap, (-float(pen[eid]), int(eid)))
        ctx = _ctx_for_lookahead(cur)

        progress = False
        while heap:
            _neg, eid = heapq.heappop(heap)
            sig = (
                int(cur.n_elements),
                frozenset(int(x) for x in cur.elements[eid]),
            )
            if sig in abandoned:
                continue
            applied = _try_lookahead_pair(
                cur, eid, ctx,
                alpha_target=alpha_target,
                min_angle_target=min_angle_target,

                max_angle_target=max_angle_target,
                op1_inventory=op1_inventory,
                op2_inventory=op2_inventory,
                coastline_projector=coastline_projector,
                gate=gate,
            )
            if applied is None:
                abandoned.add(sig)
                continue
            cur, pair_label = applied
            accepts_per_pair[pair_label] += 1
            accepts_total += 1
            progress = True
            break  # rebuild aux dicts

        if not progress:
            break

    return cur, dict(accepts_per_pair), len(abandoned)


# ---------------------------------------------------------------------------
# Pass D — cluster-scale patch re-CDT
#
# Design rationale: ``docs/patch_re_cdt_design.md``. PoC #47a-prep
# established that 51 % of the Tokyo-Bay v3 residual sits in
# face-face-adjacency clusters of size >= 3, which no 1-ring local
# edit can repair. Pass D extracts the rim polygon of a fail
# cluster, drops its interior, and re-triangulates the polygon as a
# pure Delaunay patch. Strict acceptance: every new patch element
# must pass the per-element gate, and the rim 1-ring outside the
# cluster must not gain new fails (no 2-ring drift).
# ---------------------------------------------------------------------------


#: Pass D default min cluster size. Below 3, the Pass C lookahead is
#: the right tool (PoC #44b measured a 4.6 % theoretical ceiling
#: there). Setting min=3 confines Pass D to its addressable market.
DEFAULT_PATCH_MIN_CLUSTER_SIZE: int = 3

#: Pass D default max cluster size. The PoC #47a-prep histogram
#: showed the v3 residual's largest cluster is 20; 100 is a safety
#: margin for other meshes. Larger clusters indicate a different
#: failure mode (likely Phase E under-resolved channels) and Pass D
#: should defer rather than attempt a costly re-CDT.
DEFAULT_PATCH_MAX_CLUSTER_SIZE: int = 100


def _find_fail_clusters(
    mesh: Fort14Mesh, *,
    alpha_target: float, min_angle_target: float,
    max_angle_target: float = DEFAULT_MAX_ANGLE_TARGET,
    min_cluster_size: int = DEFAULT_PATCH_MIN_CLUSTER_SIZE,
    max_cluster_size: int = DEFAULT_PATCH_MAX_CLUSTER_SIZE,
) -> list[np.ndarray]:
    """Identify connected components of fail elements under
    face-face adjacency. Returns a list of (element-id arrays),
    sorted by size (largest first). Clusters with size outside
    ``[min_cluster_size, max_cluster_size]`` are filtered out.
    Lazy-imports its scipy and diagnostics dependencies to keep the
    Phase H module's startup cost bounded.
    """
    from scipy.sparse.csgraph import (  # noqa: PLC0415
        connected_components,
    )

    from fvcom_mesh_tools.diagnostics import (  # noqa: PLC0415
        face_face_adjacency,
    )

    a, m, M = _per_element_quality(mesh.nodes, mesh.elements)
    fail = _is_fail(
        a, m, M,
        alpha_target=alpha_target, min_angle_target=min_angle_target,
        max_angle_target=max_angle_target,
    )
    if not fail.any():
        return []
    fail_eids = np.where(fail)[0]

    adj = face_face_adjacency(mesh.elements)
    fail_sub = adj[fail][:, fail]
    n_components, labels = connected_components(
        fail_sub, directed=False, return_labels=True,
    )

    clusters: list[np.ndarray] = []
    for cid in range(n_components):
        members = fail_eids[labels == cid]
        size = int(members.size)
        if min_cluster_size <= size <= max_cluster_size:
            clusters.append(members.astype(np.int64))
    clusters.sort(key=len, reverse=True)
    return clusters


def _external_rim_fail_count(
    mesh: Fort14Mesh, rim_node_ids: np.ndarray,
    exclude_eids: set[int], *,
    alpha_target: float, min_angle_target: float,
    max_angle_target: float = DEFAULT_MAX_ANGLE_TARGET,
) -> int:
    """Number of fail elements in ``mesh`` that touch any rim node
    and are NOT in ``exclude_eids``. Used by Pass D's second gate
    (no fail regression in the rim 1-ring outside the cluster).
    """
    n2e = _node_to_elements(mesh.elements, mesh.n_nodes)
    affected: set[int] = set()
    for n in rim_node_ids:
        nid = int(n)
        if 0 <= nid < mesh.n_nodes:
            for e in n2e.get(nid, ()):
                eid = int(e)
                if eid not in exclude_eids:
                    affected.add(eid)
    if not affected:
        return 0
    eids = np.fromiter(affected, dtype=np.int64)
    block = mesh.elements[eids]
    a, m, M = _per_element_quality(mesh.nodes, block)
    fail = _is_fail(
        a, m, M,
        alpha_target=alpha_target, min_angle_target=min_angle_target,
        max_angle_target=max_angle_target,
    )
    return int(fail.sum())


def _attempt_patch_recdt(
    mesh: Fort14Mesh, cluster_eids: np.ndarray, *,
    alpha_target: float, min_angle_target: float,
    max_angle_target: float = DEFAULT_MAX_ANGLE_TARGET,
    reject_boundary_clusters: bool = True,
) -> tuple[Fort14Mesh | None, dict[str, Any] | None, str]:
    """Diagnostic core of Pass D's patch re-CDT operator. Returns
    ``(new_mesh, info, reason)`` where:

    * on accept: ``(new_mesh, info, "ok")``.
    * on reject: ``(None, None, reason)`` where ``reason`` is one of
      ``empty_cluster``, ``rim_walk_failed``, ``rim_on_boundary``,
      ``retriangulate_failed``, ``gate1_alpha``, ``gate1_min_angle``,
      ``gate1_flipped``, or ``gate2_rim_regression``.

    The public ``_apply_patch_recdt`` wraps this to provide the
    ``(mesh, info) | None`` shape expected by ``_patch_recdt_round``.
    PoC #47a uses the diagnostic form directly to histogram reject
    causes.
    """
    cluster_set = {int(e) for e in cluster_eids}
    if len(cluster_set) < 1:
        return None, None, "empty_cluster"
    cluster_eids_arr = np.array(sorted(cluster_set), dtype=np.int64)

    rim_node_ids = _patch_rim_polygon(mesh.elements, cluster_eids_arr)
    if rim_node_ids is None:
        return None, None, "rim_walk_failed"

    if reject_boundary_clusters:
        bnd = _boundary_node_mask(mesh)
        if bool(bnd[rim_node_ids].any()):
            return None, None, "rim_on_boundary"

    n_rim = int(rim_node_ids.size)
    rim_xy = mesh.nodes[rim_node_ids]

    lat_centre = float(rim_xy[:, 1].mean())
    deg_per_m_lat = 1.0 / (EARTH_R_M * np.pi / 180.0)
    deg_per_m_lon = deg_per_m_lat / max(
        np.cos(np.deg2rad(lat_centre)), 1e-6,
    )
    rim_xy_m = np.column_stack([
        (rim_xy[:, 0] - rim_xy[0, 0]) / deg_per_m_lon,
        (rim_xy[:, 1] - rim_xy[0, 1]) / deg_per_m_lat,
    ])

    sx, sy = rim_xy_m[:, 0], rim_xy_m[:, 1]
    if 0.5 * float(
        np.sum(sx * np.roll(sy, -1) - np.roll(sx, -1) * sy),
    ) < 0:
        rim_node_ids = rim_node_ids[::-1].copy()
        rim_xy = mesh.nodes[rim_node_ids]
        rim_xy_m = np.column_stack([
            (rim_xy[:, 0] - rim_xy[0, 0]) / deg_per_m_lon,
            (rim_xy[:, 1] - rim_xy[0, 1]) / deg_per_m_lat,
        ])

    triangles, _reason = _retriangulate_patch(
        rim_xy_m, np.empty((0, 2), dtype=float), n_rim,
    )
    if triangles is None:
        return None, None, "retriangulate_failed"

    new_block = rim_node_ids[triangles].astype(mesh.elements.dtype)

    # Gate 1: every new patch element passes the per-element gate.
    a_new, m_new, M_new = _per_element_quality(mesh.nodes, new_block)
    if bool((a_new < alpha_target).any()):
        return None, None, "gate1_alpha"
    if bool((m_new < min_angle_target).any()):
        return None, None, "gate1_min_angle"

    # Signed-area sanity check on the new block (redundant with
    # _retriangulate_patch's own check, but cheap and defensive).
    p0 = mesh.nodes[new_block[:, 0]]
    p1 = mesh.nodes[new_block[:, 1]]
    p2 = mesh.nodes[new_block[:, 2]]
    cross = (
        (p1[:, 0] - p0[:, 0]) * (p2[:, 1] - p0[:, 1])
        - (p1[:, 1] - p0[:, 1]) * (p2[:, 0] - p0[:, 0])
    )
    if bool((cross <= 0).any()):
        return None, None, "gate1_flipped"

    # Build candidate mesh.
    keep_mask = np.ones(mesh.n_elements, dtype=bool)
    keep_mask[cluster_eids_arr] = False
    new_elements = np.vstack([mesh.elements[keep_mask], new_block])
    new_mesh = Fort14Mesh(
        title=mesh.title,
        nodes=mesh.nodes.copy(),
        depths=mesh.depths.copy(),
        elements=new_elements,
        open_boundaries=[np.asarray(s).copy() for s in mesh.open_boundaries],
        land_boundaries=[(int(ib), np.asarray(s).copy())
                         for ib, s in mesh.land_boundaries],
    )

    # Gate 2: no fail regression in the rim 1-ring outside cluster.
    fail_before = _external_rim_fail_count(
        mesh, rim_node_ids, cluster_set,
        alpha_target=alpha_target, min_angle_target=min_angle_target,
    )
    n_patch_new = int(new_block.shape[0])
    patch_eid_start = new_mesh.n_elements - n_patch_new
    new_patch_set = set(range(patch_eid_start, new_mesh.n_elements))
    fail_after = _external_rim_fail_count(
        new_mesh, rim_node_ids, new_patch_set,
        alpha_target=alpha_target, min_angle_target=min_angle_target,
    )
    if fail_after > fail_before:
        return None, None, "gate2_rim_regression"

    cluster_nodes = np.unique(mesh.elements[cluster_eids_arr].ravel())
    rim_set = {int(n) for n in rim_node_ids}
    interior_node_ids = np.array(
        [int(n) for n in cluster_nodes if int(n) not in rim_set],
        dtype=np.int64,
    )

    return new_mesh, {
        "operator": "patch_recdt",
        "cluster_size": int(cluster_eids_arr.size),
        "rim_size": int(n_rim),
        "n_new_elements": int(n_patch_new),
        "n_replaced_elements": int(cluster_eids_arr.size),
        "n_interior_orphaned": int(interior_node_ids.size),
        "external_rim_fail_before": int(fail_before),
        "external_rim_fail_after": int(fail_after),
    }, "ok"


def _apply_patch_recdt(
    mesh: Fort14Mesh, cluster_eids: np.ndarray, *,
    alpha_target: float, min_angle_target: float,
    max_angle_target: float = DEFAULT_MAX_ANGLE_TARGET,
    reject_boundary_clusters: bool = True,
) -> tuple[Fort14Mesh, dict[str, Any]] | None:
    """Pass D op: replace a fail cluster's elements with a Delaunay
    re-triangulation of the cluster's rim polygon. Interior cluster
    nodes are orphaned (kept in the node array, no longer
    referenced by any element). Accept iff:

    1. Every new patch element passes the per-element gate
       (``alpha >= alpha_target ∧ min_angle >= min_angle_target``).
    2. The rim 1-ring fail count outside the cluster does not
       increase between ``mesh`` and the post-op mesh (no 2-ring
       drift).

    v1 rejects clusters whose rim touches a boundary segment to
    avoid open / land segment book-keeping. Backwards-compatible
    wrapper over :func:`_attempt_patch_recdt`; returns ``None`` on
    every rejection. See the diagnostic helper for the reason
    histogram.
    """
    new_mesh, info, _reason = _attempt_patch_recdt(
        mesh, cluster_eids,
        alpha_target=alpha_target,
        min_angle_target=min_angle_target,

        max_angle_target=max_angle_target,
        reject_boundary_clusters=reject_boundary_clusters,
    )
    if new_mesh is None or info is None:
        return None
    return new_mesh, info


def _size_bucket(size: int) -> str:
    """Histogram bucket label for accept counters."""
    if size <= 9:
        return f"size_{size}"
    decade = (size // 10) * 10
    return f"size_{decade}_{decade + 9}"


def _patch_recdt_round(
    cur: Fort14Mesh, *,
    alpha_target: float, min_angle_target: float,
    max_angle_target: float = DEFAULT_MAX_ANGLE_TARGET,
    min_cluster_size: int,
    max_cluster_size: int,
    max_patches: int,
    reject_boundary_clusters: bool = True,
) -> tuple[Fort14Mesh, dict[str, int], int]:
    """Pass D round: process fail clusters in descending size order,
    apply ``_apply_patch_recdt`` to each, accept under the strict
    gate. Each accept rebuilds the cluster list. Returns
    ``(updated_mesh, accepts_per_size_bucket, n_rejected_by_gate)``.
    """
    accepts: dict[str, int] = defaultdict(int)
    abandoned: set = set()
    n_rejected = 0
    n_accepted = 0

    while n_accepted < max_patches:
        clusters = _find_fail_clusters(
            cur,
            alpha_target=alpha_target,
            min_angle_target=min_angle_target,

            max_angle_target=max_angle_target,
            min_cluster_size=min_cluster_size,
            max_cluster_size=max_cluster_size,
        )
        if not clusters:
            break

        progress = False
        for cluster in clusters:
            sig = frozenset(int(e) for e in cluster)
            if sig in abandoned:
                continue
            out = _apply_patch_recdt(
                cur, cluster,
                alpha_target=alpha_target,
                min_angle_target=min_angle_target,

                max_angle_target=max_angle_target,
                reject_boundary_clusters=reject_boundary_clusters,
            )
            if out is None:
                abandoned.add(sig)
                n_rejected += 1
                continue
            cur, info = out
            accepts[_size_bucket(int(info["cluster_size"]))] += 1
            n_accepted += 1
            progress = True
            break  # restart cluster enumeration on the new mesh

        if not progress:
            break

    return cur, dict(accepts), n_rejected


# ---------------------------------------------------------------------------
# Pass E — gradation refinement (FVCOM manual criterion C4,
# adjacent-element area change <= area_ratio_target)
# ---------------------------------------------------------------------------


def _apply_pass_e_swap(
    mesh: Fort14Mesh, fail_edge_uv: tuple[int, int],
    elem_pair: tuple[int, int], *,
    alpha_target: float, min_angle_target: float,
    max_angle_target: float = DEFAULT_MAX_ANGLE_TARGET,
    area_ratio_target: float = DEFAULT_AREA_RATIO_TARGET,
    max_valence: int = DEFAULT_MAX_VALENCE,
    edge_uses: dict[tuple[int, int], list[int]],
    boundary_edge_keys: set,
    c4_fail_keys: set[tuple[int, int]] | None = None,
    valence_before: np.ndarray | None = None,
) -> tuple[Fort14Mesh, dict[str, Any]] | None:
    """Pass E operator: edge_swap of the C4 fail edge.

    Topologically the quad ``L ∪ S`` is re-triangulated with the
    other diagonal ``CD`` (where ``C`` is ``L``'s apex not on the
    fail edge and ``D`` is ``S``'s apex). No vertices are added,
    so no global element-count growth and no midpoint to smooth.
    The endpoint valences (``a``, ``b``) drop by 1 each and the
    opposite valences (``C``, ``D``) rise by 1 each; the net change
    in mesh valence is zero, but ``C5`` could still regress if
    ``C`` or ``D`` was already at the cap.

    Accept gates (all required):

    1. The new diagonal ``CD``'s area_change must be at or below
       ``area_ratio_target``. If the swap merely relocates the
       violation, reject — we want a strict ``-1`` on the fail count.
    2. No edge in the local 5-edge block (``CD`` + the four edges
       ``a-C``, ``a-D``, ``b-C``, ``b-D``) may newly transition from
       non-fail to fail. Edges already in ``c4_fail_keys`` are
       counted as pre-existing and not "new".
    3. No element in ``block_after`` (= the two new triangles
       ``aCD`` and ``bDC``) may fail C1 (``min_ang >=
       min_angle_target``) or C2 (``max_ang <= max_angle_target``).
    4. ``C5`` valence cap: both ``C`` and ``D`` must have
       ``valence_before < max_valence`` (since each gains +1 from
       the swap).

    Returns ``(new_mesh, info)`` on accept, else ``None``.
    """
    fail_u = int(fail_edge_uv[0])
    fail_v = int(fail_edge_uv[1])
    fail_key = (min(fail_u, fail_v), max(fail_u, fail_v))
    if fail_key in boundary_edge_keys:
        return None  # edge_swap of a boundary edge is undefined

    e_i, e_j = int(elem_pair[0]), int(elem_pair[1])
    elem_id = e_i
    edge_local = -1
    for el in range(3):
        u = int(mesh.elements[elem_id, el])
        v = int(mesh.elements[elem_id, (el + 1) % 3])
        if (min(u, v), max(u, v)) == fail_key:
            edge_local = el
            break
    if edge_local < 0:
        # Fall back: the fail key might be ordered the other way
        # round on e_j — but our routine just needs SOMEONE that
        # contains it.
        elem_id = e_j
        for el in range(3):
            u = int(mesh.elements[elem_id, el])
            v = int(mesh.elements[elem_id, (el + 1) % 3])
            if (min(u, v), max(u, v)) == fail_key:
                edge_local = el
                break
        if edge_local < 0:
            return None
    buddy_id = e_j if elem_id == e_i else e_i

    third = int(mesh.elements[elem_id, (edge_local + 2) % 3])
    buddy_set = {int(x) for x in mesh.elements[buddy_id]}
    fourth_set = buddy_set - {fail_u, fail_v}
    if len(fourth_set) != 1:
        return None
    fourth = int(next(iter(fourth_set)))

    # Gate (4): C5 prefilter on the gainers.
    if valence_before is not None:
        if valence_before[third] >= max_valence:
            return None
        if valence_before[fourth] >= max_valence:
            return None

    # Build the swap directly. The existing ``_apply_edge_swap``
    # always emits ``[a, c, d]`` / ``[b, d, c]`` for the new pair,
    # which is CCW only when the quad walks as ``a -> c -> b -> d``
    # CCW. For "kite" geometries — typical of C4 fails, with ``c``
    # and ``d`` on opposite sides of the fail edge — the natural
    # CCW walk is the OPPOSITE rotation, ``a -> d -> b -> c``, so
    # the operator's emission would be CW and the signed-area check
    # would reject. Pass E swap tries both orientations and picks
    # whichever produces two CCW triangles.
    dtype = mesh.elements.dtype
    cand_a = np.array(
        [[fail_u, third, fourth], [fail_v, fourth, third]], dtype=dtype,
    )
    cand_b = np.array(
        [[fail_u, fourth, third], [fail_v, third, fourth]], dtype=dtype,
    )

    def _all_ccw(block):
        p0 = mesh.nodes[block[:, 0]]
        p1 = mesh.nodes[block[:, 1]]
        p2 = mesh.nodes[block[:, 2]]
        cross = (
            (p1[:, 0] - p0[:, 0]) * (p2[:, 1] - p0[:, 1])
            - (p1[:, 1] - p0[:, 1]) * (p2[:, 0] - p0[:, 0])
        )
        return bool((cross > 0).all())

    if _all_ccw(cand_a):
        block_after = cand_a
    elif _all_ccw(cand_b):
        block_after = cand_b
    else:
        return None  # concave quad — swap is geometrically invalid

    new_elements = mesh.elements.copy()
    new_elements[elem_id] = block_after[0]
    new_elements[buddy_id] = block_after[1]
    new_mesh = Fort14Mesh(
        title=mesh.title,
        nodes=mesh.nodes.copy(),
        depths=mesh.depths.copy(),
        elements=new_elements,
        open_boundaries=[
            np.asarray(s).copy() for s in mesh.open_boundaries
        ],
        land_boundaries=[
            (int(ib), np.asarray(s).copy())
            for ib, s in mesh.land_boundaries
        ],
    )

    # Gate (3): no new C1 or C2 fail in the new block.
    new_block = new_mesh.elements[[elem_id, buddy_id]]
    _a_a, m_a, M_a = _per_element_quality(new_mesh.nodes, new_block)
    if (m_a < min_angle_target).any():
        return None
    if (M_a > max_angle_target).any():
        return None

    # Gates (1) + (2): C4 status of CD and of the four surviving edges.
    new_abs_areas = np.abs(_signed_areas(new_mesh.nodes, new_mesh.elements))
    new_edge_uses = _edge_use_counts(new_mesh.elements)

    def _ac_of(k):
        bud = new_edge_uses.get(k, [])
        if len(bud) != 2:
            return -1.0  # boundary or vanished
        a1 = float(new_abs_areas[bud[0]])
        a2 = float(new_abs_areas[bud[1]])
        return abs(a1 - a2) / max(a1, a2, 1e-30)

    cd_key = (min(third, fourth), max(third, fourth))
    cd_ac = _ac_of(cd_key)
    if cd_ac < 0.0 or cd_ac > area_ratio_target:
        return None  # CD fails C4 (or vanished); swap doesn't help

    surviving_keys = [
        (min(fail_u, third), max(fail_u, third)),
        (min(fail_v, third), max(fail_v, third)),
        (min(fail_u, fourth), max(fail_u, fourth)),
        (min(fail_v, fourth), max(fail_v, fourth)),
    ]
    for k in surviving_keys:
        new_ac = _ac_of(k)
        if new_ac > area_ratio_target:
            was_fail = c4_fail_keys is not None and k in c4_fail_keys
            if not was_fail:
                return None  # new fail introduced on a surviving edge

    info_out = {
        "operator": "pass_e_swap",
        "fail_edge": [fail_key[0], fail_key[1]],
        "new_diagonal": [int(cd_key[0]), int(cd_key[1])],
        "elements_modified": [int(elem_id), int(buddy_id)],
        "area_change_before": None,
        "area_change_after": float(cd_ac),
    }
    return new_mesh, info_out


def _apply_pass_e_split(
    mesh: Fort14Mesh, fail_edge_uv: tuple[int, int],
    elem_pair: tuple[int, int], *,
    alpha_target: float, min_angle_target: float,
    max_angle_target: float = DEFAULT_MAX_ANGLE_TARGET,
    area_ratio_target: float = DEFAULT_AREA_RATIO_TARGET,
    max_valence: int = DEFAULT_MAX_VALENCE,
    edge_uses: dict[tuple[int, int], list[int]],
    boundary_edge_keys: set,
    edge_to_segment: dict[tuple[int, int], tuple[str, int, int, int]],
    c4_fail_keys: set[tuple[int, int]] | None = None,
    valence_before: np.ndarray | None = None,
    abs_areas_before: np.ndarray | None = None,
    coastline_projector: CoastlineProjector | None = None,
) -> tuple[Fort14Mesh, dict[str, Any]] | None:
    """Pass E operator: split the longest non-shared edge of the
    larger triangle adjacent to a C4-fail edge.

    Strategy: let ``L`` be the larger of the two triangles incident on
    the fail edge and ``S`` the smaller. Splitting an edge of ``L``
    (other than the shared fail edge) at its midpoint halves ``L``'s
    area, which reduces the area_change against ``S`` by roughly a
    factor of two while leaving ``S`` untouched. Splitting the shared
    edge instead halves *both* areas symmetrically and does not change
    the ratio — that is why the operator picks a non-shared edge.

    Candidate selection (cascade + indirect-regression avoidance):
    the three edges of ``L`` are filtered down to those that are
    (a) not the shared C4 fail edge itself, and (b) not themselves
    in the current C4 fail set (``c4_fail_keys``). Splitting an
    edge that is itself a C4 fail cascades the violation onto the
    two new edges (each inheriting the same area_change), so no
    net C4 reduction would result.

    A further geometric filter rejects candidates where halving
    ``L`` would introduce a *new* C4 fail on one of ``L``'s
    preserved non-shared edges (and, for interior splits, on the
    buddy triangle ``T``'s preserved non-shared edges). Concretely,
    after splitting an edge of ``L`` at its midpoint both halves
    have area ``L/2``, so any non-split edge of ``L`` previously
    matched against an external neighbour of area ``A_n`` now sits
    between an ``L/2`` triangle and ``A_n``; if ``L < A_n`` this
    asymmetry strictly worsens and may flip the edge into fail.
    The remaining candidates are sorted by length (longest first).

    Acceptance gate (all required):

    1. The original C4 fail edge's area_change must strictly drop.
       Without this we are doing topology churn for nothing.
    2. No element in the affected block may newly fail C1 or C2 —
       i.e., every post-split element must have
       ``min_ang >= min_angle_target`` AND
       ``max_ang <= max_angle_target``. The alpha quality metric is
       deliberately *not* part of this gate: it is a proxy for angle
       quality, and the FVCOM manual lists C1 and C2 directly, so
       Pass E gates against the real criteria. Alpha may degrade
       slightly because refining a coarse triangle into smaller ones
       changes its edge-length distribution.
    3. No node's valence may exceed ``max_valence`` (FVCOM manual
       criterion C5). ``edge_split_interior`` raises the valence of
       each "opposite" vertex (the third vertex of each incident
       triangle, not on the split edge) by 1, so a vertex at the
       limit before the split would regress C5 if split anyway.
       ``valence_before`` short-circuits this check without rebuilding
       the new mesh.

    Returns ``(new_mesh, info)`` on accept, else ``None``. ``info``
    matches the operator-info dict shape used by Pass B so the
    standard accept bookkeeping can chain through.
    """
    fail_u, fail_v = int(fail_edge_uv[0]), int(fail_edge_uv[1])
    fail_key = (min(fail_u, fail_v), max(fail_u, fail_v))

    a_pair = _signed_areas(mesh.nodes, mesh.elements[list(elem_pair)])
    if a_pair.size != 2:
        return None
    abs_a = np.abs(a_pair)
    if abs_a[0] >= abs_a[1]:
        larger = int(elem_pair[0])
        A_L_before, A_S_before = float(abs_a[0]), float(abs_a[1])
    else:
        larger = int(elem_pair[1])
        A_L_before, A_S_before = float(abs_a[1]), float(abs_a[0])

    # All three edges of the larger triangle, with their lengths.
    # Exclude the fail edge itself and any other C4-fail edge in L
    # (cascade avoidance — splitting a fail edge merely re-locates
    # the violation onto the two new edges).
    verts_L = mesh.elements[larger]
    candidates = []  # (edge_local, length, is_boundary)
    for el in range(3):
        u = int(verts_L[el])
        v = int(verts_L[(el + 1) % 3])
        key = (min(u, v), max(u, v))
        if key == fail_key:
            continue
        if c4_fail_keys is not None and key in c4_fail_keys:
            continue
        length = float(np.linalg.norm(mesh.nodes[u] - mesh.nodes[v]))
        candidates.append((el, length, key in boundary_edge_keys))
    if not candidates:
        return None
    candidates.sort(key=lambda x: x[1], reverse=True)

    ac_before = (
        (max(A_L_before, A_S_before) - min(A_L_before, A_S_before))
        / max(A_L_before, A_S_before, 1e-30)
    )

    for edge_local, _length, is_boundary in candidates:
        # Gate (3) prefilter — C5 valence. ``edge_split_interior``
        # raises the valence of the opposite vertex of L AND the
        # opposite vertex of the buddy triangle T (each by 1).
        # ``edge_split_boundary`` raises only L's opposite vertex.
        # If either is already at the cap, skip without paying the
        # cost of a full split + rebuild.
        c_self = int(verts_L[(edge_local + 2) % 3])
        if (
            valence_before is not None
            and valence_before[c_self] >= max_valence
        ):
            continue
        if not is_boundary:
            u = int(verts_L[edge_local])
            v = int(verts_L[(edge_local + 1) % 3])
            key = (min(u, v), max(u, v))
            incident = edge_uses.get(key, [])
            if len(incident) != 2:
                continue
            other_eid = (
                incident[0] if incident[1] == larger else incident[1]
            )
            other_verts = mesh.elements[other_eid]
            edge_set = {u, v}
            c_other = -1
            for ov in other_verts:
                ov_i = int(ov)
                if ov_i not in edge_set:
                    c_other = ov_i
                    break
            if c_other < 0:
                continue
            if (
                valence_before is not None
                and valence_before[c_other] >= max_valence
            ):
                continue

        # Indirect-regression filter — would halving L push any of L's
        # preserved non-shared edges (or, for interior splits, any of
        # T's preserved non-shared edges) over the C4 threshold?
        # The check is geometric: after a midpoint split both halves
        # have area_before / 2, so each preserved edge sits between an
        # ``area/2`` triangle and the unchanged external neighbour.
        if abs_areas_before is not None:
            split_u = int(verts_L[edge_local])
            split_v = int(verts_L[(edge_local + 1) % 3])
            split_key = (min(split_u, split_v), max(split_u, split_v))

            def _would_create_new_fail(
                eid: int, eid_area_after: float,
            ) -> bool:
                v_eid = mesh.elements[eid]
                for ell in range(3):
                    a = int(v_eid[ell])
                    b = int(v_eid[(ell + 1) % 3])
                    k = (min(a, b), max(a, b))
                    if k == fail_key:
                        continue  # gated separately by C4-drop check
                    if k == split_key:
                        continue  # being split, not "preserved"
                    if k in boundary_edge_keys:
                        continue
                    if (
                        c4_fail_keys is not None
                        and k in c4_fail_keys
                    ):
                        continue  # already a fail, not a *new* one
                    buddies = edge_uses.get(k, [])
                    if len(buddies) != 2:
                        continue
                    other_eid = (
                        buddies[0] if buddies[1] == eid else buddies[1]
                    )
                    a_other = float(abs_areas_before[other_eid])
                    denom = max(eid_area_after, a_other, 1e-30)
                    new_ac = (
                        abs(eid_area_after - a_other) / denom
                    )
                    if new_ac > area_ratio_target:
                        return True
                return False

            if _would_create_new_fail(larger, A_L_before / 2.0):
                continue
            if not is_boundary:
                t_eid = other_eid
                t_area_before = float(abs_areas_before[t_eid])
                if _would_create_new_fail(t_eid, t_area_before / 2.0):
                    continue

        # Force=True invocation of the existing split operator gives us
        # the new mesh without the C1/C2/alpha penalty gate firing; the
        # bespoke Pass E gates are checked below.
        if is_boundary:
            out = _apply_edge_split_boundary(
                mesh, larger, edge_local,
                alpha_target=alpha_target,
                min_angle_target=min_angle_target,
                max_angle_target=max_angle_target,
                edge_uses=edge_uses,
                edge_to_segment=edge_to_segment,
                coastline_projector=coastline_projector,
                force=True,
            )
        else:
            out = _apply_edge_split_interior(
                mesh, larger, edge_local,
                alpha_target=alpha_target,
                min_angle_target=min_angle_target,
                max_angle_target=max_angle_target,
                edge_uses=edge_uses,
                boundary_edge_keys=boundary_edge_keys,
                force=True,
            )
        if out is None:
            continue
        new_mesh, split_info = out

        # Gate (1): area_change at the original fail edge must drop.
        # The fail edge is still present in the new mesh (we split a
        # non-shared edge), but its two incident elements have changed.
        # Re-locate the new buddy pair for ``fail_key``.
        new_edge_uses = _edge_use_counts(new_mesh.elements)
        new_buddies = new_edge_uses.get(fail_key, [])
        if len(new_buddies) != 2:
            # The fail edge became a boundary edge or disappeared —
            # would be a topological bug here, so skip.
            continue
        areas_new = np.abs(
            _signed_areas(new_mesh.nodes, new_mesh.elements[new_buddies]),
        )
        A_L_after = float(areas_new.max())
        A_S_after = float(areas_new.min())
        ac_after = (
            (A_L_after - A_S_after) / max(A_L_after, 1e-30)
        )
        if ac_after + 1e-12 >= ac_before:
            continue  # no improvement on the target C4 fail

        # Gate (2): no new C1 or C2 fail in the affected block.
        # ``force=True`` skipped the standard penalty gate; check
        # the explicit min/max-angle thresholds here. Alpha is
        # deliberately not gated — see docstring.
        removed = split_info.get("removed_elements", [])
        if not removed:
            removed = [larger]
        n_new_tris = 2 * len(removed)
        new_eids = list(range(
            new_mesh.n_elements - n_new_tris, new_mesh.n_elements,
        ))
        block_after = new_mesh.elements[new_eids]
        _a_a, m_a, M_a = _per_element_quality(new_mesh.nodes, block_after)
        if (m_a < min_angle_target).any():
            continue  # would introduce a C1 fail
        if (M_a > max_angle_target).any():
            continue  # would introduce a C2 fail

        info_out = {
            "operator": "pass_e_split",
            "underlying": split_info["operator"],
            "fail_edge": [fail_key[0], fail_key[1]],
            "split_edge": split_info.get("edge"),
            "new_node": split_info.get("new_node"),
            "removed_elements": removed,
            "affected_elements": new_eids,
            "area_change_before": float(ac_before),
            "area_change_after": float(ac_after),
        }
        return new_mesh, info_out

    return None


def _pass_e_round(
    cur: Fort14Mesh, *,
    alpha_target: float, min_angle_target: float,
    max_angle_target: float = DEFAULT_MAX_ANGLE_TARGET,
    area_ratio_target: float = DEFAULT_AREA_RATIO_TARGET,
    max_valence: int = DEFAULT_MAX_VALENCE,
    max_splits: int,
    coastline_projector: CoastlineProjector | None = None,
) -> tuple[Fort14Mesh, int, int, int, int]:
    """Pass E round: scan internal edges for C4 failures, attempt a
    Pass-E swap (and falling back to split) on each in descending
    order of area_change. Each accept rebuilds the edge-buddy
    structures, current C4 fail set, and per-node valence vector on
    the new mesh, so cascade avoidance and the C5 prefilter both
    stay sound.

    Returns ``(updated_mesh, n_accepted, n_rejected, n_swap_accepts,
    n_split_accepts)``.
    """
    from fvcom_mesh_tools.diagnostics import (  # noqa: PLC0415
        node_valence,
    )

    n_accepted = 0
    n_rejected = 0
    n_swap_accepts = 0
    n_split_accepts = 0

    while n_accepted < max_splits:
        edge_uv, elem_pair, area_change = _per_edge_area_change(
            cur.nodes, cur.elements,
        )
        if area_change.size == 0:
            break
        fail_mask = area_change > area_ratio_target
        if not fail_mask.any():
            break
        fail_indices = np.where(fail_mask)[0]
        # Worst-first ordering — go after the largest area_change first.
        fail_indices = fail_indices[
            np.argsort(-area_change[fail_indices])
        ]

        edge_uses = _edge_use_counts(cur.elements)
        boundary_edge_keys = {
            k for k, v in edge_uses.items() if len(v) == 1
        }
        _bp, _bn, edge_to_segment = _boundary_topology(cur)
        valence_before = node_valence(cur.elements, cur.n_nodes)
        abs_areas_before = np.abs(_signed_areas(cur.nodes, cur.elements))
        c4_fail_keys = {
            (int(edge_uv[i, 0]), int(edge_uv[i, 1]))
            for i in fail_indices
        }

        progress = False
        for idx in fail_indices:
            u, v = int(edge_uv[idx, 0]), int(edge_uv[idx, 1])
            e_i, e_j = int(elem_pair[idx, 0]), int(elem_pair[idx, 1])
            # Try edge_swap first — it has no new vertices, no
            # midpoint to smooth later, and net zero valence change.
            swap_out = _apply_pass_e_swap(
                cur, (u, v), (e_i, e_j),
                alpha_target=alpha_target,
                min_angle_target=min_angle_target,
                max_angle_target=max_angle_target,
                area_ratio_target=area_ratio_target,
                max_valence=max_valence,
                edge_uses=edge_uses,
                boundary_edge_keys=boundary_edge_keys,
                c4_fail_keys=c4_fail_keys,
                valence_before=valence_before,
            )
            if swap_out is not None:
                out = swap_out
                accepted_op = "swap"
            else:
                out = _apply_pass_e_split(
                    cur, (u, v), (e_i, e_j),
                    alpha_target=alpha_target,
                    min_angle_target=min_angle_target,
                    max_angle_target=max_angle_target,
                    area_ratio_target=area_ratio_target,
                    max_valence=max_valence,
                    edge_uses=edge_uses,
                    boundary_edge_keys=boundary_edge_keys,
                    edge_to_segment=edge_to_segment,
                    c4_fail_keys=c4_fail_keys,
                    valence_before=valence_before,
                    abs_areas_before=abs_areas_before,
                    coastline_projector=coastline_projector,
                )
                accepted_op = "split"
            if out is None:
                n_rejected += 1
                continue
            cur, _info = out
            n_accepted += 1
            if accepted_op == "swap":
                n_swap_accepts += 1
            else:
                n_split_accepts += 1
            progress = True
            break  # restart enumeration on the new mesh

        if not progress:
            break

    return cur, n_accepted, n_rejected, n_swap_accepts, n_split_accepts


def phase_h_optimize(
    mesh: Fort14Mesh,
    *,
    alpha_target: float = DEFAULT_ALPHA_TARGET,
    min_angle_target: float = DEFAULT_MIN_ANGLE_TARGET,
    max_angle_target: float = DEFAULT_MAX_ANGLE_TARGET,
    max_smooth_sweeps: int = DEFAULT_MAX_SMOOTH_SWEEPS,
    max_topology_per_round: int = 10_000,
    max_outer_rounds: int = 10,
    operator_order: tuple[str, ...] = DEFAULT_OPERATOR_ORDER,
    coastline_projector: CoastlineProjector | None = None,
    lookahead_enabled: bool = False,
    max_lookahead_per_round: int = 10_000,
    lookahead_op1_inventory: tuple[str, ...] = DEFAULT_LOOKAHEAD_OP1_INVENTORY,
    lookahead_op2_inventory: tuple[str, ...] = DEFAULT_LOOKAHEAD_OP2_INVENTORY,
    lookahead_gate: str = DEFAULT_LOOKAHEAD_GATE,
    patch_recdt_enabled: bool = False,
    patch_min_cluster_size: int = DEFAULT_PATCH_MIN_CLUSTER_SIZE,
    patch_max_cluster_size: int = DEFAULT_PATCH_MAX_CLUSTER_SIZE,
    max_patches_per_round: int = 1_000,
    patch_reject_boundary_clusters: bool = True,
    pass_e_enabled: bool = False,
    pass_e_area_ratio_target: float = DEFAULT_AREA_RATIO_TARGET,
    pass_e_max_valence: int = DEFAULT_MAX_VALENCE,
    max_pass_e_splits_per_round: int = 10_000,
) -> tuple[Fort14Mesh, dict[str, Any]]:
    """Phase H driver: Pass A (batch smooth) ↔ Pass B (1-step topology)
    ↔ optional Pass C (2-step lookahead) ↔ optional Pass D
    (cluster patch re-CDT) ↔ optional Pass E (gradation refinement)
    until none make progress.

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

    **Pass C (v4 / v4.1, opt-in via ``lookahead_enabled=True``)** —
    2-step lookahead on the residual fail elements that Passes A
    and B cannot crack. For each fail element it tries op1
    (``force=True`` — validity-only) drawn from
    ``lookahead_op1_inventory``, optionally followed by op2 from
    ``lookahead_op2_inventory`` on the affected region. Acceptance
    is controlled by ``lookahead_gate``:

    * ``"target_exits_fail"`` (v4.1, default) — accept iff the
      target element ``E`` exits fail status on the post-op mesh
      (``alpha(E) >= alpha_target ∧ min_angle(E) >= min_angle_target``)
      or was removed by the op chain. This is the SMS-manual-edit
      standard: each accept guarantees the target is fixed.  An
      **op1-only** path is enabled — op1 alone is accepted if it
      already takes ``E`` out of fail, without searching op2.
    * ``"union_penalty"`` (v4 reproduction) — accept iff the
      penalty summed over elements touching ``op1.affected ∪
      op2.affected`` nodes strictly drops between the round-start
      mesh and the post-op mesh. Reproduces the PoC #45 baseline;
      use only for benchmarking against the v4 negative result.

    The PoC #44 dry-run measured 61 % "fixable" pairs under the
    ``union_penalty`` gate but PoC #45 showed the v4 driver
    delivers ≈ 0 net fail-count change due to local-gate thrashing.
    The ``target_exits_fail`` default sidesteps this by demanding
    each accept *strictly* fix the failing target. Default
    inventory restricts op1 ∈ ``{smooth_node, vertex_remove}`` and
    op2 = ``smooth_node`` from PoC #44's accepted-pair histogram.

    **Pass D (opt-in via ``patch_recdt_enabled=True``)** — cluster
    patch re-CDT. Finds connected components of fail elements under
    face-face adjacency, filters by
    ``[patch_min_cluster_size, patch_max_cluster_size]``, walks
    each cluster's rim polygon, and replaces the cluster with a
    pure-Delaunay re-triangulation of the rim. Accept iff (i) every
    new patch element passes the per-element gate, and (ii) the rim
    1-ring fail count outside the cluster does not increase
    (no 2-ring drift). Targets the 51 % of v3 residual fails that
    sit in size ≥ 3 clusters and are unreachable by 1-ring local
    edits (see ``docs/patch_re_cdt_design.md``).

    **Pass E (opt-in via ``pass_e_enabled=True``)** — gradation
    refinement targeting FVCOM manual criterion C4 (adjacent-element
    area_change <= ``pass_e_area_ratio_target``, default 0.5). Scans
    internal edges worst-area-change-first; for each C4 fail, first
    tries ``edge_swap`` on the fail edge itself (no new vertices,
    net-zero valence change). If swap fails to reduce the local
    fail count, falls back to ``edge_split`` on the longest non-
    shared edge of the larger triangle that is neither itself a
    C4 fail (cascade avoidance) nor would push a preserved non-
    shared edge of ``L`` over the C4 threshold (indirect-regression
    filter). Both operators gate against C1, C2, and C5
    (``pass_e_max_valence``, default 8 per FVCOM C5).

    The passes alternate in order A → B → C → D → E each outer round.
    Pass A runs to exhaustion, Pass B up to
    ``max_topology_per_round`` accepts, Pass C up to
    ``max_lookahead_per_round`` accepts, Pass D up to
    ``max_patches_per_round`` accepts, Pass E up to
    ``max_pass_e_splits_per_round`` accepts. After Pass D / Pass E
    accepts the Pass A loop cleans up perturbations on the next outer
    round. The loop terminates when every enabled pass contributes
    zero accepts in the same round, or when ``max_outer_rounds`` is
    hit.

    Returns ``(new_mesh, info)`` with operator histograms, sweep
    counts, lookahead pair-histogram, and pre/post quality summaries.
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
        "max_angle_target": float(max_angle_target),
        "max_smooth_sweeps": int(max_smooth_sweeps),
        "max_topology_per_round": int(max_topology_per_round),
        "max_outer_rounds": int(max_outer_rounds),
        "operator_order": list(operator_order),
        "operators_applied": defaultdict(int),
        "n_iters": 0,
        "n_smooth_sweeps": 0,
        "n_outer_rounds": 0,
        "n_abandoned": 0,
        "lookahead_enabled": bool(lookahead_enabled),
        "lookahead_op1_inventory": list(lookahead_op1_inventory),
        "lookahead_op2_inventory": list(lookahead_op2_inventory),
        "lookahead_gate": str(lookahead_gate),
        "max_lookahead_per_round": int(max_lookahead_per_round),
        "lookahead_pairs_applied": defaultdict(int),
        "n_lookahead_abandoned": 0,
        "patch_recdt_enabled": bool(patch_recdt_enabled),
        "patch_min_cluster_size": int(patch_min_cluster_size),
        "patch_max_cluster_size": int(patch_max_cluster_size),
        "max_patches_per_round": int(max_patches_per_round),
        "patch_reject_boundary_clusters": bool(patch_reject_boundary_clusters),
        "patch_recdt_accepts": defaultdict(int),
        "n_patch_recdt_rejected": 0,
        "pass_e_enabled": bool(pass_e_enabled),
        "pass_e_area_ratio_target": float(pass_e_area_ratio_target),
        "pass_e_max_valence": int(pass_e_max_valence),
        "max_pass_e_splits_per_round": int(max_pass_e_splits_per_round),
        "pass_e_accepts": 0,
        "pass_e_swap_accepts": 0,
        "pass_e_split_accepts": 0,
        "pass_e_rejected": 0,
    }
    if lookahead_gate not in LOOKAHEAD_GATES:
        raise ValueError(
            f"unknown lookahead_gate {lookahead_gate!r}; "
            f"expected one of {LOOKAHEAD_GATES}"
        )

    a0, m0, M0 = _per_element_quality(cur.nodes, cur.elements)
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

        # Pass A: batch Gauss-Seidel smooth (interior + boundary).
        if do_smooth:
            bnd_node = _boundary_node_mask(cur)
            n2e = _node_to_elements(cur.elements, cur.n_nodes)
            bnd_prev, bnd_next, _ = _boundary_topology(cur)
            for _sweep in range(max_smooth_sweeps):
                n_acc = _batch_smooth_sweep(
                    cur,
                    alpha_target=alpha_target,
                    min_angle_target=min_angle_target,

                    max_angle_target=max_angle_target,
                    boundary_node_mask=bnd_node,
                    n2e=n2e,
                    boundary_prev=bnd_prev,
                    boundary_next=bnd_next,
                    coastline_projector=coastline_projector,
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

                max_angle_target=max_angle_target,
                operator_order=operator_order,
                max_topology_accepts=max_topology_per_round,
                coastline_projector=coastline_projector,
            )
            for op_name, n in topo_acc.items():
                info["operators_applied"][op_name] += int(n)
                info["n_iters"] += int(n)
                round_accepts += int(n)
            info["n_abandoned"] = n_aband

        # Pass C (v4): 2-step lookahead, opt-in. Runs after Pass B so
        # it sees the freshest 1-step-exhausted residual.
        if lookahead_enabled:
            cur, pair_acc, n_la_aband = _lookahead_round(
                cur,
                alpha_target=alpha_target,
                min_angle_target=min_angle_target,

                max_angle_target=max_angle_target,
                op1_inventory=lookahead_op1_inventory,
                op2_inventory=lookahead_op2_inventory,
                max_lookahead_accepts=max_lookahead_per_round,
                coastline_projector=coastline_projector,
                gate=lookahead_gate,
            )
            for pair_label, n in pair_acc.items():
                info["lookahead_pairs_applied"][pair_label] += int(n)
                info["n_iters"] += int(n)
                round_accepts += int(n)
            info["n_lookahead_abandoned"] = n_la_aband

        # Pass D (v4.2): cluster patch re-CDT, opt-in. Runs last so it
        # sees the cluster residual after every cheaper pass had a go.
        if patch_recdt_enabled:
            cur, patch_accepts, n_patch_rej = _patch_recdt_round(
                cur,
                alpha_target=alpha_target,
                min_angle_target=min_angle_target,

                max_angle_target=max_angle_target,
                min_cluster_size=patch_min_cluster_size,
                max_cluster_size=patch_max_cluster_size,
                max_patches=max_patches_per_round,
                reject_boundary_clusters=patch_reject_boundary_clusters,
            )
            for bucket, n in patch_accepts.items():
                info["patch_recdt_accepts"][bucket] += int(n)
                info["n_iters"] += int(n)
                round_accepts += int(n)
            info["n_patch_recdt_rejected"] += int(n_patch_rej)

        # Pass E: gradation refinement (FVCOM C4). Targets internal
        # edges where adjacent-element area ratio exceeds the target;
        # splits the larger triangle's longest non-shared,
        # non-C4-fail edge. C5 (valence) is gated to protect the
        # FVCOM <=8 cap.
        if pass_e_enabled:
            cur, pe_acc, pe_rej, pe_swap, pe_split = _pass_e_round(
                cur,
                alpha_target=alpha_target,
                min_angle_target=min_angle_target,
                max_angle_target=max_angle_target,
                area_ratio_target=pass_e_area_ratio_target,
                max_valence=pass_e_max_valence,
                max_splits=max_pass_e_splits_per_round,
                coastline_projector=coastline_projector,
            )
            info["pass_e_accepts"] += int(pe_acc)
            info["pass_e_swap_accepts"] += int(pe_swap)
            info["pass_e_split_accepts"] += int(pe_split)
            info["pass_e_rejected"] += int(pe_rej)
            info["n_iters"] += int(pe_acc)
            round_accepts += int(pe_acc)

        if round_accepts == 0:
            break

    a_after, m_after, M_after = _per_element_quality(cur.nodes, cur.elements)
    fail_after = _is_fail(
        a_after, m_after, M_after,
        alpha_target=alpha_target, min_angle_target=min_angle_target,
        max_angle_target=max_angle_target,
    )
    info["n_output_fail"] = int(fail_after.sum())
    info["alpha_mean_after"] = (
        float(a_after.mean()) if a_after.size else float("nan")
    )
    info["n_nodes"] = int(cur.n_nodes)
    info["n_elements"] = int(cur.n_elements)
    info["operators_applied"] = dict(info["operators_applied"])
    info["lookahead_pairs_applied"] = dict(info["lookahead_pairs_applied"])
    info["patch_recdt_accepts"] = dict(info["patch_recdt_accepts"])
    return cur, info


__all__ = [
    "CoastlineProjector",
    "DEFAULT_ALPHA_TARGET",
    "DEFAULT_AREA_RATIO_TARGET",
    "DEFAULT_LOOKAHEAD_GATE",
    "DEFAULT_LOOKAHEAD_OP1_INVENTORY",
    "DEFAULT_LOOKAHEAD_OP2_INVENTORY",
    "DEFAULT_MAX_ANGLE_TARGET",
    "DEFAULT_MAX_VALENCE",
    "DEFAULT_MIN_ANGLE_TARGET",
    "DEFAULT_OPERATOR_ORDER",
    "DEFAULT_PATCH_MAX_CLUSTER_SIZE",
    "DEFAULT_PATCH_MIN_CLUSTER_SIZE",
    "LOOKAHEAD_GATES",
    "build_coastline_projector",
    "phase_h_optimize",
]
