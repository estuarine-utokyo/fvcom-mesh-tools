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
                elif op_name == "edge_split_boundary":
                    for k in range(3):
                        out = _apply_edge_split_boundary(
                            cur, int(eid), k,
                            alpha_target=alpha_target,
                            min_angle_target=min_angle_target,
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
    a, m = _per_element_quality(mesh.nodes, block)
    return float(_penalty(
        a, m,
        alpha_target=alpha_target, min_angle_target=min_angle_target,
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
            a, m = _per_element_quality(m_after.nodes, block)
            return bool(
                float(a[0]) >= alpha_target
                and float(m[0]) >= min_angle_target
            )
    return False


def _iter_op_candidates(
    mesh: Fort14Mesh, eid: int, op_name: str, *,
    force: bool, ctx: dict[str, Any],
    alpha_target: float, min_angle_target: float,
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
            ):
                return m1, f"{op1_name}+none"

            affected1_eids = _affected_elements_in_mesh(m1, info1)
            if not affected1_eids:
                continue
            block = m1.elements[affected1_eids]
            a_blk, m_blk = _per_element_quality(m1.nodes, block)
            fail_blk = _is_fail(
                a_blk, m_blk,
                alpha_target=alpha_target, min_angle_target=min_angle_target,
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
                        coastline_projector=coastline_projector,
                    ):
                        if gate == "target_exits_fail":
                            if _target_exits_fail(
                                m2, target_vset,
                                alpha_target=alpha_target,
                                min_angle_target=min_angle_target,
                            ):
                                return m2, f"{op1_name}+{op2_name}"
                        else:  # union_penalty
                            affected2_nodes = _affected_nodes_from_info(info2)
                            union = affected1_nodes | affected2_nodes
                            pen_before = _union_penalty(
                                cur, union,
                                alpha_target=alpha_target,
                                min_angle_target=min_angle_target,
                            )
                            pen_after = _union_penalty(
                                m2, union,
                                alpha_target=alpha_target,
                                min_angle_target=min_angle_target,
                            )
                            if pen_after + 1e-12 < pen_before:
                                return m2, f"{op1_name}+{op2_name}"
    return None


def _lookahead_round(
    cur: Fort14Mesh, *,
    alpha_target: float, min_angle_target: float,
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


def phase_h_optimize(
    mesh: Fort14Mesh,
    *,
    alpha_target: float = DEFAULT_ALPHA_TARGET,
    min_angle_target: float = DEFAULT_MIN_ANGLE_TARGET,
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
) -> tuple[Fort14Mesh, dict[str, Any]]:
    """Phase H driver: Pass A (batch smooth) ↔ Pass B (1-step topology)
    ↔ optional Pass C (2-step lookahead) until none make progress.

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

    The three passes alternate: Pass A runs to exhaustion, Pass B
    up to ``max_topology_per_round`` accepts, Pass C up to
    ``max_lookahead_per_round`` accepts. After Pass C accepts the
    Pass A loop cleans up perturbations on the next outer round.
    The loop terminates when every enabled pass contributes zero
    accepts in the same round, or when ``max_outer_rounds`` is hit.

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
    }
    if lookahead_gate not in LOOKAHEAD_GATES:
        raise ValueError(
            f"unknown lookahead_gate {lookahead_gate!r}; "
            f"expected one of {LOOKAHEAD_GATES}"
        )

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
    info["lookahead_pairs_applied"] = dict(info["lookahead_pairs_applied"])
    return cur, info


__all__ = [
    "CoastlineProjector",
    "DEFAULT_ALPHA_TARGET",
    "DEFAULT_LOOKAHEAD_GATE",
    "DEFAULT_LOOKAHEAD_OP1_INVENTORY",
    "DEFAULT_LOOKAHEAD_OP2_INVENTORY",
    "DEFAULT_MIN_ANGLE_TARGET",
    "DEFAULT_OPERATOR_ORDER",
    "LOOKAHEAD_GATES",
    "build_coastline_projector",
    "phase_h_optimize",
]
