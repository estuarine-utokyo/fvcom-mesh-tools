"""Exact boundary conformity: snap mesh boundary nodes onto reference
polylines (coastline) or explicit straight lines (open boundary).

Motivation (user 2026-07-04): ``om.Shoreline`` simplifies the
coastline to h0 before meshing, so the generated boundary never lies
ON the source coastline even where the coast is smooth and perfectly
representable. The SMS workflow fixes this by hand — load the
coastline, drag boundary nodes onto it. These functions automate
that edit:

* :func:`snap_boundary_to_polylines` — every (land-)boundary node
  moves to its nearest point on the reference polylines **iff** the
  move is smaller than ``max_snap_frac`` x the local boundary edge
  length and flips no incident triangle. On smooth coastline the
  nearest point is close, so nodes land EXACTLY on the data; at
  harbours/breakwaters finer than the mesh the distance cap leaves
  the simplified shape in place — precisely the requested
  "strict on smooth sections, simplified elsewhere" constraint.
* :func:`snap_nodes_to_segment` — project chosen nodes exactly onto
  a straight segment (smooth open-boundary lines), with the same
  flip guard.

Quality damage from snapping is expected and left to the standard
repair loop (``phase_h_finish`` with a coastline projector keeps the
snapped nodes ON the polylines while it repairs).
"""

from __future__ import annotations

from pathlib import Path
from typing import Any, Sequence

import numpy as np

from fvcom_mesh_tools.io import Fort14Mesh


def _signed_areas(nodes: np.ndarray, elements: np.ndarray) -> np.ndarray:
    p0, p1, p2 = nodes[elements[:, 0]], nodes[elements[:, 1]], nodes[elements[:, 2]]
    return 0.5 * (
        (p1[:, 0] - p0[:, 0]) * (p2[:, 1] - p0[:, 1])
        - (p1[:, 1] - p0[:, 1]) * (p2[:, 0] - p0[:, 0])
    )


def _boundary_topology(elements: np.ndarray, n_nodes: int):
    """(boundary_uv, node->incident elements map)."""
    raw = np.vstack([
        elements[:, [0, 1]], elements[:, [1, 2]], elements[:, [2, 0]],
    ])
    raw.sort(axis=1)
    codes = raw[:, 0].astype(np.int64) * n_nodes + raw[:, 1]
    uniq, counts = np.unique(codes, return_counts=True)
    b = uniq[counts == 1]
    boundary_uv = np.column_stack([b // n_nodes, b % n_nodes])
    rows = elements.ravel()
    cols = np.repeat(np.arange(elements.shape[0]), 3)
    order = np.argsort(rows, kind="stable")
    rows, cols = rows[order], cols[order]
    starts = np.searchsorted(rows, np.arange(n_nodes + 1))
    n2e = {v: cols[starts[v]:starts[v + 1]] for v in range(n_nodes)
           if starts[v] < starts[v + 1]}
    return boundary_uv, n2e


def _move_ok(nodes, elements, n2e, node, p_new) -> bool:
    ring = n2e.get(int(node))
    if ring is None:
        return False
    old = nodes[node].copy()
    nodes[node] = p_new
    ok = bool((_signed_areas(nodes, elements[ring]) > 0).all())
    nodes[node] = old
    return ok


def load_polylines(path: str | Path, *, to_crs: str | int | None = None) -> list:
    """Load every polyline (LineString / boundary of Polygon) from a
    vector file; optionally reproject. Files saved WITHOUT crs
    metadata are used as-is."""
    import geopandas as gpd

    gdf = gpd.read_file(Path(path))
    if to_crs is not None and gdf.crs is not None:
        gdf = gdf.to_crs(to_crs)
    lines: list = []
    for geom in gdf.geometry:
        if geom is None or geom.is_empty:
            continue
        gt = geom.geom_type
        if gt == "LineString":
            lines.append(geom)
        elif gt == "MultiLineString":
            lines.extend(list(geom.geoms))
        elif gt == "Polygon":
            lines.append(geom.boundary)
        elif gt == "MultiPolygon":
            lines.extend(sub.boundary for sub in geom.geoms)
    return lines


def snap_boundary_to_polylines(
    mesh: Fort14Mesh,
    polylines: Sequence,
    *,
    max_snap_frac: float = 0.6,
    exclude_nodes: Sequence[int] | None = None,
) -> tuple[Fort14Mesh, dict[str, Any]]:
    """Snap boundary nodes onto the nearest point of ``polylines``.

    ``polylines`` are shapely lines in the SAME coordinate system as
    ``mesh.nodes`` (use :func:`load_polylines`). Open-boundary nodes
    are always excluded; pass ``exclude_nodes`` for more. The move is
    applied only when it is shorter than ``max_snap_frac`` times the
    local boundary edge length AND flips no incident triangle.

    Returns a new mesh plus stats (n_snapped / n_far / n_flip_reverted
    and boundary-to-polyline distance percentiles before/after).
    """
    import shapely
    from shapely.strtree import STRtree

    nodes = mesh.nodes.copy()
    elements = mesh.elements
    n_nodes = mesh.n_nodes
    boundary_uv, n2e = _boundary_topology(elements, n_nodes)

    excl = set(int(v) for s in mesh.open_boundaries for v in s)
    if exclude_nodes is not None:
        excl.update(int(v) for v in exclude_nodes)

    # Local scale: mean length of the boundary edges at each node.
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
    bnodes = np.array([v for v in bnodes if int(v) not in excl], dtype=np.int64)

    tree = STRtree(list(polylines))
    geoms = np.array(list(polylines), dtype=object)
    pts = shapely.points(nodes[bnodes, 0], nodes[bnodes, 1])
    nearest_idx = tree.nearest(pts)
    d_before = shapely.distance(pts, geoms[nearest_idx])

    n_snapped = n_far = n_flip = 0
    for k, v in enumerate(bnodes):
        cap = max_snap_frac * local_h[v]
        if d_before[k] > cap:
            n_far += 1
            continue
        if d_before[k] <= 1e-9:
            continue
        line = geoms[nearest_idx[k]]
        p = line.interpolate(line.project(pts[k]))
        p_new = np.array([p.x, p.y])
        if _move_ok(nodes, elements, n2e, int(v), p_new):
            nodes[int(v)] = p_new
            n_snapped += 1
        else:
            n_flip += 1

    pts_after = shapely.points(nodes[bnodes, 0], nodes[bnodes, 1])
    d_after = shapely.distance(pts_after, geoms[tree.nearest(pts_after)])

    out = Fort14Mesh(
        title=mesh.title, nodes=nodes, depths=mesh.depths.copy(),
        elements=mesh.elements.copy(),
        open_boundaries=[np.asarray(s).copy() for s in mesh.open_boundaries],
        land_boundaries=[(ib, np.asarray(s).copy())
                         for ib, s in mesh.land_boundaries],
    )
    info = {
        "n_boundary_nodes": int(bnodes.size),
        "n_snapped": n_snapped,
        "n_far": n_far,
        "n_flip_reverted": n_flip,
        "dist_before_p50_m": float(np.percentile(d_before, 50)),
        "dist_before_p90_m": float(np.percentile(d_before, 90)),
        "dist_after_p50_m": float(np.percentile(d_after, 50)),
        "dist_after_p90_m": float(np.percentile(d_after, 90)),
    }
    return out, info


def snap_nodes_to_segment(
    mesh: Fort14Mesh,
    node_ids: Sequence[int],
    p0: tuple[float, float],
    p1: tuple[float, float],
    *,
    max_move: float | None = None,
) -> tuple[Fort14Mesh, dict[str, Any]]:
    """Project ``node_ids`` exactly onto the segment ``p0-p1`` (the
    smooth open-boundary line), clamped to the segment, with the flip
    guard. ``max_move`` (same units as coords) skips larger moves."""
    nodes = mesh.nodes.copy()
    elements = mesh.elements
    _uv, n2e = _boundary_topology(elements, mesh.n_nodes)
    a = np.asarray(p0, dtype=float)
    b = np.asarray(p1, dtype=float)
    ab = b - a
    denom = float(ab @ ab)
    n_snapped = n_skipped = n_flip = 0
    for v in node_ids:
        v = int(v)
        t = float((nodes[v] - a) @ ab) / denom
        t = min(max(t, 0.0), 1.0)
        p_new = a + t * ab
        move = float(np.hypot(*(p_new - nodes[v])))
        if move <= 1e-9:
            continue
        if max_move is not None and move > max_move:
            n_skipped += 1
            continue
        if _move_ok(nodes, elements, n2e, v, p_new):
            nodes[v] = p_new
            n_snapped += 1
        else:
            n_flip += 1
    out = Fort14Mesh(
        title=mesh.title, nodes=nodes, depths=mesh.depths.copy(),
        elements=mesh.elements.copy(),
        open_boundaries=[np.asarray(s).copy() for s in mesh.open_boundaries],
        land_boundaries=[(ib, np.asarray(s).copy())
                         for ib, s in mesh.land_boundaries],
    )
    return out, {
        "n_snapped": n_snapped,
        "n_skipped_far": n_skipped,
        "n_flip_reverted": n_flip,
    }


__all__ = [
    "load_polylines",
    "snap_boundary_to_polylines",
    "snap_nodes_to_segment",
]
