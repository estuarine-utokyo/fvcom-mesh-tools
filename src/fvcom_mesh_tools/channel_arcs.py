"""Channel ARCS: waterways as first-class polyline objects
(adopted from the SCHISM RiverMapper/OCSMesh design, owner
2026-07-12).

A channel is represented by its along-channel ARC (thalweg-like
polyline) plus a width. Meshing-side effects are derived FROM the
arc, never from polygon morphology:

* the water corridor is the arc buffered by ``width/2`` with FLAT
  caps -- erosion of land is confined to the arc's own corridor,
  transversal land barriers cannot be pierced;
* an explicit barrier check subtracts a protective buffer around
  every OTHER water surface and REFUSES (raises) if that splits
  the corridor -- a fabricated connection is structurally
  impossible;
* connectivity is explicit: the corridor either carries the arc
  end-to-end or the operation fails loudly.

Pure numpy/shapely; no oceanmesh import (license policy).
"""

from __future__ import annotations

from typing import Any

import numpy as np
import shapely
from shapely.geometry import LineString, MultiPolygon, Polygon
from shapely.ops import unary_union

__all__ = ["arc_from_points", "carve_channel_corridor"]


def _polys(geom) -> list[Polygon]:
    if geom is None or geom.is_empty:
        return []
    if isinstance(geom, Polygon):
        return [geom]
    if isinstance(geom, MultiPolygon):
        return [g for g in geom.geoms if not g.is_empty]
    return [g for g in getattr(geom, "geoms", [])
            if isinstance(g, Polygon) and not g.is_empty]


def arc_from_points(
    pts_xy: np.ndarray,
    *,
    smooth_passes: int = 2,
) -> np.ndarray:
    """Order scattered channel points (e.g. element centroids of a
    reference mesh's channel band) into an along-channel polyline:
    the diameter path of their nearest-neighbour graph, smoothed.
    """
    from scipy.sparse import coo_matrix
    from scipy.sparse.csgraph import shortest_path
    from scipy.spatial import cKDTree

    pts = np.asarray(pts_xy, dtype=float)
    n = len(pts)
    if n < 2:
        raise ValueError("arc needs at least 2 points")
    if n == 2:
        return pts.copy()
    k = min(5, n - 1)
    d, idx = cKDTree(pts).query(pts, k=k + 1)
    rows = np.repeat(np.arange(n), k)
    cols = idx[:, 1:].ravel()
    vals = d[:, 1:].ravel()
    g = coo_matrix((vals, (rows, cols)), shape=(n, n)).tocsr()
    dm, pred = shortest_path(g, directed=False,
                             return_predecessors=True)
    dm[~np.isfinite(dm)] = -1.0
    i, j = np.unravel_index(int(np.argmax(dm)), dm.shape)
    path = [int(j)]
    while path[-1] != i:
        p = int(pred[i, path[-1]])
        if p < 0:
            break
        path.append(p)
    arc = pts[path[::-1]]
    for _ in range(max(0, smooth_passes)):
        if len(arc) > 2:
            arc[1:-1] = (arc[:-2] + arc[1:-1] + arc[2:]) / 3.0
    return arc


def carve_channel_corridor(
    land_union,
    arc_ll: np.ndarray,
    width_m: float,
    *,
    min_gap_m: float,
    metric_scale: tuple[float, float],
    domain_poly=None,
    arc_on_land_tol_m: float | None = None,
) -> tuple[Any, dict[str, Any]]:
    """Cut a channel corridor of ``width_m`` along ``arc_ll`` out
    of ``land_union``, barrier-safely.

    The corridor is the flat-capped buffer of the arc. Guards
    (both fail loudly, never silently):

    1. the ARC must follow existing water -- crossing land for
       more than ``arc_on_land_tol_m`` (default ``width_m``)
       raises: a wrong arc would pierce a transversal barrier;
    2. every OTHER water surface is kept behind ``min_gap_m`` of
       land (arc-end openings exempted); if that protection stops
       the corridor from carrying the arc end to end, raise.

    Returns ``(new_land_union, info)``.
    """
    sx, sy = metric_scale
    if abs(sx - sy) / max(sx, sy) > 0.35:
        raise ValueError("metric_scale too anisotropic; project first")
    scale = 0.5 * (sx + sy)
    pts = np.asarray(arc_ll, dtype=float)
    arc = LineString(pts)
    w2 = 0.5 * width_m / scale

    # INVARIANT: the arc runs along EXISTING water (a channel
    # centreline); widening pushes the BANKS sideways into land.
    # An arc that crosses land for more than a noise tolerance is
    # a wrong arc -- it would pierce a transversal barrier, the
    # exact bug class this module exists to prevent. Fail loudly.
    tol = (arc_on_land_tol_m if arc_on_land_tol_m is not None
           else width_m)
    on_land = arc.intersection(land_union)
    on_land_m = float(on_land.length) * scale
    if on_land_m > tol:
        raise RuntimeError(
            f"channel arc crosses land for {on_land_m:.0f} m "
            f"(> tolerance {tol:.0f} m): the arc does not follow "
            "an existing waterway and would pierce a land barrier. "
            "Fix the arc; do not raise arc_on_land_tol_m unless "
            "the crossing is shoreline-data noise.")

    corridor = arc.buffer(w2, cap_style="flat")

    if domain_poly is not None:
        water = domain_poly.difference(land_union)
        other = water.difference(corridor.buffer(1e-9))
        protect = unary_union(_polys(other)).buffer(
            min_gap_m / scale)
        # the arc ENDS are where the channel must open into water:
        # exempt a disk around each end from protection, otherwise
        # the receiving water body's own buffer seals the corridor
        end_zones = unary_union(
            [shapely.Point(pts[0]).buffer(w2 + min_gap_m / scale),
             shapely.Point(pts[-1]).buffer(w2 + min_gap_m / scale)])
        protect = protect.difference(end_zones)
        carved = corridor.difference(protect)
    else:
        carved = corridor

    # the surviving corridor must still carry the arc end-to-end
    keep = None
    for gpart in _polys(carved):
        if gpart.buffer(1e-9).intersects(
                shapely.Point(arc.coords[0])) and \
           gpart.buffer(1e-9).intersects(
                shapely.Point(arc.coords[-1])):
            keep = gpart
            break
    if keep is None:
        raise RuntimeError(
            "channel corridor cannot be carved without piercing a "
            "land barrier protecting other water (min_gap "
            f"{min_gap_m:.0f} m). Options: reduce width_m, adjust "
            "the arc, or explicitly accept the connection by "
            "carving with domain_poly=None.")
    new_land = land_union.difference(keep)
    return new_land, {
        "corridor_area_m2": float(keep.area * scale * scale),
        "arc_length_m": float(arc.length * scale),
        "arc_on_land_m": on_land_m,
        "width_m": float(width_m),
    }
