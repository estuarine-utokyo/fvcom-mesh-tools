"""Mathematical waterway detection from the shoreline geometry
(owner 2026-07-12: waterways must be found from the ORIGINAL OSM
data -- no reference mesh involved).

Pipeline (all geometry, resolution-aware):

1.  Water = the sea-connected component of ``domain - land``
    (shapefile coverage gaps fabricate phantom inland water).
2.  Morphological opening with radius ``0.5 * detect_factor *
    h_mesh`` splits water into WIDE bodies and NARROW corridors.
3.  Narrow corridors and the small wide pockets between them are
    chained into waterway NETWORKS (adjacency graph).
4.  Each network becomes an ARC (+ measured width profile) via the
    channel-arc machinery -- widening happens along the arc, never
    by isotropic polygon buffering.
5.  Classification (owner rules):
    * THROUGH (both ends reach the main body, and the along-shore
      detour between the contact points is > ``shortcut_ratio`` x
      the through length) or connecting the main body to a big
      basin (>= ``min_basin_cells`` cell-equivalents): KEEP -- the
      banks are pushed into land until ``widen_rows`` standard
      rows fit. The minimum mesh size is inviolable.
    * dead ends and small-basin feeders: CLOSE (fill as land).
    * a KEEP whose corridor cannot be carved without piercing a
      protected barrier is reported as BLOCKED, loudly -- never
      silently reclassified.
"""

from __future__ import annotations

from typing import Any

import numpy as np
import shapely
from shapely.geometry import MultiPolygon, Polygon
from shapely.ops import unary_union

from fvcom_mesh_tools.channel_arcs import (
    arc_from_points,
    carve_channel_corridor,
    snap_arc_to_channel,
)

__all__ = ["detect_waterways", "apply_waterway_policy"]


def _polys(geom) -> list[Polygon]:
    if geom is None or geom.is_empty:
        return []
    if isinstance(geom, Polygon):
        return [geom]
    if isinstance(geom, MultiPolygon):
        return [g for g in geom.geoms if not g.is_empty]
    return [g for g in getattr(geom, "geoms", [])
            if isinstance(g, Polygon) and not g.is_empty]


def _interior_points(geom, spacing):
    """Point cloud inside a (multi)polygon on a regular lattice
    (fallback: representative points of its parts)."""
    x0, y0, x1, y1 = geom.bounds
    nx = max(int((x1 - x0) / spacing) + 1, 2)
    ny = max(int((y1 - y0) / spacing) + 1, 2)
    xs = np.linspace(x0, x1, nx)
    ys = np.linspace(y0, y1, ny)
    gx, gy = np.meshgrid(xs, ys)
    pts = shapely.points(gx.ravel(), gy.ravel())
    inside = shapely.covers(geom, pts)
    out = np.column_stack([gx.ravel()[inside], gy.ravel()[inside]])
    if len(out) < 2:
        out = np.array([[g.representative_point().x,
                         g.representative_point().y]
                        for g in _polys(geom)])
    return out


def detect_waterways(
    land_union,
    domain_poly,
    *,
    h_mesh_m: float,
    obc_point,
    metric_scale: tuple[float, float],
    detect_factor: float = 1.8,
    min_basin_cells: int = 6,
    shortcut_ratio: float = 2.2,
    min_extent_cells: float = 2.0,
    big_deadend_cells: float = 6.0,
    max_canal_extent_cells: float = 15.0,
    min_canal_width_frac: float = 0.5,
) -> list[dict[str, Any]]:
    """Find sub-``detect_factor*h`` waterways and decide their
    fate. Returns one record per waterway network:
    ``{arc, width_m (profile), action: keep|close, kind, geometry,
    extent_cells, basin_cells}``.
    """
    sx, sy = metric_scale
    if abs(sx - sy) / max(sx, sy) > 0.35:
        raise ValueError("metric_scale too anisotropic; project first")
    scale = 0.5 * (sx + sy)
    h = h_mesh_m / scale
    r_open = 0.5 * detect_factor * h
    a_cell = (np.sqrt(3.0) / 4.0) * h * h

    water = domain_poly.difference(land_union)
    wpolys = _polys(water)
    if not wpolys:
        raise RuntimeError("domain minus land left no water")
    obc_pt = shapely.Point(obc_point)
    water = wpolys[int(np.argmin([obc_pt.distance(g)
                                  for g in wpolys]))]

    wide = water.buffer(-r_open).buffer(
        r_open * 1.02, join_style="mitre", mitre_limit=1.2)
    wide = wide.intersection(water)
    wide_parts = _polys(wide)
    if not wide_parts:
        raise RuntimeError("opening removed ALL water -- h_mesh_m "
                           "too large for this domain?")
    main_i = int(np.argmin([obc_pt.distance(g)
                            for g in wide_parts]))
    main_poly = wide_parts[main_i]
    areas = np.array([g.area for g in wide_parts])
    big = {k for k in range(len(wide_parts))
           if k == main_i or areas[k] >= min_basin_cells * a_cell}
    small = [k for k in range(len(wide_parts)) if k not in big]

    narrow_parts = [N for N in _polys(water.difference(wide))
                    if N.area >= 0.05 * a_cell]

    # NETWORKS: corridors + small pockets that touch each other
    from scipy.sparse import coo_matrix
    from scipy.sparse.csgraph import connected_components
    eps = 0.02 * h
    items = ([("n", N) for N in narrow_parts]
             + [("s", wide_parts[k]) for k in small])
    n_it = len(items)
    rows, cols = [], []
    for i in range(n_it):
        for j in range(i + 1, n_it):
            if items[i][1].distance(items[j][1]) < eps:
                rows.append(i)
                cols.append(j)
    if n_it:
        g2 = coo_matrix((np.ones(len(rows)), (rows, cols)),
                        shape=(n_it, n_it))
        _, lab = connected_components(g2 + g2.T, directed=False)
    else:
        lab = np.zeros(0, dtype=int)

    records: list[dict[str, Any]] = []
    for net in range(lab.max() + 1 if n_it else 0):
        geoms = [items[i][1] for i in np.where(lab == net)[0]]
        union = unary_union(geoms)
        hull_len = float(union.minimum_rotated_rectangle.length) / 2
        extent_cells = hull_len / h
        # anchors on WIDE water
        touches_big = [k for k in big if k != main_i
                       and union.distance(wide_parts[k]) < eps]
        basin_cells = (max(areas[k] / a_cell for k in touches_big)
                       if touches_big else 0.0)
        inter = union.buffer(eps).intersection(main_poly)
        pieces = sorted(_polys(inter), key=lambda p: -p.area)
        through = False
        if len(pieces) >= 2:
            A = pieces[0].representative_point()
            B = pieces[1].representative_point()
            d_thru = max(A.distance(B), 0.1 * h)
            rings = [main_poly.exterior, *main_poly.interiors]
            rA = int(np.argmin([r.distance(A) for r in rings]))
            rB = int(np.argmin([r.distance(B) for r in rings]))
            if rA != rB:
                through = True
            else:
                ring = rings[rA]
                sA, sB = ring.project(A), ring.project(B)
                d_arc = abs(sA - sB)
                d_arc = min(d_arc, ring.length - d_arc)
                through = d_arc > shortcut_ratio * d_thru
        connector = through or (len(pieces) >= 1
                                and bool(touches_big))
        worthy = extent_cells >= min_extent_cells
        # PORT-CANAL rule (goto2023 design, resolution-relative):
        # a dead-end waterway of substantial-but-BOUNDED extent
        # and canal-like width (mean >= min_canal_width_frac * h)
        # is kept -- port canal systems carry tidal prism. Rivers
        # fail it on one of the two axes: below-resolution ones
        # (Hanami-gawa) are too narrow, large ones are far longer
        # than any port canal.
        mean_w_cells = (union.area / max(hull_len, 1e-9)) / h
        big_canal = (big_deadend_cells <= extent_cells
                     <= max_canal_extent_cells
                     and mean_w_cells >= min_canal_width_frac)
        keep = (connector or big_canal) and worthy

        # a narrow piece that barely touches LAND is not a
        # waterway at all -- it is an opening artifact against the
        # domain boundary or a sliver of the main body. Filling it
        # would turn open water into land: IGNORE instead.
        blen = float(union.boundary.length)
        land_frac = (float(union.boundary.intersection(
            land_union.buffer(eps)).length) / blen
            if blen > 0 else 0.0)

        narrow_members = [g for kind_, g in
                          zip([items[i][0] for i in
                               np.where(lab == net)[0]], geoms)
                          if kind_ == "n"]
        main_piece = (max(narrow_members, key=lambda g: g.area)
                      if narrow_members else None)
        rec: dict[str, Any] = {
            "kind": ("through" if through else
                     "port" if touches_big else
                     "canal" if big_canal else
                     "dead-end"),
            "action": ("keep" if keep else
                       "close" if land_frac >= 0.5 else "ignore"),
            "land_frac": round(land_frac, 2),
            "mean_width_cells": round(float(mean_w_cells), 2),
            "main_piece": main_piece,
            "extent_cells": round(float(extent_cells), 1),
            "basin_cells": round(float(basin_cells), 1),
            "geometry": union,
            "arc": None,
            "width_m": None,
        }
        if keep:
            # waterway as an ARC: interior point cloud -> ordered
            # medial guide -> snapped centreline + width profile
            pts = _interior_points(union, 0.35 * h)
            guide = arc_from_points(
                np.column_stack([pts[:, 0] * sx / scale,
                                 pts[:, 1] * sy / scale]),
                smooth_passes=2)
            guide = np.column_stack([guide[:, 0] * scale / sx,
                                     guide[:, 1] * scale / sy])
            try:
                snap = snap_arc_to_channel(
                    land_union, guide, metric_scale=metric_scale,
                    step_m=0.35 * h_mesh_m,
                    max_halfwidth_m=2.5 * h_mesh_m)
                rec["arc"] = snap["arc"]
                rec["width_m"] = snap["width_carve_m"]
            except RuntimeError as e:
                rec["action"] = "blocked"
                rec["reason"] = f"arc extraction failed: {e}"
        records.append(rec)
    return records


def apply_waterway_policy(
    land_union,
    domain_poly,
    records: list[dict[str, Any]],
    *,
    h_mesh_m: float,
    metric_scale: tuple[float, float],
    widen_rows: float = 2.0,
    min_gap_m: float = 150.0,
    h_grade_per_m: float = 0.0,
    arc_retry: str = "largest-piece",
) -> tuple[Any, dict[str, Any]]:
    """Execute the detected actions: KEEP -> carve the corridor to
    two LOCAL rows along the arc, barrier-safe; CLOSE -> fill the
    network as land.

    The per-station widen target accounts for distance-to-coast
    size growth: ``0.875 * rows * (h + h_grade_per_m * w_i / 2)``
    -- at a 350 m natural channel the mid-channel mesh edge is
    already ~h + grade*175, and a target based on bare ``h``
    realises only ~1.6 rows.

    A keep whose carve is refused is retried once with an arc from
    its LARGEST single corridor piece (``arc_retry`` =
    "largest-piece"; diameter paths over BRANCHED networks cut
    corners); the retry is recorded on the record. Still refused
    -> BLOCKED (reported, land untouched); pass ``arc_retry=None``
    to disable the retry."""
    new_land = land_union
    info = {"kept": 0, "closed": 0, "ignored": 0, "blocked": [],
            "retried": 0, "land_removed_m2": 0.0}
    fills = []
    sx, sy = metric_scale
    scale = 0.5 * (sx + sy)

    def _carve(base, arc, widths):
        w_nat = np.asarray(widths, float)
        target = 0.875 * widen_rows * (
            h_mesh_m + h_grade_per_m * w_nat / 2.0)
        w = np.maximum(w_nat, target)
        return carve_channel_corridor(
            base, arc, w, min_gap_m=min_gap_m,
            metric_scale=metric_scale, domain_poly=domain_poly,
            arc_on_land_tol_m=0.3 * h_mesh_m,
            carve_crossings=False)

    for rec in records:
        if rec["action"] == "ignore":
            info["ignored"] += 1
        elif rec["action"] == "close":
            fills.append(rec["geometry"])
            info["closed"] += 1
        elif rec["action"] == "keep":
            try:
                new_land, ci = _carve(new_land, rec["arc"],
                                      rec["width_m"])
                info["kept"] += 1
                info["land_removed_m2"] += ci["land_removed_m2"]
            except RuntimeError as e:
                retried = False
                if (arc_retry == "largest-piece"
                        and rec.get("main_piece") is not None):
                    try:
                        pts = _interior_points(
                            rec["main_piece"],
                            0.35 * h_mesh_m / scale)
                        guide = arc_from_points(
                            np.column_stack(
                                [pts[:, 0] * sx / scale,
                                 pts[:, 1] * sy / scale]),
                            smooth_passes=2)
                        guide = np.column_stack(
                            [guide[:, 0] * scale / sx,
                             guide[:, 1] * scale / sy])
                        snap = snap_arc_to_channel(
                            land_union, guide,
                            metric_scale=metric_scale,
                            step_m=0.35 * h_mesh_m,
                            max_halfwidth_m=2.5 * h_mesh_m)
                        new_land, ci = _carve(
                            new_land, snap["arc"],
                            snap["width_carve_m"])
                        rec["arc"] = snap["arc"]
                        rec["width_m"] = snap["width_carve_m"]
                        rec["retry"] = "largest-piece"
                        info["kept"] += 1
                        info["retried"] += 1
                        info["land_removed_m2"] += (
                            ci["land_removed_m2"])
                        retried = True
                    except RuntimeError as e2:
                        e = e2
                if not retried:
                    rec["action"] = "blocked"
                    rec["reason"] = str(e)
        if rec["action"] == "blocked":
            info["blocked"].append(rec)
    if fills:
        new_land = unary_union([new_land, *fills])
    return new_land, info
