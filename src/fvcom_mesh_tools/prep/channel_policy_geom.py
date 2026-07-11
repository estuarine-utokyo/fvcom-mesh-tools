"""Geometry-stage narrow-channel policy (owner 2026-07-12).

Decides the fate of narrow water BEFORE meshing, on the true
shoreline geometry — deterministic (no DistMesh realisation
dependence) and minimum-mesh-size preserving:

* a channel narrower than one cell that CONNECTS THROUGH (both
  sides reach the main water body) or leads to a big basin
  (area >= ``min_basin_cells`` equilateral cells) is **widened by
  pushing its banks into land** until two standard-size rows fit
  — never by refining the mesh (the owner's minimum size is
  inviolable);
* a dead-end channel or one leading to a small basin is **closed**
  (filled as land) together with the basin, as the goto2023 sample
  does everywhere.

Method: morphological opening of the water polygon with radius
``0.6 x h_mesh`` classifies water into WIDE bodies and NARROW
corridors; corridor components are classified by which wide bodies
they touch; widening buffers the corridor by ``widen_factor x
h_mesh / 2`` and subtracts it from land; closing unions the
corridor (and its small basins) into land.
"""

from __future__ import annotations

from typing import Any

import numpy as np
import shapely
from shapely.geometry import MultiPolygon, Polygon
from shapely.ops import unary_union

__all__ = ["apply_channel_policy_to_land"]


def _polys(geom) -> list[Polygon]:
    if geom.is_empty:
        return []
    if isinstance(geom, Polygon):
        return [geom]
    if isinstance(geom, MultiPolygon):
        return [g for g in geom.geoms if not g.is_empty]
    return [g for g in getattr(geom, "geoms", [])
            if isinstance(g, Polygon) and not g.is_empty]


def apply_channel_policy_to_land(
    land_union,
    domain_poly,
    *,
    h_mesh_m: float,
    obc_point,
    min_basin_cells: int = 6,
    detect_factor: float = 1.2,
    widen_factor: float = 2.2,
    shortcut_ratio: float = 2.2,
    metric_scale: tuple[float, float] | None = None,
) -> tuple[Any, dict[str, Any]]:
    """Apply the policy to a land-polygon union within a domain.

    Parameters
    ----------
    land_union:
        Shapely (multi)polygon of LAND, same CRS as ``domain_poly``.
    domain_poly:
        The meshing domain polygon (water = domain - land).
    h_mesh_m:
        The minimum TARGET MESH size (metres) — the inviolable
        floor. Channels are widened to ``widen_factor * h_mesh_m``.
    obc_point:
        A (x, y) point in the main water body (e.g. mid-OBC);
        identifies which wide component is "main".
    min_basin_cells:
        Basin survival threshold in equilateral-cell equivalents
        (area >= n * sqrt(3)/4 * h^2). Owner: 6 (Funabashi class).
    detect_factor / widen_factor:
        Channels narrower than ``detect_factor*h`` are policy
        targets; widened ones get ``widen_factor*h`` of width.
    metric_scale:
        ``(sx, sy)`` factors converting coordinate units to metres
        (for lon/lat inputs pass ``(111e3*cos(lat), 111e3)``);
        ``None`` means coordinates are already metres.

    Returns ``(new_land_union, info)``.
    """
    sx, sy = metric_scale if metric_scale is not None else (1.0, 1.0)
    if abs(sx - sy) / max(sx, sy) > 0.35:
        raise ValueError(
            "metric_scale is too anisotropic for isotropic "
            "buffering; project the inputs first")
    scale = 0.5 * (sx + sy)
    r_open = 0.5 * detect_factor * h_mesh_m / scale
    r_widen = 0.5 * widen_factor * h_mesh_m / scale
    a_cell = (np.sqrt(3.0) / 4.0) * (h_mesh_m / scale) ** 2
    a_basin = min_basin_cells * a_cell

    water = domain_poly.difference(land_union)
    # keep ONLY the sea-connected component: where the land layer
    # has no coverage (inland shapefile gaps), domain-minus-land
    # fabricates phantom water bodies -- a 2,261-cell phantom in
    # inland Boso made river corridors look like big-basin
    # connectors and chain-widened the Hanami-gawa (2026-07-12)
    obc_pt0 = shapely.Point(obc_point)
    wpolys = _polys(water)
    if not wpolys:
        raise RuntimeError("domain minus land left no water")
    d0 = [obc_pt0.distance(g) for g in wpolys]
    sea = wpolys[int(np.argmin(d0))]
    n_phantom = len(wpolys) - 1
    water = sea
    wide = water.buffer(-r_open).buffer(
        r_open * 1.02, join_style="mitre", mitre_limit=1.2)
    wide = wide.intersection(water)
    narrow = water.difference(wide)
    wide_parts = _polys(wide)
    if not wide_parts:
        raise RuntimeError(
            "opening removed ALL water -- h_mesh_m too large for "
            "this domain?")
    obc_pt = shapely.Point(obc_point)
    main_i = int(np.argmin([obc_pt.distance(g) for g in wide_parts]))
    areas = np.array([g.area for g in wide_parts])

    info: dict[str, Any] = {"n_narrow": 0, "widened": [],
                            "closed": [], "n_wide_parts":
                            len(wide_parts),
                            "n_phantom_water_dropped": n_phantom}
    add_water = []       # widening: subtract from land
    add_land = []        # closing: union into land
    eps = 0.02 * h_mesh_m / scale
    for N in _polys(narrow):
        if N.area < 0.05 * a_cell:
            continue     # slivers of the opening operator
        info["n_narrow"] += 1
        nb = [k for k, g in enumerate(wide_parts)
              if N.distance(g) < eps]
        c = N.representative_point()
        rec = {"center": (float(c.x), float(c.y)),
               "area_cells": float(N.area / a_cell),
               "touches_main": main_i in nb,
               "neighbor_cells": sorted(
                   round(float(areas[k] / a_cell), 1)
                   for k in nb if k != main_i)}
        big_nonmain = [k for k in nb if k != main_i
                       and areas[k] >= a_basin]
        through = False
        if main_i in nb:
            inter = N.buffer(eps).intersection(wide_parts[main_i])
            pieces = sorted(_polys(inter), key=lambda g: -g.area)
            if len(pieces) >= 2:
                # SHORTCUT significance: a genuine through-channel
                # saves real distance (Keihin: ~2 km through vs
                # 8-10 km around the island group). Nearshore
                # ribbons behind small headlands and bendy river
                # mouths also touch main twice but save nothing
                # (ratio ~1.2-1.6) -- 126 of them eroded the whole
                # coastline before this gate (2026-07-12).
                A = pieces[0].representative_point()
                B = pieces[1].representative_point()
                d_thru = max(A.distance(B),
                             0.1 * h_mesh_m / scale)
                rings = [wide_parts[main_i].exterior,
                         *wide_parts[main_i].interiors]
                rA = int(np.argmin([r.distance(A) for r in rings]))
                rB = int(np.argmin([r.distance(B) for r in rings]))
                if rA != rB:
                    through = True   # crosses to a hole boundary
                else:
                    ring = rings[rA]
                    sA = ring.project(A)
                    sB = ring.project(B)
                    d_arc = abs(sA - sB)
                    d_arc = min(d_arc, ring.length - d_arc)
                    through = d_arc > shortcut_ratio * d_thru
        # widen ONLY if the corridor CONNECTS two water bodies:
        # a significant shortcut through the main body, or a
        # main <-> big-basin port channel. A corridor touching a
        # single body -- however big -- is a river/inlet mouth and
        # must be closed (a river mouth touches the sea once).
        if through or (main_i in nb and big_nonmain):
            add_water.append(N.buffer(r_widen))
            rec["action"] = "widen"
            info["widened"].append(rec)
        else:
            add_land.append(N)
            for k in nb:
                if k != main_i and areas[k] < a_basin:
                    add_land.append(wide_parts[k])
            rec["action"] = "close"
            info["closed"].append(rec)

    new_land = land_union
    if add_land:
        new_land = unary_union([new_land, *add_land])
    if add_water:
        new_land = new_land.difference(unary_union(add_water))
    return new_land, info
