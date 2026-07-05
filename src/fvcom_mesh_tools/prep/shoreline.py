"""OSM true-land acquisition and land-opening simplification.

The land-opening operator implements the goto2023 hand-editing
policy (see docs/DESIGN_HISTORY.md): thin artificial structures
(piers, breakwaters, islets narrower than the mesh scale) are
erased; port basins and river mouths are preserved because eroding
LAND can never disconnect WATER. Opening the water instead (the
2026-07-04 mistake) deletes exactly the essential features.
"""

from __future__ import annotations

from pathlib import Path
from typing import Any

__all__ = [
    "auto_utm_epsg",
    "cut_domain_at_obc_line",
    "extend_obc_ends_perpendicular",
    "default_land_shp",
    "fetch_true_land",
    "open_land",
]


def auto_utm_epsg(lon: float, lat: float) -> int:
    """WGS84 UTM EPSG code for a coordinate (northern/southern)."""
    zone = int((lon + 180.0) // 6.0) + 1
    zone = min(max(zone, 1), 60)
    return (32600 if lat >= 0 else 32700) + zone


def default_land_shp() -> Path | None:
    """OSM land-polygons source under ``$DATA_DIR`` (GENKAI layout),
    or None when DATA_DIR is unset / the file is absent."""
    import os

    data_dir = os.environ.get("DATA_DIR")
    if not data_dir:
        return None
    cand = (Path(data_dir) / "OSM" / "land-polygons-split-4326"
            / "land_polygons.shp")
    return cand if cand.exists() else None


def fetch_true_land(
    bbox: tuple[float, float, float, float],
    *,
    land_shp_path: Path | None = None,
    min_water_area_deg2: float = 1e-5,
    cache_dir: Path | None = None,
    force: bool = False,
):
    """OSM true-land polygons (land minus rivers/lakes/docks) for
    ``bbox`` = (lon_min, lat_min, lon_max, lat_max), via xcoast.

    Returns a GeoDataFrame in EPSG:4326. Downloads and caching are
    handled by xcoast; pass ``cache_dir`` to relocate its cache
    (default: xcoast's own, typically ``~/.coastmask``).
    """
    try:
        import xcoast
    except ImportError as exc:  # pragma: no cover - env-dependent
        raise ImportError(
            "fetch_true_land requires the xcoast package "
            "(install the local clone: pip install -e ../xcoast)"
        ) from exc

    kwargs: dict[str, Any] = {"min_water_area_deg2": min_water_area_deg2}
    if land_shp_path is None:
        land_shp_path = default_land_shp()
    if land_shp_path is not None:
        kwargs["land_shp_path"] = Path(land_shp_path)
    if cache_dir is not None:
        kwargs["cache_dir"] = Path(cache_dir)
    config = xcoast.CoastmaskConfig(**kwargs)
    mask = xcoast.load(tuple(bbox), config=config, force=force)
    gdf = mask.land_gdf
    if gdf.crs is None:
        gdf = gdf.set_crs(4326)
    return gdf.to_crs(4326)


def open_land(
    land_gdf,
    *,
    r_open_m: float = 150.0,
    min_island_area_m2: float = 3.6e5,
    clip_bbox: tuple[float, float, float, float] | None = None,
    simplify_deg: float = 5e-5,
    utm_epsg: int | None = None,
):
    """Morphological opening of the LAND: erode then dilate by
    ``r_open_m`` (metres), so land features thinner than
    ``2 * r_open_m`` (piers, breakwaters, thin islets) vanish while
    water connectivity is preserved. Islands smaller than
    ``min_island_area_m2`` after opening are dropped (sub-grid at the
    target mesh scale; goto2023 keeps no interior islets).

    ``land_gdf`` is polygons in EPSG:4326; returns polygons in
    EPSG:4326, lightly simplified (``simplify_deg`` strips buffer-arc
    micro-vertices only — keep it well below the mesh scale).
    """
    import geopandas as gpd
    from shapely import make_valid, unary_union
    from shapely.geometry import MultiPolygon, Polygon, box

    gdf = land_gdf
    if gdf.crs is None:
        gdf = gdf.set_crs(4326)
    gdf = gdf.to_crs(4326)
    if clip_bbox is not None:
        gdf = gpd.clip(gdf, box(*clip_bbox))
    if len(gdf) == 0:
        return gpd.GeoDataFrame(geometry=[], crs=4326)

    land_ll = unary_union([
        make_valid(g) for g in gdf.geometry
        if g is not None and not g.is_empty
    ])
    if utm_epsg is None:
        c = land_ll.centroid
        utm_epsg = auto_utm_epsg(c.x, c.y)

    land = gpd.GeoSeries([land_ll], crs=4326).to_crs(utm_epsg).iloc[0]
    land_open = make_valid(land.buffer(-r_open_m).buffer(r_open_m))

    polys: list = []

    def _collect(g):
        if isinstance(g, Polygon):
            polys.append(g)
        elif isinstance(g, MultiPolygon) or hasattr(g, "geoms"):
            for s in g.geoms:
                _collect(s)

    _collect(land_open)
    polys = [p for p in polys
             if p.is_valid and not p.is_empty
             and p.area >= min_island_area_m2]
    out = gpd.GeoDataFrame(geometry=polys, crs=utm_epsg).to_crs(4326)
    if simplify_deg > 0:
        out.geometry = out.geometry.simplify(
            simplify_deg, preserve_topology=True,
        )
    out = out[out.geometry.is_valid & ~out.geometry.is_empty]
    return out.reset_index(drop=True)


def cut_domain_at_obc_line(
    land_gdf,
    obc_line: list[tuple[float, float]],
    bbox: tuple[float, float, float, float],
):
    """Close the domain at an artificial open-boundary LINE: all
    water on the seaward side of ``obc_line`` becomes a wall polygon
    merged into the land, so the generated mesh ends exactly at the
    line (goto2023 practice: a short smooth arc at the Uraga narrows
    instead of box edges across the Sagami trough).

    ``obc_line`` runs from its southern/eastern end to its
    northern/western end (lon, lat); the wall fills the bbox region
    south of the line.
    """
    import geopandas as gpd
    import numpy as np
    from shapely import make_valid
    from shapely.geometry import Polygon

    # Swept-band wall: the arc extruded ~30 km along its seaward
    # normal. Unlike a bbox-corner ring, this keeps coastal waters
    # BEYOND the junctions (e.g. the Boso west coast north of the
    # southern junction, which the reference mesh retains) — any
    # water remnant seaward of the band is disconnected from the bay
    # and dropped by keep_components later.
    lon_min, lat_min, lon_max, lat_max = bbox
    pts = np.asarray(obc_line, dtype=float)
    seg = np.diff(pts, axis=0)
    nrm = np.stack([-seg[:, 1], seg[:, 0]], axis=1)
    nrm = nrm / np.maximum(
        np.linalg.norm(nrm, axis=1, keepdims=True), 1e-12,
    )
    vnorm = np.vstack([nrm[:1], 0.5 * (nrm[:-1] + nrm[1:]), nrm[-1:]])
    vnorm = vnorm / np.maximum(
        np.linalg.norm(vnorm, axis=1, keepdims=True), 1e-12,
    )
    center = np.asarray([
        0.5 * (lon_min + lon_max), 0.5 * (lat_min + lat_max),
    ])
    mid = pts[len(pts) // 2]
    if np.dot(center - mid, vnorm[len(pts) // 2]) > 0:
        vnorm = -vnorm  # normals must point AWAY from the bay
    band_deg = 30000.0 / 91000.0  # ~30 km in degrees (coarse is fine)
    outer = pts + vnorm * band_deg
    ring = list(map(tuple, pts)) + list(map(tuple, outer[::-1]))
    wall = make_valid(Polygon(ring))
    gdf = land_gdf
    if gdf.crs is None:
        gdf = gdf.set_crs(4326)
    out = gpd.GeoDataFrame(
        geometry=list(gdf.to_crs(4326).geometry) + [wall], crs=4326,
    )
    return out


def extend_obc_ends_perpendicular(
    obc_line: list[tuple[float, float]],
    land_gdf,
    *,
    seaward_m: float = 1500.0,
    overshoot_m: float = 800.0,
    utm_epsg: int | None = None,
) -> list[tuple[float, float]]:
    """Replace each end of the OBC line with a straight segment along
    the local coast NORMAL, so the open boundary meets the coastline
    exactly at right angles (stability requirement: oblique OBC-coast
    junctions are a known storm-surge instability source; the sample
    mesh practice is bold simplification near the boundary).

    The terminal segment runs from ``seaward_m`` off the coast along
    the normal to ``overshoot_m`` INTO the land (the wall/domain cut
    consumes the overshoot).
    """
    import geopandas as gpd
    import numpy as np
    import shapely
    from shapely import unary_union
    from shapely.geometry import LineString

    from fvcom_mesh_tools.prep.shoreline import auto_utm_epsg

    if utm_epsg is None:
        c = obc_line[len(obc_line) // 2]
        utm_epsg = auto_utm_epsg(c[0], c[1])
    land = unary_union(list(
        land_gdf.to_crs(utm_epsg).geometry.values
    ))
    line_utm = gpd.GeoSeries(
        [LineString(obc_line)], crs=4326,
    ).to_crs(utm_epsg).iloc[0]
    pts = list(line_utm.coords)

    def _fix_end(pts, end):
        P = np.asarray(pts[-1] if end else pts[0])
        boundary = land.boundary
        q = boundary.interpolate(boundary.project(shapely.Point(P)))
        Q = np.asarray([q.x, q.y])
        # local coast tangent from a short chord around Q
        s = boundary.project(shapely.Point(P))
        q1 = boundary.interpolate(max(0.0, s - 400.0))
        q2 = boundary.interpolate(s + 400.0)
        tvec = np.asarray([q2.x - q1.x, q2.y - q1.y])
        n = np.asarray([-tvec[1], tvec[0]])
        n = n / (np.linalg.norm(n) or 1.0)
        # orient the normal seaward (away from land): the seaward
        # side is where the original end point lies
        if np.dot(P - Q, n) < 0:
            n = -n
        land_pt = Q - n * overshoot_m
        # Trim arc points closer than ~1.3*seaward_m to the coast,
        # then BLEND smoothly from the arc tangent into the coast
        # normal with a quadratic Bezier. A hard corner here (the
        # earlier sea_pt dogleg) leaves a ~1 km wedge of water
        # between the wall and the coast that feature sizing treats
        # as a narrow strait -> a fan of 200-400 m junk elements at
        # the junction.
        keep = pts[:-1] if end else pts[1:]
        if end:
            while len(keep) > 2 and np.linalg.norm(
                np.asarray(keep[-1]) - Q
            ) < 1.3 * seaward_m:
                keep = keep[:-1]
            P0 = np.asarray(keep[-1])
            t_dir = P0 - np.asarray(keep[-2])
        else:
            while len(keep) > 2 and np.linalg.norm(
                np.asarray(keep[0]) - Q
            ) < 1.3 * seaward_m:
                keep = keep[1:]
            P0 = np.asarray(keep[0])
            t_dir = P0 - np.asarray(keep[1])
        t_dir = t_dir / (np.linalg.norm(t_dir) or 1.0)
        # control point: intersection of the tangent ray with the
        # normal line through Q (falls back to the midpoint)
        A = np.column_stack([t_dir, n])
        P1 = 0.5 * (P0 + Q)
        if abs(np.linalg.det(A)) > 1e-6:
            ab = np.linalg.solve(A, Q - P0)
            if 0.0 < ab[0] < 4.0 * seaward_m:
                P1 = P0 + t_dir * ab[0]
        ts = np.linspace(0.0, 1.0, 8)[1:]
        blend = [tuple((1 - u) ** 2 * P0 + 2 * u * (1 - u) * P1
                       + u ** 2 * Q) for u in ts]
        if end:
            return keep + blend + [tuple(land_pt)]
        return [tuple(land_pt)] + blend[::-1] + keep

    pts = _fix_end(pts, end=False)
    pts = _fix_end(pts, end=True)
    out = gpd.GeoSeries(
        [LineString(pts)], crs=utm_epsg,
    ).to_crs(4326).iloc[0]
    return [(float(x), float(y)) for x, y in out.coords]
