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
    from shapely import make_valid
    from shapely.geometry import Polygon

    lon_min, lat_min, lon_max, lat_max = bbox
    pts = list(obc_line)
    s_end, n_end = pts[0], pts[-1]
    ring = (
        pts
        + [(lon_min, n_end[1]), (lon_min, lat_min),
           (lon_max, lat_min), (lon_max, s_end[1])]
        + [pts[0]]
    )
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
        sea_pt = Q + n * seaward_m
        land_pt = Q - n * overshoot_m
        if end:
            return pts[:-1] + [tuple(sea_pt), tuple(land_pt)]
        return [tuple(land_pt), tuple(sea_pt)] + pts[1:]

    pts = _fix_end(pts, end=False)
    pts = _fix_end(pts, end=True)
    out = gpd.GeoSeries(
        [LineString(pts)], crs=utm_epsg,
    ).to_crs(4326).iloc[0]
    return [(float(x), float(y)) for x, y in out.coords]
