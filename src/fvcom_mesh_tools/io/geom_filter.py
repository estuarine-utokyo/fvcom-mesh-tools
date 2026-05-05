"""Filter a shapely ``MultiPolygon`` by signed metric area.

OCSMesh's ``Geom(raster, zmax=...)`` extracts a multipolygon of every
wet-pixel patch in the input DEM. For Tokyo Bay this comes to ~500
polygons of which only the top 4 are real ocean - the rest are tiny
isolated wet pixels. The largest polygon also carries dozens of
interior holes (islands), most of which are smaller than the desired
mesh size and would be poorly resolved.

This module provides one helper that walks a multipolygon and drops:

* every outer polygon with metric area below ``min_polygon_area_m2``
* every interior hole (island) with metric area below
  ``min_island_area_m2``

The signed-metric area is computed in an auto-selected UTM zone
derived from the multipolygon centroid, so callers do not have to
plumb the zone through. The helper degrades gracefully if pyproj is
absent: it raises ``ImportError`` rather than misreporting areas.
"""

from __future__ import annotations

from typing import Any


def _auto_utm_zone(lon_deg: float) -> int:
    """Return the canonical UTM zone for a longitude in degrees."""
    z = int((lon_deg + 180.0) // 6.0) + 1
    if z < 1:
        z = 1
    if z > 60:
        z = 60
    return z


def filter_multipolygon_by_area(
    mp: Any,
    *,
    src_crs: Any = "EPSG:4326",
    metric_crs: Any | None = None,
    min_polygon_area_m2: float = 0.0,
    min_island_area_m2: float = 0.0,
) -> Any:
    """Drop polygons / holes below the supplied area thresholds.

    Parameters
    ----------
    mp:
        ``shapely.geometry.MultiPolygon`` (or ``Polygon``; promoted to
        a length-1 MultiPolygon).
    src_crs:
        CRS of ``mp`` for projection setup. Defaults to EPSG:4326.
    metric_crs:
        Metric CRS used for area calculations. If ``None``, an
        auto-selected UTM zone derived from the centroid longitude
        is used.
    min_polygon_area_m2:
        Outer polygons with metric area below this threshold are
        dropped. ``0`` (default) keeps every polygon.
    min_island_area_m2:
        Holes (islands) with metric area below this threshold are
        filled in by the surviving outer polygon. ``0`` (default)
        keeps every hole.

    Returns
    -------
    shapely.geometry.MultiPolygon
        New multipolygon with the filters applied. Empty if every
        polygon was dropped.
    """
    if min_polygon_area_m2 < 0 or min_island_area_m2 < 0:
        raise ValueError("area thresholds must be non-negative")

    try:
        from shapely.geometry import MultiPolygon, Polygon
        from shapely.ops import transform as shp_transform
    except ImportError as exc:  # pragma: no cover - environment guard
        raise ImportError(
            "filter_multipolygon_by_area requires shapely; "
            "install via `mamba install -c conda-forge shapely`."
        ) from exc

    # Promote a single Polygon to a length-1 MultiPolygon.
    if isinstance(mp, Polygon):
        mp = MultiPolygon([mp])
    if not isinstance(mp, MultiPolygon):
        raise TypeError(f"mp must be Polygon or MultiPolygon, got {type(mp)!r}")

    if min_polygon_area_m2 == 0 and min_island_area_m2 == 0:
        return mp

    try:
        from pyproj import CRS, Transformer
    except ImportError as exc:  # pragma: no cover - environment guard
        raise ImportError(
            "filter_multipolygon_by_area requires pyproj for metric area; "
            "install via `mamba install -c conda-forge pyproj`."
        ) from exc

    if metric_crs is None:
        cx, _cy = mp.centroid.x, mp.centroid.y
        zone = _auto_utm_zone(cx)
        metric_crs = f"+proj=utm +zone={zone} +ellps=WGS84"

    transformer = Transformer.from_crs(
        CRS(src_crs), CRS(metric_crs), always_xy=True,
    )
    fwd = lambda x, y: transformer.transform(x, y)  # noqa: E731

    out_polys: list[Polygon] = []
    for poly in mp.geoms:
        proj = shp_transform(fwd, poly)
        if proj.area < min_polygon_area_m2:
            continue
        if min_island_area_m2 > 0 and len(list(poly.interiors)) > 0:
            kept_holes = []
            for hole in poly.interiors:
                hole_poly = Polygon(list(hole.coords))
                hole_proj = shp_transform(fwd, hole_poly)
                if hole_proj.area >= min_island_area_m2:
                    kept_holes.append(list(hole.coords))
            new_poly = Polygon(list(poly.exterior.coords), holes=kept_holes)
        else:
            new_poly = poly
        out_polys.append(new_poly)

    return MultiPolygon(out_polys)
