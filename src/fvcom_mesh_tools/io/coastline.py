"""Load coastline shapefiles into a single ``MultiLineString`` for
OCSMesh ``Hfun.add_feature``.

The minimal pipeline only needs lon/lat polylines following the wet/dry
boundary; this module is deliberately small (no shapely buffering, no
re-segmenting) and accepts any geometry shapely understands. Polygons
are flattened into their exterior + interior rings as line strings.

The loader is independent of OCSMesh / pyproj so it stays importable
even when those optional dependencies are absent; geopandas + shapely
are required at call time and an informative error is raised otherwise.
"""

from __future__ import annotations

from pathlib import Path
from typing import Iterable


def load_coastline_as_lines(
    paths: Iterable[Path | str],
    *,
    bbox: tuple[float, float, float, float] | None = None,
):
    """Read one or more shapefiles / geo files into a ``MultiLineString``.

    Parameters
    ----------
    paths:
        Iterable of file paths recognised by ``geopandas.read_file``.
    bbox:
        ``(xmin, ymin, xmax, ymax)`` in EPSG:4326. When supplied, only
        features intersecting this rectangle are kept (massive speedup
        when the file covers a much larger region than the DEM).

    Returns
    -------
    shapely.geometry.MultiLineString
        Union of all line geometries reprojected to EPSG:4326. Empty
        ``MultiLineString`` if no inputs intersect ``bbox``.

    Raises
    ------
    ImportError:
        If geopandas / shapely are unavailable.
    ValueError:
        If a file cannot be loaded or contains an unsupported geometry
        type.
    """
    try:
        import geopandas as gpd
        from shapely.geometry import LineString, MultiLineString, box
    except ImportError as exc:  # pragma: no cover - environment guard
        raise ImportError(
            "load_coastline_as_lines requires geopandas + shapely; "
            "install via `mamba install -c conda-forge geopandas shapely`."
        ) from exc

    clip_box = box(*bbox) if bbox is not None else None

    lines: list = []
    for raw_path in paths:
        path = Path(raw_path)
        if not path.exists():
            raise ValueError(f"coastline file not found: {path}")
        gdf = gpd.read_file(path)
        # Tokyo Bay reference coastlines often ship without an explicit
        # CRS but are clearly in lon/lat. Default to EPSG:4326 in that
        # case rather than refusing to load.
        if gdf.crs is None:
            gdf = gdf.set_crs("EPSG:4326")
        else:
            gdf = gdf.to_crs("EPSG:4326")

        if clip_box is not None:
            gdf = gdf[gdf.geometry.intersects(clip_box)]

        for geom in gdf.geometry:
            if geom is None or geom.is_empty:
                continue
            gt = geom.geom_type
            if gt == "LineString":
                lines.append(geom)
            elif gt == "MultiLineString":
                lines.extend(geom.geoms)
            elif gt == "Polygon":
                _append_polygon_rings(lines, geom, LineString)
            elif gt == "MultiPolygon":
                for poly in geom.geoms:
                    _append_polygon_rings(lines, poly, LineString)
            else:
                raise ValueError(f"unsupported geometry type {gt} in {path}")

    if clip_box is not None and lines:
        clipped: list = []
        for ls in lines:
            inter = ls.intersection(clip_box)
            if inter.is_empty:
                continue
            if inter.geom_type == "LineString":
                clipped.append(inter)
            elif inter.geom_type == "MultiLineString":
                clipped.extend(inter.geoms)
            # GeometryCollection (rare): drop point parts, keep lines
            elif inter.geom_type == "GeometryCollection":
                for sub in inter.geoms:
                    if sub.geom_type == "LineString":
                        clipped.append(sub)
                    elif sub.geom_type == "MultiLineString":
                        clipped.extend(sub.geoms)
        lines = clipped

    return MultiLineString(lines)


def _append_polygon_rings(out: list, polygon, LineString) -> None:
    out.append(LineString(list(polygon.exterior.coords)))
    for ring in polygon.interiors:
        out.append(LineString(list(ring.coords)))
