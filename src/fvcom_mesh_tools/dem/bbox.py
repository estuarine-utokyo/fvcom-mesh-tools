"""Read the lon/lat extent of a DEM raster.

``fmesh-buildmesh`` needs the DEM bounding box for the open/land
boundary classification step, regardless of which mesh engine
generated the mesh. Keeping this read in the engine-agnostic ``dem``
subpackage avoids tying boundary classification to either engine
adapter.
"""

from __future__ import annotations

from pathlib import Path


def read(dem_path: str | Path) -> tuple[float, float, float, float]:
    """Return ``(minx, miny, maxx, maxy)`` of ``dem_path``.

    The bounds are reported in the raster's native CRS units. For
    EPSG:4326 inputs that is ``(minlon, minlat, maxlon, maxlat)`` in
    degrees; for projected inputs it is metres or whatever unit the
    raster carries. Callers that need lon/lat unconditionally must
    reproject themselves.
    """
    import rasterio

    with rasterio.open(dem_path) as r:
        b = r.bounds
        return float(b.left), float(b.bottom), float(b.right), float(b.top)
