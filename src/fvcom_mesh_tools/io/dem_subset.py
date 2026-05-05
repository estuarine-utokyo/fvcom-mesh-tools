"""DEM subsetting helper for ``fmesh-buildmesh``.

The ``fmesh-buildmesh`` pipeline consumes a single, CRS-tagged raster.
Global DEMs (SRTM15+, GEBCO 2024, ...) are too big to feed in directly,
and several of them ship without an embedded CRS so they cannot be
opened by ``ocsmesh.Raster`` / ``rasterio`` as-is.

:func:`subset_dem_to_geotiff` clips a source DEM to a lon/lat bounding
box and writes a CF-tagged GeoTIFF that the downstream pipeline can
ingest. Two source families are supported:

1. **Rasters that rasterio can open with an embedded CRS** (GeoTIFF,
   CF-compliant NetCDF). We do a windowed read in source pixel space
   and re-emit a GeoTIFF in the same CRS.
2. **lon/lat NetCDF without CRS** (the SRTM15+ / GEBCO style). We
   read the named scalar variable plus 1-D ``lon`` and ``lat``
   coordinates via netCDF4 directly, slice in lon/lat space, and
   tag the output with ``--src-crs`` (default ``EPSG:4326``).
"""

from __future__ import annotations

from pathlib import Path

import numpy as np
import rasterio
from rasterio.transform import from_origin
from rasterio.windows import from_bounds


def _subset_via_rasterio(
    src_path: Path,
    dst_path: Path,
    bbox: tuple[float, float, float, float],
) -> dict:
    minx, miny, maxx, maxy = bbox
    with rasterio.open(src_path) as src:
        if src.crs is None:
            raise ValueError(
                f"{src_path}: no embedded CRS; pass --src-var to use the "
                "lon/lat NetCDF path instead."
            )
        window = from_bounds(minx, miny, maxx, maxy, transform=src.transform)
        window = window.round_offsets().round_lengths()
        if window.width <= 0 or window.height <= 0:
            raise ValueError(f"bbox {bbox} does not intersect raster {src_path}.")
        data = src.read(1, window=window)
        transform = rasterio.windows.transform(window, src.transform)
        profile = {
            "driver": "GTiff",
            "height": int(window.height),
            "width": int(window.width),
            "count": 1,
            "dtype": data.dtype,
            "crs": src.crs,
            "transform": transform,
            "compress": "lzw",
        }
    with rasterio.open(dst_path, "w", **profile) as out:
        out.write(data, 1)
    return {
        "path_dependent": "rasterio",
        "shape": (int(window.height), int(window.width)),
        "crs": str(profile["crs"]),
        "bbox": bbox,
    }


def _subset_via_netcdf(
    src_path: Path,
    dst_path: Path,
    bbox: tuple[float, float, float, float],
    src_var: str,
    src_crs: str,
) -> dict:
    import netCDF4

    minlon, minlat, maxlon, maxlat = bbox
    ds = netCDF4.Dataset(src_path)
    try:
        if "lon" not in ds.variables or "lat" not in ds.variables:
            raise ValueError(
                f"{src_path}: NetCDF lon/lat path requires both 'lon' and "
                "'lat' coordinate variables."
            )
        if src_var not in ds.variables:
            raise ValueError(
                f"{src_path}: source variable '{src_var}' not found in NetCDF."
            )
        lon = ds.variables["lon"][:]
        lat = ds.variables["lat"][:]
        ix = np.where((lon >= minlon) & (lon <= maxlon))[0]
        iy = np.where((lat >= minlat) & (lat <= maxlat))[0]
        if ix.size == 0 or iy.size == 0:
            raise ValueError(f"bbox {bbox} does not intersect NetCDF {src_path}.")
        z = ds.variables[src_var][iy[0]: iy[-1] + 1, ix[0]: ix[-1] + 1]
        sub_lon = lon[ix[0]: ix[-1] + 1]
        sub_lat = lat[iy[0]: iy[-1] + 1]
    finally:
        ds.close()

    z = np.asarray(z, dtype="float32")
    z_north_down = z[::-1, :]
    dx = float(sub_lon[1] - sub_lon[0])
    dy = float(sub_lat[1] - sub_lat[0])
    transform = from_origin(
        sub_lon[0] - dx / 2.0,
        sub_lat[-1] + dy / 2.0,
        dx,
        dy,
    )
    profile = {
        "driver": "GTiff",
        "height": z_north_down.shape[0],
        "width": z_north_down.shape[1],
        "count": 1,
        "dtype": "float32",
        "crs": src_crs,
        "transform": transform,
        "compress": "lzw",
    }
    with rasterio.open(dst_path, "w", **profile) as out:
        out.write(z_north_down, 1)
    return {
        "path_dependent": "netcdf",
        "shape": tuple(z_north_down.shape),
        "crs": src_crs,
        "bbox": (
            float(sub_lon[0] - dx / 2.0),
            float(sub_lat[0] - dy / 2.0),
            float(sub_lon[-1] + dx / 2.0),
            float(sub_lat[-1] + dy / 2.0),
        ),
        "src_var": src_var,
    }


def subset_dem_to_geotiff(
    src: str | Path,
    dst: str | Path,
    bbox: tuple[float, float, float, float],
    *,
    src_var: str | None = None,
    src_crs: str = "EPSG:4326",
) -> dict:
    """Clip a DEM to ``bbox`` and write a CF-tagged GeoTIFF.

    Parameters
    ----------
    src
        Source DEM path.
    dst
        Output GeoTIFF path.
    bbox
        Lon/lat bounds ``(minlon, minlat, maxlon, maxlat)``.
    src_var
        If given, force the lon/lat-NetCDF read path and use this
        variable name (e.g. ``"z"`` for SRTM15+ / GEBCO). If ``None``,
        try rasterio first and only fall back when the source has no
        embedded CRS.
    src_crs
        CRS to tag the output with when the lon/lat-NetCDF path is
        used. Ignored on the rasterio path. Defaults to ``EPSG:4326``.

    Returns
    -------
    dict
        Metadata dict with keys ``path_dependent`` (``"rasterio"`` or
        ``"netcdf"``), ``shape`` (rows, cols), ``crs``, and ``bbox``
        (the actual bbox written, which may be slightly larger than
        the requested one due to pixel-edge alignment).
    """
    src = Path(src).resolve()
    dst = Path(dst).resolve()
    dst.parent.mkdir(parents=True, exist_ok=True)
    if not src.exists():
        raise FileNotFoundError(src)

    if src_var is None:
        try:
            return _subset_via_rasterio(src, dst, bbox)
        except (rasterio.errors.RasterioIOError, ValueError):
            return _subset_via_netcdf(src, dst, bbox, src_var="z", src_crs=src_crs)
    return _subset_via_netcdf(src, dst, bbox, src_var=src_var, src_crs=src_crs)
