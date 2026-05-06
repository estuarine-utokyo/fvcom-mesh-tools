"""DEM sampling at mesh-node coordinates.

Mesh engines emit ``(points, cells)`` only; depth interpolation is
the caller's job (``fmesh-buildmesh``). :func:`at_points` samples a
DEM raster at the lon/lat rows of a points array and returns
*depth* in metres, positive-down (``fort.14`` convention): water
columns where the DEM elevation is ``-z`` get depth ``+z``; land
(elevation ``> 0``) becomes ``0`` m.
"""

from __future__ import annotations

from pathlib import Path

import numpy as np
import rasterio
from rasterio.warp import transform as warp_transform


def at_points(
    dem_path: Path,
    points: np.ndarray,
    *,
    method: str = "linear",
    fill_value: float = 0.0,
) -> np.ndarray:
    """Sample ``dem_path`` at ``(lon, lat)`` rows of ``points``.

    Parameters
    ----------
    dem_path
        DEM raster (must have a CRS).
    points
        ``(NP, 2)`` array of ``(lon, lat)`` coordinates in EPSG:4326.
    method
        ``"linear"`` (bilinear) or ``"nearest"``.
    fill_value
        Returned for points outside the raster footprint or where the
        DEM has nodata. Default 0 m (sea level).

    Returns
    -------
    depths
        ``(NP,)`` array of *depth* in metres (positive down). For points
        where the DEM is land (elevation > 0), depth is set to 0.
    """
    if method not in ("linear", "nearest"):
        raise ValueError(f"unknown method: {method!r}")

    with rasterio.open(dem_path) as r:
        # Read the band into memory (Tokyo Bay-class DEMs are ~1-20 MB).
        band = r.read(1).astype(np.float64)
        nodata = r.nodata
        if nodata is not None:
            band = np.where(band == nodata, np.nan, band)

        src_crs = r.crs
        if src_crs is None:
            raise ValueError(f"{dem_path}: DEM has no CRS")

        # Reproject query points if needed
        lons = np.asarray(points[:, 0], dtype=np.float64)
        lats = np.asarray(points[:, 1], dtype=np.float64)
        if src_crs.to_string() != "EPSG:4326":
            xs, ys = warp_transform(
                "EPSG:4326", src_crs,
                lons.tolist(), lats.tolist(),
            )
            xs = np.asarray(xs)
            ys = np.asarray(ys)
        else:
            xs, ys = lons, lats

        # Convert (x, y) world coords to (col, row) raster coords using
        # the inverse affine transform. col + 0.5 / row + 0.5 represent
        # pixel centres.
        inv = ~r.transform
        cols, rows = inv * (xs, ys)

    H, W = band.shape

    if method == "nearest":
        ci = np.clip(np.round(cols - 0.5).astype(int), 0, W - 1)
        ri = np.clip(np.round(rows - 0.5).astype(int), 0, H - 1)
        z = band[ri, ci]
    else:
        # Bilinear on pixel-centre grid: shift cols/rows by -0.5 so
        # that integer indices align with cell centres.
        cf = cols - 0.5
        rf = rows - 0.5
        c0 = np.clip(np.floor(cf).astype(int), 0, W - 2)
        r0 = np.clip(np.floor(rf).astype(int), 0, H - 2)
        dc = np.clip(cf - c0, 0.0, 1.0)
        dr = np.clip(rf - r0, 0.0, 1.0)
        z00 = band[r0, c0]
        z01 = band[r0, c0 + 1]
        z10 = band[r0 + 1, c0]
        z11 = band[r0 + 1, c0 + 1]
        z = (
            (1 - dr) * ((1 - dc) * z00 + dc * z01)
            + dr * ((1 - dc) * z10 + dc * z11)
        )

    out_of_bounds = (cols < 0) | (cols > W) | (rows < 0) | (rows > H)
    z = np.where(out_of_bounds | np.isnan(z), fill_value, z)

    # DEM convention: positive up (elevation). fort.14 wants positive
    # down (depth). Land (z > 0) becomes 0 m depth.
    depth = np.where(z >= 0.0, 0.0, -z)
    return depth.astype(np.float64)
