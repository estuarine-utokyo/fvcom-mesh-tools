"""Tests for ``dem.subset.to_geotiff`` and the ``fmesh-subset-dem`` CLI."""

from __future__ import annotations

from pathlib import Path

import numpy as np
import pytest

rasterio = pytest.importorskip("rasterio")
netCDF4 = pytest.importorskip("netCDF4")

from fvcom_mesh_tools.cli.subset_dem import main as cli_main  # noqa: E402
from fvcom_mesh_tools.dem.subset import to_geotiff  # noqa: E402


def _write_nc_lonlat_z(path: Path, lon, lat, z, var_name: str = "z") -> None:
    """Write a SRTM15+/GEBCO-style NetCDF: lon/lat 1-D, z 2-D, no CRS."""
    ds = netCDF4.Dataset(path, "w")
    try:
        ds.createDimension("lon", lon.size)
        ds.createDimension("lat", lat.size)
        v_lon = ds.createVariable("lon", "f8", ("lon",))
        v_lat = ds.createVariable("lat", "f8", ("lat",))
        v_z = ds.createVariable(var_name, "f4", ("lat", "lon"))
        v_lon[:] = lon
        v_lat[:] = lat
        v_z[:] = z
    finally:
        ds.close()


def _write_geotiff(path: Path, data, transform, crs="EPSG:4326") -> None:
    profile = {
        "driver": "GTiff",
        "height": data.shape[0],
        "width": data.shape[1],
        "count": 1,
        "dtype": data.dtype,
        "crs": crs,
        "transform": transform,
    }
    with rasterio.open(path, "w", **profile) as out:
        out.write(data, 1)


def test_subset_via_netcdf_roundtrip(tmp_path: Path) -> None:
    lon = np.linspace(135.0, 136.0, 21)  # 0.05 deg step
    lat = np.linspace(34.0, 35.0, 21)
    z = (
        np.arange(lat.size).reshape(-1, 1) * 100.0
        + np.arange(lon.size).reshape(1, -1)
    ).astype("f4")

    src = tmp_path / "src.nc"
    dst = tmp_path / "out.tif"
    _write_nc_lonlat_z(src, lon, lat, z)

    info = to_geotiff(
        src, dst, (135.20, 34.30, 135.60, 34.70), src_var="z"
    )

    assert info["path_dependent"] == "netcdf"
    assert info["crs"] == "EPSG:4326"
    assert dst.exists()

    with rasterio.open(dst) as out:
        assert out.crs.to_string() == "EPSG:4326"
        # Re-emitted GeoTIFF is north-down, so row 0 is the *highest* lat.
        # Sanity-check: first row's max value > last row's max value.
        arr = out.read(1)
        assert arr[0].max() > arr[-1].max()
        # Bounds cover at least the requested bbox.
        b = out.bounds
        assert b.left <= 135.20 + 1e-6
        assert b.right >= 135.60 - 1e-6
        assert b.bottom <= 34.30 + 1e-6
        assert b.top >= 34.70 - 1e-6


def test_subset_via_rasterio_roundtrip(tmp_path: Path) -> None:
    from rasterio.transform import from_origin

    # 0.05 deg pixels, north-down convention.
    width, height = 21, 21
    pixel = 0.05
    lon0 = 135.0  # left edge
    lat0 = 35.0   # top edge
    transform = from_origin(lon0, lat0, pixel, pixel)
    data = (
        np.arange(height).reshape(-1, 1) * 100.0
        + np.arange(width).reshape(1, -1)
    ).astype("f4")

    src = tmp_path / "src.tif"
    dst = tmp_path / "out.tif"
    _write_geotiff(src, data, transform)

    info = to_geotiff(src, dst, (135.20, 34.30, 135.60, 34.70))

    assert info["path_dependent"] == "rasterio"
    with rasterio.open(dst) as out:
        assert out.crs.to_string() == "EPSG:4326"
        b = out.bounds
        assert b.left <= 135.20 + 1e-6
        assert b.right >= 135.60 - 1e-6
        assert b.bottom <= 34.30 + 1e-6
        assert b.top >= 34.70 - 1e-6


def test_bbox_outside_raster_raises(tmp_path: Path) -> None:
    lon = np.linspace(135.0, 136.0, 11)
    lat = np.linspace(34.0, 35.0, 11)
    z = np.zeros((11, 11), dtype="f4")
    src = tmp_path / "src.nc"
    _write_nc_lonlat_z(src, lon, lat, z)

    with pytest.raises(ValueError, match="does not intersect"):
        to_geotiff(
            src, tmp_path / "x.tif", (140.0, 30.0, 141.0, 31.0), src_var="z"
        )


def test_missing_variable_raises(tmp_path: Path) -> None:
    lon = np.linspace(135.0, 136.0, 11)
    lat = np.linspace(34.0, 35.0, 11)
    z = np.zeros((11, 11), dtype="f4")
    src = tmp_path / "src.nc"
    _write_nc_lonlat_z(src, lon, lat, z, var_name="z")

    with pytest.raises(ValueError, match="not found"):
        to_geotiff(
            src, tmp_path / "x.tif", (135.2, 34.2, 135.8, 34.8),
            src_var="elevation",
        )


def test_cli_main_smoke(tmp_path: Path, capsys) -> None:
    lon = np.linspace(135.0, 136.0, 11)
    lat = np.linspace(34.0, 35.0, 11)
    z = np.zeros((11, 11), dtype="f4")
    src = tmp_path / "src.nc"
    dst = tmp_path / "out.tif"
    _write_nc_lonlat_z(src, lon, lat, z)

    rc = cli_main([
        str(src), str(dst),
        "--bbox", "135.2", "34.2", "135.8", "34.8",
        "--src-var", "z",
    ])
    assert rc == 0
    assert dst.exists()
    out = capsys.readouterr().out
    assert "path:   netcdf" in out
    assert "EPSG:4326" in out


def test_cli_rejects_bad_bbox(tmp_path: Path, capsys) -> None:
    src = tmp_path / "src.nc"
    src.write_bytes(b"")  # exists but unreadable
    rc = cli_main([
        str(src), str(tmp_path / "x.tif"),
        "--bbox", "136", "35", "135", "34",  # min > max
    ])
    assert rc == 2
    assert "minlon < maxlon" in capsys.readouterr().err
