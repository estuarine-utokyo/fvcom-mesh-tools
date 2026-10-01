"""Named bathymetry sources and the priority stack (dem/sources.py)."""

from pathlib import Path

import numpy as np
import pytest

from fvcom_mesh_tools.dem import sources
from fvcom_mesh_tools.dem.sources import CaoNested, Grid, M7001Points, sample


def _grid_nc(path: Path, lon, lat, z, var="z"):
    import netCDF4

    path.parent.mkdir(parents=True, exist_ok=True)
    with netCDF4.Dataset(path, "w") as ds:
        ds.createDimension("lon", len(lon))
        ds.createDimension("lat", len(lat))
        ds.createVariable("lon", "f8", ("lon",))[:] = lon
        ds.createVariable("lat", "f8", ("lat",))[:] = lat
        ds.createVariable(var, "f4", ("lat", "lon"))[:] = z


@pytest.fixture
def fake_sources(tmp_path, monkeypatch):
    """Two grids: a fine one over a small window, a coarse one everywhere."""
    lon = np.linspace(139.0, 141.0, 21)
    lat = np.linspace(34.0, 36.0, 21)
    _grid_nc(tmp_path / "coarse.nc", lon, lat, -100.0 * np.ones((21, 21)))
    flon = np.linspace(139.5, 139.7, 5)
    flat = np.linspace(35.0, 35.2, 5)
    _grid_nc(tmp_path / "fine.nc", flon, flat, -np.add.outer(flat, flon) * 0 - 7.0,
             var="elevation")
    monkeypatch.setitem(sources.SOURCES, "fine", Grid("fine.nc", "elevation"))
    monkeypatch.setitem(sources.SOURCES, "coarse", Grid("coarse.nc", "z"))
    monkeypatch.setitem(sources.DATUM, "fine", "T.P.")
    monkeypatch.setitem(sources.DATUM, "coarse", "MSL")
    return tmp_path


def test_the_first_source_covering_a_point_wins(fake_sources):
    lon = np.array([139.6, 140.5])
    lat = np.array([35.1, 35.5])
    d, w = sample(["fine", "coarse"], lon, lat, data_dir=fake_sources)
    assert np.allclose(d, [7.0, 100.0]) and w.tolist() == [0, 1]
    d, w = sample(["coarse", "fine"], lon, lat, data_dir=fake_sources)
    assert np.allclose(d, [100.0, 100.0]) and w.tolist() == [0, 0]


def test_a_point_no_source_covers_is_nan(fake_sources):
    d, w = sample(["fine"], np.array([140.5]), np.array([35.5]), data_dir=fake_sources)
    assert np.isnan(d[0]) and w[0] == -1


def test_unknown_or_repeated_names_are_refused(fake_sources):
    with pytest.raises(ValueError, match="unknown"):
        sample(["nope"], np.array([140.0]), np.array([35.0]), data_dir=fake_sources)
    with pytest.raises(ValueError, match="twice"):
        sample(["fine", "fine"], np.array([140.0]), np.array([35.0]), data_dir=fake_sources)
    with pytest.raises(ValueError, match="at least one"):
        sample([], np.array([140.0]), np.array([35.0]), data_dir=fake_sources)


def test_a_zipped_grid_is_refused_with_a_reason(tmp_path, monkeypatch):
    (tmp_path / "g.nc").write_bytes(b"PK\x03\x04 not netcdf")
    monkeypatch.setitem(sources.SOURCES, "zipped", Grid("g.nc", "z"))
    with pytest.raises(ValueError, match="zip archive"):
        sample(["zipped"], np.array([140.0]), np.array([35.0]), data_dir=tmp_path)


def test_m7001_points_interpolate_on_tp_and_ignore_the_coastline(tmp_path, monkeypatch):
    import pandas as pd

    pts = pd.DataFrame({
        "mark": ["N", "N", "N", "M", "L"],
        "lon": [139.0, 139.2, 139.0, 139.2, 139.1],
        "lat": [35.0, 35.0, 35.2, 35.2, 35.1],
        "z_tp": [-10.0, -20.0, -10.0, -20.0, np.nan],
    })
    pts.to_parquet(tmp_path / "m.parquet")
    monkeypatch.setitem(sources.SOURCES, "pts", M7001Points("m.parquet"))
    d, w = sample(["pts"], np.array([139.1, 139.5]), np.array([35.1, 35.1]),
                  data_dir=tmp_path)
    assert d[0] == pytest.approx(15.0) and np.isnan(d[1])


def test_cao_nested_takes_the_finest_area_and_reads_the_fixed_width_file(tmp_path):
    """Two areas in zone IX: the fine one inside the coarse one."""
    from pyproj import Transformer

    cao = CaoNested(rel="cao", zones={"09": 2451})
    to_xy = Transformer.from_crs(4326, 2451, always_xy=True)
    x, y = to_xy.transform(139.8, 35.5)
    coarse = dict(zone="09", area="0090-01", h=90.0, x0=x - 900, y0=y - 900, nx=20, ny=20)
    fine = dict(zone="09", area="0030-01", h=30.0, x0=x - 150, y0=y - 150, nx=10, ny=10)
    cao.areas = lambda root: [coarse, fine]
    d = tmp_path / "cao" / "地形データ" / "地形データ_第09系"
    d.mkdir(parents=True)

    def write(area, value, nx, ny):
        vals = [f"{value:8.2f}"] * (nx * ny)
        lines = ["".join(vals[i:i + 10]) for i in range(0, len(vals), 10)]
        (d / f"depth_{area}.dat").write_text("\n".join(lines) + "\n")

    write("0090-01", 12.5, 20, 20)
    write("0030-01", -3.25, 10, 10)           # land, negative
    far_lon, far_lat = Transformer.from_crs(2451, 4326, always_xy=True).transform(x + 600, y)
    out = cao.depth(np.array([139.8, far_lon]), np.array([35.5, far_lat]), tmp_path)
    assert out.tolist() == pytest.approx([-3.25, 12.5])


def test_cao_refuses_a_file_of_the_wrong_size(tmp_path):
    cao = CaoNested(rel="cao", zones={"09": 2451})
    d = tmp_path / "cao" / "地形データ" / "地形データ_第09系"
    d.mkdir(parents=True)
    (d / "depth_0010-01.dat").write_text("    1.00    2.00\n")
    with pytest.raises(ValueError, match="bytes, expected"):
        cao.grid(tmp_path, "09", "0010-01", 3, 1)


def test_points_from_a_source_not_on_tp_are_named(fake_sources):
    lon = np.array([139.6, 140.5, 140.6])
    lat = np.array([35.1, 35.5, 35.6])
    with pytest.warns(UserWarning, match=r"2 point\(s\).*not on T\.P\..*coarse"):
        d, w = sample(["fine", "coarse"], lon, lat, data_dir=fake_sources)
    assert sources.non_tp_count(["fine", "coarse"], w) == (2, ["coarse"])
    assert sources.non_tp_count(["fine", "coarse"], np.array([0, 0])) == (0, [])


def test_every_registered_source_has_a_datum():
    assert set(sources.DATUM) >= set(sources.SOURCES)
    assert sources.DATUM["m7001"] == "T.P." and sources.DATUM["srtm15plus"] == "MSL"
