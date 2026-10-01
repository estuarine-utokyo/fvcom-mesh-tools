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


def test_masked_cells_are_no_data_and_the_next_source_fills_them(tmp_path, monkeypatch):
    """Review F1: a masked -9999 must not come back as a 9999 m depth."""
    import netCDF4

    lon = np.linspace(139.0, 141.0, 21)
    lat = np.linspace(34.0, 36.0, 21)
    with netCDF4.Dataset(tmp_path / "masked.nc", "w") as ds:
        ds.createDimension("lon", 21)
        ds.createDimension("lat", 21)
        ds.createVariable("lon", "f8", ("lon",))[:] = lon
        ds.createVariable("lat", "f8", ("lat",))[:] = lat
        v = ds.createVariable("z", "f4", ("lat", "lon"), fill_value=-9999.0)
        v[:] = np.ma.masked_all((21, 21))
    _grid_nc(tmp_path / "coarse.nc", lon, lat, -100.0 * np.ones((21, 21)))
    monkeypatch.setitem(sources.SOURCES, "masked", Grid("masked.nc", "z"))
    monkeypatch.setitem(sources.SOURCES, "coarse", Grid("coarse.nc", "z"))
    monkeypatch.setitem(sources.DATUM, "masked", "T.P.")
    monkeypatch.setitem(sources.DATUM, "coarse", "T.P.")
    d, w = sample(["masked", "coarse"], np.array([140.0]), np.array([35.0]), data_dir=tmp_path)
    assert d[0] == pytest.approx(100.0) and w[0] == 1


def test_m7001_degenerate_window_is_uncovered_not_an_error(tmp_path, monkeypatch):
    """Review F15: two distinct points (or collinear ones) cannot be triangulated."""
    import pandas as pd

    pd.DataFrame({"mark": ["N", "N", "N"], "lon": [139.0, 139.0, 139.2],
                  "lat": [35.0, 35.0, 35.2], "z_tp": [-10.0, -10.0, -20.0]}
                 ).to_parquet(tmp_path / "two.parquet")
    pd.DataFrame({"mark": ["N"] * 3, "lon": [139.0, 139.1, 139.2],
                  "lat": [35.0, 35.1, 35.2], "z_tp": [-10.0, -15.0, -20.0]}
                 ).to_parquet(tmp_path / "line.parquet")
    for f in ("two.parquet", "line.parquet"):
        out = M7001Points(f).depth(np.array([139.1]), np.array([35.1]), tmp_path)
        assert np.isnan(out[0])


def _cao_area(tmp_path, root_name, value):
    from pyproj import Transformer

    x, y = Transformer.from_crs(4326, 2451, always_xy=True).transform(139.8, 35.5)
    d = tmp_path / root_name / "cao" / "地形データ" / "地形データ_第09系"
    d.mkdir(parents=True)
    vals = [f"{value:8.2f}"] * 100
    (d / "depth_0030-01.dat").write_text(
        "\n".join("".join(vals[i:i + 10]) for i in range(0, 100, 10)) + "\n")
    return dict(zone="09", area="0030-01", h=30.0, x0=x - 150, y0=y - 150, nx=10, ny=10)


def test_cao_cache_is_per_data_root(tmp_path):
    """Review F13: a second data root must be read, not served from the first."""
    cao = CaoNested(rel="cao", zones={"09": 2451})
    a = _cao_area(tmp_path, "r1", 10.0)
    _cao_area(tmp_path, "r2", 20.0)
    cao.areas = lambda root: [a]
    assert cao.depth(np.array([139.8]), np.array([35.5]), tmp_path / "r1")[0] == 10.0
    assert cao.depth(np.array([139.8]), np.array([35.5]), tmp_path / "r2")[0] == 20.0


def test_cao_finer_grid_without_data_keeps_the_coarser_value(tmp_path):
    """Review F14: NaN in the finer grid must not erase the coarser depth."""
    cao = CaoNested(rel="cao", zones={"09": 2451})
    fine = _cao_area(tmp_path, "r", 0.0)
    coarse = dict(fine, area="0090-01", h=90.0, x0=fine["x0"] - 600, y0=fine["y0"] - 600,
                  nx=20, ny=20)
    d = tmp_path / "r" / "cao" / "地形データ" / "地形データ_第09系"
    vals = [f"{12.5:8.2f}"] * 400
    (d / "depth_0090-01.dat").write_text(
        "\n".join("".join(vals[i:i + 10]) for i in range(0, 400, 10)) + "\n")
    cao.areas = lambda root: [coarse, fine]
    f = d / "depth_0030-01.dat"
    st = f.stat()
    cao._cache[(str(f.resolve()), st.st_size, st.st_mtime_ns, 10, 10)] = np.full((10, 10), np.nan)
    assert cao.depth(np.array([139.8]), np.array([35.5]), tmp_path / "r")[0] == 12.5


def test_cao_provenance_lists_the_depth_files(tmp_path):
    """Review F16: the depth files, not only the area tables, are provenance."""
    _cao_area(tmp_path, "r", 5.0)
    (tmp_path / "r" / "cao" / "計算範囲設定").mkdir()
    (tmp_path / "r" / "cao" / "計算範囲設定" / "計算範囲設定_第09系.xls").write_bytes(b"x")
    names = [p.name for p in CaoNested(rel="cao").files(tmp_path / "r")]
    assert "depth_0030-01.dat" in names and "計算範囲設定_第09系.xls" in names


def test_a_masked_neighbour_with_zero_weight_keeps_the_sample():
    """Review round 3 F7: a sample on a valid node beside a masked cell was NaN."""
    from fvcom_mesh_tools.dem.sources import _bilinear

    gx, gy = np.array([0.0, 1.0]), np.array([0.0, 1.0])
    z = np.array([[-10.0, -20.0], [-30.0, np.nan]])
    out = _bilinear(gx, gy, z, np.array([0.0, 0.5, 0.25]), np.array([0.0, 0.0, 0.75]))
    assert out[0] == -10.0 and out[1] == -15.0 and np.isnan(out[2])
    assert np.isnan(_bilinear(gx, gy, z, np.array([2.0]), np.array([0.5]))[0])


def test_bilinear_keeps_the_query_shape():
    """Review round 4 F3: a (2, 2) query came back flat."""
    from fvcom_mesh_tools.dem.sources import _bilinear

    gx, gy = np.array([0.0, 1.0]), np.array([0.0, 1.0])
    q = np.array([[0.0, 0.5], [0.5, 1.0]])
    out = _bilinear(gx, gy, np.full((2, 2), 3.0), q, q)
    assert out.shape == (2, 2) and np.all(out == 3.0)


def test_cao_keeps_valid_centres_beside_no_data_and_on_the_last_column():
    """Review round 4 F9: a valid centre beside a NaN cell, and the last
    column of centres, came back NaN."""
    from pyproj import Transformer

    h, x0, y0 = 30.0, -30000.0, -90000.0
    area = dict(zone="09", area="0030-01", h=h, x0=x0, y0=y0, nx=2, ny=2)
    to_ll = Transformer.from_crs(2451, 4326, always_xy=True)

    def centre(i, j):   # row 0 is north
        return to_ll.transform(x0 + (i + 0.5) * h, y0 + 2 * h - (j + 0.5) * h)

    for grid, (i, j), want in (([[10.0, 20.0], [30.0, np.nan]], (0, 0), 10.0),
                               ([[7.0, 7.0], [7.0, 7.0]], (1, 0), 7.0)):
        cao = CaoNested(rel="cao", zones={"09": 2451})
        cao.areas = lambda root: [area]
        cao.grid = lambda *a, g=np.array(grid): g
        lon, lat = centre(i, j)
        assert cao.depth(np.array([lon]), np.array([lat]), Path("."))[0] == pytest.approx(want)


def test_cao_keeps_the_query_shape():
    """Review round 5 F9."""
    cao = CaoNested(rel="cao", zones={"09": 2451})
    cao.areas = lambda root: []
    q = np.full((2, 2), 139.8)
    assert cao.depth(q, q - 4.3, Path(".")).shape == (2, 2)


def test_m7001_depth_does_not_depend_on_the_other_query_points(tmp_path, monkeypatch):
    """Review round 7 F2: a centre inside four soundings 0.3 deg away was
    uncovered when asked alone and covered when asked with its neighbours."""
    import pandas as pd

    pd.DataFrame({"mark": ["N"] * 4, "lon": [139.7, 140.3, 139.7, 140.3],
                  "lat": [34.7, 34.7, 35.3, 35.3], "z_tp": [-10.0] * 4}
                 ).to_parquet(tmp_path / "m.parquet")
    monkeypatch.setitem(sources.SOURCES, "pts", M7001Points("m.parquet"))
    alone, _ = sample(["pts"], np.array([140.0]), np.array([35.0]), data_dir=tmp_path)
    batch, _ = sample(["pts"], np.array([140.0, 139.7, 140.3]), np.array([35.0, 34.7, 35.3]),
                      data_dir=tmp_path)
    assert alone[0] == pytest.approx(10.0) and batch[0] == alone[0]


def test_cao_reads_a_depth_file_again_after_it_changes(tmp_path):
    """Review round 8 F9: the cache served the old depths."""
    import os

    cao = CaoNested(rel="cao", zones={"09": 2451})
    a = _cao_area(tmp_path, "r", 10.0)
    cao.areas = lambda root: [a]
    assert cao.depth(np.array([139.8]), np.array([35.5]), tmp_path / "r")[0] == 10.0
    f = next((tmp_path / "r").rglob("depth_0030-01.dat"))
    st = f.stat()
    f.write_text(f.read_text().replace("   10.00", "   20.00"))
    os.utime(f, ns=(st.st_atime_ns, st.st_mtime_ns + 10**9))
    assert cao.depth(np.array([139.8]), np.array([35.5]), tmp_path / "r")[0] == 20.0


def test_a_grid_stored_lon_by_lat_is_not_transposed(tmp_path, monkeypatch):
    """Review round 8 F10."""
    import netCDF4

    with netCDF4.Dataset(tmp_path / "t.nc", "w") as ds:
        ds.createDimension("lon", 2)
        ds.createDimension("lat", 2)
        ds.createVariable("lon", "f8", ("lon",))[:] = [139.0, 139.1]
        ds.createVariable("lat", "f8", ("lat",))[:] = [35.0, 35.1]
        ds.createVariable("z", "f4", ("lon", "lat"))[:] = [[-10.0, -20.0], [-30.0, -40.0]]
    monkeypatch.setitem(sources.SOURCES, "t", Grid("t.nc", "z"))
    monkeypatch.setitem(sources.DATUM, "t", "T.P.")
    d, _ = sample(["t"], np.array([139.0]), np.array([35.1]), data_dir=tmp_path)
    assert d[0] == pytest.approx(20.0)
