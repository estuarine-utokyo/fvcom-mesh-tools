"""Named bathymetry sources, and a priority stack over them.

A recipe names the sources it wants, in order; every point takes its depth
from the first source that covers it.  Depth is **positive down**, metres;
land comes out negative.  Every source reads from ``$DATA_DIR`` only.

========================  ==================================================
name                      what it is
========================  ==================================================
``cao_shutochokka_2025``  Cabinet Office tsunami topography (首都直下地震モデル
                          検討会, 令和7年～): nested 10-2430 m grids, T.P.
                          The finest grid covering a point is used.
``m7001``                 JHA M7001 (Southern Kanto) soundings and low-tide
                          line on T.P., interpolated linearly between the
                          points (``M7001/TP/M7001_TP.parquet``).
``m7001_tokyobay``        M7001 gridded at ~180 m, Tokyo Bay only, T.P.
``srtm15_kanto``          SRTM15+ V2.6, Kanto window (15 arc-seconds).
``srtm15plus``            SRTM15+ V2.6, global (15 arc-seconds).
``gebco_2024``            GEBCO_2024, global (15 arc-seconds).
========================  ==================================================

The global grids are referred to mean sea level, which differs from T.P. by
well under the grids' own error (tens of metres near the coast); they are
there to fill what the survey products do not cover.

The Cabinet Office grids may not be redistributed as they are (their licence,
article 4): this module reads them in place and hands out interpolated
values only.
"""

from __future__ import annotations

import glob
import os
from dataclasses import dataclass, field
from pathlib import Path
from typing import Any

import numpy as np

__all__ = ["DATUM", "SOURCES", "non_tp_count", "sample", "source_files"]


def _data_dir(data_dir=None) -> Path:
    root = data_dir or os.environ.get("DATA_DIR")
    if not root:
        raise RuntimeError("DATA_DIR is not set -- the bathymetry sources live there")
    return Path(root)


def _window(lon: np.ndarray, lat: np.ndarray, pad: float):
    return (float(np.nanmin(lon)) - pad, float(np.nanmax(lon)) + pad,
            float(np.nanmin(lat)) - pad, float(np.nanmax(lat)) + pad)


def _snap(f) -> np.ndarray:
    """Fractional grid indices, put on a grid line when within 1e-9 of it."""
    f = np.asarray(f, float)
    with np.errstate(invalid="ignore"):
        return np.where(np.abs(f - np.round(f)) < 1e-9, np.round(f), f)


def _bilinear_index(z, fi, fj) -> np.ndarray:
    """Bilinear interpolation of ``z[j, i]`` at fractional indices.

    A corner with no data (NaN) spoils a sample only when its weight is not
    zero: a point on a valid grid node or edge keeps its value beside a
    masked cell (review of the extend tools, round 3 F7). Indices outside
    ``[0, n-1]`` or not finite give NaN; the last row and column are inside.
    An index within 1e-9 of a grid line is put on it, so that a projection's
    round-off does not give a no-data corner a weight of 1e-20 (round 4 F9).
    """
    fi = _snap(fi)
    fj = _snap(fj)
    ny, nx = z.shape
    out = np.full(fi.shape, np.nan)
    ok = (np.isfinite(fi) & np.isfinite(fj) & (fi >= 0) & (fi <= nx - 1)
          & (fj >= 0) & (fj <= ny - 1))
    if nx < 2 or ny < 2 or not ok.any():
        return out
    i0 = np.clip(np.floor(fi[ok]).astype(int), 0, nx - 2)
    j0 = np.clip(np.floor(fj[ok]).astype(int), 0, ny - 2)
    ti, tj = fi[ok] - i0, fj[ok] - j0
    acc = np.zeros(ok.sum())
    bad = np.zeros(ok.sum(), bool)
    for di, dj, w in ((0, 0, (1 - ti) * (1 - tj)), (1, 0, ti * (1 - tj)),
                      (0, 1, (1 - ti) * tj), (1, 1, ti * tj)):
        v = z[j0 + dj, i0 + di]
        use = w > 0
        bad |= use & ~np.isfinite(v)
        acc += np.where(use & np.isfinite(v), w * np.nan_to_num(v), 0.0)
    out[ok] = np.where(bad, np.nan, acc)
    return out


def _bilinear(gx, gy, z, x, y) -> np.ndarray:
    """Bilinear interpolation on a regular grid; outside the grid is NaN.

    The result has the shape of the query ``x`` (review round 4 F3).
    """
    x = np.asarray(x, float)
    y = np.asarray(y, float)
    if x.shape != y.shape:
        raise ValueError(f"query shapes differ: {x.shape} and {y.shape}")
    fi = np.interp(np.ravel(x), gx, np.arange(len(gx)), left=np.nan, right=np.nan)
    fj = np.interp(np.ravel(y), gy, np.arange(len(gy)), left=np.nan, right=np.nan)
    return _bilinear_index(z, fi, fj).reshape(x.shape)


@dataclass
class Grid:
    """A regular lon/lat grid; elevation positive up in ``var``."""

    rel: str
    var: str

    def files(self, root: Path) -> list[Path]:
        return [root / self.rel]

    def depth(self, lon, lat, root: Path) -> np.ndarray:
        import netCDF4

        path = root / self.rel
        if not path.exists():
            raise FileNotFoundError(path)
        with open(path, "rb") as fh:
            if fh.read(2) == b"PK":
                raise ValueError(f"{path} is a zip archive, not netCDF -- unpack it first")
        x0, x1, y0, y1 = _window(lon, lat, 0.02)
        with netCDF4.Dataset(path) as ds:
            glon = np.asarray(ds["lon"][:], float)
            glat = np.asarray(ds["lat"][:], float)
            if not (np.all(np.diff(glon) > 0) and np.all(np.diff(glat) > 0)):
                raise ValueError(f"{path}: lon and lat must increase")
            i0 = max(int(np.searchsorted(glon, x0)) - 1, 0)
            i1 = min(int(np.searchsorted(glon, x1)) + 1, glon.size)
            j0 = max(int(np.searchsorted(glat, y0)) - 1, 0)
            j1 = min(int(np.searchsorted(glat, y1)) + 1, glat.size)
            out = np.full(np.shape(lon), np.nan)
            if i1 - i0 < 2 or j1 - j0 < 2:
                return out
            # the variable's dimensions are those of lat and lon, in either
            # order; equal axis lengths would hide a transposed grid (review
            # round 8 F10)
            lon_dim, lat_dim = ds["lon"].dimensions, ds["lat"].dimensions
            var_dims = ds[self.var].dimensions
            if len(lon_dim) != 1 or len(lat_dim) != 1:
                raise ValueError(f"{path}: lon and lat must be one-dimensional")
            if var_dims == lat_dim + lon_dim:
                block = ds[self.var][j0:j1, i0:i1]
            elif var_dims == lon_dim + lat_dim:
                block = ds[self.var][i0:i1, j0:j1].T
            else:
                raise ValueError(f"{path}: {self.var} has dimensions {var_dims}, "
                                 f"not ({lat_dim[0]}, {lon_dim[0]})")
            # masked cells (the file's _FillValue) are no data, not depths (review F1)
            z = np.ma.filled(np.ma.asarray(block, float), np.nan)
        return -_bilinear(glon[i0:i1], glat[j0:j1], z, lon, lat)


@dataclass
class M7001Points:
    """M7001 soundings (N) and low-tide line (M) on T.P., linear between points."""

    rel: str = "geodata/bathymetry/M7001/TP/M7001_TP.parquet"
    _cache: dict = field(default_factory=dict, repr=False)

    def files(self, root: Path) -> list[Path]:
        return [root / self.rel]

    def interpolator(self, root: Path):
        """One triangulation of the whole dataset, built once per file.

        A triangulation of the soundings near the query points made a depth,
        and whether a point was covered at all, depend on which other points
        were asked in the same call (review of the extend tools, round 7 F2).
        None when the data cannot be triangulated.
        """
        import pandas as pd
        from scipy.interpolate import LinearNDInterpolator

        path = (root / self.rel).resolve()
        st = path.stat()
        key = (str(path), st.st_size, st.st_mtime_ns)
        if key not in self._cache:
            df = pd.read_parquet(path, columns=["mark", "lon", "lat", "z_tp"])
            df = df[df["mark"].isin(["N", "M"]) & np.isfinite(df["z_tp"])
                    & np.isfinite(df["lon"]) & np.isfinite(df["lat"])]
            # coincident points (a sounding repeated on two contour lines) make
            # the triangulation degenerate; keep one value per location
            df = df.drop_duplicates(["lon", "lat"])
            pts = df[["lon", "lat"]].to_numpy()
            # fewer than three distinct points, or all on a line, cannot be
            # triangulated: nothing is covered, and the next source may cover
            # it (review F15)
            if len(pts) < 3 or np.linalg.matrix_rank(pts[1:] - pts[0], tol=1e-12) < 2:
                f = None
            else:
                f = LinearNDInterpolator(pts, -df["z_tp"].to_numpy())
            self._cache.clear()
            self._cache[key] = f
        return self._cache[key]

    def depth(self, lon, lat, root: Path) -> np.ndarray:
        lon = np.asarray(lon, float)
        lat = np.asarray(lat, float)
        f = self.interpolator(root)
        if f is None or lon.size == 0:
            return np.full(lon.shape, np.nan)
        return f(np.ravel(lon), np.ravel(lat)).reshape(lon.shape)


@dataclass
class CaoNested:
    """Cabinet Office nested grids: the finest one covering each point wins."""

    rel: str = "geodata/bathymetry/CAO_shutochokka_2025/extracted"
    zones: dict = field(default_factory=lambda: {"08": 2450, "09": 2451, "14": 2456})
    _cache: dict = field(default_factory=dict, repr=False)

    def files(self, root: Path) -> list[Path]:
        """The area tables and every depth file (review F16: the depths were missing)."""
        base = root / self.rel
        tables = glob.glob(str(base / "計算範囲設定" / "*.xls"))
        depths = glob.glob(str(base / "地形データ" / "**" / "depth_*.dat"), recursive=True)
        return sorted(Path(p) for p in tables + depths)

    def areas(self, root: Path):
        """Every nested area: zone, name, cell size, SW corner (m) and counts.

        In the Cabinet Office tables X is the west-east coordinate and Y the
        south-north one (the north-west corner shares X with the south-west).
        """
        import pandas as pd

        rows = []
        for zone in self.zones:
            path = root / self.rel / "計算範囲設定" / f"計算範囲設定_第{zone}系.xls"
            if not path.exists():
                raise FileNotFoundError(path)
            book = pd.ExcelFile(path)
            for sheet in book.sheet_names:
                if not sheet.endswith("m"):
                    continue
                table = book.parse(sheet, header=None)
                for _, r in table.iterrows():
                    name = r[1]
                    if not (isinstance(name, str) and "-" in name and name[:4].isdigit()):
                        continue
                    rows.append(dict(zone=zone, area=name, h=float(r[2]), x0=float(r[3]),
                                     y0=float(r[4]), nx=int(r[9]), ny=int(r[10])))
        return rows

    def grid(self, root: Path, zone: str, area: str, nx: int, ny: int) -> np.ndarray:
        """One area's depths, row 0 = north; parsed once per run."""
        # the file itself (resolved path, size, modification time) and the
        # shape are the key: one process may read two roots (review F13), and
        # a file changed since must be read again (review round 8 F9)
        hits = glob.glob(str(root / self.rel / "地形データ" / f"地形データ_第{zone}系"
                             / "**" / f"depth_{area}.dat"), recursive=True)
        if len(hits) != 1:
            raise FileNotFoundError(f"depth_{area}.dat for zone {zone}: {hits}")
        st = Path(hits[0]).stat()
        key = (str(Path(hits[0]).resolve()), st.st_size, st.st_mtime_ns, nx, ny)
        if key not in self._cache:
            # Fortran (10f8.2): every value is exactly eight characters, so the
            # file minus its line breaks is a flat array of 8-byte fields
            raw = Path(hits[0]).read_bytes().replace(b"\r", b"").replace(b"\n", b"")
            if len(raw) != 8 * nx * ny:
                raise ValueError(f"depth_{area}.dat: {len(raw)} bytes, expected {8 * nx * ny}")
            self._cache[key] = np.frombuffer(raw, dtype="S8").astype(float).reshape(ny, nx)
        return self._cache[key]

    def depth(self, lon, lat, root: Path) -> np.ndarray:
        from pyproj import Transformer

        lon = np.asarray(lon, float)
        lat = np.asarray(lat, float)
        if lon.shape != lat.shape:
            raise ValueError(f"query shapes differ: {lon.shape} and {lat.shape}")
        shape = lon.shape                  # given back (review round 5 F9)
        lon, lat = np.ravel(lon), np.ravel(lat)
        out = np.full(lon.shape, np.nan)
        best = np.full(lon.shape, np.inf)
        for a in sorted(self.areas(root), key=lambda r: -r["h"]):   # coarse first
            tf = Transformer.from_crs(4326, self.zones[a["zone"]], always_xy=True)
            x, y = tf.transform(lon, lat)
            fi = _snap((x - a["x0"]) / a["h"] - 0.5)           # cell centres
            fj = _snap((a["y0"] + a["ny"] * a["h"] - y) / a["h"] - 0.5)
            # the last row and column of centres are inside (review round 4 F9)
            ok = ((fi >= 0) & (fi <= a["nx"] - 1) & (fj >= 0) & (fj <= a["ny"] - 1)
                  & (a["h"] <= best))
            if not ok.any():
                continue
            g = self.grid(root, a["zone"], a["area"], a["nx"], a["ny"])
            # zero-weight NaN corners do not spoil a valid centre (round 4 F9)
            val = _bilinear_index(g, fi[ok], fj[ok])
            # a finer grid replaces a coarser one only where it has a value
            # (review F14)
            idx = np.flatnonzero(ok)[np.isfinite(val)]
            out[idx] = val[np.isfinite(val)]
            best[idx] = a["h"]
        return out.reshape(shape)


SOURCES: dict[str, Any] = {
    "cao_shutochokka_2025": CaoNested(),
    "m7001": M7001Points(),
    "m7001_tokyobay": Grid("geodata/bathymetry/M7001/TP/M7001_dem_tokyobay.nc", "elevation"),
    "srtm15_kanto": Grid("geodata/bathymetry/tokyo_bay/SRTM15_kanto_15s.nc", "z"),
    "srtm15plus": Grid("geodata/bathymetry/SRTM15plus/SRTM15+.nc", "z"),
    "gebco_2024": Grid("geodata/bathymetry/GEBCO/GEBCO_2024.nc", "elevation"),
}

#: Vertical datum of each source. Depths for a model on T.P. should come from
#: T.P. sources; the global grids are on mean sea level (tens of cm off T.P.,
#: below their own error offshore but not on tidal flats).
DATUM: dict[str, str] = {
    "cao_shutochokka_2025": "T.P.", "m7001": "T.P.", "m7001_tokyobay": "T.P.",
    "srtm15_kanto": "MSL", "srtm15plus": "MSL", "gebco_2024": "MSL",
}
if set(DATUM) != set(SOURCES):
    raise RuntimeError("every bathymetry source needs a datum in DATUM")


def _check(names) -> list[str]:
    names = list(names)
    if not names:
        raise ValueError("name at least one bathymetry source")
    unknown = [n for n in names if n not in SOURCES]
    if unknown:
        raise ValueError(f"unknown bathymetry source(s) {unknown}; known: {sorted(SOURCES)}")
    if len(set(names)) != len(names):
        raise ValueError(f"a bathymetry source is named twice: {names}")
    return names


def source_files(names, data_dir=None) -> dict[str, list[Path]]:
    """The files each named source reads, for provenance."""
    root = _data_dir(data_dir)
    return {n: SOURCES[n].files(root) for n in _check(names)}


def sample(names, lon, lat, data_dir=None) -> tuple[np.ndarray, np.ndarray]:
    """Depth (positive down, m) at the points, from the first source covering each.

    Returns ``(depth, which)``: ``which`` is the index into ``names`` of the
    source used, -1 where none covers the point (depth NaN there). A warning
    names the points taken from a source not on T.P. (:data:`DATUM`).
    """
    names = _check(names)
    root = _data_dir(data_dir)
    lon = np.asarray(lon, float)
    lat = np.asarray(lat, float)
    depth = np.full(lon.shape, np.nan)
    which = np.full(lon.shape, -1, dtype=np.int16)
    for k, name in enumerate(names):
        todo = np.isnan(depth)
        if not todo.any():
            break
        d = SOURCES[name].depth(lon[todo], lat[todo], root)
        got = np.isfinite(d)
        idx = np.flatnonzero(todo.ravel())[got.ravel()]
        depth.ravel()[idx] = d.ravel()[got.ravel()]
        which.ravel()[idx] = k
    n_off, off = non_tp_count(names, which)
    if n_off:
        import warnings

        warnings.warn(f"{n_off} point(s) take their depth from a source not on T.P. "
                      f"({', '.join(off)})", stacklevel=2)
    return depth, which


def non_tp_count(names, which) -> tuple[int, list[str]]:
    """How many points (``which`` from :func:`sample`) came from non-T.P. sources, and which."""
    names = _check(names)
    which = np.asarray(which)
    off = [n for k, n in enumerate(names)
           if DATUM.get(n, "unknown") != "T.P." and (which == k).any()]
    n = int(sum((which == names.index(o)).sum() for o in off))
    return n, off
