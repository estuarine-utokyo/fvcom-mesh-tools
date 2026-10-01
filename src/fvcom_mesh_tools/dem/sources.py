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


def _bilinear(gx, gy, z, x, y) -> np.ndarray:
    """Bilinear interpolation on a regular grid; NaN only where it matters.

    A corner with no data (NaN) spoils a sample only when its weight is not
    zero: a point on a valid grid node or edge keeps its value beside a
    masked cell (review of the extend tools, round 3 F7). Outside the grid
    is NaN.
    """
    x = np.ravel(np.asarray(x, float))
    y = np.ravel(np.asarray(y, float))
    shape_out = np.shape(x)
    fi = np.interp(x, gx, np.arange(len(gx)), left=np.nan, right=np.nan)
    fj = np.interp(y, gy, np.arange(len(gy)), left=np.nan, right=np.nan)
    out = np.full(len(x), np.nan)
    ok = np.isfinite(fi) & np.isfinite(fj)
    i0 = np.clip(np.floor(fi[ok]).astype(int), 0, len(gx) - 2)
    j0 = np.clip(np.floor(fj[ok]).astype(int), 0, len(gy) - 2)
    ti, tj = fi[ok] - i0, fj[ok] - j0
    acc = np.zeros(ok.sum())
    bad = np.zeros(ok.sum(), bool)
    for di, dj, w in ((0, 0, (1 - ti) * (1 - tj)), (1, 0, ti * (1 - tj)),
                      (0, 1, (1 - ti) * tj), (1, 1, ti * tj)):
        v = z[j0 + dj, i0 + di]
        use = w > 0
        bad |= use & ~np.isfinite(v)
        acc += np.where(use & np.isfinite(v), w * np.nan_to_num(v), 0.0)
    out[np.flatnonzero(ok)] = np.where(bad, np.nan, acc)
    return out.reshape(shape_out)


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
            # masked cells (the file's _FillValue) are no data, not depths (review F1)
            z = np.ma.filled(np.ma.asarray(ds[self.var][j0:j1, i0:i1], float), np.nan)
        return -_bilinear(glon[i0:i1], glat[j0:j1], z, lon, lat)


@dataclass
class M7001Points:
    """M7001 soundings (N) and low-tide line (M) on T.P., linear between points."""

    rel: str = "geodata/bathymetry/M7001/TP/M7001_TP.parquet"

    def files(self, root: Path) -> list[Path]:
        return [root / self.rel]

    def depth(self, lon, lat, root: Path) -> np.ndarray:
        import pandas as pd
        from scipy.interpolate import LinearNDInterpolator

        # the window must reach past the query points to the soundings that
        # surround them, or the triangulation stops short of its own hull
        x0, x1, y0, y1 = _window(lon, lat, 0.25)
        df = pd.read_parquet(root / self.rel, columns=["mark", "lon", "lat", "z_tp"])
        df = df[df["mark"].isin(["N", "M"]) & df["lon"].between(x0, x1)
                & df["lat"].between(y0, y1) & np.isfinite(df["z_tp"])]
        out = np.full(np.shape(lon), np.nan)
        if len(df) < 3:
            return out
        # coincident points (a sounding repeated on two contour lines) make the
        # triangulation degenerate; keep one value per location
        df = df.drop_duplicates(["lon", "lat"])
        pts = df[["lon", "lat"]].to_numpy()
        # fewer than three distinct points, or all on a line, cannot be
        # triangulated: the window is uncovered, and the next source may cover
        # it (review F15)
        if len(pts) < 3 or np.linalg.matrix_rank(pts[1:] - pts[0], tol=1e-12) < 2:
            return out
        f = LinearNDInterpolator(pts, -df["z_tp"].to_numpy())
        return f(np.ravel(lon), np.ravel(lat)).reshape(np.shape(lon))


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
        # the data root and the shape are part of the key: one process may read
        # two roots (review F13)
        key = (str(Path(root).resolve()), zone, area, nx, ny)
        if key not in self._cache:
            hits = glob.glob(str(root / self.rel / "地形データ" / f"地形データ_第{zone}系"
                                 / "**" / f"depth_{area}.dat"), recursive=True)
            if len(hits) != 1:
                raise FileNotFoundError(f"depth_{area}.dat for zone {zone}: {hits}")
            # Fortran (10f8.2): every value is exactly eight characters, so the
            # file minus its line breaks is a flat array of 8-byte fields
            raw = Path(hits[0]).read_bytes().replace(b"\r", b"").replace(b"\n", b"")
            if len(raw) != 8 * nx * ny:
                raise ValueError(f"depth_{area}.dat: {len(raw)} bytes, expected {8 * nx * ny}")
            self._cache[key] = np.frombuffer(raw, dtype="S8").astype(float).reshape(ny, nx)
        return self._cache[key]

    def depth(self, lon, lat, root: Path) -> np.ndarray:
        from pyproj import Transformer

        lon = np.ravel(np.asarray(lon, float))
        lat = np.ravel(np.asarray(lat, float))
        out = np.full(lon.shape, np.nan)
        best = np.full(lon.shape, np.inf)
        for a in sorted(self.areas(root), key=lambda r: -r["h"]):   # coarse first
            tf = Transformer.from_crs(4326, self.zones[a["zone"]], always_xy=True)
            x, y = tf.transform(lon, lat)
            fi = (x - a["x0"]) / a["h"] - 0.5                  # cell centres
            fj = (a["y0"] + a["ny"] * a["h"] - y) / a["h"] - 0.5
            ok = ((fi >= 0) & (fi < a["nx"] - 1) & (fj >= 0) & (fj < a["ny"] - 1)
                  & (a["h"] <= best))
            if not ok.any():
                continue
            g = self.grid(root, a["zone"], a["area"], a["nx"], a["ny"])
            i0 = np.floor(fi[ok]).astype(int)
            j0 = np.floor(fj[ok]).astype(int)
            ti, tj = fi[ok] - i0, fj[ok] - j0
            val = (g[j0, i0] * (1 - ti) * (1 - tj) + g[j0, i0 + 1] * ti * (1 - tj)
                   + g[j0 + 1, i0] * (1 - ti) * tj + g[j0 + 1, i0 + 1] * ti * tj)
            # a finer grid replaces a coarser one only where it has a value
            # (review F14)
            idx = np.flatnonzero(ok)[np.isfinite(val)]
            out[idx] = val[np.isfinite(val)]
            best[idx] = a["h"]
        return out


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
