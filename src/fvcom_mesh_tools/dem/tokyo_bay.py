"""The Tokyo Bay depth ladder: which product answers at a given point.

`DATA_INVENTORY.md` resolved the precedence before this module existed --
M7001 is the survey authority, the `tokyo_bay` 30 m grid fills the intertidal
gap and anywhere M7001 has no data, and the Kanto blend covers the bay mouth
and the shelf. This applies that order per node and says which rung answered,
because a depth whose provenance is unknown cannot be argued with.

Three things it is careful about, each because the alternative is a plausible
mistake:

**A finite grid value is not a sounding.** The M7001 product is gridded at
181 x 222 m and interpolates across its own gaps, so "M7001 covers this
point" says nothing about whether a survey point is near it. Measured inside
the four fishery polygons this repository refines, the grid is 100 % finite
(Banzu 99.9 %) while the soundings are 128-188 m apart. :func:`sounding_distance`
is therefore separate from :func:`sample`, and the driver reports both.

**Extrapolation is a rung, not a silence.** Where no product answers, the
nearest finite cell of whichever product has one is used and the distance is
returned, so the area that was invented can be drawn rather than counted
(owner, 2026-09-23).

**`mark == 'L'` rows are not soundings.** All 145,372 of them have
`depth_cd == 0` and no `z_tp`: they are the chart-datum zero line, and
counting them as observations would make a coastline look like bathymetry.

Depth here is positive down; every product stores elevation positive up.
"""

from __future__ import annotations

import os
from functools import lru_cache
from pathlib import Path
from typing import Any

import numpy as np

__all__ = [
    "LADDER",
    "EXTRAPOLATED",
    "sample",
    "sounding_distance",
]

#: The rungs, in the order they are consulted.  The index into this tuple is
#: what :func:`sample` returns as the provenance of each node.
LADDER = ("m7001", "grid30", "kanto")

#: The provenance of a node no product covered.
EXTRAPOLATED = -1

_PRODUCTS: dict[str, tuple[str, str]] = {
    # rung: (path relative to DATA_DIR, variable name)
    "m7001": ("geodata/bathymetry/M7001/TP/M7001_dem_tokyobay.nc", "elevation"),
    "grid30": ("geodata/bathymetry/tokyo_bay/depth_0030-11+12+13+14+15.nc", "Band1"),
    "kanto": ("geodata/bathymetry/tokyo_bay/kanto_M7001_srtm_15s.nc", "z"),
}

_SOUNDINGS = "geodata/bathymetry/M7001/TP/M7001_TP.parquet"

#: Metres per degree of latitude, and the longitude scale is this times
#: cos(lat).  Distances here are reported, never used to place a node, so a
#: spherical approximation about the query's own latitude is enough.
_M_PER_DEG = 111_000.0


def _data_dir() -> Path:
    root = os.environ.get("DATA_DIR")
    if not root:
        raise RuntimeError("DATA_DIR is not set -- required for the depth ladder")
    return Path(root)


@lru_cache(maxsize=None)
def _grid(rung: str) -> tuple[np.ndarray, np.ndarray, np.ndarray]:
    """``(lat, lon, elevation)`` for one rung, latitude ascending."""
    import netCDF4

    rel, var = _PRODUCTS[rung]
    path = _data_dir() / rel
    if not path.exists():
        raise FileNotFoundError(f"{rung}: {path}")
    with netCDF4.Dataset(path) as ds:
        lat = np.asarray(ds["lat"][:], dtype=float)
        lon = np.asarray(ds["lon"][:], dtype=float)
        z = np.ma.filled(np.asarray(ds[var][:], dtype=float), np.nan)
    if z.shape != (lat.size, lon.size):
        raise ValueError(f"{path.name}: expected {(lat.size, lon.size)}, got {z.shape}")
    # RegularGridInterpolator needs ascending axes and will not say so clearly.
    if lat[0] > lat[-1]:
        lat, z = lat[::-1], z[::-1, :]
    if lon[0] > lon[-1]:
        lon, z = lon[::-1], z[:, ::-1]
    return lat, lon, z


def _interp(rung: str, lon: np.ndarray, lat: np.ndarray) -> np.ndarray:
    from scipy.interpolate import RegularGridInterpolator

    glat, glon, z = _grid(rung)
    f = RegularGridInterpolator((glat, glon), z, bounds_error=False, fill_value=np.nan)
    return np.asarray(f(np.column_stack([lat, lon])), dtype=float)


def _nearest_finite(rung: str, lon: np.ndarray, lat: np.ndarray,
                    lat0: float) -> tuple[np.ndarray, np.ndarray]:
    """Nearest finite cell of one rung: its elevation and the distance in metres.

    The candidate set is restricted to a box around the queries and widened
    until it holds something, because the 30 m grid has 4.9 million cells and
    building a tree over all of them to answer a handful of points would cost
    more than the answer is worth.
    """
    from scipy.spatial import cKDTree

    glat, glon, z = _grid(rung)
    kx = _M_PER_DEG * float(np.cos(np.radians(lat0)))
    query = np.column_stack([lon * kx, lat * _M_PER_DEG])
    pad = 0.02
    while True:
        # Past a few degrees the box has stopped being an optimisation, so
        # take the whole product rather than keep doubling: a query that far
        # from the data is going to be answered from far away whatever the
        # search did, and the distance it comes back with says so.
        whole = pad > 5.0
        ilat = (np.arange(glat.size) if whole else
                np.flatnonzero((glat >= lat.min() - pad) & (glat <= lat.max() + pad)))
        ilon = (np.arange(glon.size) if whole else
                np.flatnonzero((glon >= lon.min() - pad) & (glon <= lon.max() + pad)))
        if ilat.size and ilon.size:
            sub = z[np.ix_(ilat, ilon)]
            ok = np.isfinite(sub)
            if ok.any():
                jj, ii = np.nonzero(ok)
                pts = np.column_stack([glon[ilon][ii] * kx, glat[ilat][jj] * _M_PER_DEG])
                d, k = cKDTree(pts).query(query)
                return sub[jj[k], ii[k]], d
        if whole:                       # this product has no finite cell at all
            return np.full(lon.shape, np.nan), np.full(lon.shape, np.inf)
        pad *= 2.0


def sample(lon, lat, *, extrapolate: bool = True
           ) -> tuple[np.ndarray, np.ndarray, np.ndarray]:
    """Depth (positive down, metres, T.P.) from the ladder.

    Returns ``(depth, rung, distance_m)``.  ``rung`` indexes :data:`LADDER`,
    or is :data:`EXTRAPOLATED` for a point no product covered.  ``distance_m``
    is zero for a point a product answered and the distance to the nearest
    real cell for an extrapolated one, so the caller can report the area that
    was invented rather than merely count it.

    ``extrapolate=False`` leaves uncovered points as NaN, which is what a test
    wants and what a caller who would rather refuse than guess wants.
    """
    lon = np.asarray(lon, dtype=float).ravel()
    lat = np.asarray(lat, dtype=float).ravel()
    if lon.shape != lat.shape:
        raise ValueError("lon and lat must have the same shape")
    if not lon.size:
        return (np.zeros(0), np.zeros(0, dtype=np.int8), np.zeros(0))
    if not (np.isfinite(lon).all() and np.isfinite(lat).all()):
        raise ValueError("the ladder was given non-finite coordinates")

    elev = np.full(lon.shape, np.nan)
    rung = np.full(lon.shape, EXTRAPOLATED, dtype=np.int8)
    for k, name in enumerate(LADDER):
        todo = ~np.isfinite(elev)
        if not todo.any():
            break
        v = _interp(name, lon[todo], lat[todo])
        idx = np.flatnonzero(todo)[np.isfinite(v)]
        elev[idx] = v[np.isfinite(v)]
        rung[idx] = k

    distance = np.zeros(lon.shape)
    missing = ~np.isfinite(elev)
    if missing.any() and extrapolate:
        lat0 = float(np.mean(lat[missing]))
        best_d = np.full(int(missing.sum()), np.inf)
        best_z = np.full(int(missing.sum()), np.nan)
        for name in LADDER:
            z, d = _nearest_finite(name, lon[missing], lat[missing], lat0)
            take = d < best_d
            best_d[take], best_z[take] = d[take], z[take]
        elev[missing] = best_z
        distance[missing] = best_d
    return -elev, rung, distance


def sounding_distance(lon, lat) -> np.ndarray:
    """Distance in metres to the nearest real M7001 sounding.

    This is what "where M7001 has nothing" means, measured.  A node whose
    nearest sounding is further away than the gridded product's own 181 m
    interval is carrying that product's interpolation, not a survey value,
    and no threshold has to be chosen in advance for the number to say so.

    ``mark == 'L'`` rows are excluded: they are the chart-datum zero line and
    carry no depth at all.
    """
    import pandas as pd
    from scipy.spatial import cKDTree

    lon = np.asarray(lon, dtype=float).ravel()
    lat = np.asarray(lat, dtype=float).ravel()
    if not lon.size:
        return np.zeros(0)
    path = _data_dir() / _SOUNDINGS
    if not path.exists():
        raise FileNotFoundError(path)
    df = pd.read_parquet(path, columns=["lon", "lat", "z_tp"])
    pad = 0.05
    while True:
        m = (np.isfinite(df["z_tp"].to_numpy())
             & (df["lon"].to_numpy() >= lon.min() - pad)
             & (df["lon"].to_numpy() <= lon.max() + pad)
             & (df["lat"].to_numpy() >= lat.min() - pad)
             & (df["lat"].to_numpy() <= lat.max() + pad))
        if m.any() or pad > 20.0:
            break
        pad *= 2.0
    if not m.any():
        return np.full(lon.shape, np.inf)
    kx = _M_PER_DEG * float(np.cos(np.radians(float(np.mean(lat)))))
    pts = np.column_stack([df["lon"].to_numpy()[m] * kx,
                           df["lat"].to_numpy()[m] * _M_PER_DEG])
    d, _ = cKDTree(pts).query(np.column_stack([lon * kx, lat * _M_PER_DEG]))
    return np.asarray(d, dtype=float)


def provenance_report(rung: np.ndarray, distance: np.ndarray) -> dict[str, Any]:
    """Which rung answered, as counts and fractions."""
    rung = np.asarray(rung)
    n = int(rung.size)
    out: dict[str, Any] = {"n_nodes": n}
    for k, name in enumerate(LADDER):
        out[f"n_from_{name}"] = int((rung == k).sum())
    ex = rung == EXTRAPOLATED
    out["n_extrapolated"] = int(ex.sum())
    out["extrapolated_fraction"] = float(ex.mean()) if n else 0.0
    out["extrapolated_distance_max_m"] = (
        float(np.asarray(distance)[ex].max()) if ex.any() else 0.0)
    out["extrapolated_distance_median_m"] = (
        float(np.median(np.asarray(distance)[ex])) if ex.any() else 0.0)
    return out
