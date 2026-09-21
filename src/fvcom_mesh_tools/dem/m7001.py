"""Build TB-FVCOM production bathymetry for an arbitrary Tokyo Bay mesh.

The goto2023 production depth file is
``TokyoBay_dep_m7001tp_rfac0p2_cap300.dat``: the JHA/JCG **M7001** survey on
the **Tokyo Peil** datum, floored at 3 m, r-factor smoothed to ``r <= 0.2``,
then capped at 300 m (TB-FVCOM ``input/goto2023/MESH.md``; the smoother is
``hydro/analysis/build_bathy_variants.py``).

Comparing two meshes only makes sense if their depths come from the same
source through the same steps, so this module applies that recipe to any mesh
instead of interpolating one mesh's finished node depths onto another.

Sources, in the order they are consulted:

1. ``M7001/TP/M7001_dem_tokyobay.nc`` -- the M7001 survey gridded at ~180 m on
   the T.P. datum.  Covers Tokyo Bay proper, which is where the survey is.
2. ``tokyo_bay/kanto_M7001_srtm_15s.nc`` -- M7001 blended into SRTM15 over the
   whole Kanto region at ~380 m, for the bay mouth and the shelf beyond, which
   the first product does not reach.

Both carry elevation positive up; depth here is positive down.
"""

from __future__ import annotations

import os
from pathlib import Path
from typing import Any

import numpy as np

__all__ = [
    "MIN_DEPTH_M",
    "MAX_DEPTH_M",
    "RMAX",
    "interpolate_m7001_tp",
    "node_edges",
    "production_depths",
    "rfactor_smooth",
]

#: The production recipe's three numbers (TB-FVCOM build_bathy_variants.py).
MIN_DEPTH_M = 3.0
MAX_DEPTH_M = 300.0
RMAX = 0.2

_FINE = "geodata/bathymetry/M7001/TP/M7001_dem_tokyobay.nc"
_WIDE = "geodata/bathymetry/tokyo_bay/kanto_M7001_srtm_15s.nc"


def _data_dir() -> Path:
    root = os.environ.get("DATA_DIR")
    if not root:
        raise RuntimeError("DATA_DIR is not set -- required for the M7001 products")
    return Path(root)


def _grid(path: Path, var: str) -> tuple[np.ndarray, np.ndarray, np.ndarray]:
    import netCDF4

    with netCDF4.Dataset(path) as ds:
        lat = np.asarray(ds["lat"][:], float)
        lon = np.asarray(ds["lon"][:], float)
        z = np.asarray(ds[var][:], float)
    if z.shape != (lat.size, lon.size):
        raise ValueError(f"{path.name}: expected (lat, lon) = "
                         f"{(lat.size, lon.size)}, got {z.shape}")
    return lat, lon, z


def interpolate_m7001_tp(lon: np.ndarray, lat: np.ndarray) -> tuple[np.ndarray, np.ndarray]:
    """Depth (positive down, metres, T.P. datum) at the given geographic points.

    Returns ``(depth, source)``, where ``source`` is 0 for the fine M7001 grid
    and 1 for the Kanto-wide blend.  Raises if any point is covered by neither.
    """
    from scipy.interpolate import RegularGridInterpolator

    lon = np.asarray(lon, float)
    lat = np.asarray(lat, float)
    pts = np.column_stack([lat, lon])
    depth = np.full(len(pts), np.nan)
    source = np.full(len(pts), -1, dtype=np.int8)
    root = _data_dir()
    for tag, (rel, var) in enumerate(((_FINE, "elevation"), (_WIDE, "z"))):
        path = root / rel
        if not path.exists():
            raise FileNotFoundError(path)
        glat, glon, z = _grid(path, var)
        g = RegularGridInterpolator((glat, glon), z, bounds_error=False, fill_value=np.nan)
        zi = g(pts)
        take = np.isnan(depth) & np.isfinite(zi)
        depth[take] = -zi[take]
        source[take] = tag
    if np.isnan(depth).any():
        raise ValueError(f"{int(np.isnan(depth).sum())} nodes outside every M7001 product")
    return depth, source


def node_edges(elements: np.ndarray) -> tuple[np.ndarray, np.ndarray]:
    """Unique undirected node-node edges of a triangulation (0-based)."""
    tri = np.asarray(elements, dtype=np.int64)
    e = np.sort(np.vstack([tri[:, [0, 1]], tri[:, [1, 2]], tri[:, [2, 0]]]), axis=1)
    e = np.unique(e, axis=0)
    return e[:, 0], e[:, 1]


def rfactor_smooth(
    h0: np.ndarray,
    ei: np.ndarray,
    ej: np.ndarray,
    *,
    rmax: float = RMAX,
    hmin: float = MIN_DEPTH_M,
    max_iter: int = 2000,
) -> tuple[np.ndarray, int, float]:
    """Beckmann-Haidvogel r-factor limiter, as TB-FVCOM applies it.

    For every edge with ``r = |hi - hj| / (hi + hj) > rmax``, move the deeper
    node up and the shallower down by equal amounts so ``r -> rmax``; average
    the moves each node receives and iterate.  Deterministic for a given
    ``rmax`` and edge set.  Returns ``(depth, iterations, final max r)``.
    """
    h = np.asarray(h0, float).copy()
    for it in range(int(max_iter)):
        hi, hj = h[ei], h[ej]
        r = np.abs(hi - hj) / (hi + hj)
        bad = r > rmax + 1e-9
        if not bad.any():
            return h, it, float(r.max())
        delta = np.where(bad, np.maximum((np.abs(hi - hj) - rmax * (hi + hj)) / 2.0, 0.0), 0.0)
        sgn = np.sign(hi - hj)
        add = np.zeros_like(h)
        cnt = np.zeros_like(h)
        np.add.at(add, ei, -sgn * delta)
        np.add.at(cnt, ei, bad.astype(float))
        np.add.at(add, ej, +sgn * delta)
        np.add.at(cnt, ej, bad.astype(float))
        h += add / np.where(cnt > 0, cnt, 1.0)
        h = np.maximum(h, hmin)
    hi, hj = h[ei], h[ej]
    return h, int(max_iter), float((np.abs(hi - hj) / (hi + hj)).max())


def production_depths(
    lon: np.ndarray,
    lat: np.ndarray,
    elements: np.ndarray,
    *,
    hmin: float = MIN_DEPTH_M,
    hmax: float = MAX_DEPTH_M,
    rmax: float = RMAX,
) -> tuple[np.ndarray, dict[str, Any]]:
    """The full ``m7001tp_rfac0p2_cap300`` recipe for one mesh.

    ``lon``/``lat`` are the node coordinates in EPSG:4326; ``elements`` is the
    0-based triangle table.  Returns ``(depth, report)``.
    """
    raw, source = interpolate_m7001_tp(lon, lat)
    floored = np.maximum(raw, hmin)
    ei, ej = node_edges(elements)
    r_before = float((np.abs(floored[ei] - floored[ej]) / (floored[ei] + floored[ej])).max())
    smoothed, iters, r_after = rfactor_smooth(floored, ei, ej, rmax=rmax, hmin=hmin)
    depth = np.minimum(smoothed, hmax)
    report = {
        "n_nodes": int(len(depth)),
        "n_edges": int(len(ei)),
        "from_fine_grid": int((source == 0).sum()),
        "from_kanto_blend": int((source == 1).sum()),
        "raw_min_m": float(raw.min()),
        "raw_max_m": float(raw.max()),
        "r_before": r_before,
        "r_after": r_after,
        "rfactor_iterations": iters,
        "depth_min_m": float(depth.min()),
        "depth_max_m": float(depth.max()),
        "n_capped": int((smoothed > hmax).sum()),
    }
    return depth, report
