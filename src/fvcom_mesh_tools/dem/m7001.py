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


def interpolate_m7001_tp(lon: np.ndarray, lat: np.ndarray) -> tuple[np.ndarray, np.ndarray]:
    """Depth (positive down, metres, T.P. datum) at the given geographic points.

    Returns ``(depth, source)``, where ``source`` is 0 for the fine M7001 grid
    and 1 for the Kanto-wide blend.  Raises if any point is covered by neither.

    Each grid is read by ``dem.sources.Grid``, the loader reviewed for masked
    cells, dimension order, axes and zero-weight missing corners (review rounds
    29-31): a missing cell is NaN and the point goes to the next source.
    """
    from fvcom_mesh_tools.dem.sources import Grid

    if np.ma.is_masked(lon) or np.ma.is_masked(lat):
        raise ValueError("the query coordinates have masked values")
    lon = np.asarray(lon, float)
    lat = np.asarray(lat, float)
    if lon.ndim != 1 or lon.shape != lat.shape:
        raise ValueError(f"lon and lat must be matching (N,) arrays, not {lon.shape} and "
                         f"{lat.shape}")
    if not (np.isfinite(lon).all() and np.isfinite(lat).all()):
        raise ValueError("the query coordinates must be finite")
    depth = np.full(len(lon), np.nan)
    source = np.full(len(lon), -1, dtype=np.int8)
    root = _data_dir()
    for tag, (rel, var) in enumerate(((_FINE, "elevation"), (_WIDE, "z"))):
        if not lon.size:
            break
        d = Grid(rel, var).depth(lon, lat, root)
        take = np.isnan(depth) & np.isfinite(d)
        depth[take] = d[take]
        source[take] = tag
    if np.isnan(depth).any():
        raise ValueError(f"{int(np.isnan(depth).sum())} nodes outside every M7001 product")
    return depth, source


def node_edges(elements: np.ndarray) -> tuple[np.ndarray, np.ndarray]:
    """Unique undirected node-node edges of a triangulation (0-based)."""
    from fvcom_mesh_tools.io.fvcom_native import _indices

    # whole, unmasked, non-negative indices, checked before the cast (review
    # round 31 F8)
    raw = np.asarray(elements)
    n = int(np.nanmax(raw)) + 1 if raw.size and raw.dtype.kind in "iuf" \
        and np.isfinite(raw).all() else 1
    tri = _indices(elements, n, "elements", ndim=2)
    if tri.size and tri.shape[1] != 3:
        raise ValueError(f"elements must be (M, 3), not {tri.shape}")
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
    if np.ma.is_masked(ei) or np.ma.is_masked(ej):
        raise ValueError("ei or ej has masked values")
    ei, ej = np.asarray(ei), np.asarray(ej)
    if ei.size == 0 and ej.size == 0:
        ei = ej = np.empty(0, np.int64)
    # known, finite, positive depths, edge ends that are nodes, and sane
    # controls; the floor holds from the start, and depths whose sums would
    # overflow are refused (review round 30 F8, F9)
    if np.ma.is_masked(h0) or h.ndim != 1 or not (np.isfinite(h).all() and (h > 0).all()):
        raise ValueError("depths must be a known, finite, positive (N,) array")
    if ei.shape != ej.shape or ei.dtype.kind not in "iu" or ej.dtype.kind not in "iu" \
            or (ei.size and (min(ei.min(), ej.min()) < 0 or max(ei.max(), ej.max()) >= len(h))):
        raise ValueError("ei and ej must be matching whole node indices")
    if not (0 < rmax < 1 and np.isfinite(hmin) and hmin > 0 and int(max_iter) >= 1):
        raise ValueError(f"bad controls: rmax {rmax}, hmin {hmin}, max_iter {max_iter}")
    deg = np.bincount(np.r_[ei, ej], minlength=len(h)).max(initial=0)
    big = np.finfo(float).max / (4 * (int(deg) + 1))
    if max(h.max(initial=0.0), hmin) > big:
        raise ValueError(f"depths beyond {big:.1e} m cannot be limited without overflow")
    h = np.maximum(h, hmin)
    for it in range(int(max_iter)):
        hi, hj = h[ei], h[ej]
        r = np.abs(hi - hj) / (hi + hj)
        bad = r > rmax + 1e-9
        if not bad.any():
            return h, it, float(r.max(initial=0.0))
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
    r = float((np.abs(hi - hj) / (hi + hj)).max(initial=0.0))
    if not r <= rmax + 1e-9:         # not reached: a result that breaks the limit is no result
        raise ValueError(f"r-factor limit {rmax} not reached in {max_iter} iterations "
                         f"(r = {r:.6f})")
    return h, int(max_iter), r


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
    if not (np.isfinite([hmin, hmax]).all() and 0 < hmin <= hmax):
        raise ValueError(f"the depth bounds must be finite with 0 < hmin <= hmax, not "
                         f"[{hmin}, {hmax}]")
    raw, source = interpolate_m7001_tp(lon, lat)
    floored = np.maximum(raw, hmin)
    ei, ej = node_edges(elements)
    r_before = float((np.abs(floored[ei] - floored[ej]) / (floored[ei] + floored[ej])).max())
    smoothed, iters, r_after = rfactor_smooth(floored, ei, ej, rmax=rmax, hmin=hmin)
    depth = np.minimum(smoothed, hmax)
    if len(elements) and int(np.max(elements)) >= len(depth):
        raise ValueError("the elements name a node the coordinates do not have")
    if not (np.isfinite(depth).all() and (depth >= hmin).all() and (depth <= hmax).all()):
        raise ValueError("the production depths are not finite and within their bounds")
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
