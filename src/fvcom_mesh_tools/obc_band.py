"""Open-boundary band construction (generalised from the goto2023
sample anatomy, owner policy 2026-07-11).

Principles
----------
* Boundaries are instability-prone: the mesh near an open boundary
  must be LARGE, and interior gradation must never shrink it.
* The first interior row must form orderly "ladder" rows: a SMOOTH
  guide line parallel to the OBC, constrained into the mesh
  (pfix+egfix), gives every OBC node a near-perpendicular partner.

Measured constants (Tokyo Bay sample): the inner-line offset tracks
``K_OFFSET x local target mesh size`` (coast end ~1050 m where the
local size is ~840 m; deep end ~2100 m where it is ~1680 m, i.e.
K ~ 1.25), and the artificial-closure spacing tapers from the
junction value to the local coastal size (1901/1314/912/.../628 m).

This module is intentionally free of any oceanmesh import (license
policy): callers evaluate their sizing field themselves and pass
plain arrays.  Distances here are geodesic-approximate (local
``cos(lat)`` metric), adequate for band construction at bay scale.
"""

from __future__ import annotations

from typing import Any

import numpy as np

__all__ = ["build_obc_band", "corridor_targets", "apply_corridor"]

DEG_PER_M = 1.0 / 111e3


def _metric(ll: np.ndarray, cosw: float) -> np.ndarray:
    return np.column_stack([ll[:, 0] * cosw, ll[:, 1]]) * 111e3


def _unmetric(xy: np.ndarray, cosw: float) -> np.ndarray:
    return np.column_stack([xy[:, 0] / (111e3 * cosw),
                            xy[:, 1] / 111e3])


def build_obc_band(
    arc_ll: np.ndarray,
    h_arc_m: np.ndarray,
    *,
    k_offset: float = 1.25,
    skip_ends: int = 1,
    smooth_passes: int = 2,
    taper: str = "linear-ends",
) -> dict[str, Any]:
    """Constrained inner guide line parallel to an OBC arc.

    Parameters
    ----------
    arc_ll:
        ``(N, 2)`` lon/lat of the user-drawn open-boundary line, in
        along-boundary order. ``N >= 3``.
    h_arc_m:
        ``(N,)`` local TARGET MESH size (metres) at each arc node —
        the caller evaluates its sizing field there (including any
        field->mesh rendering factor, e.g. the DistMesh ~1.2x).
    k_offset:
        Inner-line offset as a multiple of ``h_arc_m``
        (sample-calibrated 1.25).
    skip_ends:
        Number of arc nodes at EACH end that get no inner partner
        (at a corner with an artificial closure the inward normal
        runs along the land line; end wedges mesh fine freely).
    smooth_passes:
        Moving-average passes over the offsets so a noisy sizing
        field cannot wobble the guide line.
    taper:
        ``"linear-ends"`` (default, the shape measured on the
        goto2023 sample): offsets interpolate linearly along the
        arc between ``k_offset*h`` at the two ENDS — a smooth
        monotone band. ``"local"``: offsets follow
        ``k_offset*h_arc_m`` per node (smoothed) — use when the
        size along the boundary is genuinely non-monotone.

    Returns
    -------
    dict with ``pfix`` ((N+M, 2) lonlat: arc then inner points),
    ``egfix`` ((K, 2) int indices into pfix), ``inner_ll``
    ((M, 2)), ``offsets_m`` ((N,)).
    """
    arc_ll = np.asarray(arc_ll, dtype=float)
    h_arc_m = np.asarray(h_arc_m, dtype=float)
    n = len(arc_ll)
    if n < 3:
        raise ValueError(
            f"OBC arc needs >= 3 nodes, got {n}. A 2-point line has "
            "no interior to ladder; densify the input line first.")
    if h_arc_m.shape != (n,):
        raise ValueError(
            f"h_arc_m must have shape ({n},), got {h_arc_m.shape}")
    if (h_arc_m <= 0).any():
        raise ValueError("h_arc_m must be positive")
    if 2 * skip_ends >= n - 1:
        raise ValueError(
            f"skip_ends={skip_ends} leaves no inner nodes for an "
            f"arc of {n} nodes")

    if taper == "linear-ends":
        off = np.linspace(k_offset * float(h_arc_m[0]),
                          k_offset * float(h_arc_m[-1]), n)
    elif taper == "local":
        off = k_offset * h_arc_m.copy()
        for _ in range(max(0, smooth_passes)):
            off[1:-1] = (off[:-2] + off[1:-1] + off[2:]) / 3.0
    else:
        raise ValueError(
            f"taper={taper!r} is not valid; choose 'linear-ends' "
            "(sample-measured default) or 'local'.")

    cosw = float(np.cos(np.deg2rad(arc_ll[:, 1].mean())))
    xy = _metric(arc_ll, cosw)
    t = np.gradient(xy, axis=0)
    t /= np.linalg.norm(t, axis=1)[:, None]
    # inward = left of the walk direction; the caller orients the
    # arc so the domain lies to its LEFT
    nrm = np.column_stack([-t[:, 1], t[:, 0]])
    inner_xy = xy + nrm * off[:, None]
    inner_ll = _unmetric(inner_xy, cosw)
    lo = skip_ends
    hi = n - skip_ends
    inner_ll = inner_ll[lo:hi]
    m = len(inner_ll)
    arc_seg = np.column_stack([np.arange(n - 1), np.arange(1, n)])
    inner_seg = (np.column_stack([np.arange(m - 1),
                                  np.arange(1, m)]) + n)
    return {
        "pfix": np.vstack([arc_ll, inner_ll]),
        "egfix": np.vstack([arc_seg, inner_seg]),
        "inner_ll": inner_ll,
        "offsets_m": off,
    }


def corridor_targets(
    arc_ll: np.ndarray,
    h_arc_m: np.ndarray,
    *,
    closure_ll: np.ndarray | None = None,
    h_closure_end_m: float | None = None,
    mesh_factor: float = 1.2,
    step_m: float = 100.0,
) -> tuple[np.ndarray, np.ndarray]:
    """Densified boundary-crossing polyline with per-point FIELD
    size targets (metres): the arc carries ``h_arc_m/mesh_factor``;
    an optional artificial closure line tapers from the junction
    value to ``h_closure_end_m/mesh_factor`` at its far (coastal)
    end. Feed the result to :func:`apply_corridor`.

    Returns ``(points_m, targets_m)`` in the local metric frame of
    ``arc_ll`` (same frame ``apply_corridor`` uses).
    """
    arc_ll = np.asarray(arc_ll, dtype=float)
    h_arc_m = np.asarray(h_arc_m, dtype=float)
    cosw = float(np.cos(np.deg2rad(arc_ll[:, 1].mean())))
    xy = _metric(arc_ll, cosw)
    pts: list[np.ndarray] = []
    tgt: list[float] = []
    for i in range(len(xy) - 1):
        a, b = xy[i], xy[i + 1]
        length = float(np.linalg.norm(b - a))
        for f in np.arange(0.0, 1.0, step_m / max(length, step_m)):
            pts.append(a * (1 - f) + b * f)
            tgt.append((h_arc_m[i] * (1 - f) + h_arc_m[i + 1] * f)
                       / mesh_factor)
    if closure_ll is not None:
        closure_ll = np.asarray(closure_ll, dtype=float)
        if h_closure_end_m is None:
            raise ValueError(
                "closure_ll given without h_closure_end_m: the "
                "coastal-end target size must be explicit")
        cxy = _metric(closure_ll, cosw)
        seglen = np.r_[0.0, np.cumsum(
            np.linalg.norm(np.diff(cxy, axis=0), axis=1))]
        total = float(seglen[-1])
        h0 = float(h_arc_m[-1]) / mesh_factor
        h1 = float(h_closure_end_m) / mesh_factor
        for i in range(len(cxy) - 1):
            a, b = cxy[i], cxy[i + 1]
            length = float(np.linalg.norm(b - a))
            for f in np.arange(0.0, 1.0001,
                               step_m / max(length, step_m)):
                s = (seglen[i] + f * length) / max(total, 1.0)
                pts.append(a * (1 - f) + b * f)
                tgt.append(h0 * (1 - s) + h1 * s)
    return np.asarray(pts), np.asarray(tgt)


def apply_corridor(
    lon_g: np.ndarray,
    lat_g: np.ndarray,
    values_deg: np.ndarray,
    points_m: np.ndarray,
    targets_m: np.ndarray,
    *,
    grade: float,
    arc_mean_lat: float,
) -> tuple[np.ndarray, int]:
    """Boundary-priority size override: raise the sizing lattice to
    the corridor target near the boundary crossing, tapering
    outward at ``grade``. Apply AFTER gradation so interior
    smoothing can never shrink the boundary band.

    Returns ``(new_values_deg, n_raised)``.
    """
    from scipy.spatial import cKDTree

    cosw = float(np.cos(np.deg2rad(arc_mean_lat)))
    q = np.column_stack([lon_g.ravel() * cosw * 111e3,
                         lat_g.ravel() * 111e3])
    d, i = cKDTree(points_m).query(q, workers=-1)
    t = targets_m[i]
    corr = np.maximum(t - grade * np.maximum(0.0, d - t), 0.0)
    corr = corr.reshape(values_deg.shape) * DEG_PER_M
    out = np.maximum(np.asarray(values_deg, dtype=float), corr)
    return out, int((corr > values_deg).sum())
