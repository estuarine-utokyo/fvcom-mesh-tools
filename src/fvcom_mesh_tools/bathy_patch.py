"""Depths for a refined patch, taken from a survey product rather than the base.

`refine.depths_from_base` is the other option and stays the default: the base
mesh's own field, evaluated at the new nodes. This module is what the ``hires``
branch uses instead -- the Tokyo Bay ladder (:mod:`fvcom_mesh_tools.dem.tokyo_bay`)
sampled at the delivered coordinates, ramped into the base field across the
transition so that the frozen depths are met exactly.

**The ramp is not an r-factor argument.** An earlier revision of this design
claimed that a continuous blend makes the seam safe by construction. It does
not: continuity bounds no edge difference. A source-base difference ``D``
spread over a width ``W`` has slope ``D/W``, so an edge of length ``h`` at
depth ``H`` carries about ``r = (D/W) * h / (2H)`` -- and a 297 m difference
over 330 m puts ``r = 0.818`` on the first edge inside the interface.
:func:`feasibility_margin` computes that before any meshing happens, and the
minimum depth and the r-factor smoothing are a later step.

**The weight is geometric, not the size field.** ``patch.patch_sizing``
divides by ``distmesh_scale`` and ramps to the base mesh's own local size, so
there is no ``(ambient, target)`` pair that inverts it: reading a weight out of
it gave 1.014 in the core and 0.180 at the nominal outer edge. And the weight
has to vanish at the interface the cut ACTUALLY has -- ``select_patch`` takes
whole faces by centroid and then grows the selection to repair pinches, so the
analytic buffer is not it.
"""

from __future__ import annotations

from typing import Any

import numpy as np

__all__ = [
    "blend_weights",
    "edge_slopes",
    "feasibility_margin",
    "interface_lines",
    "patch_depths",
]


def interface_lines(nodes, elements, retained, *, coast_nodes=()):
    """The retained interface of a cut, as lines the weight can measure from.

    The rim of the retained region has two kinds of edge: the ones that face
    the hole -- the interface, where the patch must meet frozen depths -- and
    the ones on the physical coastline. Only the first kind belongs here.

    Forcing ``w = 0`` on the coast is the mistake this argument exists to
    prevent: the free coastline inside the cut is hole boundary too, and
    weighting it to zero would give the newly resolved shore its BASE depths,
    which is the opposite of what the branch is for.
    """
    import shapely

    xy = np.asarray(nodes, dtype=float)[:, :2]
    tri = np.asarray(elements, dtype=np.int64)
    keep = np.asarray(retained, dtype=np.int64)
    if not keep.size:
        return shapely.MultiLineString([])
    sel = tri[keep]
    e = np.sort(np.vstack([sel[:, [0, 1]], sel[:, [1, 2]], sel[:, [2, 0]]]), axis=1)
    u, c = np.unique(e, axis=0, return_counts=True)
    rim = u[c == 1]
    # The mesh's own boundary is not an interface: an edge on the outer
    # boundary of the BASE mesh borders no hole.
    a = np.sort(np.vstack([tri[:, [0, 1]], tri[:, [1, 2]], tri[:, [2, 0]]]), axis=1)
    ua, ca = np.unique(a, axis=0, return_counts=True)
    outer = set(map(tuple, ua[ca == 1].tolist()))
    coast = set(np.asarray(coast_nodes, dtype=np.int64).ravel().tolist())
    lines = [xy[[i, j]] for i, j in rim.tolist()
             if (i, j) not in outer and not (i in coast and j in coast)]
    return shapely.MultiLineString([shapely.LineString(ln) for ln in lines]) \
        if lines else shapely.MultiLineString([])


def blend_weights(points, regions, interface) -> np.ndarray:
    """``w = d_int / (d_int + d_core)``: 1 on the core, 0 at the interface.

    ``points`` are ``(N, 2)`` in the mesh CRS, ``regions`` the declared
    geometries in the same CRS (their union is used, so overlapping regions
    need no separate rule), and ``interface`` the lines of
    :func:`interface_lines`.

    Both distances are zero where a core touches the interface -- a transition
    clipped by land does this -- and the weight is 1 there, because a node on
    the core is in the core whatever else is true of it.
    """
    import shapely

    p = np.asarray(points, dtype=float)[:, :2]
    if not p.size:
        return np.zeros(0)
    geoms = [regions] if hasattr(regions, "geom_type") else list(regions)
    core = shapely.union_all(geoms) if geoms else None
    q = shapely.points(p[:, 0], p[:, 1])
    d_core = np.asarray(shapely.distance(q, core)) if core is not None \
        else np.full(len(p), np.inf)
    d_int = np.asarray(shapely.distance(q, interface)) \
        if interface is not None and not interface.is_empty \
        else np.full(len(p), np.inf)
    total = d_int + d_core
    # Both zero means a core node that is also ON the interface -- a
    # transition clipped by land does it -- and a node on the core is in the
    # core whatever else is true of it, so w = 1.  An infinite d_int is a
    # patch with no interface at all (every rim edge was coast), and there the
    # source wins everywhere; dividing would give inf/inf.
    ok = np.isfinite(total) & (total > 0)
    w = np.ones(len(p))
    np.divide(d_int, total, out=w, where=ok)
    return np.clip(w, 0.0, 1.0)


def feasibility_margin(target_h_m: float, transition_m: float,
                       depth_m: float, difference_m: float,
                       rmax: float = 0.2) -> dict[str, Any]:
    """Will the later r-factor step be able to do its job here?

    ``r ~ (D / W) * h / (2H)``, so ``r <= rmax`` needs
    ``h <= 2 * H * W * rmax / D``.  Returns the allowed edge length beside the
    one asked for, so a region whose seam cannot be smoothed says so before
    anything is meshed rather than after.
    """
    d = float(abs(difference_m))
    allowed = (2.0 * float(depth_m) * float(transition_m) * float(rmax) / d
               if d > 0 else float("inf"))
    return {
        "target_h_m": float(target_h_m),
        "transition_m": float(transition_m),
        "depth_m": float(depth_m),
        "source_minus_base_m": float(difference_m),
        "rmax": float(rmax),
        "allowed_edge_m": float(allowed),
        "expected_r_at_target": float(
            d * target_h_m / (2.0 * float(depth_m) * float(transition_m)))
        if depth_m > 0 and transition_m > 0 else float("inf"),
        "ok": bool(allowed >= target_h_m),
    }


def patch_depths(lonlat, xy, base_depths, *, regions, interface,
                 scope: str = "hole", blend: str = "ramp",
                 extrapolate: bool = True) -> tuple[np.ndarray, dict[str, Any]]:
    """Source depths for the new nodes, ramped into the base field.

    ``lonlat`` and ``xy`` are the SAME nodes in geographic and mesh
    coordinates -- the new nodes only, at their final positions, because the
    repair moves them and a depth sampled before the move belongs to a
    coordinate the mesh no longer has.  ``base_depths`` is the base field
    already evaluated at those positions, which is what the ramp blends
    towards and what a ``w = 0`` node receives exactly.

    Returns ``(depth, report)``.  Nothing is floored, capped or smoothed: the
    depths are the source's, and the minimum depth and the r-factor smoothing
    are the next step (owner, 2026-09-23).  So the depths may be zero or
    negative on a tidal flat, which is a tidal flat and not a defect -- FVCOM
    integrates it with wetting and drying.
    """
    from fvcom_mesh_tools.dem import tokyo_bay

    if scope not in ("hole", "core"):
        raise ValueError(f"scope must be 'hole' or 'core', got {scope!r}")
    if blend not in ("ramp", "none"):
        raise ValueError(f"blend must be 'ramp' or 'none', got {blend!r}")
    if scope == "core" and blend == "ramp":
        raise ValueError("scope: core with blend: ramp is contradictory -- a ramp "
                         "that is not evaluated over the transition is not a ramp")

    ll = np.asarray(lonlat, dtype=float)
    p = np.asarray(xy, dtype=float)[:, :2]
    base = np.asarray(base_depths, dtype=float)
    if not (ll.shape[0] == p.shape[0] == base.shape[0]):
        raise ValueError("lonlat, xy and base_depths must describe the same nodes")
    if not base.size:
        return base.copy(), {"n_nodes": 0}
    if not np.isfinite(base).all():
        raise ValueError("the base depths handed to the blend are not all finite")

    src, rung, dist = tokyo_bay.sample(ll[:, 0], ll[:, 1], extrapolate=extrapolate)
    if not np.isfinite(src).all():
        raise ValueError(
            f"{int((~np.isfinite(src)).sum())} nodes were covered by no product "
            "and extrapolation was refused")

    import shapely

    if blend == "ramp":
        w = blend_weights(p, regions, interface)
    else:
        geoms = [regions] if hasattr(regions, "geom_type") else list(regions)
        core = shapely.union_all(geoms) if geoms else None
        inside = np.asarray(shapely.contains(
            core, shapely.points(p[:, 0], p[:, 1]))) if core is not None \
            else np.zeros(len(p), dtype=bool)
        w = np.ones(len(p)) if scope == "hole" else inside.astype(float)

    depth = w * src + (1.0 - w) * base
    report = {
        **tokyo_bay.provenance_report(rung, dist),
        "scope": scope,
        "blend": blend,
        "w_min": float(w.min()),
        "w_max": float(w.max()),
        "n_at_w1": int((w >= 1.0 - 1e-12).sum()),
        "n_at_w0": int((w <= 1e-12).sum()),
        "source_min_m": float(src.min()),
        "source_max_m": float(src.max()),
        "base_min_m": float(base.min()),
        "base_max_m": float(base.max()),
        "depth_min_m": float(depth.min()),
        "depth_max_m": float(depth.max()),
        "n_at_or_below_zero": int((depth <= 0.0).sum()),
        "max_source_minus_base_m": float(np.abs(src - base).max()),
        "median_source_minus_base_m": float(np.median(src - base)),
    }
    return depth, report


def edge_slopes(nodes, elements, depths, changed) -> dict[str, Any]:
    """Seabed gradient across the edges the patch created or moved.

    Two measures, because one of them is not always defined.  ``|dh| / L`` in
    m/m is defined on every edge including an intertidal one; the classic
    ``r = |dh| / (h_i + h_j)`` is not, and is reported only where
    ``h_i + h_j > 0``.  What is undefined on a tidal flat is that particular
    expression, not the gradient (owner, 2026-09-23).

    ``changed`` marks the new nodes, so the retained-to-new edges -- the seam
    the next step inherits -- are separated from the wholly frozen ones, which
    are unchanged by construction and are not the interesting set.
    """
    xy = np.asarray(nodes, dtype=float)[:, :2]
    tri = np.asarray(elements, dtype=np.int64)
    h = np.asarray(depths, dtype=float)
    new = np.asarray(changed, dtype=bool)
    e = np.unique(np.sort(np.vstack([tri[:, [0, 1]], tri[:, [1, 2]],
                                     tri[:, [2, 0]]]), axis=1), axis=0)
    i, j = e[:, 0], e[:, 1]
    length = np.linalg.norm(xy[i] - xy[j], axis=1)
    slope = np.where(length > 0, np.abs(h[i] - h[j]) / np.where(length > 0, length, 1.0), 0.0)
    total = h[i] + h[j]
    defined = total > 0
    r = np.where(defined, np.abs(h[i] - h[j]) / np.where(defined, total, 1.0), np.nan)
    seam = new[i] ^ new[j]                  # retained on one side, new on the other
    touched = new[i] | new[j]

    def stats(mask):
        if not mask.any():
            return {"n": 0}
        out = {"n": int(mask.sum()),
               "slope_max": float(slope[mask].max()),
               "slope_p99": float(np.percentile(slope[mask], 99))}
        d = mask & defined
        out["n_r_defined"] = int(d.sum())
        out["r_max"] = float(r[d].max()) if d.any() else None
        return out

    return {"all_changed": stats(touched), "seam": stats(seam),
            "n_edges": int(len(e)),
            "n_r_undefined": int((~defined).sum())}
