"""Channel ARCS: waterways as first-class polyline objects
(adopted from the SCHISM RiverMapper/OCSMesh design, owner
2026-07-12).

A channel is represented by its along-channel ARC (thalweg-like
polyline) plus a width. Meshing-side effects are derived FROM the
arc, never from polygon morphology:

* the water corridor is the arc buffered by ``width/2`` with FLAT
  caps -- erosion of land is confined to the arc's own corridor,
  transversal land barriers cannot be pierced;
* an explicit barrier check subtracts a protective buffer around
  every OTHER water surface and REFUSES (raises) if that splits
  the corridor -- a fabricated connection is structurally
  impossible;
* connectivity is explicit: the corridor either carries the arc
  end-to-end or the operation fails loudly.

Pure numpy/shapely; no oceanmesh import (license policy).
"""

from __future__ import annotations

from typing import Any

import numpy as np
import shapely
from shapely.geometry import LineString, MultiPolygon, Polygon
from shapely.ops import unary_union

__all__ = ["arc_from_points", "carve_channel_corridor",
           "snap_arc_to_channel", "bank_chains",
           "skeleton_branches"]


def _polys(geom) -> list[Polygon]:
    if geom is None or geom.is_empty:
        return []
    if isinstance(geom, Polygon):
        return [geom]
    if isinstance(geom, MultiPolygon):
        return [g for g in geom.geoms if not g.is_empty]
    return [g for g in getattr(geom, "geoms", [])
            if isinstance(g, Polygon) and not g.is_empty]


def arc_from_points(
    pts_xy: np.ndarray,
    *,
    smooth_passes: int = 2,
) -> np.ndarray:
    """Order scattered channel points (e.g. element centroids of a
    reference mesh's channel band) into an along-channel polyline:
    the diameter path of their nearest-neighbour graph, smoothed.
    """
    from scipy.sparse import coo_matrix
    from scipy.sparse.csgraph import shortest_path
    from scipy.spatial import cKDTree

    pts = np.asarray(pts_xy, dtype=float)
    n = len(pts)
    if n < 2:
        raise ValueError("arc needs at least 2 points")
    if n == 2:
        return pts.copy()
    k = min(5, n - 1)
    d, idx = cKDTree(pts).query(pts, k=k + 1)
    rows = np.repeat(np.arange(n), k)
    cols = idx[:, 1:].ravel()
    vals = d[:, 1:].ravel()
    g = coo_matrix((vals, (rows, cols)), shape=(n, n)).tocsr()
    dm, pred = shortest_path(g, directed=False,
                             return_predecessors=True)
    dm[~np.isfinite(dm)] = -1.0
    i, j = np.unravel_index(int(np.argmax(dm)), dm.shape)
    path = [int(j)]
    while path[-1] != i:
        p = int(pred[i, path[-1]])
        if p < 0:
            break
        path.append(p)
    arc = pts[path[::-1]]
    for _ in range(max(0, smooth_passes)):
        if len(arc) > 2:
            arc[1:-1] = (arc[:-2] + arc[1:-1] + arc[2:]) / 3.0
    return arc


def skeleton_branches(
    poly,
    *,
    metric_scale: tuple[float, float],
    density_m: float = 60.0,
    prune_m: float = 200.0,
) -> list[np.ndarray]:
    """Medial-axis BRANCH decomposition of a waterway polygon
    (owner 2026-07-12: a diameter-path arc over a large branched
    canal network cuts corners at every junction and the carve is
    refused -- the Keihin system stayed unwidened and got cut).

    Voronoi diagram of the densified boundary -> keep edges inside
    the polygon -> prune leaf spurs shorter than ``prune_m`` ->
    concatenate degree-2 chains into one polyline per BRANCH.
    Each branch is a simple centreline the corridor carve can
    follow without shortcutting.
    """
    from collections import defaultdict

    sx, sy = metric_scale
    if abs(sx - sy) / max(sx, sy) > 0.35:
        raise ValueError("metric_scale too anisotropic; project first")
    scale = 0.5 * (sx + sy)
    dens = density_m / scale
    boundary = poly.boundary
    n_pts = max(int(boundary.length / dens), 8)
    samples = [boundary.interpolate(t, normalized=True)
               for t in np.linspace(0.0, 1.0, n_pts,
                                    endpoint=False)]
    from shapely.ops import voronoi_diagram as _voronoi
    vor = _voronoi(shapely.MultiPoint(samples), edges=True)
    inner = poly.buffer(-0.02 * dens)
    segs = []
    for g in getattr(vor, "geoms", [vor]):
        for ls in getattr(g, "geoms", [g]):
            c = np.asarray(ls.coords)
            for k in range(len(c) - 1):
                seg = shapely.LineString(c[k:k + 2])
                if inner.covers(seg):
                    segs.append((tuple(np.round(c[k], 8)),
                                 tuple(np.round(c[k + 1], 8))))
    if not segs:
        return []
    adj = defaultdict(set)
    for u, v in segs:
        if u != v:
            adj[u].add(v)
            adj[v].add(u)

    def _dist(u, v):
        return float(np.hypot((u[0] - v[0]) * sx,
                              (u[1] - v[1]) * sy))

    # prune short leaf spurs (Voronoi noise toward the banks)
    changed = True
    while changed:
        changed = False
        for u in [u for u in list(adj) if len(adj[u]) == 1]:
            if len(adj.get(u, ())) != 1:   # removed this sweep
                continue
            path = [u]
            v = next(iter(adj[u]))
            plen = _dist(u, v)
            while len(adj[v]) == 2 and plen < prune_m / scale * 1.0:
                nxt = [w for w in adj[v] if w != path[-1]]
                if not nxt:
                    break
                path.append(v)
                plen += _dist(v, nxt[0])
                v = nxt[0]
            if plen < prune_m / scale:
                for k in range(len(path)):
                    a = path[k]
                    for w in list(adj[a]):
                        adj[w].discard(a)
                    adj.pop(a, None)
                changed = True

    # branches = simple chains between nodes of degree != 2
    ends = [u for u in adj if len(adj[u]) != 2]
    seen_edges = set()
    branches = []

    def _walk(u, v):
        path = [u, v]
        while len(adj[v]) == 2:
            nxt = [w for w in adj[v] if w != path[-2]]
            if not nxt or nxt[0] == path[-1]:
                break
            v = nxt[0]
            path.append(v)
            if len(adj[v]) != 2:
                break
        return path

    for u in ends:
        for v in adj[u]:
            key = (u, v)
            if key in seen_edges:
                continue
            path = _walk(u, v)
            for k in range(len(path) - 1):
                seen_edges.add((path[k], path[k + 1]))
                seen_edges.add((path[k + 1], path[k]))
            branches.append(np.asarray(path, dtype=float))
    # isolated loops (all degree 2): walk them once
    for u in list(adj):
        if len(adj[u]) == 2:
            v = next(iter(adj[u]))
            if (u, v) not in seen_edges:
                path = _walk(u, v)
                for k in range(len(path) - 1):
                    seen_edges.add((path[k], path[k + 1]))
                    seen_edges.add((path[k + 1], path[k]))
                branches.append(np.asarray(path, dtype=float))
    return [b for b in branches
            if len(b) >= 2
            and shapely.LineString(b).length * scale >= prune_m]


def snap_arc_to_channel(
    land_union,
    arc_ll: np.ndarray,
    *,
    metric_scale: tuple[float, float],
    step_m: float = 120.0,
    max_halfwidth_m: float = 600.0,
    smooth_passes: int = 1,
) -> dict[str, np.ndarray]:
    """Snap a rough guide arc onto the REAL channel: at stations
    every ``step_m``, cast a perpendicular cross-section, take the
    water gap nearest the station, and move the station to the gap
    centre; record the gap width and the two bank points (exactly
    on the shoreline). Blocked stations (piers/bridges drawn as
    land, or no gap within ``max_halfwidth_m``) get centre, width
    and banks interpolated from the clean neighbours.

    Returns ``{"arc", "width_m", "bank_left", "bank_right",
    "snapped"}`` -- arrays of length K (``snapped`` marks stations
    measured on real water, not interpolated).
    """
    sx, sy = metric_scale
    if abs(sx - sy) / max(sx, sy) > 0.35:
        raise ValueError("metric_scale too anisotropic; project first")
    scale = 0.5 * (sx + sy)
    arc = LineString(np.asarray(arc_ll, dtype=float))
    n_st = max(int(round(arc.length * scale / step_m)) + 1, 5)
    ss = np.linspace(0.0, 1.0, n_st)
    centers = np.zeros((n_st, 2))
    banks_l = np.zeros((n_st, 2))
    banks_r = np.zeros((n_st, 2))
    widths = np.full(n_st, np.nan)
    half = max_halfwidth_m / scale
    for k, s in enumerate(ss):
        p = arc.interpolate(s, normalized=True)
        pa = arc.interpolate(max(s - 0.02, 0.0), normalized=True)
        pb = arc.interpolate(min(s + 0.02, 1.0), normalized=True)
        t = np.array([pb.x - pa.x, pb.y - pa.y])
        t /= np.hypot(*t) + 1e-15
        nv = np.array([-t[1], t[0]])
        sec = LineString([(p.x - nv[0] * half, p.y - nv[1] * half),
                          (p.x + nv[0] * half, p.y + nv[1] * half)])
        wat = sec.difference(land_union)
        segs = ([g for g in getattr(wat, "geoms", [wat])
                 if not g.is_empty and g.length > 0]
                if not wat.is_empty else [])
        centers[k] = (p.x, p.y)
        if not segs:
            continue
        best = min(segs, key=lambda g: g.distance(p))
        if best.distance(p) * scale > 0.5 * max_halfwidth_m:
            continue
        c = best.interpolate(0.5, normalized=True)
        ends = np.asarray(best.coords)
        centers[k] = (c.x, c.y)
        banks_l[k], banks_r[k] = ends[0], ends[-1]
        widths[k] = best.length * scale
    span = 2.0 * max_halfwidth_m
    saturated = np.isfinite(widths) & (widths >= 0.95 * span)
    clean = np.isfinite(widths)
    if clean.sum() < 2:
        raise RuntimeError(
            "snap_arc_to_channel: fewer than 2 stations found real "
            "water -- the guide arc does not follow a channel")
    # pinch detection: reference = median CHANNEL width. Stations
    # whose section saturates (open-mouth water spanning the whole
    # section) would inflate the median, so exclude them from the
    # reference; a pinch is a sharp drop below that reference
    # (piers/bridges drawn as land across a real passage).
    ref = widths[clean & ~saturated]
    med = float(np.median(ref)) if len(ref) else np.nan
    if np.isfinite(med):
        clean &= ~(widths < 0.4 * med)    # pinches -> interpolate
    idx = np.arange(n_st, dtype=float)
    cl = np.where(clean)[0]
    for arr in (centers, banks_l, banks_r):
        for col in (0, 1):
            arr[~clean, col] = np.interp(idx[~clean], idx[cl],
                                         arr[cl, col])
    widths[~clean] = np.interp(idx[~clean], idx[cl], widths[cl])
    # TRIM leading/trailing open-mouth stations (section saturated
    # = wide-open water): their "centres" are mouth centres, not
    # the channel axis, and smoothing would drag the corridor onto
    # the real banks at the channel entrance. The corridor's round
    # end still overlaps the mouth water, so connectivity holds.
    nz = np.where(~saturated)[0]
    if len(nz) >= 2:
        s0, s1 = int(nz[0]), int(nz[-1]) + 1
        centers = centers[s0:s1]
        banks_l, banks_r = banks_l[s0:s1], banks_r[s0:s1]
        widths, clean = widths[s0:s1], clean[s0:s1]
        n_st = len(centers)
    for _ in range(max(0, smooth_passes)):
        if n_st > 2:
            centers[1:-1] = (centers[:-2] + centers[1:-1]
                             + centers[2:]) / 3.0
    # carving width: at open mouths the section saturates and a
    # corridor that wide would bulge onto the REAL banks at the
    # open->channel joints (round segment joins). All water there
    # anyway -- cap the carve width near the channel's own scale.
    w_carve = (np.minimum(widths, 1.6 * med)
               if np.isfinite(med) else widths.copy())
    return {"arc": centers, "width_m": widths,
            "width_carve_m": w_carve,
            "bank_left": banks_l, "bank_right": banks_r,
            "snapped": clean}


def bank_chains(
    snap: dict[str, np.ndarray],
    *,
    spacing_m: float,
    metric_scale: tuple[float, float],
    max_width_factor: float = 1.8,
    inset: float = 0.48,
) -> tuple[np.ndarray, np.ndarray]:
    """Build pfix points + egfix segments constraining BOTH banks
    of the sub-cell-width portion of a snapped channel (the proven
    OBC-ladder primitive applied to channel banks).

    Stations are kept where the measured width is below
    ``max_width_factor * spacing_m`` (wider water meshes fine
    unconstrained), subsampled to ~``spacing_m``, and each bank
    point is pulled ``inset`` (fraction of local width, 0.48 =
    ~2 % inside the water) toward the centre so a fixed node never
    lands outside the meshing domain when the carve differs from
    the measured bank by a metre or two.

    Returns ``(pfix (P,2), egfix (E,2) int indices into pfix)``.
    """
    sx, sy = metric_scale
    scale = 0.5 * (sx + sy)
    arc = snap["arc"]
    w = snap["width_m"]
    seg = np.hypot(*(np.diff(arc, axis=0).T)) * scale
    cum = np.concatenate([[0.0], np.cumsum(seg)])
    narrow = w < max_width_factor * spacing_m
    pf: list[np.ndarray] = []
    eg: list[list[int]] = []
    k = 0
    while k < len(arc):
        if not narrow[k]:
            k += 1
            continue
        run = [k]
        while k + 1 < len(arc) and narrow[k + 1]:
            k += 1
            run.append(k)
        k += 1
        if cum[run[-1]] - cum[run[0]] < 0.8 * spacing_m:
            continue
        keep = [run[0]]
        for i in run[1:]:
            if cum[i] - cum[keep[-1]] >= 0.8 * spacing_m:
                keep.append(i)
        if len(keep) < 2:
            continue
        for side in ("bank_left", "bank_right"):
            base = len(pf)
            for i in keep:
                b = snap[side][i]
                pt = arc[i] + (b - arc[i]) * 2.0 * inset
                pf.append(pt)
            eg.extend([[base + j, base + j + 1]
                       for j in range(len(keep) - 1)])
    if not pf:
        return np.zeros((0, 2)), np.zeros((0, 2), dtype=int)
    return np.asarray(pf), np.asarray(eg, dtype=int)


def carve_channel_corridor(
    land_union,
    arc_ll: np.ndarray,
    width_m,
    *,
    min_gap_m: float,
    metric_scale: tuple[float, float],
    domain_poly=None,
    arc_on_land_tol_m: float | None = None,
    carve_crossings: bool = True,
) -> tuple[Any, dict[str, Any]]:
    """Cut a channel corridor of ``width_m`` along ``arc_ll`` out
    of ``land_union``, barrier-safely.

    The corridor is the flat-capped buffer of the arc. Guards
    (both fail loudly, never silently):

    1. the ARC must follow existing water -- crossing land for
       more than ``arc_on_land_tol_m`` (default ``width_m``)
       raises: a wrong arc would pierce a transversal barrier;
    2. every OTHER water surface is kept behind ``min_gap_m`` of
       land (arc-end openings exempted);
    3. with ``carve_crossings=False`` (DETECTED waterways) land
       around within-tolerance arc/land crossings is preserved --
       an arc clipping a bank corner must not carve the corner
       into a passage. ``True`` (explicit manual edits, e.g. a
       pier drawn as land) carves through;
    4. after carving, the two arc ENDS must lie in the same
       connected component of the resulting water, else raise.

    Returns ``(new_land_union, info)``.
    """
    sx, sy = metric_scale
    if abs(sx - sy) / max(sx, sy) > 0.35:
        raise ValueError("metric_scale too anisotropic; project first")
    scale = 0.5 * (sx + sy)
    pts = np.asarray(arc_ll, dtype=float)
    arc = LineString(pts)
    w = np.asarray(width_m, dtype=float)
    w_max = float(w.max()) if w.ndim else float(w)

    # INVARIANT: the arc runs along EXISTING water (a channel
    # centreline); widening pushes the BANKS sideways into land.
    # An arc that crosses land for more than a noise tolerance is
    # a wrong arc -- it would pierce a transversal barrier, the
    # exact bug class this module exists to prevent. Fail loudly.
    tol = (arc_on_land_tol_m if arc_on_land_tol_m is not None
           else w_max)
    on_land = arc.intersection(land_union)
    on_land_m = float(on_land.length) * scale
    if on_land_m > tol:
        raise RuntimeError(
            f"channel arc crosses land for {on_land_m:.0f} m "
            f"(> tolerance {tol:.0f} m): the arc does not follow "
            "an existing waterway and would pierce a land barrier. "
            "Fix the arc; do not raise arc_on_land_tol_m unless "
            "the crossing is shoreline-data noise.")

    if w.ndim == 0:
        corridor = arc.buffer(0.5 * float(w) / scale,
                              cap_style="flat")
    else:
        # per-vertex widths (measured channel profile): union of
        # per-segment buffers, round joints covering the seams.
        # Segment width = MIN of its endpoints: carving stays on
        # the conservative side, so a wide->narrow transition can
        # never bulge onto the narrow stretch's real banks
        if len(w) != len(pts):
            raise ValueError("width_m array must match arc points")
        corridor = unary_union([
            LineString(pts[i:i + 2]).buffer(
                0.5 * float(min(w[i], w[i + 1])) / scale)
            for i in range(len(pts) - 1)])
        # trim the round END overshoot with flat cuts: at a
        # dead-end head the w/2 bulb beyond the last station would
        # be carved but never meshed, leaving one flat bank-to-
        # bank wedge cell (C1/C4 tail, element 3894 run 6184675)
        big = 2.0 * w_max / scale
        for i_end, i_prev in ((0, 1), (len(pts) - 1,
                                       len(pts) - 2)):
            t = pts[i_end] - pts[i_prev]
            t = t / (np.hypot(*t) + 1e-15)
            nvec = np.array([-t[1], t[0]])
            p0 = pts[i_end]
            half = shapely.Polygon([
                p0 - nvec * big, p0 + nvec * big,
                p0 + nvec * big + t * big,
                p0 - nvec * big + t * big])
            corridor = corridor.difference(half)

    if domain_poly is not None:
        water = domain_poly.difference(land_union)
        # protect every LOCAL water component the arc does NOT run
        # through. Protecting only water OUTSIDE the corridor was
        # blind to a thin barrier lying INSIDE the corridor width:
        # the water beyond it fell inside the corridor too, lost
        # its protection, and the barrier was carved into a
        # fabricated passage (breach at 139.8991/35.3703, run
        # 6184643).
        local = water.intersection(
            corridor.buffer(2.0 * min_gap_m / scale))
        arcb = arc.buffer(1e-9)
        other_parts = [g for g in _polys(local)
                       if not g.intersects(arcb)]
        protect = (unary_union(other_parts).buffer(
            min_gap_m / scale) if other_parts
            else shapely.Polygon())
        # the arc ENDS are where the channel must open into water:
        # exempt a disk around each end from protection, otherwise
        # a receiving body just missed by the arc seals the mouth
        w2e = 0.5 * w_max / scale
        end_zones = unary_union(
            [shapely.Point(pts[0]).buffer(w2e + min_gap_m / scale),
             shapely.Point(pts[-1]).buffer(w2e + min_gap_m / scale)])
        protect = protect.difference(end_zones)
        carved = corridor.difference(protect)
        if not carve_crossings and on_land_m > 1.0:
            # DETECTED arcs: a within-tolerance land crossing is
            # bank-corner clipping noise, NOT a licence to carve
            # the corner into a passage (the corridor covers both
            # sides there, so min_gap protection is blind to it --
            # 5 fabricated passages in comparator run 6184563).
            cross = [g for g in getattr(on_land, "geoms",
                                        [on_land])
                     if g.length * scale > 20.0]
            if cross:
                carved = carved.difference(unary_union(
                    [g.buffer(min_gap_m / scale) for g in cross]))
        # the guard: after carving, the two arc ENDS must lie in
        # the SAME connected component of the resulting water --
        # fails exactly when a land barrier still separates them.
        merged = water.union(carved)
        pA = shapely.Point(pts[0]).buffer(2.0 / scale)
        pB = shapely.Point(pts[-1]).buffer(2.0 / scale)
        connected = any(
            gpart.intersects(pA) and gpart.intersects(pB)
            for gpart in _polys(merged))
        if not connected:
            raise RuntimeError(
                "channel cannot be carved without piercing a land "
                "barrier protecting other water (the arc ends stay "
                f"in separate water bodies at min_gap "
                f"{min_gap_m:.0f} m). Options: reduce width_m, "
                "adjust the arc, or explicitly accept the "
                "connection with carve_crossings=True and an "
                "explicit arc_on_land_tol_m.")
    else:
        carved = corridor
    new_land = land_union.difference(carved)
    land_removed = land_union.intersection(carved)
    return new_land, {
        "corridor_area_m2": float(carved.area * scale * scale),
        "land_removed_m2": float(land_removed.area * scale * scale),
        "arc_length_m": float(arc.length * scale),
        "arc_on_land_m": on_land_m,
        "width_max_m": w_max,
    }
