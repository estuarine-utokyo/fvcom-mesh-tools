"""Mathematical waterway detection from the shoreline geometry
(owner 2026-07-12: waterways must be found from the ORIGINAL OSM
data -- no reference mesh involved).

Pipeline (all geometry, resolution-aware):

1.  Water = the sea-connected component of ``domain - land``
    (shapefile coverage gaps fabricate phantom inland water).
2.  Morphological opening with radius ``0.5 * detect_factor *
    h_mesh`` splits water into WIDE bodies and NARROW corridors.
3.  Narrow corridors and the small wide pockets between them are
    chained into waterway NETWORKS (adjacency graph).
4.  Each network becomes an ARC (+ measured width profile) via the
    channel-arc machinery -- widening happens along the arc, never
    by isotropic polygon buffering.
5.  Classification (owner rules):
    * THROUGH (both ends reach the main body, and the along-shore
      detour between the contact points is > ``shortcut_ratio`` x
      the through length) or connecting the main body to a big
      basin (>= ``min_basin_cells`` cell-equivalents): KEEP -- the
      banks are pushed into land until ``widen_rows`` standard
      rows fit. The minimum mesh size is inviolable.
    * dead ends and small-basin feeders: CLOSE (fill as land).
    * a KEEP whose corridor cannot be carved without piercing a
      protected barrier is reported as BLOCKED, loudly -- never
      silently reclassified.
"""

from __future__ import annotations

from typing import Any

import numpy as np
import shapely
from shapely.geometry import MultiPolygon, Polygon
from shapely.ops import unary_union

from fvcom_mesh_tools.channel_arcs import (
    arc_from_points,
    carve_channel_corridor,
    skeleton_branches,
    snap_arc_to_channel,
)

__all__ = ["detect_waterways", "apply_waterway_policy"]


def _polys(geom) -> list[Polygon]:
    if geom is None or geom.is_empty:
        return []
    if isinstance(geom, Polygon):
        return [geom]
    if isinstance(geom, MultiPolygon):
        return [g for g in geom.geoms if not g.is_empty]
    return [g for g in getattr(geom, "geoms", [])
            if isinstance(g, Polygon) and not g.is_empty]


def _interior_points(geom, spacing):
    """Point cloud inside a (multi)polygon on a regular lattice
    (fallback: representative points of its parts)."""
    x0, y0, x1, y1 = geom.bounds
    nx = max(int((x1 - x0) / spacing) + 1, 2)
    ny = max(int((y1 - y0) / spacing) + 1, 2)
    xs = np.linspace(x0, x1, nx)
    ys = np.linspace(y0, y1, ny)
    gx, gy = np.meshgrid(xs, ys)
    pts = shapely.points(gx.ravel(), gy.ravel())
    inside = shapely.covers(geom, pts)
    out = np.column_stack([gx.ravel()[inside], gy.ravel()[inside]])
    if len(out) < 2:
        out = np.array([[g.representative_point().x,
                         g.representative_point().y]
                        for g in _polys(geom)])
    return out


def _ladder_constraints(
    snapd: dict,
    h_mesh_m: float,
    metric_scale: tuple[float, float],
    inset: float = 0.47,
) -> tuple[np.ndarray, np.ndarray]:
    """Two-row LADDER constraints for a marginal channel band
    (owner 2026-07-13): fix three node rows -- both banks (on the
    carved shoreline, pulled ``inset`` inward) and a staggered
    centreline row -- with egfix chains along each row. The CDT
    then forms an orderly two-row band, the same primitive that
    keeps the OBC ladder tidy.

    Station spacing ~= the local two-row edge length; the first
    and last stations are trimmed so junctions and open mouths
    stay unconstrained."""
    arc = np.asarray(snapd["arc"], float)
    bl = np.asarray(snapd["bank_left"], float)
    br = np.asarray(snapd["bank_right"], float)
    wm = np.asarray(snapd["width_m"], float)
    ok = np.asarray(snapd["snapped"], bool)
    sx, sy = metric_scale
    scale = 0.5 * (sx + sy)
    seg = np.hypot(*(np.diff(arc, axis=0).T)) * scale
    cum = np.concatenate([[0.0], np.cumsum(seg)])
    L = float(cum[-1])
    sat = wm >= 0.95 * 2.0 * 2.5 * h_mesh_m
    wn = wm[~sat] if bool((~sat).any()) else wm
    w_med = float(np.median(wn))
    a = max(0.8 * h_mesh_m, w_med / 1.8)
    # refinements after runs 6188830/6188855 (27 -> 17 C1):
    # ladder stations must sit on MEASURED cross-sections whose
    # width is close to the band median; interpolated (pinch/
    # bridge) stretches and width outliers get no fixed nodes,
    # and generous end margins let the band taper into the
    # unconstrained mesh.
    usable = (ok & ~sat
              & (np.abs(wm - w_med) <= 0.25 * w_med))
    if L < 6.0 * a:
        return np.zeros((0, 2)), np.zeros((0, 2), dtype=int)

    def _interp(rowpts, s_at):
        out = np.empty((len(s_at), 2))
        for c2 in (0, 1):
            out[:, c2] = np.interp(s_at, cum, rowpts[:, c2])
        return out

    def _usable_at(s_at):
        j = np.searchsorted(cum, s_at)
        j0 = np.clip(j - 1, 0, len(usable) - 1)
        j1 = np.clip(j, 0, len(usable) - 1)
        return usable[j0] & usable[j1]

    s_all = np.arange(1.5 * a, L - 1.5 * a, a)
    if len(s_all) < 3:
        return np.zeros((0, 2)), np.zeros((0, 2), dtype=int)
    keepm = _usable_at(s_all)
    pf_list = []
    eg = []
    # contiguous usable runs become independent ladder segments
    runs = []
    start = None
    for k, good in enumerate(keepm):
        if good and start is None:
            start = k
        elif not good and start is not None:
            runs.append((start, k))
            start = None
    if start is not None:
        runs.append((start, len(keepm)))
    for r0, r1 in runs:
        if r1 - r0 < 3:      # too short to constrain
            continue
        s_nodes = s_all[r0:r1]
        cL = _interp(arc, s_nodes)
        pL = cL + (_interp(bl, s_nodes) - cL) * 2.0 * inset
        pR = cL + (_interp(br, s_nodes) - cL) * 2.0 * inset
        s_mid = s_nodes[:-1] + 0.5 * a
        pC = _interp(arc, s_mid)
        base = sum(len(x) for x in pf_list)
        nL = len(s_nodes)
        pf_list += [pL, pR, pC]
        for k in range(nL - 1):
            eg.append([base + k, base + k + 1])
            eg.append([base + nL + k, base + nL + k + 1])
        # centre row: pfix only -- an egfix centre chain
        # over-constrains the diagonals (C2 tail, run 6188855)
    if not pf_list:
        return np.zeros((0, 2)), np.zeros((0, 2), dtype=int)
    return np.vstack(pf_list), np.asarray(eg, dtype=int)


def _ladder_size_targets(snapd, a_m):
    """Size-field override points for a ladder band: the fixed
    rows only work when the ambient sizing matches their spacing
    (the OBC-band lesson) -- return (points_ll, field_targets_m
    = a/1.2) sampled along the band arc."""
    arc = np.asarray(snapd["arc"], float)
    pts = arc[::1]
    tgt = np.full(len(pts), a_m / 1.2)
    return pts, tgt


def detect_waterways(
    land_union,
    domain_poly,
    *,
    h_mesh_m: float,
    obc_point,
    metric_scale: tuple[float, float],
    detect_factor: float = 1.8,
    min_basin_cells: int = 6,
    shortcut_ratio: float = 2.2,
    min_extent_cells: float = 2.0,
    big_deadend_cells: float = 6.0,
    max_canal_extent_cells: float = 15.0,
    min_canal_width_frac: float = 0.5,
    min_resolve_width_frac: float = 0.2,
) -> list[dict[str, Any]]:
    """Find sub-``detect_factor*h`` waterways and decide their
    fate. Returns one record per waterway network:
    ``{arc, width_m (profile), action: keep|close, kind, geometry,
    extent_cells, basin_cells}``.
    """
    sx, sy = metric_scale
    if abs(sx - sy) / max(sx, sy) > 0.35:
        raise ValueError("metric_scale too anisotropic; project first")
    scale = 0.5 * (sx + sy)
    h = h_mesh_m / scale
    r_open = 0.5 * detect_factor * h
    a_cell = (np.sqrt(3.0) / 4.0) * h * h

    water = domain_poly.difference(land_union)
    wpolys = _polys(water)
    if not wpolys:
        raise RuntimeError("domain minus land left no water")
    obc_pt = shapely.Point(obc_point)
    water = wpolys[int(np.argmin([obc_pt.distance(g)
                                  for g in wpolys]))]

    wide = water.buffer(-r_open).buffer(
        r_open * 1.02, join_style="mitre", mitre_limit=1.2)
    wide = wide.intersection(water)
    wide_parts = _polys(wide)
    if not wide_parts:
        raise RuntimeError("opening removed ALL water -- h_mesh_m "
                           "too large for this domain?")
    main_i = int(np.argmin([obc_pt.distance(g)
                            for g in wide_parts]))
    main_poly = wide_parts[main_i]
    areas = np.array([g.area for g in wide_parts])
    big = {k for k in range(len(wide_parts))
           if k == main_i or areas[k] >= min_basin_cells * a_cell}
    small = [k for k in range(len(wide_parts)) if k not in big]

    narrow_parts = [N for N in _polys(water.difference(wide))
                    if N.area >= 0.05 * a_cell]

    # NETWORKS: corridors + small pockets that touch each other,
    # PLUS bridge-gap chaining (owner 2026-07-12, Daishi canal
    # severance): OSM draws bridges/roads as land strips that chop
    # a canal into disjoint pieces, each of which then classifies
    # as a dead-end stub. Pieces whose gap is a SHORT land neck
    # (<= 0.7 h) with end-to-end facing (small shared frontage --
    # two PARALLEL canals across a land strip share a long
    # frontage and must NOT chain) belong to one waterway; the
    # gap is recorded and carved open for kept networks.
    from scipy.sparse import coo_matrix
    from scipy.sparse.csgraph import connected_components
    eps = 0.02 * h
    bridge_gap = 0.7 * h
    items = ([("n", N) for N in narrow_parts]
             + [("s", wide_parts[k]) for k in small])
    n_it = len(items)
    rows, cols = [], []
    bridge_lines: dict[tuple[int, int], Any] = {}
    for i in range(n_it):
        for j in range(i + 1, n_it):
            d = items[i][1].distance(items[j][1])
            if d < eps:
                rows.append(i)
                cols.append(j)
            elif (d < bridge_gap
                  and items[i][0] == "n" and items[j][0] == "n"):
                gi, gj = items[i][1], items[j][1]
                front = gi.boundary.intersection(
                    gj.buffer(d + 0.1 * h)).length
                w_i = 2.0 * gi.area / max(gi.boundary.length,
                                          1e-9)
                w_j = 2.0 * gj.area / max(gj.boundary.length,
                                          1e-9)
                if front <= 3.0 * max(w_i, w_j, 0.3 * h):
                    rows.append(i)
                    cols.append(j)
                    bridge_lines[(i, j)] = (
                        shapely.shortest_line(gi, gj), gi, gj)
    if n_it:
        g2 = coo_matrix((np.ones(len(rows)), (rows, cols)),
                        shape=(n_it, n_it))
        _, lab = connected_components(g2 + g2.T, directed=False)
    else:
        lab = np.zeros(0, dtype=int)

    records: list[dict[str, Any]] = []
    for net in range(lab.max() + 1 if n_it else 0):
        members = np.where(lab == net)[0]
        geoms = [items[i][1] for i in members]
        union = unary_union(geoms)
        mset = set(int(v) for v in members)
        net_bridges = [tup for (i, j), tup in
                       bridge_lines.items()
                       if i in mset and j in mset]
        hull_len = float(union.minimum_rotated_rectangle.length) / 2
        extent_cells = hull_len / h
        # anchors on WIDE water
        touches_big = [k for k in big if k != main_i
                       and union.distance(wide_parts[k]) < eps]
        basin_cells = (max(areas[k] / a_cell for k in touches_big)
                       if touches_big else 0.0)
        inter = union.buffer(eps).intersection(main_poly)
        pieces = sorted(_polys(inter), key=lambda p: -p.area)
        through = False
        if len(pieces) >= 2:
            A = pieces[0].representative_point()
            B = pieces[1].representative_point()
            d_thru = max(A.distance(B), 0.1 * h)
            rings = [main_poly.exterior, *main_poly.interiors]
            rA = int(np.argmin([r.distance(A) for r in rings]))
            rB = int(np.argmin([r.distance(B) for r in rings]))
            if rA != rB:
                through = True
            else:
                ring = rings[rA]
                sA, sB = ring.project(A), ring.project(B)
                d_arc = abs(sA - sB)
                d_arc = min(d_arc, ring.length - d_arc)
                through = d_arc > shortcut_ratio * d_thru
        # ANCHOR rule (owner 2026-07-12, Keihin canal severance):
        # canal systems are chains main - corridor - big pocket -
        # corridor - big pocket ...; a link between two big
        # anchors never touches the main body directly, yet
        # closing it cuts the chain. Any network touching >= 2
        # anchors (main and/or big pockets) is a CONNECTOR.
        n_anchors = (1 if len(pieces) >= 1 else 0) \
            + len(touches_big)
        connector = through or n_anchors >= 2
        worthy = extent_cells >= min_extent_cells
        # PORT-CANAL rule (goto2023 design, resolution-relative):
        # a dead-end waterway of substantial-but-BOUNDED extent
        # and canal-like width (mean >= min_canal_width_frac * h)
        # is kept -- port canal systems carry tidal prism. Rivers
        # fail it on one of the two axes: below-resolution ones
        # (Hanami-gawa) are too narrow, large ones are far longer
        # than any port canal.
        # width of the DOMINANT corridor piece (2A/P): micro
        # opening-artifacts chained into the union dilute a
        # union-level estimate below the resolve floor
        narrow_members = [items[i][1] for i in members
                          if items[i][0] == "n"]
        main_piece = (max(narrow_members, key=lambda g: g.area)
                      if narrow_members else None)
        wsrc = main_piece if main_piece is not None else union
        mean_w_cells = (2.0 * wsrc.area
                        / max(wsrc.boundary.length, 1e-9)) / h
        big_canal = (big_deadend_cells <= extent_cells
                     <= max_canal_extent_cells
                     and mean_w_cells >= min_canal_width_frac)
        # RESOLVE-WIDTH floor (owner 2026-07-12): a channel whose
        # NATURAL width is far below the minimum mesh size is not
        # a resolve target AT ALL -- the sample leaves such
        # ditches untouched, and keeping them forced 609 m
        # corridors across land with one-wide remnants. Applies
        # to EVERY keep path (through/anchor rules had no width
        # condition).
        resolvable = mean_w_cells >= min_resolve_width_frac
        keep = (connector or big_canal) and worthy and resolvable

        # a narrow piece that barely touches LAND is not a
        # waterway at all -- it is an opening artifact against the
        # domain boundary or a sliver of the main body. Filling it
        # would turn open water into land: IGNORE instead.
        # EXCEPTION: a strip pinned between the ARTIFICIAL domain
        # edge and land (land_frac >= 0.25, touching the domain
        # exterior) is a data CRACK, not water -- left open it
        # meshes as huge sliver cells at the edge (16.8 deg C1,
        # run 6186561): close it.
        blen = float(union.boundary.length)
        land_frac = (float(union.boundary.intersection(
            land_union.buffer(eps)).length) / blen
            if blen > 0 else 0.0)

        rec: dict[str, Any] = {
            "kind": ("through" if through else
                     "port" if touches_big else
                     "canal" if big_canal else
                     "dead-end"),
            "action": ("keep" if keep else
                       "close" if (land_frac >= 0.5
                                   or (land_frac >= 0.25
                                       and union.distance(
                                           domain_poly.exterior)
                                       < eps))
                       else "ignore"),
            "land_frac": round(land_frac, 2),
            "mean_width_cells": round(float(mean_w_cells), 2),
            "main_piece": main_piece,
            "bridges": net_bridges,
            "extent_cells": round(float(extent_cells), 1),
            "basin_cells": round(float(basin_cells), 1),
            "geometry": union,
            "arc": None,
            "width_m": None,
        }
        if keep:
            # waterway as an ARC: interior point cloud -> ordered
            # medial guide -> snapped centreline + width profile
            pts = _interior_points(union, 0.35 * h)
            guide = arc_from_points(
                np.column_stack([pts[:, 0] * sx / scale,
                                 pts[:, 1] * sy / scale]),
                smooth_passes=2)
            guide = np.column_stack([guide[:, 0] * scale / sx,
                                     guide[:, 1] * scale / sy])
            try:
                snap = snap_arc_to_channel(
                    land_union, guide, metric_scale=metric_scale,
                    step_m=0.35 * h_mesh_m,
                    max_halfwidth_m=2.5 * h_mesh_m)
                rec["arc"] = snap["arc"]
                rec["width_m"] = snap["width_carve_m"]
            except RuntimeError as e:
                rec["action"] = "blocked"
                rec["reason"] = f"arc extraction failed: {e}"
        records.append(rec)
    return records


def apply_waterway_policy(
    land_union,
    domain_poly,
    records: list[dict[str, Any]],
    *,
    h_mesh_m: float,
    metric_scale: tuple[float, float],
    widen_rows: float = 2.0,
    min_gap_m: float = 150.0,
    h_grade_per_m: float = 0.0,
    branch_floor_frac: float = 0.45,
    close_blocked: bool = True,
    open_bridges: bool | str = "auto",
    waterway_lines=None,
    force_two_rows: bool = False,
    widen_factor: float = 1.0,
    attain_bar_h: float = 1.7,
) -> tuple[Any, dict[str, Any]]:
    """Execute the detected actions: KEEP -> carve the corridor to
    two LOCAL rows along the arc, barrier-safe; CLOSE -> fill the
    network as land.

    The per-station widen target accounts for distance-to-coast
    size growth: ``widen_factor * rows * (h + h_grade_per_m *
    w_i / 2)`` -- at a 350 m natural channel the mid-channel mesh
    edge is already ~h + grade*175, and a target based on bare
    ``h`` realises only ~1.6 rows.  ``widen_factor`` was 0.875
    through run 6188868; the realized widths came out 1.83-1.96 h
    on width-limited stretches and 5 of the 13 one-wide chokes in
    the ledger sat exactly there (DistMesh needs ~2 h to place
    two rows reliably), so the factor is now 1.0.

    A keep whose carve is refused is retried with a MEDIAL-AXIS
    branch decomposition (``arc_retry`` = "skeleton"): each simple
    branch is snapped and carved on its own, so junction
    shortcuts cannot occur; the retry is recorded. Still refused
    -> BLOCKED, and with ``close_blocked=True`` (owner rule
    2026-07-12: a channel that cannot be made two rows wide must
    not be meshed at all) the network is FILLED like a close --
    never left as sub-cell water that meshes one cell wide. The
    record keeps action "blocked" plus ``closed=True`` so the
    decision stays visible."""
    new_land = land_union
    info = {"kept": 0, "closed": 0, "ignored": 0, "blocked": [],
            "retried": 0, "bridges_opened": 0, "stub_fills": 0,
            "marginal_kept": 0, "dup_skipped": 0,
            "line_branches": 0, "skel_branches": 0,
            "refine_arcs": [],
            "band_pfix": [], "band_egfix": [], "band_n": 0,
            "band_size_pts": [], "band_size_tgt": [],
            "land_removed_m2": 0.0}
    fills = []
    sx, sy = metric_scale
    scale = 0.5 * (sx + sy)
    wlist = ([g for g in waterway_lines]
             if waterway_lines is not None else [])
    from shapely.strtree import STRtree
    wtree_g = STRtree(wlist) if wlist else None

    def _carve(base, arc, widths, tol_extra_m=0.0,
               attain_bar_h=attain_bar_h):
        w_nat = np.asarray(widths, float)
        target = widen_factor * widen_rows * (
            h_mesh_m + h_grade_per_m * w_nat / 2.0)
        w = np.maximum(w_nat, target)
        cand, ci = carve_channel_corridor(
            base, arc, w, min_gap_m=min_gap_m,
            metric_scale=metric_scale, domain_poly=domain_poly,
            arc_on_land_tol_m=0.3 * h_mesh_m + tol_extra_m,
            carve_crossings=False)
        # ACHIEVED-width check (owner 2026-07-12: a channel that
        # cannot actually be made ~two rows wide must NOT be
        # meshed): barrier protection can legally trim the
        # corridor back to the natural width, in which case the
        # carve "succeeds" but the channel would mesh one cell
        # wide. Measure the result; too narrow -> refuse.
        ach = snap_arc_to_channel(
            cand, np.asarray(arc, float),
            metric_scale=metric_scale, step_m=0.5 * h_mesh_m,
            max_halfwidth_m=2.5 * h_mesh_m)
        wa = ach["width_m"]
        # 1.7 h attainability bar (was 1.5 h under widen_factor
        # 0.875): DistMesh places two rows reliably only from
        # ~1.9-2.0 h realized, and 1.5-1.7 h is exactly the
        # marginal regime the ledger chokes came from.
        # attain_bar_h=0 disables the bar (connectivity-critical
        # branches, see _closure_severs).
        if attain_bar_h > 0.0:
            narrow_frac = float(np.mean(
                wa < attain_bar_h * h_mesh_m))
            if narrow_frac > 0.35:
                raise RuntimeError(
                    f"achieved width < {attain_bar_h:.1f} h over "
                    f"{narrow_frac:.0%} of the arc (min "
                    f"{np.min(wa):.0f} m): two rows are not "
                    "attainable under the barrier constraints")
        return cand, ci, w, ach

    def _closure_severs(base, arc, widths):
        """Would refusing this branch (leaving it sub-cell, hence
        unmeshed) sever local water connectivity?  CONNECTIVITY IS
        SACRED (owner): a real passage must never be cut, even
        when two rows are not attainable -- the marginal cells go
        to the one-wide ledger instead.  Test: remove the branch
        tube from the water in a local window; if the water just
        beyond the two arc ends no longer falls in one piece, the
        branch is the only local connection (run 6190466: the N7
        Chiba-canal branch, barrier-capped below the 1.7 h bar,
        severed a 24-element sample passage when refused)."""
        a = np.asarray(arc, float)
        if len(a) < 2:
            return False
        ln = shapely.LineString(a)
        rad = max(0.9 * float(np.median(np.asarray(widths,
                                                   float))),
                  0.5 * h_mesh_m) / scale
        win = shapely.box(*ln.buffer(10.0 * h_mesh_m
                                     / scale).bounds)
        water = win.intersection(domain_poly).difference(base)
        rest = water.difference(ln.buffer(rad))
        pieces = [p for p in getattr(rest, "geoms", [rest])
                  if not p.is_empty and p.area > 1e-12]
        if not pieces:
            return False
        probes = []
        for k, other in ((0, 1), (-1, -2)):
            u = a[k] - a[other]
            n_u = float(np.hypot(*u))
            if n_u < 1e-12:
                return False
            probes.append(shapely.Point(
                a[k] + u / n_u * (rad + 0.5 * h_mesh_m / scale)))
        hit = []
        for pr in probes:
            d = [pr.distance(p) for p in pieces]
            j = int(np.argmin(d))
            if d[j] * scale > 0.7 * h_mesh_m:
                return False          # ends in land: no passage
            hit.append(j)
        return hit[0] != hit[1]

    for rec in records:
        if rec["action"] == "ignore":
            info["ignored"] += 1
        elif rec["action"] == "close":
            fills.append(rec["geometry"])
            info["closed"] += 1
        elif rec["action"] == "keep":
            # PER-BRANCH keep/close (owner 2026-07-12, G8-d4/e4
            # fabrication + H8-c5 closure): a chained network
            # mixes wide canal reaches with sub-resolution
            # ditches, and any network-level scalar misjudges
            # both. Skeletonise ALWAYS, measure EACH branch, carve
            # the resolvable ones, leave the rest to the stub
            # fill.
            br_m = sum(float(tup[0].length) for tup in
                       rec.get("bridges") or []) * scale
            done = []
            fails = []
            # ARC SOURCE 1 -- OSM waterway centrelines (owner
            # 2026-07-13, RiverMapper-style): real centreline
            # data beats a Voronoi skeleton wherever it exists.
            # Clip lines to the network, merge, keep pieces
            # longer than 0.7 h as branches.
            brs = []
            n_line_br = 0
            geom_b = rec["geometry"].buffer(0.3 * h_mesh_m / scale)
            if wtree_g is not None:
                from shapely.ops import linemerge
                segs = []
                for k in wtree_g.query(geom_b):
                    cut = wlist[int(k)].intersection(geom_b)
                    for ls in getattr(cut, "geoms", [cut]):
                        if (not ls.is_empty
                                and ls.geom_type == "LineString"
                                and ls.length > 0):
                            segs.append(ls)
                if segs:
                    u_segs = unary_union(segs)
                    merged = (linemerge(u_segs)
                              if u_segs.geom_type
                              == "MultiLineString" else u_segs)
                    for ls in getattr(merged, "geoms", [merged]):
                        if (ls.geom_type == "LineString"
                                and ls.length * scale
                                >= 0.7 * h_mesh_m):
                            brs.append(np.asarray(ls.coords,
                                                  dtype=float))
                    n_line_br = len(brs)
            # ARC SOURCE 2 -- Voronoi medial-axis skeleton over
            # the FULL network geometry, always. Pre-subtracting
            # the centreline tubes left a hole whenever a line
            # branch failed to carve (Chiba severance, run
            # 6188791); a skeleton branch over an already-widened
            # channel is a harmless near-no-op, so overlap is the
            # safe choice.
            try:
                brs += skeleton_branches(
                    rec["geometry"].buffer(0),
                    metric_scale=metric_scale,
                    density_m=0.15 * h_mesh_m,
                    prune_m=0.7 * h_mesh_m)
            except (RuntimeError, ValueError) as e2:
                fails.append(str(e2)[:90])
            info["line_branches"] += n_line_br
            info["skel_branches"] += len(brs) - n_line_br
            rec["n_line_branches"] = n_line_br
            snaps = []
            for br in brs:
                try:
                    snaps.append(snap_arc_to_channel(
                        land_union, br,
                        metric_scale=metric_scale,
                        step_m=0.35 * h_mesh_m,
                        max_halfwidth_m=2.5 * h_mesh_m))
                except (RuntimeError, ValueError) as e2:
                    snaps.append(None)
                    fails.append(str(e2)[:90])
            wmeds = []
            for sp in snaps:
                if sp is None:
                    wmeds.append(-1.0)
                    continue
                wm = sp["width_m"]
                sat = wm >= 0.95 * 2.0 * 2.5 * h_mesh_m
                wn = wm[~sat] if bool((~sat).any()) else wm
                wmeds.append(float(np.median(wn)))
            floor_w = branch_floor_frac * h_mesh_m
            wide_ix = {i for i, wv in enumerate(wmeds)
                       if wv >= floor_w}

            # a SHORT sub-floor branch joining two wide branches
            # is a pinch OF the wide waterway: widening it is
            # mandated (a wide channel must never be split)
            def _endpts(i):
                return (tuple(np.round(brs[i][0], 6)),
                        tuple(np.round(brs[i][-1], 6)))

            jmap: dict[Any, set] = {}
            for i in range(len(brs)):
                for pnt in _endpts(i):
                    jmap.setdefault(pnt, set()).add(i)
            for i in range(len(brs)):
                if i in wide_ix or snaps[i] is None:
                    continue
                if (shapely.LineString(brs[i]).length
                        * scale > 4.0 * h_mesh_m):
                    continue
                touches = set()
                for pnt in _endpts(i):
                    touches |= jmap.get(pnt, set()) & wide_ix
                if len(touches - {i}) >= 2:
                    wide_ix.add(i)
            def _do_branch(i):
                nonlocal new_land
                try:
                    new_land, ci, w_used, ach = _carve(
                        new_land, snaps[i]["arc"],
                        snaps[i]["width_carve_m"],
                        tol_extra_m=br_m + 30.0)
                    done.append((snaps[i]["arc"], w_used))
                    info["land_removed_m2"] += (
                        ci["land_removed_m2"])
                    # achieved-width profile feeds the channel
                    # REFINEMENT corridors (owner 2026-07-13:
                    # finer cells allowed where CFL is untouched)
                    info["refine_arcs"].append(
                        (np.asarray(ach["arc"], float),
                         np.asarray(ach["width_m"], float)))
                    # FORCED two-row ladder (owner 2026-07-13):
                    # EXPERIMENTAL, default OFF. Bare rows made
                    # 27 C1 violations (run 6188830); coupling
                    # the size field helped but 17 remained
                    # (6188855) -- curved bands need end tapering
                    # and junction-aware trimming before this can
                    # meet the quality gates. The machinery is
                    # kept for that follow-up.
                    if force_two_rows:
                        wa2 = ach["width_m"]
                        sat2 = wa2 >= 0.95 * 2.0 * 2.5 * h_mesh_m
                        wn2 = wa2[~sat2] if bool(
                            (~sat2).any()) else wa2
                        wmed2 = float(np.median(wn2))
                        if 1.4 * h_mesh_m <= wmed2 \
                                <= 2.3 * h_mesh_m:
                            pf2, eg2 = _ladder_constraints(
                                ach, h_mesh_m, metric_scale)
                            if len(pf2):
                                info["band_pfix"].append(pf2)
                                info["band_egfix"].append(
                                    eg2 + info["band_n"])
                                info["band_n"] += len(pf2)
                                a2m = max(0.8 * h_mesh_m,
                                          wmed2 / 1.8)
                                bp2, bt2 = _ladder_size_targets(
                                    ach, a2m)
                                info["band_size_pts"].append(bp2)
                                info["band_size_tgt"].append(bt2)
                except (RuntimeError, ValueError) as e2:
                    # CONNECTIVITY OVERRIDE: a branch that fails
                    # the attainability bar but is the only local
                    # connection must still be carved -- as wide
                    # as the barriers allow. Its residual narrow
                    # cells surface in the one-wide ledger.
                    kept_marginal = False
                    if "not attainable" in str(e2):
                        try:
                            if _closure_severs(
                                    new_land, snaps[i]["arc"],
                                    snaps[i]["width_carve_m"]):
                                new_land, ci, w_used, ach = \
                                    _carve(
                                        new_land,
                                        snaps[i]["arc"],
                                        snaps[i]["width_carve_m"],
                                        tol_extra_m=br_m + 30.0,
                                        attain_bar_h=0.0)
                                done.append((snaps[i]["arc"],
                                             w_used))
                                info["land_removed_m2"] += (
                                    ci["land_removed_m2"])
                                info["refine_arcs"].append((
                                    np.asarray(ach["arc"],
                                               float),
                                    np.asarray(ach["width_m"],
                                               float)))
                                info["marginal_kept"] += 1
                                rec.setdefault(
                                    "marginal_branches",
                                    []).append(
                                    "connectivity-critical, "
                                    + str(e2)[:80])
                                kept_marginal = True
                        except (RuntimeError, ValueError) as e3:
                            fails.append(
                                "critical-but-uncarvable: "
                                + str(e3)[:90])
                            kept_marginal = True
                    if not kept_marginal:
                        fails.append(str(e2)[:110])

            # NEAR-DUPLICATE BRANCH SKIP with VERIFICATION (runs
            # 6190495/6190543): the overlap hybrid keeps
            # centreline AND skeleton branches of the same
            # channel; snapped tens of metres apart they carve
            # two offset corridors whose union has jagged banks
            # (all 6 G8-d1/e4 violations sat on doubled rec#110
            # arcs). A branch mostly inside corridors ALREADY
            # CARVED here is deferred -- then MEASURED: if the
            # water along it is genuinely ~two rows wide it was a
            # true duplicate (skip); if not, the coverage test
            # caught a real parallel channel (bare 0.35w/70%
            # skip made 91 skips and new chokes OW19-OW23, run
            # 6190543) and it is carved after all. Coverage is
            # tested against DONE arcs only, so a failed
            # centreline branch never suppresses its skeleton
            # twin (the pre-subtraction Chiba severance,
            # run 6188791).
            deferred_dup = []
            for i in sorted(wide_ix):
                if snaps[i] is None:
                    continue
                if done:
                    ln_i = shapely.LineString(
                        np.asarray(snaps[i]["arc"], float))
                    tube_d = unary_union([
                        shapely.LineString(
                            np.asarray(a2, float)).buffer(
                            0.35 * float(np.median(w2)) / scale)
                        for a2, w2 in done])
                    if (ln_i.intersection(tube_d).length
                            > 0.7 * ln_i.length):
                        deferred_dup.append(i)
                        continue
                _do_branch(i)
            for i in deferred_dup:
                dup_ok = False
                try:
                    achv = snap_arc_to_channel(
                        new_land,
                        np.asarray(snaps[i]["arc"], float),
                        metric_scale=metric_scale,
                        step_m=0.5 * h_mesh_m,
                        max_halfwidth_m=2.5 * h_mesh_m)
                    wv = achv["width_m"]
                    satv = wv >= 0.95 * 2.0 * 2.5 * h_mesh_m
                    wnv = wv[~satv] if bool((~satv).any()) else wv
                    dup_ok = (float(np.median(wnv))
                              >= 1.8 * h_mesh_m)
                except (RuntimeError, ValueError):
                    dup_ok = False    # unmeasurable: carve it
                if dup_ok:
                    info["dup_skipped"] += 1
                else:
                    _do_branch(i)
            if done:
                li = int(np.argmax(
                    [shapely.LineString(a2).length
                     for a2, _ in done]))
                rec["arc"] = done[li][0]
                rec["width_m"] = done[li][1]
                rec["w_used"] = done[li][1]
                rec["arcs_done"] = done
                rec["branches"] = (
                    f"{len(done)}/{len(brs)} carved, "
                    f"{len(brs) - len(wide_ix)} below floor")
                if fails:
                    rec["branch_failures"] = fails
                info["kept"] += 1
            else:
                rec["action"] = "blocked"
                rec["reason"] = (
                    ("carve failed on every eligible branch: "
                     + "; ".join(fails[:2])) if fails else
                    "all branches below the resolve floor "
                    "(0.45h) -- do-not-mesh, closed")
                if close_blocked:
                    fills.append(rec["geometry"])
                    rec["closed"] = True
        if (open_bridges and rec["action"] == "keep"
                and rec.get("bridges")):
            # BRIDGE vs LEVEE (owner 2026-07-13): geometry alone
            # cannot tell a road bridge over one canal from a
            # levee between two SEPARATE waters (G8-d4/e4
            # fabrication). "auto" opens a gap ONLY when an OSM
            # WATERWAY CENTRELINE (river/canal/stream line data)
            # passes THROUGH the strip -- flow continuity is then
            # attested by data, not guessed from shape. True
            # forces all open (tests); False disables.
            # open the OSM bridge strips (roads drawn as land
            # across the canal, Daishi-canal severance): a short
            # TRANSVERSAL carve with crossings allowed. Only
            # bridges BETWEEN CARVED branch corridors are opened
            # -- a bridge toward an unresolved ditch would
            # fabricate a connection (G8-d4/e4).
            tub = unary_union([
                shapely.LineString(np.asarray(a2, float)).buffer(
                    0.6 * float(np.max(w2)) / scale)
                for a2, w2 in rec.get("arcs_done") or []])
            wtree = wtree_g if open_bridges == "auto" else None

            def _authorized(ln, gi, gj):
                """A waterway centreline must run through BOTH
                flanking water pieces AND pass near this gap --
                data-attested flow continuity."""
                if open_bridges is True:
                    return True
                if wtree is None:
                    return False
                mid = ln.interpolate(0.5, normalized=True)
                near = mid.buffer(2.0 * h_mesh_m / scale)
                for k in wtree.query(near):
                    L = wlist[int(k)]
                    if (L.intersects(near)
                            and L.intersects(gi)
                            and L.intersects(gj)):
                        return True
                return False

            opened = 0
            for ln, gi_b, gj_b in rec["bridges"]:
                if not _authorized(ln, gi_b, gj_b):
                    continue
                if not (tub.intersects(shapely.Point(
                            ln.coords[0]).buffer(
                            0.3 * h_mesh_m / scale))
                        and tub.intersects(shapely.Point(
                            ln.coords[-1]).buffer(
                            0.3 * h_mesh_m / scale))):
                    continue
                pts_b = np.asarray(ln.coords, float)
                d = pts_b[-1] - pts_b[0]
                if np.hypot(*d) < 1e-12:
                    continue
                u = d / np.hypot(*d)
                ext = 0.4 * h_mesh_m / scale
                arc_b = np.vstack([pts_b[0] - u * ext,
                                   pts_b[-1] + u * ext])
                try:
                    new_land, cb = carve_channel_corridor(
                        new_land, arc_b,
                        0.875 * widen_rows * h_mesh_m,
                        min_gap_m=min_gap_m,
                        metric_scale=metric_scale,
                        domain_poly=domain_poly,
                        arc_on_land_tol_m=float(ln.length)
                        * scale + h_mesh_m,
                        carve_crossings=True)
                    opened += 1
                    info["land_removed_m2"] += \
                        cb["land_removed_m2"]
                except RuntimeError as e3:
                    rec.setdefault("bridge_failures",
                                   []).append(str(e3))
            rec["bridges_opened"] = opened
            info["bridges_opened"] += opened
        if rec["action"] == "keep" and rec.get("w_used") is not None:
            # STUB HEADS (owner do-not-mesh rule applied to kept
            # networks): residual network water OUTSIDE the arc
            # corridor that stays narrower than 0.5 h cannot be
            # widened (the arc does not reach it) and would mesh
            # one cell wide or terminate the mesh in a jagged
            # wedge (element 4797, run 6185775). Fill it.
            pairs_aw = rec.get("arcs_done") or [
                (rec["arc"], rec["w_used"])]
            tube = unary_union([
                shapely.LineString(np.asarray(a2, float)).buffer(
                    0.55 * float(np.max(w2)) / scale)
                for a2, w2 in pairs_aw])
            n_stub = 0
            for pz in _polys(rec["geometry"].difference(tube)):
                wz = 2.0 * pz.area / max(pz.boundary.length,
                                         1e-9) * scale
                if wz < 0.5 * h_mesh_m:
                    fills.append(pz)
                    n_stub += 1
            if n_stub:
                rec["stub_fills"] = n_stub
                info["stub_fills"] += n_stub
        if rec["action"] == "blocked":
            info["blocked"].append(rec)
    if fills:
        # NOTE: do NOT smooth/dilate the fills. A closing pass
        # (buffer +0.15h/-0.1h) interacted catastrophically with
        # the Shoreline stage (623 cells meshed over Tokyo city,
        # one real severance -- run 6186580). Raw fills leave a
        # 1-2 cell realization-sensitive quality tail at the
        # artificial west edge instead, tracked in the ledger.
        new_land = unary_union([new_land, *fills])
    # LAND-CRUMB CLEANUP (run 6190495): corridor carving can
    # shave sub-cell land fragments off a bigger polygon (a
    # 290 m2 speck at G8-e4 fed C1 slivers). Drop pieces far
    # below meshable size -- UNLESS the piece already existed as
    # its own island in the ORIGINAL land (real islets are data,
    # never deleted).
    crumb_max_m2 = (0.25 * h_mesh_m) ** 2
    pieces = _polys(new_land)
    small = [p for p in pieces if p.area * sx * sy < crumb_max_m2]
    if small:
        from shapely.strtree import STRtree
        orig = _polys(land_union)
        otree = STRtree(orig)
        kept_pieces = [p for p in pieces
                       if p.area * sx * sy >= crumb_max_m2]
        n_crumb = 0
        for p in small:
            islet = False
            for k in otree.query(p):
                o = orig[k]
                if (o.area <= 2.0 * p.area and
                        o.intersection(p).area > 0.6 * p.area):
                    islet = True
                    break
            if islet:
                kept_pieces.append(p)
            else:
                n_crumb += 1
        if n_crumb:
            info["crumbs_dropped"] = n_crumb
            new_land = unary_union(kept_pieces)
    return new_land, info
