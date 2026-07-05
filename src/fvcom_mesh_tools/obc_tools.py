"""Open-boundary construction on the outer ring (package version of
the PoC 6x-series helpers).

The OBC is an ARTIFICIAL smooth line (user requirement): nodes
assigned to it are snapped exactly onto straight segments along the
domain's open-sea edges, junctions trimmed, and the FVCOM
land/open boundary lists rebuilt. Perpendicularity is then enforced
by the local first-ring fixer (interior nodes only).
"""

from __future__ import annotations

from typing import Any

import numpy as np
from pyproj import Transformer

from fvcom_mesh_tools.algorithms import (
    align_open_boundary_local,
    snap_nodes_to_segment,
)
from fvcom_mesh_tools.algorithms.boundary import (
    boundary_edges_from_tris,
    chain_edges_to_loops,
    outer_loop,
)
from fvcom_mesh_tools.io import Fort14Mesh

__all__ = ["assign_west_south_obc"]


def assign_west_south_obc(
    mesh: Fort14Mesh,
    *,
    utm_epsg: int = 32654,
    band_deg: float = 0.012,
    shoreline_shp=None,
    coast_tol_m: float = 500.0,
    trim: int = 1,
    max_move_m: float = 600.0,
    land_ibtype: int = 20,
    perp_seed: int = 9500,
    min_depth_m: float | None = None,
    snap: bool = True,
    obc_line_lonlat=None,
    log=print,
) -> tuple[Fort14Mesh, dict[str, Any]]:
    """Detect the west+south open-sea edges of the outer ring, snap
    each run exactly onto its straight chord, set the FVCOM open/land
    boundary lists, run the local perpendicularity fixer, and
    (optionally) clip depths at ``min_depth_m``.

    ``mesh`` is in UTM metres; returns a new mesh in the same CRS.

    ``min_depth_m`` is OFF by default: depth-field construction is
    out of scope for mesh generation (the depths carried here are
    buildmesh's provisional DEM interpolation, not a designed
    product). When set, it clamps VALUES only — node positions and
    the coastline are untouched.
    """
    to_ll = Transformer.from_crs(f"EPSG:{utm_epsg}", "EPSG:4326",
                                 always_xy=True)
    lon, lat = to_ll.transform(mesh.nodes[:, 0], mesh.nodes[:, 1])

    loops = chain_edges_to_loops(boundary_edges_from_tris(mesh.elements))
    outer = outer_loop(loops, mesh.nodes)
    ring = outer[:-1]
    rlon, rlat = lon[ring], lat[ring]
    if obc_line_lonlat is not None:
        # The OBC line is explicit: candidates are simply the ring
        # nodes NEAR the arc (the wall polygon that cuts the domain
        # at the arc is itself part of the engineered shoreline, so
        # distance-to-shoreline cannot see the arc: PoC #105 v2
        # found no open-sea run at all).
        import shapely
        from shapely.geometry import LineString

        to_m2 = Transformer.from_crs("EPSG:4326", f"EPSG:{utm_epsg}",
                                     always_xy=True)
        ax, ay = to_m2.transform(
            [q[0] for q in obc_line_lonlat],
            [q[1] for q in obc_line_lonlat],
        )
        arc0 = LineString(list(zip(ax, ay)))
        pts0 = shapely.points(mesh.nodes[ring, 0], mesh.nodes[ring, 1])
        d_arc = shapely.distance(pts0, arc0)
        # The wall's arc edge is CDT-constrained, so true arc nodes
        # sit within lattice noise of the line; a wide corridor
        # (2*max_move) swallowed coastline nodes near the junctions
        # and drew a staircase OBC up the Miura shore. Require BOTH
        # a tight corridor and an interior projection (clamped-end
        # projections = coast beyond the junction).
        s_arc = np.array([
            arc0.project(shapely.Point(mesh.nodes[v, 0],
                                       mesh.nodes[v, 1]))
            for v in ring
        ])
        # Corridor: pre-snap boundary nodes sit up to several
        # hundred metres off the line (DistMesh lattice + cleanup),
        # so use max_move_m; the staircase came from clamped-END
        # projections, which the interior filter removes.
        eps = 30.0
        mask = ((d_arc < max_move_m)
                & (s_arc > eps) & (s_arc < arc0.length - eps))
        south_b = mask
        west_b = np.zeros_like(mask)
    elif shoreline_shp is not None:
        # Open-sea edge = outer-ring nodes FAR from the engineered
        # shoreline (the lat/lon band heuristic mixes in coastline
        # nodes on non-rectangular domains — PoC #96: west run 0/107
        # snapped, perpendicularity worst 88 deg).
        import shapely
        from shapely.strtree import STRtree

        from fvcom_mesh_tools.algorithms.boundary_snap import (
            load_polylines,
        )

        lines_utm = load_polylines(shoreline_shp, to_crs=utm_epsg)
        tree = STRtree(lines_utm)
        arr = np.array(lines_utm, dtype=object)
        pts = shapely.points(mesh.nodes[ring, 0], mesh.nodes[ring, 1])
        d_coast = shapely.distance(pts, arr[tree.nearest(pts)])
        mask = d_coast > coast_tol_m
        # west/south attribution decided later by chord splitting.
        south_b = mask
        west_b = np.zeros_like(mask)
    else:
        south_b = rlat <= rlat.min() + band_deg
        west_b = rlon <= rlon.min() + band_deg
        mask = south_b | west_b
    idx = np.where(mask)[0]
    if idx.size < 4:
        raise ValueError("no open-sea run found on the outer ring")
    runs = []
    s = p = int(idx[0])
    for q in idx[1:]:
        q = int(q)
        if q == p + 1:
            p = q
        else:
            runs.append((s, p))
            s = p = q
    runs.append((s, p))
    # Merge the wrap-around run (ring is cyclic).
    if len(runs) > 1 and runs[0][0] == 0 and runs[-1][1] == len(ring) - 1:
        s0, e0 = runs.pop(0)
        s1, e1 = runs.pop(-1)
        runs.append((s1 - len(ring), e0))
    a, b = max(runs, key=lambda r: r[1] - r[0])
    arc_pos = np.arange(a, b + 1) % len(ring)
    arc_nodes = ring[arc_pos]
    if obc_line_lonlat is not None:
        parts = [arc_nodes]
        labels = ["arc"]
    elif shoreline_shp is not None:
        # Split the arc at its corner: the node farthest from the
        # end-to-end chord (if its deviation is significant) — the
        # west/south rectangle corner on box-cut domains.
        P = mesh.nodes[arc_nodes]
        p0v, p1v = P[0], P[-1]
        chord = p1v - p0v
        nrm = np.linalg.norm(chord) or 1.0
        dev = np.abs(
            (P[:, 0] - p0v[0]) * chord[1]
            - (P[:, 1] - p0v[1]) * chord[0]
        ) / nrm
        k_corner = int(np.argmax(dev))
        if dev[k_corner] > 2.0 * max_move_m and 2 <= k_corner \
                <= len(arc_nodes) - 3:
            parts = [arc_nodes[: k_corner + 1],
                     arc_nodes[k_corner:]]
        else:
            parts = [arc_nodes]
        labels = ["seg1", "seg2"][: len(parts)]
    else:
        parts = [arc_nodes[west_b[arc_pos]],
                 arc_nodes[south_b[arc_pos] & ~west_b[arc_pos]]]
        labels = ["west", "south"]
    open_len = (b - a)

    info: dict[str, Any] = {"n_arc": int(open_len + 1)}

    # --- SNAP ONCE: nodes move only here. Membership is tracked by
    # COORDINATES: element deletion + keep_components renumber the
    # nodes, so id sets go stale (the n_obc=1 collapse). ---------------
    def _key(xy):
        return (round(float(xy[0]), 3), round(float(xy[1]), 3))

    snapped_keys: set = set()
    snapped_ids: set[int] = set()
    if obc_line_lonlat is not None and snap:
        # Project every open-sea node onto the artificial ARC line
        # (goto2023-style curved OBC) instead of straight chords.
        import shapely
        from shapely.geometry import LineString

        to_m = Transformer.from_crs("EPSG:4326", f"EPSG:{utm_epsg}",
                                    always_xy=True)
        lx, ly = to_m.transform(
            [q[0] for q in obc_line_lonlat],
            [q[1] for q in obc_line_lonlat],
        )
        arc = LineString(list(zip(lx, ly)))
        all_arc = np.concatenate(parts) if len(parts) else np.array([])
        n_arc_snap = 0
        els_a = mesh.elements

        def _ring_ok(v2):
            we = np.where((els_a == v2).any(axis=1))[0]
            tri = els_a[we]
            p0 = mesh.nodes[tri[:, 0]]
            p1 = mesh.nodes[tri[:, 1]]
            p2 = mesh.nodes[tri[:, 2]]
            a2 = ((p1[:, 0] - p0[:, 0]) * (p2[:, 1] - p0[:, 1])
                  - (p1[:, 1] - p0[:, 1]) * (p2[:, 0] - p0[:, 0]))
            return bool((a2 > 0).all())

        for v in all_arc:
            v = int(v)
            q = arc.interpolate(arc.project(
                shapely.Point(mesh.nodes[v, 0], mesh.nodes[v, 1])
            ))
            move = float(np.hypot(q.x - mesh.nodes[v, 0],
                                  q.y - mesh.nodes[v, 1]))
            if move <= max_move_m:
                old_pos = mesh.nodes[v].copy()
                mesh.nodes[v] = (q.x, q.y)
                if not _ring_ok(v):
                    mesh.nodes[v] = old_pos
                else:
                    n_arc_snap += 1
            snapped_keys.add(_key(mesh.nodes[v]))
            snapped_ids.add(v)
        info["arc_snap"] = {"n": int(n_arc_snap),
                            "of": int(all_arc.size)}
        log(f"[obc] arc snap: {n_arc_snap}/{all_arc.size}")
        parts = []
    for run_nodes, label in zip(parts, labels):
        if run_nodes.size >= 2:
            if snap:
                p0 = tuple(mesh.nodes[run_nodes[0]])
                p1 = tuple(mesh.nodes[run_nodes[-1]])
                mesh, li = snap_nodes_to_segment(
                    mesh, [int(v) for v in run_nodes], p0, p1,
                    max_move=max_move_m,
                )
                info[f"line_{label}"] = li
                log(f"[obc] {label} line: {li}")
            for v in run_nodes:
                snapped_keys.add(_key(mesh.nodes[int(v)]))
                snapped_ids.add(int(v))

    # Gated relax of interior neighbours of the snapped runs: the
    # chord snap shears its 1-ring exactly like coastline snapping
    # did; absorb it HERE because no stage after obc may move nodes.
    from fvcom_mesh_tools.algorithms.perp_local import _tri_quality

    if not snap:
        # finalize mode: lists + perp only, zero node motion.
        info["n_relaxed"] = 0
        els = mesh.elements
        touch = np.zeros(els.shape[0], dtype=bool)
    else:
        els = mesh.elements
        touch = np.isin(els, list(snapped_ids)).any(axis=1)
    ring_e = np.where(touch)[0]
    bnd_uv = boundary_edges_from_tris(els)
    bnd_nodes = set(int(x) for e in bnd_uv for x in e)
    interior = [int(w) for w in np.unique(els[ring_e].ravel())
                if int(w) not in bnd_nodes]
    n_relax = 0
    for _sweep in range(2):
        for w in interior:
            we = np.where((els == w).any(axis=1))[0]
            nbrs = sorted({int(x) for e in we for x in els[e]
                           if int(x) != w})
            if len(nbrs) < 3:
                continue
            tri_w = els[we]
            mn0, mx0, tw0 = _tri_quality(mesh.nodes[tri_w])
            bad0 = int((mn0 < 30).sum() + (mx0 > 130).sum()
                       + (tw0 <= 0).sum())
            old_pos = mesh.nodes[w].copy()
            mesh.nodes[w] = mesh.nodes[nbrs].mean(axis=0)
            mn1, mx1, tw1 = _tri_quality(mesh.nodes[tri_w])
            bad1 = int((mn1 < 30).sum() + (mx1 > 130).sum()
                       + (tw1 <= 0).sum())
            if bad1 > bad0 or (tw1 <= 0).any():
                mesh.nodes[w] = old_pos
            else:
                n_relax += 1
    info["n_relaxed"] = n_relax
    log(f"[obc] gated relax around snapped runs: {n_relax} accepts")

    # --- Assignment + bounded structural aftercare: deletions and
    # boundary-list rebuilds ONLY — no further node motion (a
    # re-snapping recursion here dragged nodes between shifting
    # chords and wrecked quality; see the #98 v3 log). -----------------
    from fvcom_mesh_tools.mesh_clean import keep_components
    from fvcom_mesh_tools.mesh_clean import remove_elements as _rm
    from fvcom_mesh_tools.qa import fvcom_boundary_element_flags

    def _rebuild(mesh_in: Fort14Mesh) -> Fort14Mesh:
        loops2 = chain_edges_to_loops(
            boundary_edges_from_tris(mesh_in.elements)
        )
        outer2 = outer_loop(loops2, mesh_in.nodes)
        ring2 = outer2[:-1]
        member = np.array([
            _key(mesh_in.nodes[int(v)]) in snapped_keys for v in ring2
        ])
        idx2 = np.where(member)[0]
        if idx2.size < 4:
            raise ValueError("open-sea run lost during cleanup")
        runs2 = []
        s2 = p2 = int(idx2[0])
        for q2 in idx2[1:]:
            q2 = int(q2)
            if q2 == p2 + 1:
                p2 = q2
            else:
                runs2.append((s2, p2))
                s2 = p2 = q2
        runs2.append((s2, p2))
        if len(runs2) > 1 and runs2[0][0] == 0 \
                and runs2[-1][1] == len(ring2) - 1:
            s0, e0 = runs2.pop(0)
            s1, e1 = runs2.pop(-1)
            runs2.append((s1 - len(ring2), e0))
        a2, b2 = max(runs2, key=lambda r: r[1] - r[0])
        ring2 = np.roll(ring2, -(a2 % len(ring2)))
        open_len2 = b2 - a2
        lo2, hi2 = trim, open_len2 - trim
        open_seg2 = ring2[lo2:hi2 + 1].copy()
        land_seg2 = np.concatenate([ring2[hi2:], ring2[:lo2 + 1]])
        islands2 = [lp[:-1].copy() for lp in loops2
                    if lp is not outer2]
        return Fort14Mesh(
            title=mesh_in.title, nodes=mesh_in.nodes,
            depths=mesh_in.depths, elements=mesh_in.elements,
            open_boundaries=[open_seg2],
            land_boundaries=[(land_ibtype, land_seg2)]
            + [(land_ibtype, i2) for i2 in islands2],
        )

    mesh = _rebuild(mesh)
    n_deleted = 0
    from fvcom_mesh_tools.qa import _edge_topology

    def _pinch_elements(mesh_in):
        topo = _edge_topology(mesh_in.elements, mesh_in.n_nodes)
        buv = topo.uv[topo.counts == 1]
        cnt = np.zeros(mesh_in.n_nodes, dtype=np.int64)
        if buv.size:
            np.add.at(cnt, buv.ravel(), 1)
        pinch = np.where(cnt > 2)[0]
        bad_p = np.zeros(mesh_in.n_elements, dtype=bool)
        if pinch.size:
            bad_p |= np.isin(mesh_in.elements, pinch).any(axis=1)
        for u2, v2 in topo.uv[topo.counts > 2]:
            bad_p |= ((mesh_in.elements == u2).any(axis=1)
                      & (mesh_in.elements == v2).any(axis=1))
        return bad_p

    def _junction_closure(mesh_in):
        """Extend OBC membership along the outer ring to the coast
        (aftercare deletions retreat the tail; unmarked boundary
        edges across the water strip act as a WALL in FVCOM)."""
        if obc_line_lonlat is None:
            return mesh_in, 0
        import shapely
        from shapely.geometry import LineString

        to_m3 = Transformer.from_crs("EPSG:4326", f"EPSG:{utm_epsg}",
                                     always_xy=True)
        cx, cy = to_m3.transform(
            [q[0] for q in obc_line_lonlat],
            [q[1] for q in obc_line_lonlat],
        )
        arc3 = LineString(list(zip(cx, cy)))
        topo3 = _edge_topology(mesh_in.elements, mesh_in.n_nodes)
        buv = topo3.uv[topo3.counts == 1]
        nbr: dict[int, list[int]] = {}
        for u4, v4 in buv:
            nbr.setdefault(int(u4), []).append(int(v4))
            nbr.setdefault(int(v4), []).append(int(u4))

        def _ring_ok3(v2):
            we = np.where((mesh_in.elements == v2).any(axis=1))[0]
            tri = mesh_in.elements[we]
            p0 = mesh_in.nodes[tri[:, 0]]
            p1 = mesh_in.nodes[tri[:, 1]]
            p2 = mesh_in.nodes[tri[:, 2]]
            a2 = ((p1[:, 0] - p0[:, 0]) * (p2[:, 1] - p0[:, 1])
                  - (p1[:, 1] - p0[:, 1]) * (p2[:, 0] - p0[:, 0]))
            return bool((a2 > 0).all())

        obc_list = [int(v) for v in mesh_in.open_boundaries[0]]
        n_closed = 0
        for end_i in (0, -1):
            if len(obc_list) < 2:
                break
            cur = obc_list[end_i]
            prev = obc_list[end_i + 1 if end_i == 0 else end_i - 1]
            for _ in range(12):
                nxt = [w for w in nbr.get(cur, [])
                       if w != prev and w not in obc_list]
                if len(nxt) != 1:
                    break
                w = nxt[0]
                pt = shapely.Point(mesh_in.nodes[w, 0],
                                   mesh_in.nodes[w, 1])
                d = float(shapely.distance(pt, arc3))
                s = float(arc3.project(pt))
                interior = 30.0 < s < arc3.length - 30.0
                if not ((interior and d <= max_move_m)
                        or (not interior and d <= 120.0)):
                    break
                q = arc3.interpolate(s)
                old_pos = mesh_in.nodes[w].copy()
                mesh_in.nodes[w] = (q.x, q.y)
                if not _ring_ok3(w):
                    mesh_in.nodes[w] = old_pos
                if end_i == 0:
                    obc_list.insert(0, w)
                else:
                    obc_list.append(w)
                snapped_keys.add(_key(mesh_in.nodes[w]))
                snapped_ids.add(int(w))
                n_closed += 1
                prev, cur = cur, w
        if n_closed:
            bnd = mesh_in.open_boundaries
            bnd[0] = np.array(obc_list, dtype=bnd[0].dtype)
        return mesh_in, n_closed

    def _derive_arc_membership(mesh_in):
        """Re-derive the OBC node list PURELY from geometry: outer-
        ring nodes in the line corridor with interior projections
        (plus near-line clamped ends), ordered by arc length.
        Incremental closure walks + coordinate-keyed rebuilds go
        stale as soon as siteops/polish move nodes (review23: list
        collapsed 21 -> 8 and the tails to both coasts were left
        unmarked = artificial walls)."""
        if obc_line_lonlat is None:
            return mesh_in, 0
        import shapely
        from shapely.geometry import LineString

        to_m5 = Transformer.from_crs("EPSG:4326", f"EPSG:{utm_epsg}",
                                     always_xy=True)
        ax5, ay5 = to_m5.transform(
            [q[0] for q in obc_line_lonlat],
            [q[1] for q in obc_line_lonlat],
        )
        arc5 = LineString(list(zip(ax5, ay5)))
        loops5 = chain_edges_to_loops(
            boundary_edges_from_tris(mesh_in.elements)
        )
        ring5 = outer_loop(loops5, mesh_in.nodes)[:-1]
        pts5 = shapely.points(mesh_in.nodes[ring5, 0],
                              mesh_in.nodes[ring5, 1])
        d5 = shapely.distance(pts5, arc5)
        s5 = np.array([
            arc5.project(shapely.Point(*mesh_in.nodes[v]))
            for v in ring5
        ])
        interior5 = (s5 > 30.0) & (s5 < arc5.length - 30.0)
        keep5 = ((interior5 & (d5 <= max_move_m))
                 | (~interior5 & (d5 <= 120.0)))
        kidx = np.where(keep5)[0]
        if kidx.size < 2:
            return mesh_in, 0
        # OBC = the contiguous RING PATH between the two arc-length
        # extreme members (guarantees boundary-edge adjacency; s-
        # sorted member lists violated obc_ordering where corridor
        # gaps left off-chain stragglers, review24). Direction: the
        # ring path that contains the most members.
        i0 = int(kidx[np.argmin(s5[kidx])])
        i1 = int(kidx[np.argmax(s5[kidx])])
        n5 = len(ring5)
        fwd = np.arange(i0, i0 + (i1 - i0) % n5 + 1) % n5
        bwd = np.arange(i1, i1 + (i0 - i1) % n5 + 1) % n5
        kset = set(kidx.tolist())
        n_fwd = sum(1 for j in fwd if j in kset)
        n_bwd = sum(1 for j in bwd if j in kset)
        path = fwd if (n_fwd, -len(fwd)) >= (n_bwd, -len(bwd)) else bwd
        if len(path) > 3 * kidx.size + 8:
            # degenerate direction pick: fall back to the other side
            path = bwd if path is fwd else fwd
        members = ring5[path]
        if members.size < 2:
            return mesh_in, 0
        # orient by increasing arc length for a stable list
        if s5[path[0]] > s5[path[-1]]:
            members = members[::-1]
        bnd5 = mesh_in.open_boundaries
        bnd5[0] = members.astype(bnd5[0].dtype)
        for v in members:
            snapped_keys.add(_key(mesh_in.nodes[int(v)]))
            snapped_ids.add(int(v))
        return mesh_in, int(members.size)

    mesh, n_m1 = _derive_arc_membership(mesh)
    if n_m1:
        log(f"[obc] arc membership (pre-aftercare): {n_m1} nodes")

    n_split = 0
    for _round in range(8):
        flags = fvcom_boundary_element_flags(mesh)
        r4 = flags["r4_mask"]
        if obc_line_lonlat is not None and r4.any():
            # Arc mode: SPLIT R4 elements (1-to-3 centroid insert)
            # instead of deleting them — deletion retreats the
            # boundary off the OBC line (review25: an ~800 m step
            # with an 851 m off-line node the straightener cannot
            # recover). The centroid becomes the interior node that
            # de-R4s the element while the boundary stays put.
            ids = np.where(r4)[0]
            new_tris = []
            keep_mask = np.ones(mesh.n_elements, dtype=bool)
            for e in ids:
                tri = mesh.elements[e]
                cxy = mesh.nodes[tri].mean(axis=0)
                vc = mesh.n_nodes
                mesh.nodes = np.vstack([mesh.nodes, cxy[None, :]])
                mesh.depths = np.append(
                    mesh.depths, float(mesh.depths[tri].mean())
                )
                keep_mask[e] = False
                a1, b1, c1 = (int(tri[0]), int(tri[1]), int(tri[2]))
                new_tris += [[a1, b1, vc], [b1, c1, vc], [c1, a1, vc]]
            mesh.elements = np.vstack(
                [mesh.elements[keep_mask],
                 np.asarray(new_tris, dtype=mesh.elements.dtype)]
            )
            n_split += int(ids.size)
            mesh = _rebuild(mesh)
            continue
        # Pinch nodes break FVCOM's own NBE pairing at setup
        # ("ELEMENT ... HAS NO NEIGHBORS" despite valid adjacency):
        # delete them with the fatal classes.
        bad = r4 | flags["fake_open_mask"] | _pinch_elements(mesh)
        if not bad.any():
            break
        n_deleted += int(bad.sum())
        mesh = _rm(mesh, ~bad)
        mesh, _ = keep_components(mesh)
        mesh = _rebuild(mesh)
    if n_split:
        log(f"[obc] R4 centroid splits: {n_split}")
    info["n_structural_deleted"] = n_deleted
    info["n_obc"] = int(mesh.open_boundaries[0].size)
    log(f"[obc] structural aftercare: deleted {n_deleted}, "
        f"n_obc={info['n_obc']}")

    mesh, n_m2 = _derive_arc_membership(mesh)
    if n_m2:
        log(f"[obc] arc membership (post-aftercare): {n_m2} nodes")
        info["n_obc"] = int(mesh.open_boundaries[0].size)

    if obc_line_lonlat is not None:
        # OBC line straightener (AI manual-editing pass, single
        # sweep): any OBC node still off the effective line (e.g.
        # the Boso touchdown step, user review 2026-07-05) is
        # projected onto it; whitelist-compliant (motion along/onto
        # the OBC line only), flip-guarded, no iteration.
        import shapely as _shp
        from shapely.geometry import LineString as _LS

        to_m4 = Transformer.from_crs("EPSG:4326", f"EPSG:{utm_epsg}",
                                     always_xy=True)
        sx, sy = to_m4.transform(
            [q[0] for q in obc_line_lonlat],
            [q[1] for q in obc_line_lonlat],
        )
        arc4 = _LS(list(zip(sx, sy)))
        els4 = mesh.elements

        def _ring_ok4(v2):
            we = np.where((els4 == v2).any(axis=1))[0]
            tri = els4[we]
            p0 = mesh.nodes[tri[:, 0]]
            p1 = mesh.nodes[tri[:, 1]]
            p2 = mesh.nodes[tri[:, 2]]
            a2 = ((p1[:, 0] - p0[:, 0]) * (p2[:, 1] - p0[:, 1])
                  - (p1[:, 1] - p0[:, 1]) * (p2[:, 0] - p0[:, 0]))
            return bool((a2 > 0).all())

        n_str = 0
        worst = 0.0
        for v in [int(q) for q in mesh.open_boundaries[0]]:
            pt = _shp.Point(mesh.nodes[v, 0], mesh.nodes[v, 1])
            d = float(_shp.distance(pt, arc4))
            if d <= 20.0 or d > max_move_m:
                worst = max(worst, 0.0 if d <= 20.0 else d)
                continue
            q4 = arc4.interpolate(arc4.project(pt))
            old_pos = mesh.nodes[v].copy()
            mesh.nodes[v] = (q4.x, q4.y)
            if _ring_ok4(v):
                n_str += 1
            else:
                mesh.nodes[v] = old_pos
                worst = max(worst, d)
        if n_str or worst:
            log(f"[obc] line straightener: moved {n_str}, "
                f"worst residual {worst:.0f} m")
            info["straightener"] = {"moved": int(n_str),
                                    "worst_residual_m": float(worst)}

    mesh, pinfo = align_open_boundary_local(
        mesh, seed=perp_seed, max_outer=1,
    )
    info["perp"] = {"accepted": pinfo["accepted_total"],
                    "remaining": len(pinfo["remaining"])}
    log(f"[obc] perp: {info['perp']}")

    if min_depth_m is not None:
        n_clip = int((mesh.depths < min_depth_m).sum())
        mesh.depths[:] = np.maximum(mesh.depths, min_depth_m)
        info["depth_clipped"] = n_clip
        log(f"[obc] depth clip @{min_depth_m:g} m: {n_clip} nodes")
    return mesh, info
