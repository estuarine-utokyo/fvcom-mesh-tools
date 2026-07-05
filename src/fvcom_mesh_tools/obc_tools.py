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
    if shoreline_shp is not None:
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
    if shoreline_shp is not None:
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

    for _round in range(8):
        flags = fvcom_boundary_element_flags(mesh)
        # Pinch nodes break FVCOM's own NBE pairing at setup
        # ("ELEMENT ... HAS NO NEIGHBORS" despite valid adjacency):
        # delete them with the fatal classes.
        bad = flags["r4_mask"] | flags["fake_open_mask"] \
            | _pinch_elements(mesh)
        if not bad.any():
            break
        n_deleted += int(bad.sum())
        mesh = _rm(mesh, ~bad)
        mesh, _ = keep_components(mesh)
        mesh = _rebuild(mesh)
    info["n_structural_deleted"] = n_deleted
    info["n_obc"] = int(mesh.open_boundaries[0].size)
    log(f"[obc] structural aftercare: deleted {n_deleted}, "
        f"n_obc={info['n_obc']}")

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
