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
    ring = np.roll(ring, -(a % len(ring)))
    open_len = (b - a)

    info: dict[str, Any] = {"n_arc": int(open_len + 1)}
    for run_nodes, label in zip(parts, labels):
        if run_nodes.size >= 2:
            p0 = tuple(mesh.nodes[run_nodes[0]])
            p1 = tuple(mesh.nodes[run_nodes[-1]])
            mesh, li = snap_nodes_to_segment(
                mesh, [int(v) for v in run_nodes], p0, p1,
                max_move=max_move_m,
            )
            info[f"line_{label}"] = li
            log(f"[obc] {label} line: {li}")

    lo, hi = trim, open_len - trim
    open_seg = ring[lo:hi + 1].copy()
    land_seg = np.concatenate([ring[hi:], ring[:lo + 1]])
    islands = [lp[:-1].copy() for lp in loops if lp is not outer]
    mesh = Fort14Mesh(
        title=mesh.title, nodes=mesh.nodes, depths=mesh.depths,
        elements=mesh.elements,
        open_boundaries=[open_seg],
        land_boundaries=[(land_ibtype, land_seg)]
        + [(land_ibtype, i) for i in islands],
    )
    info["n_obc"] = int(open_seg.size)

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
