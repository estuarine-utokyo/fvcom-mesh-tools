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
    trim: int = 1,
    max_move_m: float = 600.0,
    land_ibtype: int = 20,
    perp_seed: int = 9500,
    min_depth_m: float | None = 2.0,
    log=print,
) -> tuple[Fort14Mesh, dict[str, Any]]:
    """Detect the west+south open-sea edges of the outer ring, snap
    each run exactly onto its straight chord, set the FVCOM open/land
    boundary lists, run the local perpendicularity fixer, and
    (optionally) clip depths at ``min_depth_m``.

    ``mesh`` is in UTM metres; returns a new mesh in the same CRS.
    """
    to_ll = Transformer.from_crs(f"EPSG:{utm_epsg}", "EPSG:4326",
                                 always_xy=True)
    lon, lat = to_ll.transform(mesh.nodes[:, 0], mesh.nodes[:, 1])

    loops = chain_edges_to_loops(boundary_edges_from_tris(mesh.elements))
    outer = outer_loop(loops, mesh.nodes)
    ring = outer[:-1]
    rlon, rlat = lon[ring], lat[ring]
    south_b = rlat <= rlat.min() + band_deg
    west_b = rlon <= rlon.min() + band_deg
    mask = south_b | west_b
    idx = np.where(mask)[0]
    if idx.size < 4:
        raise ValueError("no west/south open-sea band found")
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
    west_nodes = arc_nodes[west_b[arc_pos]]
    south_nodes = arc_nodes[south_b[arc_pos] & ~west_b[arc_pos]]
    ring = np.roll(ring, -(a % len(ring)))
    open_len = (b - a)

    info: dict[str, Any] = {"n_arc": int(open_len + 1)}
    for run_nodes, label in ((west_nodes, "west"), (south_nodes, "south")):
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
