"""Constraint-respecting finishing of CDT-boundary meshes.

Extraction of the PoC #93 stages into a reusable function: the mesh
arrives with its boundary ON the engineered shoreline (CDT egfix)
plus the sliver population the disabled boundary-deleters would have
removed by retreating. Finishing removes the slivers WITHOUT giving
up the boundary:

1. cap-triangle deletion — slivers whose three vertices all lie on
   the mesh boundary; deleting them exposes the constrained edges,
   so conformity is untouched;
2. weld of near-coincident nodes (chain corners);
3. ONE budgeted `phase_h_optimize` round with the shoreline
   projector (boundary nodes may slide ALONG the line only).

No convergence loops: each stage runs once (user discipline,
docs/DESIGN_HISTORY.md).
"""

from __future__ import annotations

from pathlib import Path
from typing import Any, Callable

import numpy as np

from fvcom_mesh_tools.algorithms.perp_local import _tri_quality
from fvcom_mesh_tools.io import Fort14Mesh
from fvcom_mesh_tools.mesh_clean import (
    compact_nodes,
    keep_components,
    remove_elements,
    weld_close_nodes,
)
from fvcom_mesh_tools.mesh_clean_phase_h import (
    build_coastline_projector,
    phase_h_optimize,
)

# jsasaki 2026-11-16: new helpers for channel protection, OBC carry-over and lone nodes
__all__ = ["finish_constrained_mesh", "keep_components_protected", "obc_perpendicularity",
           "quality_counts",
           "rebuild_land_boundaries", "remove_elements_keep_boundaries",
           "remove_lone_nodes"]

_DEG_PER_M = 1.0 / 111_194.92664455873


def quality_counts(mesh: Fort14Mesh) -> dict[str, Any]:
    """C1/C2/flip counts and the minimum interior angle."""
    mn, mx, tw = _tri_quality(mesh.nodes[mesh.elements])
    return {
        "c1": int((mn < 30).sum()),
        "c2": int((mx > 130).sum()),
        "flipped": int((tw <= 0).sum()),
        "min_angle_deg": float(mn.min()) if mn.size else 90.0,
    }


def _boundary_node_set(mesh: Fort14Mesh) -> set[int]:
    els = mesh.elements
    raw = np.vstack([els[:, [0, 1]], els[:, [1, 2]], els[:, [2, 0]]])
    raw.sort(axis=1)
    codes = raw[:, 0].astype(np.int64) * mesh.n_nodes + raw[:, 1]
    uniq, counts = np.unique(codes, return_counts=True)
    b = uniq[counts == 1]
    return set(np.column_stack([b // mesh.n_nodes,
                                b % mesh.n_nodes]).ravel().tolist())


def finish_constrained_mesh(
    mesh: Fort14Mesh,
    shoreline_shp: Path,
    work_dir: Path,
    *,
    min_angle: float = 30.0,
    max_angle: float = 130.0,
    cap_min_angle: float = 20.0,
    cap_max_angle: float = 150.0,
    cap_rounds: int = 3,
    weld_tol_m: float = 2.0,
    optimize_budget_s: float = 600.0,
    projector_snap_m: float = 200.0,
    # jsasaki 2026-11-16: no Tokyo default (was 32654): the caller names the CRS
    utm_epsg: int | None = None,
    #---> jsasaki 2026-11-16: channel protection, OBC carry-over, valence, lone nodes
    protect_polygons=None,
    on_dropped_protected: str = "fail",
    max_valence: int | None = None,
    valence_min_angle_deg: float = 0.0,
    remove_lone: bool = False,
    freeze_obc_band: bool = False,
    #<--- jsasaki 2026-11-16
    log: Callable[[str], None] = print,
) -> tuple[Fort14Mesh, dict[str, Any]]:
    """Finish a CDT-boundary mesh (UTM coordinates) as in PoC #93.

    ``shoreline_shp`` is the engineered shoreline (EPSG:4326
    polygons) the boundary is constrained to; a no-CRS UTM copy for
    the projector is cached under ``work_dir``.

    ``utm_epsg`` (required) is the EPSG code of the mesh coordinates.
    ``protect_polygons`` (shapely geometry in mesh coordinates): inside
    them a cap triangle is deleted only when it is an ear (at most one
    neighbouring element), so cap deletion can never cut a narrow
    channel; a component dropped afterwards that holds protected water
    raises (``on_dropped_protected="fail"``) or is reported
    (``"report"``). An open boundary on the input is carried through
    (element removal keeps it; the components touching it are kept).
    ``max_valence`` adds a valence-reducing edge-swap pass after the
    optimizer; ``remove_lone`` deletes the elements of single-element
    nodes off the open boundary (an FVCOM QA gate). ``freeze_obc_band``
    restores, after the optimizer, the open-boundary nodes and their
    neighbours to their positions before it where that flips nothing and
    worsens no angle below 30 deg (the optimizer otherwise undoes a
    perpendicular open-boundary ladder: VR v3 job 7590960, 26 nodes over
    20 deg). With none of these and no open boundary the function is the
    PoC #93 sequence unchanged.
    """
    import geopandas as gpd

    # jsasaki 2026-11-16: a silent Tokyo CRS default is a wrong answer elsewhere
    if utm_epsg is None:
        raise ValueError("finish_constrained_mesh: utm_epsg is required (EPSG of the mesh)")
    info: dict[str, Any] = {"input": quality_counts(mesh)}
    log(f"[finish] input NP={mesh.n_nodes:,} quality={info['input']}")

    #---> jsasaki 2026-11-16: guarded path when an OBC or protected polygons are present
    has_obc = any(len(s) for s in mesh.open_boundaries)
    guarded = has_obc or protect_polygons is not None
    if guarded:
        prot_mask = _protected_node_mask(mesh, protect_polygons)
    #<--- jsasaki 2026-11-16
    # 1. cap-triangle deletion (bounded rounds, monotone).
    removed = 0
    for _ in range(cap_rounds):
        bset = _boundary_node_set(mesh)
        mn, mx, _tw = _tri_quality(mesh.nodes[mesh.elements])
        bad = (mn < cap_min_angle) | (mx > cap_max_angle)
        allb = np.isin(mesh.elements, list(bset)).all(axis=1)
        kill = bad & allb
        #---> jsasaki 2026-11-16: never cut a protected channel or the open boundary
        if guarded:
            kill &= _cap_deletion_allowed(mesh, prot_mask)
        #<--- jsasaki 2026-11-16
        if not kill.any():
            break
        #---> jsasaki 2026-11-16: the guarded path keeps boundaries and reports drops
        if guarded:
            keep_tag = np.zeros(mesh.n_nodes, bool) if prot_mask is None else prot_mask
            mesh, prot_mask = remove_elements_keep_boundaries(mesh, ~kill, keep_tag)
            mesh, kinfo, prot_mask = keep_components_protected(
                mesh, prot_mask, policy=on_dropped_protected, log=log)
            info.setdefault("components_dropped", []).extend(kinfo["dropped"])
        else:
            mesh = remove_elements(mesh, ~kill)
            mesh, _ = keep_components(mesh)
        #<--- jsasaki 2026-11-16
        removed += int(kill.sum())
    mesh, _ = compact_nodes(mesh)
    info["caps_removed"] = removed
    log(f"[finish] caps removed: {removed} quality={quality_counts(mesh)}")

    # 2. weld.
    mesh, winfo = weld_close_nodes(mesh, tol=weld_tol_m)
    info["weld"] = winfo
    log(f"[finish] weld: {winfo}")

    # 3. one budgeted optimize round, boundary sliding on-line only.
    work_dir.mkdir(parents=True, exist_ok=True)
    shp_utm = work_dir / "finish_shoreline_utm_nocrs.shp"
    if not shp_utm.exists():
        gdf = gpd.read_file(shoreline_shp).to_crs(utm_epsg)
        gdf = gdf.set_crs(None, allow_override=True)
        gdf.to_file(shp_utm)
    projector = build_coastline_projector(
        [shp_utm],
        max_snap_distance_m=projector_snap_m / _DEG_PER_M,
        mean_latitude_deg=0.0,
    )
    # jsasaki 2026-11-16: tag nodes so the open-boundary band can be restored
    band_state = _tag_obc_band(mesh) if freeze_obc_band else None
    mesh, oinfo = phase_h_optimize(
        mesh,
        min_angle_target=min_angle,
        max_angle_target=max_angle,
        pass_f_enabled=True,
        pass_g_enabled=True,
        pass_g_min_angle_target=min_angle,
        max_outer_rounds=1,
        coastline_projector=projector,
        time_budget_s=optimize_budget_s,
    )
    mesh, cinfo = compact_nodes(mesh)
    # jsasaki 2026-11-16: restore the open-boundary band moved by the optimizer
    if band_state is not None:
        mesh, info["obc_band"] = _restore_obc_band(mesh, band_state, log=log)
    info["optimize"] = {
        "budget_exhausted": bool(oinfo.get("budget_exhausted")),
        "compact": cinfo,
    }
    #---> jsasaki 2026-11-16: optional valence reduction and lone-node removal
    if max_valence is not None:
        from fvcom_mesh_tools.mesh_clean import repair_overconnected_nodes

        v0 = int(np.bincount(mesh.elements.ravel()).max())
        mesh, vinfo = repair_overconnected_nodes(
            mesh, max_nbr_elem=int(max_valence),
            min_angle_floor_deg=float(valence_min_angle_deg))
        v1 = int(np.bincount(mesh.elements.ravel()).max())
        info["valence"] = {"max_before": v0, "max_after": v1,
                           "swaps": int(vinfo.get("total_swaps", vinfo.get("n_swaps", 0)))}
        log(f"[finish] valence <= {max_valence}: max {v0} -> {v1} {info['valence']}")
    if remove_lone:
        mesh, linfo = remove_lone_nodes(mesh, protect_polygons=protect_polygons,
                                        log=log)
        info["lone_nodes"] = linfo
    if guarded or remove_lone:
        mesh = rebuild_land_boundaries(mesh)
    #<--- jsasaki 2026-11-16
    info["output"] = quality_counts(mesh)
    log(f"[finish] done NP={mesh.n_nodes:,} quality={info['output']}")
    return mesh, info


#---> jsasaki 2026-11-16: new-region helpers (channel protection, boundary bookkeeping,
# open-boundary band restore)
# ---------------------------------------------------------------------------
# Channel protection and boundary bookkeeping (2026-11-16, new region support)
# ---------------------------------------------------------------------------


def obc_perpendicularity(mesh: Fort14Mesh) -> np.ndarray:
    """Best-edge deviation from perpendicular (deg) of every open-boundary
    node, in the order of the concatenated open strings (the QA gate's
    metric, :func:`fvcom_mesh_tools.qa.run_qa` ``obc_perpendicularity``)."""
    from fvcom_mesh_tools.algorithms.perpendicularity import boundary_tangents

    out = []
    els = mesh.elements
    e = np.unique(np.sort(np.vstack([els[:, [0, 1]], els[:, [1, 2]], els[:, [2, 0]]]),
                          axis=1), axis=0)
    for seg in mesh.open_boundaries:
        seg = np.asarray(seg, np.int64)
        dev = np.full(seg.size, np.inf)
        if seg.size >= 2:
            tang = boundary_tangents(mesh.nodes[seg])
            inv = np.full(mesh.n_nodes, -1, np.int64)
            inv[seg] = np.arange(seg.size)
            a_in, b_in = inv[e[:, 0]] >= 0, inv[e[:, 1]] >= 0
            inc = e[a_in ^ b_in]
            bnode = np.where(inv[inc[:, 0]] >= 0, inc[:, 0], inc[:, 1])
            onode = np.where(inv[inc[:, 0]] >= 0, inc[:, 1], inc[:, 0])
            vec = mesh.nodes[onode] - mesh.nodes[bnode]
            vec = vec / np.linalg.norm(vec, axis=1, keepdims=True)
            cos = np.clip(np.abs((vec * tang[inv[bnode]]).sum(axis=1)), 0.0, 1.0)
            np.minimum.at(dev, inv[bnode], 90.0 - np.degrees(np.arccos(cos)))
        out.append(dev)
    return np.concatenate(out) if out else np.empty(0)


def _tag_obc_band(mesh: Fort14Mesh):
    """Replace the depths by random node tags (carried through the optimizer:
    kept per node, averaged for split nodes) and record the open-boundary
    band (open nodes and their element neighbours) before the optimizer."""
    rng = np.random.default_rng(20261116)
    tags = 10.0 + rng.random(mesh.n_nodes)
    obc = np.zeros(mesh.n_nodes, bool)
    for s in mesh.open_boundaries:
        obc[np.asarray(s, int)] = True
    band = np.zeros(mesh.n_nodes, bool)
    band[np.unique(mesh.elements[obc[mesh.elements].any(axis=1)])] = True
    state = {"tags": tags, "depths": mesh.depths.copy(), "xy": mesh.nodes.copy(),
             "band": band}
    mesh.depths = tags.copy()
    return state


def _restore_obc_band(mesh: Fort14Mesh, state, *, min_angle: float = 30.0,
                      max_area_change: float = 0.5,
                      log: Callable[[str], None] = print):
    """Undo :func:`_tag_obc_band`: depths back (new nodes take the depth of
    the nearest original node) and band nodes back to their positions where
    no incident element flips, the incident minimum angle stays >= the
    smaller of its current value and ``min_angle``, and the adjacent-area
    change around the node stays <= the larger of its current value and
    ``max_area_change``."""
    from scipy.spatial import cKDTree

    idx = {float(t): i for i, t in enumerate(state["tags"])}
    orig = np.array([idx.get(float(t), -1) for t in mesh.depths], int)
    known = orig >= 0
    dep = np.empty(mesh.n_nodes)
    dep[known] = state["depths"][orig[known]]
    if (~known).any():
        _d, j = cKDTree(state["xy"]).query(mesh.nodes[~known])
        dep[~known] = state["depths"][j]
    mesh.depths = dep
    n2e: dict[int, list[int]] = {}
    for e, tri in enumerate(mesh.elements):
        for v in tri:
            n2e.setdefault(int(v), []).append(e)
    # element pairs across each interior edge (for the area-change check)
    raw = np.sort(np.stack([mesh.elements[:, [0, 1]], mesh.elements[:, [1, 2]],
                            mesh.elements[:, [2, 0]]], axis=1), axis=2).reshape(-1, 2)
    code = raw[:, 0].astype(np.int64) * mesh.n_nodes + raw[:, 1]
    eid = np.repeat(np.arange(mesh.n_elements), 3)
    order = np.argsort(code, kind="stable")
    cs, es = code[order], eid[order]
    same = np.nonzero(cs[1:] == cs[:-1])[0]
    pairs = np.column_stack([es[same], es[same + 1]])
    e2p: dict[int, list[int]] = {}
    for k, (a, b) in enumerate(pairs):
        e2p.setdefault(int(a), []).append(k)
        e2p.setdefault(int(b), []).append(k)

    def area_change(pk):
        _mn, _mx, tw = _tri_quality(mesh.nodes[mesh.elements[pairs[pk].ravel()]])
        ar = 0.5 * np.abs(tw).reshape(-1, 2)
        return float(((ar.max(axis=1) - ar.min(axis=1)) / ar.max(axis=1)).max())

    restored = refused = 0
    for v in np.nonzero(known)[0]:
        o = orig[v]
        if not state["band"][o] or np.array_equal(mesh.nodes[v], state["xy"][o]):
            continue
        inc = n2e.get(int(v), [])
        els = mesh.elements[inc]
        pk = sorted({k for e in inc for k in e2p.get(e, [])})
        mn0, _mx0, _tw0 = _tri_quality(mesh.nodes[els])
        c40 = area_change(pk) if pk else 0.0
        old = mesh.nodes[v].copy()
        mesh.nodes[v] = state["xy"][o]
        mn1, _mx1, tw1 = _tri_quality(mesh.nodes[els])
        c41 = area_change(pk) if pk else 0.0
        if (tw1 <= 0).any() or mn1.min() < min(mn0.min(), min_angle) \
                or c41 > max(c40, max_area_change):
            mesh.nodes[v] = old
            refused += 1
        else:
            restored += 1
    info = {"restored": restored, "refused": refused}
    log(f"[finish] open-boundary band restored: {info}")
    return mesh, info


def _protected_node_mask(mesh: Fort14Mesh, protect_polygons):
    """Nodes inside ``protect_polygons`` (None when there is no protection)."""
    if protect_polygons is None:
        return None
    import shapely

    return shapely.contains_xy(protect_polygons, mesh.nodes[:, 0], mesh.nodes[:, 1])


def _element_neighbour_count(elements: np.ndarray) -> np.ndarray:
    raw = np.sort(np.stack([elements[:, [0, 1]], elements[:, [1, 2]],
                            elements[:, [2, 0]]], axis=1), axis=2)
    n = int(elements.max()) + 1
    code = raw[..., 0].astype(np.int64) * n + raw[..., 1]
    uniq, inv, counts = np.unique(code.ravel(), return_inverse=True, return_counts=True)
    shared = (counts[inv] == 2).reshape(code.shape)
    return shared.sum(axis=1)


def _cap_deletion_allowed(mesh: Fort14Mesh, prot_mask) -> np.ndarray:
    """Per element: may a cap triangle be deleted here?

    Never when it holds two or more open-boundary nodes (it would break
    the open string). Inside protected polygons only an ear (at most one
    neighbouring element) may go: deleting an ear cannot disconnect the
    element graph, deleting a bridge can (the two-element-wide channel
    corridor that VR v2 lost).
    """
    ok = np.ones(mesh.n_elements, bool)
    if any(len(s) for s in mesh.open_boundaries):
        obc = np.zeros(mesh.n_nodes, bool)
        obc[np.concatenate([np.asarray(s, int) for s in mesh.open_boundaries])] = True
        ok &= obc[mesh.elements].sum(axis=1) < 2
    if prot_mask is not None:
        in_prot = prot_mask[mesh.elements].any(axis=1)
        ear = _element_neighbour_count(mesh.elements) <= 1
        ok &= ~in_prot | ear
    return ok


def remove_elements_keep_boundaries(mesh: Fort14Mesh, keep_elem_mask, node_tag=None):
    """:func:`mesh_clean.remove_elements` that remaps the open boundaries
    (dropped nodes leave the string) instead of clearing them; land
    boundaries are cleared (use :func:`rebuild_land_boundaries`).
    ``node_tag`` (bool per node) is carried along and returned."""
    keep_elem_mask = np.asarray(keep_elem_mask, bool)
    els = mesh.elements[keep_elem_mask]
    used = np.unique(els)
    remap = np.full(mesh.n_nodes, -1, np.int64)
    remap[used] = np.arange(len(used))
    opens = []
    for seg in mesh.open_boundaries:
        s2 = remap[np.asarray(seg, np.int64)]
        s2 = s2[s2 >= 0]
        if s2.size:
            opens.append(s2)
    out = Fort14Mesh(title=mesh.title, obc_type=mesh.obc_type,
                     nodes=mesh.nodes[used].copy(), depths=mesh.depths[used].copy(),
                     elements=remap[els].astype(mesh.elements.dtype),
                     open_boundaries=opens, land_boundaries=[])
    tag = None if node_tag is None else np.asarray(node_tag, bool)[used]
    return out, tag


def keep_components_protected(mesh: Fort14Mesh, node_tag=None, *, policy: str = "fail",
                              log: Callable[[str], None] = print):
    """Keep the element components that touch the open boundary (the
    largest one when there is none); a dropped component holding a node
    tagged in ``node_tag`` (protected water) raises (``policy="fail"``)
    or is logged (``"report"``). Every dropped component is listed in
    ``info["dropped"]`` (elements, area, centroid). Open boundaries are
    remapped, not cleared. Returns ``(mesh, info, node_tag)``."""
    from scipy.sparse.csgraph import connected_components

    from fvcom_mesh_tools.diagnostics import face_face_adjacency

    if policy not in ("fail", "report"):
        raise ValueError(f"policy must be 'fail' or 'report', not {policy!r}")
    info: dict[str, Any] = {"dropped": []}
    if mesh.n_elements == 0:
        return mesh, info, node_tag
    n_comp, lab = connected_components(face_face_adjacency(mesh.elements), directed=False)
    if n_comp == 1:
        return mesh, info, node_tag
    sizes = np.bincount(lab, minlength=n_comp)
    obc = [int(v) for s in mesh.open_boundaries for v in s]
    if obc:
        m = np.zeros(mesh.n_nodes, bool)
        m[obc] = True
        keep = set(np.unique(lab[m[mesh.elements].any(axis=1)]).tolist())
    else:
        keep = {int(sizes.argmax())}
    _mn, _mx, tw = _tri_quality(mesh.nodes[mesh.elements])
    area = 0.5 * np.abs(tw)
    cen = mesh.nodes[mesh.elements].mean(axis=1)
    bad = []
    for c in range(n_comp):
        if c in keep:
            continue
        sel = lab == c
        rec = {"n_elements": int(sel.sum()), "area_m2": float(area[sel].sum()),
               "x": float(cen[sel, 0].mean()), "y": float(cen[sel, 1].mean()),
               "protected": bool(node_tag is not None
                                 and node_tag[np.unique(mesh.elements[sel])].any())}
        info["dropped"].append(rec)
        if rec["protected"]:
            bad.append(rec)
    msg = "; ".join(f"{r['n_elements']} el / {r['area_m2'] / 1e6:.3f} km2 at "
                    f"({r['x']:.0f}, {r['y']:.0f})" for r in bad)
    if bad and policy == "fail":
        raise RuntimeError(f"components holding protected water would be dropped: {msg}")
    if info["dropped"]:
        log(f"[components] dropped {len(info['dropped'])} component(s), "
            f"{sum(r['n_elements'] for r in info['dropped'])} elements"
            + (f"; PROTECTED: {msg}" if bad else ""))
    out, tag = remove_elements_keep_boundaries(mesh, np.isin(lab, list(keep)), node_tag)
    return out, info, tag


def remove_lone_nodes(mesh: Fort14Mesh, *, protect_polygons=None, max_rounds: int = 10,
                      log: Callable[[str], None] = print):
    """Delete the elements of nodes used by one element only (off the open
    boundary), until none is left. Such an element is an ear (it has at most
    one neighbour), so its deletion never disconnects the mesh. Returns
    ``(mesh, info)`` with the deleted count and area."""
    n_del, area = 0, 0.0
    for _ in range(max_rounds):
        cnt = np.bincount(mesh.elements.ravel(), minlength=mesh.n_nodes)
        lone = cnt == 1
        for s in mesh.open_boundaries:
            lone[np.asarray(s, int)] = False
        kill = lone[mesh.elements].any(axis=1)
        if not kill.any():
            break
        _mn, _mx, tw = _tri_quality(mesh.nodes[mesh.elements[kill]])
        area += float(0.5 * np.abs(tw).sum())
        n_del += int(kill.sum())
        mesh, _ = remove_elements_keep_boundaries(mesh, ~kill)
    info = {"elements_deleted": n_del, "area_m2": area}
    log(f"[finish] single-element nodes: {n_del} elements deleted ({area:.0f} m2)")
    return mesh, info


def rebuild_land_boundaries(mesh: Fort14Mesh, *, land_ibtype: int = 20,
                            island_ibtype: int = 21) -> Fort14Mesh:
    """Land boundary lists from the boundary loops: the outer loop minus
    the open string (one land string, ``land_ibtype``) and every other
    loop as a closed island (``island_ibtype``). The open boundaries are
    kept as they are (at most one open string is supported)."""
    from fvcom_mesh_tools.algorithms.boundary import (
        boundary_edges_from_tris,
        chain_edges_to_loops,
        outer_loop,
    )

    loops = chain_edges_to_loops(boundary_edges_from_tris(mesh.elements))
    if not loops:
        return mesh
    outer = outer_loop(loops, mesh.nodes)
    land: list[tuple[int, np.ndarray]] = []
    opens = [np.asarray(s, np.int64) for s in mesh.open_boundaries if len(s)]
    ring = outer[:-1]
    if len(opens) == 1 and len(opens[0]) >= 2 and set(opens[0].tolist()) <= set(ring.tolist()):
        o = opens[0]
        r = np.roll(ring, -int(np.nonzero(ring == o[0])[0][0]))
        if r[1] != o[1]:
            r = np.roll(r[::-1], 1)
        stop = int(np.nonzero(r == o[-1])[0][0])
        land.append((land_ibtype, np.append(r[stop:], r[0])))
    elif not opens:
        land.append((land_ibtype, outer.copy()))
    else:
        return Fort14Mesh(title=mesh.title, obc_type=mesh.obc_type, nodes=mesh.nodes,
                          depths=mesh.depths, elements=mesh.elements,
                          open_boundaries=mesh.open_boundaries, land_boundaries=[])
    for lp in loops:
        if lp is outer:
            continue
        land.append((island_ibtype, lp.copy()))
    return Fort14Mesh(title=mesh.title, obc_type=mesh.obc_type, nodes=mesh.nodes,
                      depths=mesh.depths, elements=mesh.elements,
                      open_boundaries=mesh.open_boundaries, land_boundaries=land)
#<--- jsasaki 2026-11-16
