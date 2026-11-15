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
# jsasaki 2026-11-16: review newregion r1 #1: protected_node_mask exported
# jsasaki 2026-11-16: review newregion r2 #1: protected_element_mask, protected_overlap_m2 exported
__all__ = ["finish_constrained_mesh", "keep_components_protected", "obc_perpendicularity",
           "protected_element_mask", "protected_node_mask", "protected_overlap_m2",
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
    # jsasaki 2026-11-16: review newregion r3 #2: protected-water loss allowed per guarded step
    protected_loss_tol_m2: float = 0.0,
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
    # jsasaki 2026-11-16: review newregion r3 #2: the input, for the invariants checked at the end
    mesh_in = _snapshot(mesh) if guarded else None
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
            # jsasaki 2026-11-16: review newregion r2 #1: element overlap, recomputed per round
            prot_el = (None if protect_polygons is None else
                       protected_element_mask(mesh.nodes, mesh.elements, protect_polygons))
            kill &= _cap_deletion_allowed(mesh, prot_el)
        #<--- jsasaki 2026-11-16
        if not kill.any():
            break
        #---> jsasaki 2026-11-16: the guarded path keeps boundaries and reports drops
        if guarded:
            # jsasaki 2026-11-16: review newregion r2 #1: drops judged by element overlap
            mesh, _ = remove_elements_keep_boundaries(mesh, ~kill)
            mesh, kinfo, _ = keep_components_protected(
                mesh, protect_polygons=protect_polygons, policy=on_dropped_protected, log=log)
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
    #---> jsasaki 2026-11-16: review newregion r3 #2: guarded weld: open-boundary nodes are never
    # welded, and a weld that loses protected water is refused
    if guarded:
        # jsasaki 2026-11-16: review newregion r3 #2: an open-boundary node is welded only onto a
        # COINCIDENT node (< 1 mm; the OBC assignment can leave such a pair, VR job 7591832): that
        # repairs the string without moving its geometry
        from scipy.spatial import cKDTree

        coinc = {int(v) for pr in cKDTree(mesh.nodes).query_pairs(1e-3) for v in pr}
        obc_ids = [int(v) for s_ in mesh.open_boundaries for v in s_ if int(v) not in coinc]
        pre = _snapshot(mesh)
        mesh, winfo = weld_close_nodes(mesh, tol=weld_tol_m, exclude=obc_ids)
        if protect_polygons is not None and winfo.get("n_welded"):
            lost = protected_loss_m2(pre, mesh, protect_polygons)
            if lost > protected_loss_tol_m2:
                mesh = pre
                winfo = {"n_welded": 0, "refused_loss_m2": lost}
    else:
        mesh, winfo = weld_close_nodes(mesh, tol=weld_tol_m)
    #<--- jsasaki 2026-11-16
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
    # jsasaki 2026-11-16: review newregion r3 #2: also with protected water (its bank nodes are
    # restored when the optimizer loses protected water)
    pre_opt = _snapshot(mesh) if guarded else None     # before the depth tags
    band_state = (_tag_obc_band(mesh) if freeze_obc_band or protect_polygons is not None
                  else None)
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
        if not freeze_obc_band:
            band_state["band"][:] = False     # review newregion r3 #2: only the depth tags
            band_state["obc"][:] = False      # review newregion r3 #2: OBC restored if frozen
        mesh, info["obc_band"] = _restore_obc_band(mesh, band_state, log=log)
    #---> jsasaki 2026-11-16: review newregion r3 #2: the optimizer may not lose protected water
    # or break the open boundary: bank nodes go back where it lost water; if that is not enough,
    # the whole optimizer round is refused
    if guarded:
        mesh, info["optimizer_guard"] = _guard_optimizer(
            pre_opt, mesh, band_state, protect_polygons, protected_loss_tol_m2,
            freeze_obc_band, log)
    #<--- jsasaki 2026-11-16
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
        # jsasaki 2026-11-16: review newregion r1 #2: protection policy applies here too
        mesh, linfo = remove_lone_nodes(mesh, protect_polygons=protect_polygons,
                                        policy=on_dropped_protected, log=log)
        info["lone_nodes"] = linfo
    if guarded or remove_lone:
        mesh = rebuild_land_boundaries(mesh)
    #<--- jsasaki 2026-11-16
    #---> jsasaki 2026-11-16: review newregion r3 #2: invariants before returning
    if guarded:
        # jsasaki 2026-11-16: review newregion r3 #2: string intact here; its GEOMETRY is re-derived
        # and validated against the line by the obc/obcfinal stages the recipe must have (r1 #3)
        bad = validate_open_boundary(mesh_in, mesh, same_positions=False)
        if bad:
            raise RuntimeError(f"finish broke the open boundary: {bad}")
        if protect_polygons is not None:
            lost = protected_loss_m2(mesh_in, mesh, protect_polygons)
            info["protected_loss_m2"] = lost
            log(f"[finish] protected water lost: {lost:.1f} m2 (tolerance "
                f"{protected_loss_tol_m2:g} m2 per step)")
    #<--- jsasaki 2026-11-16
    info["output"] = quality_counts(mesh)
    log(f"[finish] done NP={mesh.n_nodes:,} quality={info['output']}")
    return mesh, info


#---> jsasaki 2026-11-16: new-region helpers (channel protection, boundary bookkeeping,
# open-boundary band restore)
# ---------------------------------------------------------------------------
# Channel protection and boundary bookkeeping (2026-11-16, new region support)
# ---------------------------------------------------------------------------


#---> jsasaki 2026-11-16: review newregion r3 #2/#3/#4: invariants shared by finish, siteops and
# the finishing stage
def _snapshot(mesh: Fort14Mesh) -> Fort14Mesh:
    return Fort14Mesh(title=mesh.title, obc_type=mesh.obc_type, nodes=mesh.nodes.copy(),
                      depths=mesh.depths.copy(), elements=mesh.elements.copy(),
                      open_boundaries=[np.asarray(s).copy() for s in mesh.open_boundaries],
                      land_boundaries=[(ib, np.asarray(s).copy())
                                       for ib, s in mesh.land_boundaries])


def _tri_keys(nodes, elements) -> list:
    pts = np.round(np.asarray(nodes, float)[np.asarray(elements, np.int64)], 6)
    return [tuple(sorted(map(tuple, t))) for t in pts]


# jsasaki 2026-11-16: review newregion r3 #3: numerical floor of the geometric loss (unions of
# coinciding triangles leave slivers of ~1e-6 m2); NOT a tolerance
LOSS_EPS_M2 = 1e-3


def protected_loss_m2(before: Fort14Mesh, after: Fort14Mesh, protect_polygons) -> float:
    """Protected water covered by ``before`` and not by ``after`` (m2): the
    geometric loss, NOT netted against water gained elsewhere. Only the
    triangles that differ between the two meshes are compared (identical
    triangles cover the same water)."""
    import shapely

    kb = _tri_keys(before.nodes, before.elements)
    ka = _tri_keys(after.nodes, after.elements)
    sa, sb = set(ka), set(kb)
    rm = np.array([k not in sa for k in kb])
    add = np.array([k not in sb for k in ka])
    if not rm.any():
        return 0.0
    pieces, tree = _protect_tiles(protect_polygons)
    t_rm = shapely.polygons(before.nodes[before.elements[rm]])
    hit = np.unique(tree.query(t_rm, predicate="intersects")[0])
    if not hit.size:
        return 0.0
    gone = shapely.intersection(shapely.union_all(t_rm[hit]), protect_polygons)
    if add.any():
        t_add = shapely.polygons(after.nodes[after.elements[add]])
        gone = shapely.difference(gone, shapely.union_all(t_add))
    a = float(shapely.area(gone))
    return 0.0 if a <= LOSS_EPS_M2 else a


def validate_open_boundary(before: Fort14Mesh, after: Fort14Mesh, *,
                           same_positions: bool = True) -> list[str]:
    """Problems with the open boundary of ``after`` relative to ``before``:
    lost strings, changed node count, nodes off the original open line or out
    of order (``same_positions``), consecutive nodes not joined by a boundary
    edge."""
    probs: list[str] = []
    ob = [np.asarray(s, np.int64) for s in before.open_boundaries if len(s)]
    oa = [np.asarray(s, np.int64) for s in after.open_boundaries if len(s)]
    if len(ob) != len(oa):
        return [f"{len(ob)} open string(s) before, {len(oa)} after"]
    els = after.elements
    raw = np.sort(np.vstack([els[:, [0, 1]], els[:, [1, 2]], els[:, [2, 0]]]), axis=1)
    code = raw[:, 0].astype(np.int64) * after.n_nodes + raw[:, 1]
    u, c = np.unique(code, return_counts=True)
    # jsasaki 2026-11-16: review newregion r3 #2: a coincident (< 1 mm) repeat in the string is one
    # node (a weld of an exact duplicate repairs the string, it does not shorten it)
    def _dedup(nodes, seg):
        keep = np.ones(len(seg), bool)
        keep[1:] = np.hypot(*(nodes[seg[1:]] - nodes[seg[:-1]]).T) > 1e-3
        return seg[keep]

    ob = [_dedup(before.nodes, b) for b in ob]
    oa = [_dedup(after.nodes, a) for a in oa]
    for k, (b, a) in enumerate(zip(ob, oa)):
        if len(a) != len(b):
            probs.append(f"string {k}: {len(b)} nodes before, {len(a)} after")
            continue
        # jsasaki 2026-11-16: review newregion r3 #2: the open-boundary GEOMETRY must stay: every
        # node on the original string's polyline (1 mm), in the original order (sliding along the
        # line is allowed; a weld that moved (1,-1) to (0,0) is not)
        if same_positions:
            import shapely

            line = shapely.LineString(before.nodes[b]) if len(b) >= 2 else None
            if line is None:
                if not np.allclose(before.nodes[b], after.nodes[a], atol=1e-3):
                    probs.append(f"string {k}: node moved")
            else:
                pts = shapely.points(after.nodes[a])
                off = shapely.distance(pts, line)
                pos = np.array([line.project(q) for q in pts])
                if (off > 1e-3).any():
                    probs.append(f"string {k}: {int((off > 1e-3).sum())} node(s) off the open line "
                                 f"(max {off.max():.3f} m)")
                elif (np.diff(pos) <= 0).any():
                    probs.append(f"string {k}: node order changed")
        # jsasaki 2026-11-16: review newregion r4 #4: coincident (< 1 mm) nodes count as one node
        # (the OBC assignment can leave such a pair; VR job 7592472): edges between their classes
        from scipy.spatial import cKDTree

        rep_ = np.arange(after.n_nodes)
        for i, j in cKDTree(after.nodes).query_pairs(1e-3):
            rep_[max(i, j)] = rep_[min(i, j)]
        bu = np.sort(rep_[np.column_stack([u[c == 1] // after.n_nodes,
                                           u[c == 1] % after.n_nodes])], axis=1)
        bset = {(int(x), int(y)) for x, y in bu}
        ee = np.sort(rep_[np.column_stack([a[:-1], a[1:]])], axis=1)
        nb = sum((int(x), int(y)) not in bset for x, y in ee)
        if nb:
            probs.append(f"string {k}: {nb} consecutive pair(s) not on a boundary edge")
    return probs


def _guard_optimizer(pre, mesh, band_state, protect_polygons, tol, freeze_obc_band, log):
    """After the optimizer: where protected water was lost, the mesh-boundary
    nodes of protected elements go back to their pre-optimizer positions
    (flip-free moves only); a loss still above ``tol``, or a broken open
    boundary, refuses the whole optimizer round (the pre-optimizer mesh)."""
    info: dict[str, Any] = {}
    if protect_polygons is not None:
        lost = protected_loss_m2(pre, mesh, protect_polygons)
        info["loss_m2"] = lost
        orig = None if band_state is None else band_state.get("orig_of")
        if lost > tol and orig is not None:
            pe = protected_element_mask(mesh.nodes, mesh.elements, protect_polygons)
            raw = np.sort(np.vstack([mesh.elements[:, [0, 1]], mesh.elements[:, [1, 2]],
                                     mesh.elements[:, [2, 0]]]), axis=1)
            u, c = np.unique(raw, axis=0, return_counts=True)
            bnd = np.zeros(mesh.n_nodes, bool)
            bnd[u[c == 1].ravel()] = True
            cand = np.unique(mesh.elements[pe])
            cand = cand[bnd[cand]]
            n2e: dict[int, list[int]] = {}
            for e, tri in enumerate(mesh.elements):
                for v in tri:
                    n2e.setdefault(int(v), []).append(e)
            moved = 0
            for v in cand:
                o = int(orig[v])
                if o < 0 or np.array_equal(mesh.nodes[v], pre.nodes[o]):
                    continue
                old = mesh.nodes[v].copy()
                mesh.nodes[v] = pre.nodes[o]
                _mn, _mx, tw = _tri_quality(mesh.nodes[mesh.elements[n2e.get(int(v), [])]])
                if (tw <= 0).any():
                    mesh.nodes[v] = old
                else:
                    moved += 1
            lost = protected_loss_m2(pre, mesh, protect_polygons)
            info.update(bank_nodes_restored=moved, loss_after_restore_m2=lost)
        if lost > tol:
            log(f"[finish] optimizer refused: it lost {lost:.1f} m2 of protected water")
            info["refused"] = True
            return pre, info
    # jsasaki 2026-11-16: review newregion r3 #2: string intact (geometry: obc/obcfinal stages)
    bad = validate_open_boundary(pre, mesh, same_positions=False)
    if bad:
        log(f"[finish] optimizer refused: open boundary {bad}")
        info.update(refused=True, obc=bad)
        return pre, info
    log(f"[finish] optimizer guard: {info}")
    return mesh, info
#<--- jsasaki 2026-11-16


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
    # jsasaki 2026-11-16: review newregion r1 #4: triangulation kept for new-node depths
    # jsasaki 2026-11-16: review newregion r3 #2: the open-boundary nodes (restored first, flip
    # test only: the optimizer slid an end node 295 m round the junction, VR job 7591939)
    state = {"tags": tags, "depths": mesh.depths.copy(), "xy": mesh.nodes.copy(),
             "elements": mesh.elements.copy(), "band": band, "obc": obc}
    mesh.depths = tags.copy()
    return state


def _restore_obc_band(mesh: Fort14Mesh, state, *, min_angle: float = 30.0,
                      max_area_change: float = 0.5,
                      log: Callable[[str], None] = print):
    """Undo :func:`_tag_obc_band`: depths back (a new node takes the depth
    interpolated linearly on the original triangulation, the nearest original
    node's depth only outside it) and band nodes back to their positions where
    no incident element flips, the incident minimum angle stays >= the
    smaller of its current value and ``min_angle``, and the adjacent-area
    change of EVERY pair around the node stays <= the larger of its current
    value and ``max_area_change``."""
    from scipy.spatial import cKDTree

    idx = {float(t): i for i, t in enumerate(state["tags"])}
    orig = np.array([idx.get(float(t), -1) for t in mesh.depths], int)
    known = orig >= 0
    # jsasaki 2026-11-16: review newregion r3 #2: node -> pre-optimizer node, for the bank restore
    state["orig_of"] = orig
    dep = np.empty(mesh.n_nodes)
    dep[known] = state["depths"][orig[known]]
    if (~known).any():
        # jsasaki 2026-11-16: review newregion r1 #4: interpolate (a split midpoint keeps the
        # optimizer's mean depth) instead of copying the nearest node
        from matplotlib.tri import LinearTriInterpolator, Triangulation

        tri = Triangulation(state["xy"][:, 0], state["xy"][:, 1], state["elements"])
        q = mesh.nodes[~known]
        val = np.ma.filled(LinearTriInterpolator(tri, state["depths"])(q[:, 0], q[:, 1]),
                           np.nan)
        miss = ~np.isfinite(val)
        if miss.any():
            _d, j = cKDTree(state["xy"]).query(q[miss])
            val[miss] = state["depths"][j]
        dep[~known] = val
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

    # jsasaki 2026-11-16: review newregion r2 #6: per pair (the maximum hid newly failing pairs)
    def area_change(pk):
        _mn, _mx, tw = _tri_quality(mesh.nodes[mesh.elements[pairs[pk].ravel()]])
        ar = 0.5 * np.abs(tw).reshape(-1, 2)
        return (ar.max(axis=1) - ar.min(axis=1)) / ar.max(axis=1)

    restored = refused = 0
    # jsasaki 2026-11-16: review newregion r3 #2: open-boundary nodes first, flip test only
    obc_o = state.get("obc", np.zeros(len(state["xy"]), bool))
    cand = np.nonzero(known)[0]
    cand = np.concatenate([cand[obc_o[orig[cand]]], cand[~obc_o[orig[cand]]]])
    for v in cand:
        o = orig[v]
        is_obc = bool(obc_o[o])
        if not (state["band"][o] or is_obc) or np.array_equal(mesh.nodes[v], state["xy"][o]):
            continue
        inc = n2e.get(int(v), [])
        els = mesh.elements[inc]
        pk = sorted({k for e in inc for k in e2p.get(e, [])})
        mn0, _mx0, _tw0 = _tri_quality(mesh.nodes[els])
        c40 = area_change(pk) if pk else np.zeros(0)
        old = mesh.nodes[v].copy()
        mesh.nodes[v] = state["xy"][o]
        mn1, _mx1, tw1 = _tri_quality(mesh.nodes[els])
        c41 = area_change(pk) if pk else np.zeros(0)
        # jsasaki 2026-11-16: review newregion r3 #2: OBC nodes get the same quality gate (a forced
        # restore kept the west junction at 25.9 deg, VR job 7592228); the obc stages re-derive them
        if (tw1 <= 0).any() or mn1.min() < min(mn0.min(), min_angle) \
                or (c41 > np.maximum(c40, max_area_change) + 1e-12).any():
            mesh.nodes[v] = old
            refused += 1
        else:
            restored += 1
    info = {"restored": restored, "refused": refused}
    log(f"[finish] open-boundary band restored: {info}")
    return mesh, info


# jsasaki 2026-11-16: review newregion r1 #1: nodes ON the protected outline (channel banks)
# are protected too: strict containment left a bank-to-bank channel unprotected
PROTECT_TOL_M = 1.0


def protected_node_mask(nodes: np.ndarray, protect_polygons, tol_m: float = PROTECT_TOL_M):
    """Nodes inside ``protect_polygons`` or within ``tol_m`` of its outline
    (a channel meshed bank to bank has every node ON the outline)."""
    import shapely

    shapely.prepare(protect_polygons)
    nodes = np.asarray(nodes, float)
    inside = shapely.contains_xy(protect_polygons, nodes[:, 0], nodes[:, 1])
    out = ~inside
    if out.any():
        pts = shapely.points(nodes[out, 0], nodes[out, 1])
        inside[out] = shapely.dwithin(protect_polygons, pts, tol_m)
    return inside


# jsasaki 2026-11-16: review newregion r2 #1: elements are protected by positive overlap with the
# protected water (a corridor narrower than the mesh spacing has no vertex in it)
def protected_element_mask(nodes: np.ndarray, elements: np.ndarray, protect_polygons,
                           tol_m: float = PROTECT_TOL_M, min_overlap_m2: float = 1e-6):
    """Elements with a protected vertex (:func:`protected_node_mask`) or a
    positive overlap area (> ``min_overlap_m2``) with ``protect_polygons``."""

    nodes = np.asarray(nodes, float)
    elements = np.asarray(elements, np.int64)
    out = protected_node_mask(nodes, protect_polygons, tol_m)[elements].any(axis=1)
    rest = np.nonzero(~out)[0]
    if rest.size:
        out[rest[_tri_overlap(nodes, elements[rest], protect_polygons) > min_overlap_m2]] = True
    return out


_TILE_CACHE: dict = {}


def _protect_tiles(protect_polygons, tile_m: float = 2000.0):
    """The protected geometry cut into ``tile_m`` squares (intersections with a
    large many-vertex polygon are slow; a 49k-element mesh took 67 s untiled)."""
    import shapely
    from shapely.geometry import box

    key = id(protect_polygons)
    hit = _TILE_CACHE.get(key)
    if hit is not None and hit[0] is protect_polygons:
        return hit[1], hit[2]
    x0, y0, x1, y1 = protect_polygons.bounds
    cells = [box(x, y, x + tile_m, y + tile_m)
             for x in np.arange(x0, x1, tile_m) for y in np.arange(y0, y1, tile_m)]
    pieces = shapely.intersection(np.array(cells, dtype=object), protect_polygons)
    pieces = pieces[~shapely.is_empty(pieces)]
    tree = shapely.STRtree(pieces)
    _TILE_CACHE.clear()
    _TILE_CACHE[key] = (protect_polygons, pieces, tree)
    return pieces, tree


def _tri_overlap(nodes, elements, protect_polygons) -> np.ndarray:
    """Per-element intersection area with ``protect_polygons`` (tiled)."""
    import shapely

    tris = shapely.polygons(np.asarray(nodes, float)[np.asarray(elements, np.int64)])
    pieces, tree = _protect_tiles(protect_polygons)
    ti, pj = tree.query(tris, predicate="intersects")
    out = np.zeros(len(tris))
    if ti.size:
        np.add.at(out, ti, shapely.area(shapely.intersection(tris[ti], pieces[pj])))
    return out


def protected_overlap_m2(nodes: np.ndarray, elements: np.ndarray, protect_polygons) -> float:
    """Total intersection area of the elements with ``protect_polygons``."""
    if len(elements) == 0:
        return 0.0
    return float(_tri_overlap(nodes, elements, protect_polygons).sum())


def _protected_node_mask(mesh: Fort14Mesh, protect_polygons):
    """Nodes in ``protect_polygons`` (None when there is no protection)."""
    if protect_polygons is None:
        return None
    # jsasaki 2026-11-16: review newregion r1 #1: boundary-inclusive
    return protected_node_mask(mesh.nodes, protect_polygons)


def _element_neighbour_count(elements: np.ndarray) -> np.ndarray:
    raw = np.sort(np.stack([elements[:, [0, 1]], elements[:, [1, 2]],
                            elements[:, [2, 0]]], axis=1), axis=2)
    n = int(elements.max()) + 1
    code = raw[..., 0].astype(np.int64) * n + raw[..., 1]
    uniq, inv, counts = np.unique(code.ravel(), return_inverse=True, return_counts=True)
    shared = (counts[inv] == 2).reshape(code.shape)
    return shared.sum(axis=1)


# jsasaki 2026-11-16: review newregion r2 #1: takes the protected ELEMENT mask
def _cap_deletion_allowed(mesh: Fort14Mesh, prot_elem) -> np.ndarray:
    """Per element: may a cap triangle be deleted here?

    Never when it holds an open-boundary node (it could break the open
    string). Never inside protected water: deleting a bridge can
    disconnect the element graph (the two-element-wide channel corridor that
    VR v2 lost), and deleting ears round after round eats a channel from its
    end (review newregion r1 #2); the optimizer fixes the angles there.
    """
    ok = np.ones(mesh.n_elements, bool)
    if any(len(s) for s in mesh.open_boundaries):
        obc = np.zeros(mesh.n_nodes, bool)
        obc[np.concatenate([np.asarray(s, int) for s in mesh.open_boundaries])] = True
        # jsasaki 2026-11-16: review newregion r3 #2: no element with an open-boundary node (its
        # deletion could take the node out of the string; finish validates the string at the end)
        ok &= ~obc[mesh.elements].any(axis=1)
    if prot_elem is not None:
        # jsasaki 2026-11-16: review newregion r1 #2: no cap deletion at all in protected water
        # jsasaki 2026-11-16: review newregion r2 #1: per element (overlap), not per vertex
        ok &= ~np.asarray(prot_elem, bool)
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


# jsasaki 2026-11-16: review newregion r2 #1: protect_polygons (element overlap) next to node_tag
def keep_components_protected(mesh: Fort14Mesh, node_tag=None, *, protect_polygons=None,
                              policy: str = "fail", log: Callable[[str], None] = print):
    """Keep the element components that touch the open boundary (the
    largest one when there is none); a dropped component holding protected
    water raises (``policy="fail"``) or is logged (``"report"``). Protected
    water: elements overlapping ``protect_polygons`` (outline included), or
    nodes tagged in ``node_tag``. Every dropped component is listed in
    ``info["dropped"]`` (elements, area, centroid, protected overlap). Open
    boundaries are remapped, not cleared. Returns ``(mesh, info, node_tag)``."""
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
    # jsasaki 2026-11-16: review newregion r2 #1: element overlap with the protected water
    prot_el = (protected_element_mask(mesh.nodes, mesh.elements, protect_polygons)
               if protect_polygons is not None else np.zeros(mesh.n_elements, bool))
    bad = []
    for c in range(n_comp):
        if c in keep:
            continue
        sel = lab == c
        rec = {"n_elements": int(sel.sum()), "area_m2": float(area[sel].sum()),
               "x": float(cen[sel, 0].mean()), "y": float(cen[sel, 1].mean()),
               # jsasaki 2026-11-16: review newregion r2 #1: element overlap or node tag
               "protected": bool(prot_el[sel].any() or (
                   node_tag is not None and node_tag[np.unique(mesh.elements[sel])].any()))}
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


# jsasaki 2026-11-16: review newregion r3 #2: lone node by edge flip or split (coverage unchanged)
def _flip_lone_nodes(mesh: Fort14Mesh) -> tuple[Fort14Mesh, int]:
    """A single-element node v of element (v, a, b) whose edge a-b is shared
    with (a, b, c): flip a-b to v-c, giving (v, a, c) and (v, c, b), when the
    quad is convex and a and b keep two elements each; otherwise split a-b at
    its midpoint m into (v, a, m), (v, m, b), (c, m, a), (c, b, m). Either way
    the covered water is unchanged and no single-element node is created.
    Returns ``(mesh, n_fixed)``."""
    els = [list(map(int, t)) for t in mesh.elements]
    nodes = [tuple(map(float, q)) for q in mesh.nodes]
    depths = list(map(float, mesh.depths))

    def area2(t):
        (x0, y0), (x1, y1), (x2, y2) = (nodes[t[0]], nodes[t[1]], nodes[t[2]])
        return (x1 - x0) * (y2 - y0) - (y1 - y0) * (x2 - x0)

    def oriented(t, sign):
        return t if np.sign(area2(t)) == sign else [t[0], t[2], t[1]]

    n_fix = 0
    n0 = int((np.bincount(mesh.elements.ravel(), minlength=mesh.n_nodes) == 1).sum())
    for _ in range(2 * n0 + 1):
        arr = np.asarray(els, np.int64)
        cnt = np.bincount(arr.ravel(), minlength=len(nodes))
        lone = [int(v) for v in np.nonzero(cnt == 1)[0]]
        fixed = False
        for v in lone:
            e0 = int(np.nonzero((arr == v).any(axis=1))[0][0])
            a, b = [x for x in els[e0] if x != v]
            nb = [e for e in np.nonzero((arr == a).any(axis=1) & (arr == b).any(axis=1))[0]
                  if e != e0]
            if len(nb) != 1:
                continue
            e1 = int(nb[0])
            c = next(x for x in els[e1] if x not in (a, b))
            s0 = np.sign(area2(els[e0]))
            t1, t2 = oriented([v, a, c], s0), oriented([v, c, b], s0)
            convex = (np.sign(area2([v, a, c])) == np.sign(area2([v, a, b]))
                      == np.sign(area2([v, c, b])))
            if convex and cnt[a] >= 3 and cnt[b] >= 3 and abs(area2(t1)) > 1e-9 \
                    and abs(area2(t2)) > 1e-9:
                els[e0], els[e1] = t1, t2
            else:
                m = len(nodes)
                nodes.append(tuple(0.5 * (np.asarray(nodes[a]) + np.asarray(nodes[b]))))
                depths.append(0.5 * (depths[a] + depths[b]))
                s1 = np.sign(area2(els[e1]))
                els[e0] = oriented([v, a, m], s0)
                els.append(oriented([v, m, b], s0))
                els[e1] = oriented([c, m, a], s1)
                els.append(oriented([c, b, m], s1))
            n_fix += 1
            fixed = True
            break       # topology changed: recount
        if not fixed:
            break
    if not n_fix:
        return mesh, 0
    out = _snapshot(mesh)
    out.nodes = np.asarray(nodes, float)
    out.depths = np.asarray(depths, float)
    out.elements = np.asarray(els, mesh.elements.dtype)
    return out, n_fix


#---> jsasaki 2026-11-16: review newregion r1 #2: protected water and open-boundary elements are
# never deleted; lone nodes left there fail (policy "fail") or are reported
def remove_lone_nodes(mesh: Fort14Mesh, *, protect_polygons=None, max_rounds: int = 10,
                      policy: str = "fail", log: Callable[[str], None] = print):
    """Delete the elements of nodes used by one element only (off the open
    boundary), until none is left. Such an element is an ear (it has at most
    one neighbour), so its deletion never disconnects the mesh. An element
    holding a protected node (``protect_polygons``, outline included) or an
    open-boundary node is never deleted; single-element nodes left for that
    reason raise (``policy="fail"``) or are reported (``"report"``). Returns
    ``(mesh, info)`` with the deleted count and area and the nodes kept."""
    if policy not in ("fail", "report"):
        raise ValueError(f"policy must be 'fail' or 'report', not {policy!r}")
    # jsasaki 2026-11-16: review newregion r3 #2: an edge flip first (no water lost; VR job 7592034
    # left two single-element nodes in protected water that deletion may not touch)
    mesh, n_flip = _flip_lone_nodes(mesh)
    n_del, area = 0, 0.0
    for _ in range(max_rounds):
        cnt = np.bincount(mesh.elements.ravel(), minlength=mesh.n_nodes)
        lone = cnt == 1
        obc = np.zeros(mesh.n_nodes, bool)
        for s in mesh.open_boundaries:
            obc[np.asarray(s, int)] = True
        lone &= ~obc
        kill = lone[mesh.elements].any(axis=1) & ~obc[mesh.elements].any(axis=1)
        # jsasaki 2026-11-16: review newregion r2 #1: protected by element overlap
        if protect_polygons is not None:
            kill &= ~protected_element_mask(mesh.nodes, mesh.elements, protect_polygons)
        if not kill.any():
            break
        _mn, _mx, tw = _tri_quality(mesh.nodes[mesh.elements[kill]])
        area += float(0.5 * np.abs(tw).sum())
        n_del += int(kill.sum())
        mesh, _ = remove_elements_keep_boundaries(mesh, ~kill)
    cnt = np.bincount(mesh.elements.ravel(), minlength=mesh.n_nodes)
    left = cnt == 1
    for s in mesh.open_boundaries:
        left[np.asarray(s, int)] = False
    info = {"elements_deleted": n_del, "area_m2": area, "lone_kept": int(left.sum()),
            "flipped": n_flip}  # jsasaki 2026-11-16: review newregion r3 #2
    log(f"[finish] single-element nodes: {n_flip} fixed by an edge flip / split, "
        f"{n_del} elements deleted ({area:.0f} m2), "
        f"{int(left.sum())} kept (protected water / open boundary)")
    if left.any():
        where = "; ".join(f"({x:.0f}, {y:.0f})" for x, y in mesh.nodes[left][:5])
        msg = (f"{int(left.sum())} single-element node(s) in protected water or at the open "
               f"boundary were not removed: {where}")
        if policy == "fail":
            raise RuntimeError(msg)
        log(f"[finish] WARNING {msg}")
    return mesh, info
#<--- jsasaki 2026-11-16


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
