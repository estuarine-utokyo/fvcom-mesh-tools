"""Mesh-clean operations: prune disjoint pools, trim dead-end elements,
repair 1-cell-wide channels, balance over-connected node valence,
widen under-resolved channels, remove skewed elements, and smooth
interior nodes.

Seven phases, applied in order:

    Phase A — :func:`keep_components`
        Drop dual-graph connected components by size and / or whether
        they touch the open boundary. Default keeps only the largest
        component (the standard "remove disjoint wet pools" repair).

    Phase B — :func:`trim_dead_ends`
        Iteratively delete degree-1 elements that have no
        open-boundary edge. Each round can expose new dead-ends;
        ``max_iters`` caps the loop.

    Phase C — :func:`repair_thin_chains`
        Detect chains of "thin" elements (triangles whose three
        vertices are all on a boundary) of at least
        ``min_thin_chain`` length — the 1-cell-wide-channel
        signature — and either widen each member by inserting a
        centroid (default), turning every thin triangle into 3
        sub-triangles with one interior vertex, or delete the whole
        chain.

    Phase D — :func:`repair_overconnected_nodes`
        Greedy Lawson-style edge swap that strictly reduces the
        per-edge "valence excess"
        ``sum_{n in {i,j,k,m}} max(0, valence(n) - max_nbr_elem)``.
        Off by default; PoC #27 found the FVCOM-safe 20° quality
        floor rejects every candidate on real meshes, so the floor
        defaults to 0° (only triangle inversion forbidden) when the
        phase is enabled.

    Phase E — :func:`repair_under_resolved_channels`
        Widen elements flagged by the medial-axis-style channel-width
        detector (detector 6, ``w/h < min_w_h``) using the same
        centroid-insertion mechanism as Phase C-widen. Catches
        2- and 3-cell-wide narrowed channels that Phase C does not
        flag. Off by default — detector 6 typically flags thousands
        of elements on real meshes, so enable it deliberately.

    Phase F — :func:`repair_skewed_elements`
        Delete triangles whose minimum interior angle is below
        ``min_angle_deg`` or whose maximum is at or above
        ``max_angle_deg``. Wraps ``ocsmesh.utils.cleanup_skewed_el``;
        no element is improved, only removed. Off by default — slivers
        are sometimes load-bearing near the boundary, and removing
        them silently can introduce holes the user didn't expect.

    Phase G — :func:`smooth_mesh_laplacian`
        Move every non-boundary node to the average of its connected
        neighbours (Laplacian smoothing). Wraps ``oceanmesh.laplacian2``,
        which automatically pins all topological boundary nodes; the
        connectivity, depths, and boundary lists are unchanged. Off by
        default — smoothing improves angle distribution but moves
        nodes away from the depth values they were originally
        interpolated at, so re-interpolating depths from the source
        DEM is an out-of-band step.

After every Phase A / B / C-delete / D / F deletion or topology
change the boundaries are re-derived via
:func:`fvcom_mesh_tools.algorithms.classify_boundaries_by_bbox`, so
the output mesh has a consistent open / land segmentation in the
same style as ``fmesh-buildmesh``. Phase C-widen, Phase D, and
Phase E-widen do not change the boundary edge set, so their
re-derive is a no-op.
"""

from __future__ import annotations

from typing import Any, Literal

import numpy as np
from scipy.sparse.csgraph import connected_components

from fvcom_mesh_tools.algorithms import (
    classify_boundaries_by_bbox,
    swap_edges_for_valence,
)
from fvcom_mesh_tools.diagnostics import (
    DEFAULT_ARC_SEPARATION_FACTOR,
    DEFAULT_CHANNEL_SAMPLE_DS_M,
    DEFAULT_MIN_W_H,
    DEFAULT_OPPOSITE_BANK_COS_MAX,
    dead_end_elements_flag,
    face_face_adjacency,
    open_boundary_node_mask,
    thin_chain_elements_flag,
    thin_elements_flag,
    under_resolved_channels_flag,
)
from fvcom_mesh_tools.io import Fort14Mesh

ThinChainMode = Literal["widen", "delete", "none"]
UnderResolvedMode = Literal["widen", "delete", "none"]

EARTH_R_M: float = 6_371_000.0

# Default bbox-classification tolerance in metres. 150 m matches
# ``fmesh-buildmesh``'s ``0.75 * default-hmin`` (hmin defaults to 200 m).
DEFAULT_BBOX_TOL_M: float = 150.0


def _deg_per_metre(lat_deg: float) -> float:
    """Conservative degrees-per-metre at ``lat_deg`` (longitude direction)."""
    return 1.0 / (EARTH_R_M * np.cos(np.deg2rad(lat_deg)) * np.pi / 180.0)


# ---------------------------------------------------------------------------
# Primitives
# ---------------------------------------------------------------------------


def remove_elements(
    mesh: Fort14Mesh, keep_elem_mask: np.ndarray,
) -> Fort14Mesh:
    """Drop elements where ``keep_elem_mask`` is False, drop unused nodes,
    renumber surviving elements / depths, and reset the boundary lists.

    Boundary lists are intentionally cleared; callers should re-derive
    them via :func:`rebuild_boundaries` once the deletion sequence is
    complete.
    """
    keep_elem_mask = np.asarray(keep_elem_mask, dtype=bool)
    if keep_elem_mask.shape != (mesh.n_elements,):
        raise ValueError(
            f"keep_elem_mask shape {keep_elem_mask.shape} != ({mesh.n_elements},)"
        )
    new_elements = mesh.elements[keep_elem_mask]
    if new_elements.size == 0:
        return Fort14Mesh(
            title=mesh.title,
            nodes=np.empty((0, 2), dtype=mesh.nodes.dtype),
            depths=np.empty((0,), dtype=mesh.depths.dtype),
            elements=np.empty((0, 3), dtype=mesh.elements.dtype),
            open_boundaries=[],
            land_boundaries=[],
        )
    used = np.unique(new_elements)
    remap = np.full(mesh.n_nodes, -1, dtype=np.int64)
    remap[used] = np.arange(len(used), dtype=np.int64)
    new_elements = remap[new_elements].astype(mesh.elements.dtype)
    return Fort14Mesh(
        title=mesh.title,
        nodes=mesh.nodes[used].copy(),
        depths=mesh.depths[used].copy(),
        elements=new_elements,
        open_boundaries=[],
        land_boundaries=[],
    )


def rebuild_boundaries(
    mesh: Fort14Mesh,
    *,
    bbox: tuple[float, float, float, float],
    tol_deg: float,
    land_ibtype: int = 0,
    open_merge_coast_gap: int = 0,
) -> Fort14Mesh:
    """Return ``mesh`` with open / land boundaries re-derived from the
    surviving outer ring via bbox proximity. Nodes / elements / depths
    are unchanged.
    """
    open_segs, land_bnds = classify_boundaries_by_bbox(
        mesh,
        bbox=bbox,
        tol=tol_deg,
        land_ibtype=land_ibtype,
        open_merge_coast_gap=open_merge_coast_gap,
    )
    return Fort14Mesh(
        title=mesh.title,
        nodes=mesh.nodes,
        depths=mesh.depths,
        elements=mesh.elements,
        open_boundaries=open_segs,
        land_boundaries=land_bnds,
    )


# ---------------------------------------------------------------------------
# Phase A: keep components
# ---------------------------------------------------------------------------


def keep_components(
    mesh: Fort14Mesh,
    *,
    min_elements: int | None = None,
    require_open_boundary: bool = False,
) -> tuple[Fort14Mesh, dict[str, Any]]:
    """Drop dual-graph connected components.

    Selection policy (in order of precedence):

        * ``require_open_boundary=True`` — keep components that contain
          at least one open-boundary node.
        * ``min_elements=N``           — keep components whose element
          count is >= ``N``.
        * default                       — keep only the single largest
          component.

    The two filters compose: ``require_open_boundary=True,
    min_elements=N`` keeps components that both touch the open boundary
    and have >= N elements. As a safety net, if the policy would empty
    the mesh, the largest component is kept instead.

    The returned mesh has its boundaries reset to empty; pair this with
    :func:`rebuild_boundaries` to get a usable Fort14Mesh.
    """
    if mesh.n_elements == 0:
        return mesh, {
            "n_components_before": 0,
            "n_components_kept": 0,
            "kept_component_sizes": [],
            "all_component_sizes": [],
            "n_elements_removed": 0,
            "n_nodes_removed": 0,
        }
    adj = face_face_adjacency(mesh.elements)
    n_comp, labels = connected_components(adj, directed=False, return_labels=True)
    sizes = np.bincount(labels, minlength=n_comp)

    if require_open_boundary:
        ob_mask = open_boundary_node_mask(mesh)
        elem_has_ob = ob_mask[mesh.elements].any(axis=1)
        candidates = {int(c) for c in np.unique(labels[elem_has_ob])}
        if min_elements is not None:
            candidates = {c for c in candidates if sizes[c] >= min_elements}
        if not candidates:
            candidates = {int(sizes.argmax())}
    elif min_elements is not None:
        candidates = {int(c) for c in range(n_comp) if sizes[c] >= min_elements}
        if not candidates:
            candidates = {int(sizes.argmax())}
    else:
        candidates = {int(sizes.argmax())}

    keep_mask = np.isin(labels, list(candidates))
    new_mesh = remove_elements(mesh, keep_mask)
    info: dict[str, Any] = {
        "n_components_before": int(n_comp),
        "n_components_kept": int(len(candidates)),
        "kept_component_sizes": sorted(
            [int(sizes[c]) for c in candidates], reverse=True,
        ),
        "all_component_sizes": [int(s) for s in sizes.tolist()],
        "n_elements_removed": int((~keep_mask).sum()),
        "n_nodes_removed": int(mesh.n_nodes - new_mesh.n_nodes),
    }
    return new_mesh, info


# ---------------------------------------------------------------------------
# Phase B: trim dead-end elements
# ---------------------------------------------------------------------------


def trim_dead_ends(
    mesh: Fort14Mesh,
    *,
    max_iters: int = 10,
    bbox: tuple[float, float, float, float],
    tol_deg: float,
    land_ibtype: int = 0,
    open_merge_coast_gap: int = 0,
) -> tuple[Fort14Mesh, dict[str, Any]]:
    """Iteratively delete degree-1 elements that have no open-boundary edge.

    Each iteration recomputes the dual-graph adjacency, applies
    :func:`fvcom_mesh_tools.diagnostics.dead_end_elements_flag` (which
    uses ``mesh.open_boundaries`` to identify OB edges), deletes the
    flagged elements, and re-derives boundaries via bbox proximity so
    the next iteration's OB edge set is up to date. The loop stops when
    a round flags zero dead-ends or ``max_iters`` is reached.

    Pre-condition: ``mesh`` should have its open boundary populated. If
    it does not, every degree-1 element looks like a dead-end (no OB
    edges to filter against), so the caller should run
    :func:`rebuild_boundaries` first.
    """
    history: list[int] = []
    cur = mesh
    for _ in range(max_iters):
        if cur.n_elements == 0:
            break
        adj = face_face_adjacency(cur.elements)
        flag = dead_end_elements_flag(adj, cur)
        n = int(flag.sum())
        history.append(n)
        if n == 0:
            break
        cur = remove_elements(cur, ~flag)
        if cur.n_elements == 0:
            break
        cur = rebuild_boundaries(
            cur,
            bbox=bbox,
            tol_deg=tol_deg,
            land_ibtype=land_ibtype,
            open_merge_coast_gap=open_merge_coast_gap,
        )
    converged = (len(history) == 0) or (history[-1] == 0)
    info: dict[str, Any] = {
        "per_iter_dead_end_count": history,
        "iterations_run": len(history),
        "converged": bool(converged),
        "total_elements_removed": int(mesh.n_elements - cur.n_elements),
    }
    return cur, info


# ---------------------------------------------------------------------------
# Phase C: repair thin chains (widen via centroid insert, or delete)
# ---------------------------------------------------------------------------


def _detect_thin_chain_flag(
    mesh: Fort14Mesh, *, min_chain_length: int,
) -> np.ndarray:
    """Element-level mask: thin elements that are part of a connected
    run of at least ``min_chain_length`` adjacent thin elements.

    Same definition as :func:`fvcom_mesh_tools.diagnostics.thin_chain_elements_flag`.
    """
    thin = thin_elements_flag(mesh)
    if not thin.any() or min_chain_length <= 1:
        return thin
    adj = face_face_adjacency(mesh.elements)
    return thin_chain_elements_flag(adj, thin, min_chain_length=min_chain_length)


def widen_thin_elements_at_centroid(
    mesh: Fort14Mesh, target_flag: np.ndarray,
) -> tuple[Fort14Mesh, dict[str, Any]]:
    """Replace each flagged element with three sub-triangles fanning out
    from its centroid.

    For each flagged triangle ``(P0, P1, P2)`` we add a new node ``C``
    at the triangle centroid (depth = mean of the 3 vertex depths) and
    emit three sub-triangles ``(P0, P1, C)``, ``(P1, P2, C)``,
    ``(P2, P0, C)``. This widens any "all-3-vertices-on-boundary"
    element into a fan with a strictly interior centroid, giving
    cross-channel resolution of two cells where there was previously
    one.

    Boundary lists are preserved unchanged: centroid insertion only
    adds interior nodes and edges; the existing boundary edge set is
    not modified, and previously-stored boundary node IDs remain
    valid because the new node IDs are appended.
    """
    target_flag = np.asarray(target_flag, dtype=bool)
    if target_flag.shape != (mesh.n_elements,):
        raise ValueError(
            f"target_flag shape {target_flag.shape} != ({mesh.n_elements},)"
        )
    n_widen = int(target_flag.sum())
    info: dict[str, Any] = {
        "n_widened": n_widen,
        "n_new_nodes": 0,
        "n_new_elements": 0,
    }
    if n_widen == 0:
        return mesh, info

    target_idx = np.where(target_flag)[0]
    target_elems = mesh.elements[target_idx]
    centroids = mesh.nodes[target_elems].mean(axis=1)
    centroid_depths = mesh.depths[target_elems].mean(axis=1)

    new_node_ids = np.arange(
        mesh.n_nodes, mesh.n_nodes + n_widen, dtype=mesh.elements.dtype,
    )
    e0 = target_elems[:, 0]
    e1 = target_elems[:, 1]
    e2 = target_elems[:, 2]
    sub_a = np.column_stack([e0, e1, new_node_ids])
    sub_b = np.column_stack([e1, e2, new_node_ids])
    sub_c = np.column_stack([e2, e0, new_node_ids])
    new_subtris = np.vstack([sub_a, sub_b, sub_c]).astype(mesh.elements.dtype)

    keep_mask = ~target_flag
    new_elements = np.vstack([mesh.elements[keep_mask], new_subtris])
    new_nodes = np.vstack([mesh.nodes, centroids])
    new_depths = np.concatenate([mesh.depths, centroid_depths])

    info["n_new_nodes"] = n_widen
    info["n_new_elements"] = int(new_elements.shape[0] - mesh.n_elements)

    return Fort14Mesh(
        title=mesh.title,
        nodes=new_nodes,
        depths=new_depths,
        elements=new_elements,
        open_boundaries=[np.asarray(s).copy() for s in mesh.open_boundaries],
        land_boundaries=[(int(ib), np.asarray(s).copy())
                         for ib, s in mesh.land_boundaries],
    ), info


def repair_thin_chains(
    mesh: Fort14Mesh,
    *,
    mode: ThinChainMode = "widen",
    min_chain_length: int = 3,
    bbox: tuple[float, float, float, float] | None = None,
    tol_deg: float | None = None,
    land_ibtype: int = 0,
    open_merge_coast_gap: int = 0,
) -> tuple[Fort14Mesh, dict[str, Any]]:
    """Apply Phase C: widen or delete 1-cell-wide channels.

    ``mode='widen'`` (default) inserts a centroid in every thin-chain
    element; boundaries are preserved. ``mode='delete'`` removes every
    thin-chain element and re-derives boundaries via bbox proximity
    (``bbox`` and ``tol_deg`` must be provided). ``mode='none'`` is a
    no-op apart from emitting a small info dict.
    """
    if mode == "none":
        return mesh, {"mode": "none", "n_chain_elements": 0, "skipped": True}

    chain_flag = _detect_thin_chain_flag(mesh, min_chain_length=min_chain_length)
    n_chain = int(chain_flag.sum())
    info: dict[str, Any] = {
        "mode": mode,
        "min_chain_length": int(min_chain_length),
        "n_chain_elements": n_chain,
    }
    if n_chain == 0:
        info["skipped"] = True
        return mesh, info

    if mode == "widen":
        new_mesh, w_info = widen_thin_elements_at_centroid(mesh, chain_flag)
        info.update(w_info)
        return new_mesh, info

    if mode == "delete":
        if bbox is None or tol_deg is None:
            raise ValueError(
                "delete mode requires bbox and tol_deg for boundary rebuild"
            )
        new_mesh = remove_elements(mesh, ~chain_flag)
        new_mesh = rebuild_boundaries(
            new_mesh, bbox=bbox, tol_deg=tol_deg,
            land_ibtype=land_ibtype,
            open_merge_coast_gap=open_merge_coast_gap,
        )
        info["n_elements_removed"] = n_chain
        info["n_nodes_removed"] = int(mesh.n_nodes - new_mesh.n_nodes)
        return new_mesh, info

    raise ValueError(f"unknown thin_chain_mode: {mode!r}")


# ---------------------------------------------------------------------------
# Phase D: repair over-connected nodes via valence-balancing edge swap
# ---------------------------------------------------------------------------


def repair_overconnected_nodes(
    mesh: Fort14Mesh,
    *,
    max_nbr_elem: int = 8,
    max_iters: int = 50,
    min_angle_floor_deg: float = 0.0,
) -> tuple[Fort14Mesh, dict[str, Any]]:
    """Phase D: drive every node valence to at most ``max_nbr_elem`` via
    iterative Lawson edge flips that strictly decrease total valence
    excess.

    Boundary edges are excluded from the candidate list, so this never
    changes the open / land segmentation. ``min_angle_floor_deg``
    defaults to 0 (only triangle inversion forbidden); raise to e.g.
    20 to forbid sliver creation, at the cost of typically rejecting
    every candidate on fan-like local topology.
    """
    out, info = swap_edges_for_valence(
        mesh,
        max_nbr_elem=max_nbr_elem,
        max_iters=max_iters,
        min_angle_floor_deg=min_angle_floor_deg,
    )
    return out, info


# ---------------------------------------------------------------------------
# Phase E: widen under-resolved channels (detector 6 repair)
# ---------------------------------------------------------------------------


def repair_under_resolved_channels(
    mesh: Fort14Mesh,
    *,
    mode: UnderResolvedMode = "widen",
    min_w_h: float = DEFAULT_MIN_W_H,
    sample_ds_m: float = DEFAULT_CHANNEL_SAMPLE_DS_M,
    arc_separation_factor: float = DEFAULT_ARC_SEPARATION_FACTOR,
    opposite_bank_cos_max: float = DEFAULT_OPPOSITE_BANK_COS_MAX,
    bbox: tuple[float, float, float, float] | None = None,
    tol_deg: float | None = None,
    land_ibtype: int = 0,
    open_merge_coast_gap: int = 0,
) -> tuple[Fort14Mesh, dict[str, Any]]:
    """Apply Phase E: widen or delete under-resolved channel elements.

    The target flag is :func:`under_resolved_channels_flag` (detector 6)
    with the supplied ``min_w_h`` threshold and same medial-axis
    parameters used by ``fmesh-mesh-check``. ``mode='widen'`` (default)
    inserts a centroid in every flagged element; boundaries are
    preserved. ``mode='delete'`` removes the flagged elements and
    re-derives boundaries via bbox proximity (``bbox`` and ``tol_deg``
    must be provided). ``mode='none'`` is a no-op.

    Centroid insertion shrinks each new sub-triangle's median edge
    length to ~0.577 × the parent's, while the geometric channel
    width (set by the bank polylines) is unchanged. So the post-widen
    w/h ratio is ~1.73 × the pre-widen ratio: borderline-flagged
    elements (ratio just below ``min_w_h``) cross the threshold, but
    very narrow channels (ratio well below) stay flagged. Phase E is
    therefore best read as "lift local resolution one step" rather
    than "guarantee 3 cells across every narrow channel"; the latter
    needs medial-axis insertion which is outside this driver.
    Empirically on the PoC #19 cleaned Tokyo-Bay mesh: 3,178 → 3,032
    flagged (4.6 % reduction); 6.7 % → 5.6 % per-mesh flagged
    fraction.

    Detector 6 typically flags thousands of elements on real meshes,
    so ``widen`` adds one interior node and two net elements per
    flagged element — the output mesh size grows proportionally.
    ``delete`` is useful for stripping cosmetic narrow inlets but is
    destructive on meshes where the inlets matter.
    """
    if mode == "none":
        return mesh, {"mode": "none", "n_flagged": 0, "skipped": True}

    flag, _metric = under_resolved_channels_flag(
        mesh,
        min_w_h=min_w_h,
        sample_ds_m=sample_ds_m,
        arc_separation_factor=arc_separation_factor,
        opposite_bank_cos_max=opposite_bank_cos_max,
    )
    n_flagged = int(flag.sum())
    info: dict[str, Any] = {
        "mode": mode,
        "min_w_h": float(min_w_h),
        "sample_ds_m": float(sample_ds_m),
        "n_flagged": n_flagged,
    }
    if n_flagged == 0:
        info["skipped"] = True
        return mesh, info

    if mode == "widen":
        new_mesh, w_info = widen_thin_elements_at_centroid(mesh, flag)
        info.update(w_info)
        return new_mesh, info

    if mode == "delete":
        if bbox is None or tol_deg is None:
            raise ValueError(
                "delete mode requires bbox and tol_deg for boundary rebuild"
            )
        new_mesh = remove_elements(mesh, ~flag)
        new_mesh = rebuild_boundaries(
            new_mesh, bbox=bbox, tol_deg=tol_deg,
            land_ibtype=land_ibtype,
            open_merge_coast_gap=open_merge_coast_gap,
        )
        info["n_elements_removed"] = n_flagged
        info["n_nodes_removed"] = int(mesh.n_nodes - new_mesh.n_nodes)
        return new_mesh, info

    raise ValueError(f"unknown under_resolved_mode: {mode!r}")


# ---------------------------------------------------------------------------
# Phase F: angle-based skewed-element removal (wraps ocsmesh.utils)
# ---------------------------------------------------------------------------


# Defaults match ``ocsmesh.utils.cleanup_skewed_el``'s triangle thresholds.
DEFAULT_SKEWED_MIN_ANGLE_DEG: float = 1.0
DEFAULT_SKEWED_MAX_ANGLE_DEG: float = 175.0


def repair_skewed_elements(
    mesh: Fort14Mesh,
    *,
    min_angle_deg: float = DEFAULT_SKEWED_MIN_ANGLE_DEG,
    max_angle_deg: float = DEFAULT_SKEWED_MAX_ANGLE_DEG,
    bbox: tuple[float, float, float, float] | None = None,
    tol_deg: float | None = None,
    land_ibtype: int = 0,
    open_merge_coast_gap: int = 0,
) -> tuple[Fort14Mesh, dict[str, Any]]:
    """Phase F: delete triangles whose minimum interior angle is below
    ``min_angle_deg`` or whose maximum is at or above ``max_angle_deg``.

    Wraps :func:`ocsmesh.utils.cleanup_skewed_el`. ocsmesh is a
    library-only dependency here (no gmsh involved); see
    ``docs/engine_complementarity.md`` for the rationale.

    Boundaries are re-derived via bbox proximity if any element is
    removed; the caller therefore must provide ``bbox`` and
    ``tol_deg`` whenever Phase F may delete. (If you know in advance
    that no triangle will be skewed, you can pass ``bbox=None`` and
    skip the rebuild — the function detects the no-op case.)

    Returns ``(new_mesh, info)`` where ``info`` includes the angle
    thresholds used, the count of deleted elements, and the count of
    nodes dropped.
    """
    if min_angle_deg < 0 or max_angle_deg > 180 or min_angle_deg >= max_angle_deg:
        raise ValueError(
            "expected 0 <= min_angle_deg < max_angle_deg <= 180, "
            f"got [{min_angle_deg}, {max_angle_deg}]"
        )
    if mesh.n_elements == 0:
        return mesh, {
            "min_angle_deg": float(min_angle_deg),
            "max_angle_deg": float(max_angle_deg),
            "n_elements_removed": 0,
            "n_nodes_removed": 0,
            "skipped": True,
        }

    from ocsmesh import utils as ocsmesh_utils

    from fvcom_mesh_tools.mesh_compose.convert import (
        fort14_to_meshdata,
        meshdata_to_fort14,
    )

    md_in = fort14_to_meshdata(mesh)
    md_out = ocsmesh_utils.cleanup_skewed_el(
        md_in,
        lw_bound_tri=float(min_angle_deg),
        up_bound_tri=float(max_angle_deg),
    )
    cleaned = meshdata_to_fort14(md_out, title=mesh.title)

    n_removed = int(mesh.n_elements - cleaned.n_elements)
    info: dict[str, Any] = {
        "min_angle_deg": float(min_angle_deg),
        "max_angle_deg": float(max_angle_deg),
        "n_elements_removed": n_removed,
        "n_nodes_removed": int(mesh.n_nodes - cleaned.n_nodes),
    }
    if n_removed == 0:
        # Topology unchanged — preserve the input's boundary lists so
        # the no-op stays a no-op even when bbox isn't supplied.
        info["skipped"] = True
        return Fort14Mesh(
            title=mesh.title,
            nodes=cleaned.nodes,
            depths=cleaned.depths,
            elements=cleaned.elements,
            open_boundaries=[np.asarray(s).copy() for s in mesh.open_boundaries],
            land_boundaries=[(int(ib), np.asarray(s).copy())
                             for ib, s in mesh.land_boundaries],
        ), info

    if bbox is None or tol_deg is None:
        raise ValueError(
            "Phase F removed elements but bbox / tol_deg are missing for "
            "the boundary rebuild. Pass them when enabling Phase F."
        )
    cleaned = rebuild_boundaries(
        cleaned, bbox=bbox, tol_deg=tol_deg,
        land_ibtype=land_ibtype,
        open_merge_coast_gap=open_merge_coast_gap,
    )
    return cleaned, info


# ---------------------------------------------------------------------------
# Phase G: Laplacian smoothing (wraps oceanmesh.laplacian2)
# ---------------------------------------------------------------------------


# Defaults match ``oceanmesh.laplacian2`` exactly.
DEFAULT_SMOOTH_LAPLACIAN_ITERS: int = 20
DEFAULT_SMOOTH_LAPLACIAN_TOL: float = 0.01


def smooth_mesh_laplacian(
    mesh: Fort14Mesh,
    *,
    max_iter: int = DEFAULT_SMOOTH_LAPLACIAN_ITERS,
    tol: float = DEFAULT_SMOOTH_LAPLACIAN_TOL,
    pfix: np.ndarray | None = None,
) -> tuple[Fort14Mesh, dict[str, Any]]:
    """Phase G: move every non-boundary node to the mean position of its
    connected neighbours (Laplacian smoothing).

    Wraps :func:`oceanmesh.laplacian2`. oceanmesh derives the boundary
    set from the mesh topology and pins it automatically, so the
    output's open / land boundary node coordinates and the boundary
    lists themselves are unchanged. Connectivity and depth indices
    are preserved (``laplacian2`` only moves vertices). ``pfix`` adds
    extra fixed coordinates beyond the topological boundary.

    .. note::
       Importing :mod:`oceanmesh` (GPL-3.0-or-later) propagates GPL
       into the redistributed combined work. See
       ``THIRD_PARTY_NOTICES.md`` and
       ``docs/engine_complementarity.md``. Callers without oceanmesh
       installed will see :class:`ImportError` at call time.

    .. note::
       Smoothing moves interior nodes away from the coordinates at
       which their depths were interpolated. For high-fidelity FVCOM
       runs you should re-interpolate depths from the source DEM
       after Phase G; the toolkit does not currently do this
       automatically.
    """
    if max_iter < 1:
        raise ValueError(f"max_iter must be >= 1, got {max_iter}")
    if tol <= 0:
        raise ValueError(f"tol must be > 0, got {tol}")
    if mesh.n_elements == 0:
        return mesh, {
            "max_iter": int(max_iter),
            "tol": float(tol),
            "n_pfix": 0,
            "n_nodes_moved": 0,
            "displacement_max": 0.0,
            "displacement_mean": 0.0,
            "skipped": True,
        }

    from oceanmesh import laplacian2  # GPL-3.0-or-later, lazy import

    pre = np.asarray(mesh.nodes, dtype=float).copy()
    new_vertices, _ = laplacian2(
        pre.copy(),
        np.asarray(mesh.elements, dtype=int),
        max_iter=int(max_iter),
        tol=float(tol),
        pfix=pfix,
    )
    new_vertices = np.asarray(new_vertices, dtype=float)

    disp = new_vertices - pre
    disp_norm = np.linalg.norm(disp, axis=1)
    eps = 1e-12

    info: dict[str, Any] = {
        "max_iter": int(max_iter),
        "tol": float(tol),
        "n_pfix": 0 if pfix is None else int(np.asarray(pfix).reshape(-1, 2).shape[0]),
        "n_nodes_moved": int((disp_norm > eps).sum()),
        "displacement_max": float(disp_norm.max()) if disp_norm.size else 0.0,
        "displacement_mean": float(disp_norm.mean()) if disp_norm.size else 0.0,
    }

    out = Fort14Mesh(
        title=mesh.title,
        nodes=new_vertices,
        depths=mesh.depths.copy(),
        elements=mesh.elements.copy(),
        open_boundaries=[np.asarray(s).copy() for s in mesh.open_boundaries],
        land_boundaries=[(int(ib), np.asarray(s).copy())
                         for ib, s in mesh.land_boundaries],
    )
    return out, info


# ---------------------------------------------------------------------------
# Driver: clean_mesh
# ---------------------------------------------------------------------------


def clean_mesh(
    mesh: Fort14Mesh,
    *,
    bbox: tuple[float, float, float, float],
    bbox_tol_m: float = DEFAULT_BBOX_TOL_M,
    land_ibtype: int = 0,
    open_merge_coast_gap: int = 0,
    remove_disjoint: bool = True,
    min_component_elements: int | None = None,
    require_open_boundary: bool = False,
    trim_dead_ends_iters: int = 10,
    thin_chain_mode: ThinChainMode = "widen",
    min_thin_chain: int = 3,
    repair_overconnected_iters: int = 0,
    max_nbr_elem: int = 8,
    overconn_min_angle_floor_deg: float = 0.0,
    under_resolved_mode: UnderResolvedMode = "none",
    under_resolved_min_w_h: float = DEFAULT_MIN_W_H,
    under_resolved_sample_ds_m: float = DEFAULT_CHANNEL_SAMPLE_DS_M,
    under_resolved_arc_separation_factor: float = DEFAULT_ARC_SEPARATION_FACTOR,
    under_resolved_opposite_bank_cos_max: float = DEFAULT_OPPOSITE_BANK_COS_MAX,
    repair_skewed: bool = False,
    repair_skewed_min_angle_deg: float = DEFAULT_SKEWED_MIN_ANGLE_DEG,
    repair_skewed_max_angle_deg: float = DEFAULT_SKEWED_MAX_ANGLE_DEG,
    smooth_laplacian: bool = False,
    smooth_laplacian_iters: int = DEFAULT_SMOOTH_LAPLACIAN_ITERS,
    smooth_laplacian_tol: float = DEFAULT_SMOOTH_LAPLACIAN_TOL,
) -> tuple[Fort14Mesh, dict[str, Any]]:
    """Run Phase A (component pruning) and Phase B (dead-end trimming).

    Boundaries are always re-derived via bbox proximity so the output
    mesh has a consistent open / land segmentation. If both phases are
    disabled the input mesh is still passed through bbox classification
    (an explicit no-op cleanup is still a useful normaliser).

    Parameters
    ----------
    bbox:
        ``(xmin, ymin, xmax, ymax)`` in the same coordinate system as
        ``mesh.nodes``. Boundary nodes within the latitude-converted
        ``bbox_tol_m`` of this rectangle are flagged as open.
    bbox_tol_m:
        Tolerance in metres; converted to degrees at the surviving
        mesh's mid-latitude before being passed into the bbox
        classifier. Default ``150`` m matches ``fmesh-buildmesh``.
    land_ibtype:
        ``ibtype`` written for re-derived land segments (default 0).
    open_merge_coast_gap:
        See :func:`classify_boundaries_by_bbox`.
    remove_disjoint:
        Run Phase A.
    min_component_elements, require_open_boundary:
        Phase A selection policy. Default keeps only the largest
        component.
    trim_dead_ends_iters:
        Phase B iteration cap. Set to 0 to skip Phase B.
    thin_chain_mode:
        Phase C policy. ``"widen"`` (default) inserts centroids into
        every thin-chain element so 1-cell channels become 2-cell;
        ``"delete"`` removes the chain entirely and re-derives
        boundaries; ``"none"`` skips Phase C.
    min_thin_chain:
        Minimum length of a connected thin-element run that we treat
        as a 1-cell-wide channel run for Phase C. Default 3, matching
        ``fmesh-mesh-check``.
    repair_overconnected_iters:
        Phase D iteration cap. ``0`` (default) disables Phase D; set
        to e.g. 50 to enable.
    max_nbr_elem:
        FVCOM ``MAX_NBR_ELEM`` cap that Phase D drives every valence
        to. Default 8 (matches ``fmesh-mesh-check``).
    overconn_min_angle_floor_deg:
        Phase D quality floor in degrees. Default 0 (only triangle
        inversion forbidden) — raise to e.g. 20 to forbid sliver
        creation, at the cost of typically rejecting every candidate
        on fan-like local topology (PoC #27).
    under_resolved_mode:
        Phase E policy for under-resolved channel elements (detector
        6). ``"widen"`` inserts a centroid in every flagged element
        so a 2-cell channel becomes 3-cell; ``"delete"`` removes the
        flagged elements and re-derives boundaries; ``"none"``
        (default) skips Phase E. Detector 6 typically flags thousands
        of elements on real meshes — enable deliberately.
    under_resolved_min_w_h:
        Threshold for Phase E. An element is flagged when its local
        channel width divided by the median edge length is below
        this value. Default ``DEFAULT_MIN_W_H`` matches
        ``fmesh-mesh-check``.
    under_resolved_sample_ds_m, under_resolved_arc_separation_factor,
    under_resolved_opposite_bank_cos_max:
        Medial-axis detector parameters. See
        :func:`channel_width_metric`.
    repair_skewed:
        Phase F switch (off by default). When True, delete triangles
        whose minimum interior angle is below
        ``repair_skewed_min_angle_deg`` or whose maximum is at or
        above ``repair_skewed_max_angle_deg`` via
        :func:`repair_skewed_elements` (a wrapper around
        :func:`ocsmesh.utils.cleanup_skewed_el`).
    repair_skewed_min_angle_deg, repair_skewed_max_angle_deg:
        Phase F thresholds in degrees. Defaults match
        ``ocsmesh.utils.cleanup_skewed_el``: ``[1.0, 175.0]``. Tighten
        (e.g. raise to 5° / lower to 170°) to remove milder skews at
        the cost of more aggressive deletion.
    smooth_laplacian:
        Phase G switch (off by default). When True, run Laplacian
        smoothing of all interior nodes via
        :func:`smooth_mesh_laplacian` (wraps
        :func:`oceanmesh.laplacian2`). Boundary nodes are auto-pinned
        by oceanmesh; connectivity and depth indices are preserved.
    smooth_laplacian_iters, smooth_laplacian_tol:
        Phase G iteration cap and convergence tolerance. Defaults
        match ``oceanmesh.laplacian2``: ``max_iter=20, tol=0.01``.
    """
    if thin_chain_mode not in ("widen", "delete", "none"):
        raise ValueError(
            f"thin_chain_mode must be one of widen / delete / none, got {thin_chain_mode!r}"
        )
    if under_resolved_mode not in ("widen", "delete", "none"):
        raise ValueError(
            "under_resolved_mode must be one of widen / delete / none, "
            f"got {under_resolved_mode!r}"
        )
    info: dict[str, Any] = {
        "input": {
            "n_nodes": int(mesh.n_nodes),
            "n_elements": int(mesh.n_elements),
            "n_open_boundaries": len(mesh.open_boundaries),
            "n_land_boundaries": len(mesh.land_boundaries),
        },
        "config": {
            "bbox": list(bbox),
            "bbox_tol_m": float(bbox_tol_m),
            "land_ibtype": int(land_ibtype),
            "open_merge_coast_gap": int(open_merge_coast_gap),
            "remove_disjoint": bool(remove_disjoint),
            "min_component_elements": min_component_elements,
            "require_open_boundary": bool(require_open_boundary),
            "trim_dead_ends_iters": int(trim_dead_ends_iters),
            "thin_chain_mode": thin_chain_mode,
            "min_thin_chain": int(min_thin_chain),
            "repair_overconnected_iters": int(repair_overconnected_iters),
            "max_nbr_elem": int(max_nbr_elem),
            "overconn_min_angle_floor_deg": float(overconn_min_angle_floor_deg),
            "under_resolved_mode": under_resolved_mode,
            "under_resolved_min_w_h": float(under_resolved_min_w_h),
            "under_resolved_sample_ds_m": float(under_resolved_sample_ds_m),
            "under_resolved_arc_separation_factor":
                float(under_resolved_arc_separation_factor),
            "under_resolved_opposite_bank_cos_max":
                float(under_resolved_opposite_bank_cos_max),
            "repair_skewed": bool(repair_skewed),
            "repair_skewed_min_angle_deg": float(repair_skewed_min_angle_deg),
            "repair_skewed_max_angle_deg": float(repair_skewed_max_angle_deg),
            "smooth_laplacian": bool(smooth_laplacian),
            "smooth_laplacian_iters": int(smooth_laplacian_iters),
            "smooth_laplacian_tol": float(smooth_laplacian_tol),
        },
        "phases": [],
    }
    cur = mesh

    def _tol_deg(m: Fort14Mesh) -> float:
        if m.n_nodes == 0:
            return 0.0
        return bbox_tol_m * _deg_per_metre(float(m.nodes[:, 1].mean()))

    if remove_disjoint:
        cur, p_info = keep_components(
            cur,
            min_elements=min_component_elements,
            require_open_boundary=require_open_boundary,
        )
        info["phases"].append({"name": "keep_components", **p_info})
        cur = rebuild_boundaries(
            cur,
            bbox=bbox,
            tol_deg=_tol_deg(cur),
            land_ibtype=land_ibtype,
            open_merge_coast_gap=open_merge_coast_gap,
        )

    if trim_dead_ends_iters > 0:
        if not cur.open_boundaries:
            cur = rebuild_boundaries(
                cur,
                bbox=bbox,
                tol_deg=_tol_deg(cur),
                land_ibtype=land_ibtype,
                open_merge_coast_gap=open_merge_coast_gap,
            )
        cur, t_info = trim_dead_ends(
            cur,
            max_iters=trim_dead_ends_iters,
            bbox=bbox,
            tol_deg=_tol_deg(cur),
            land_ibtype=land_ibtype,
            open_merge_coast_gap=open_merge_coast_gap,
        )
        info["phases"].append({"name": "trim_dead_ends", **t_info})

    if thin_chain_mode != "none":
        cur, c_info = repair_thin_chains(
            cur,
            mode=thin_chain_mode,
            min_chain_length=min_thin_chain,
            bbox=bbox,
            tol_deg=_tol_deg(cur),
            land_ibtype=land_ibtype,
            open_merge_coast_gap=open_merge_coast_gap,
        )
        info["phases"].append({"name": "repair_thin_chains", **c_info})

    if repair_overconnected_iters > 0:
        cur, d_info = repair_overconnected_nodes(
            cur,
            max_nbr_elem=max_nbr_elem,
            max_iters=repair_overconnected_iters,
            min_angle_floor_deg=overconn_min_angle_floor_deg,
        )
        info["phases"].append({"name": "repair_overconnected_nodes", **d_info})

    if under_resolved_mode != "none":
        cur, e_info = repair_under_resolved_channels(
            cur,
            mode=under_resolved_mode,
            min_w_h=under_resolved_min_w_h,
            sample_ds_m=under_resolved_sample_ds_m,
            arc_separation_factor=under_resolved_arc_separation_factor,
            opposite_bank_cos_max=under_resolved_opposite_bank_cos_max,
            bbox=bbox,
            tol_deg=_tol_deg(cur),
            land_ibtype=land_ibtype,
            open_merge_coast_gap=open_merge_coast_gap,
        )
        info["phases"].append({"name": "repair_under_resolved_channels", **e_info})

    if repair_skewed:
        cur, f_info = repair_skewed_elements(
            cur,
            min_angle_deg=repair_skewed_min_angle_deg,
            max_angle_deg=repair_skewed_max_angle_deg,
            bbox=bbox,
            tol_deg=_tol_deg(cur),
            land_ibtype=land_ibtype,
            open_merge_coast_gap=open_merge_coast_gap,
        )
        info["phases"].append({"name": "repair_skewed_elements", **f_info})

    if smooth_laplacian:
        cur, g_info = smooth_mesh_laplacian(
            cur,
            max_iter=smooth_laplacian_iters,
            tol=smooth_laplacian_tol,
        )
        info["phases"].append({"name": "smooth_mesh_laplacian", **g_info})

    if not info["phases"]:
        # Every phase disabled — still re-derive boundaries so the
        # output is a normalised pass-through.
        cur = rebuild_boundaries(
            cur,
            bbox=bbox,
            tol_deg=_tol_deg(cur),
            land_ibtype=land_ibtype,
            open_merge_coast_gap=open_merge_coast_gap,
        )

    info["output"] = {
        "n_nodes": int(cur.n_nodes),
        "n_elements": int(cur.n_elements),
        "n_open_boundaries": len(cur.open_boundaries),
        "n_land_boundaries": len(cur.land_boundaries),
    }
    return cur, info


__all__ = [
    "DEFAULT_BBOX_TOL_M",
    "DEFAULT_SKEWED_MAX_ANGLE_DEG",
    "DEFAULT_SKEWED_MIN_ANGLE_DEG",
    "DEFAULT_SMOOTH_LAPLACIAN_ITERS",
    "DEFAULT_SMOOTH_LAPLACIAN_TOL",
    "ThinChainMode",
    "UnderResolvedMode",
    "clean_mesh",
    "keep_components",
    "rebuild_boundaries",
    "remove_elements",
    "repair_overconnected_nodes",
    "repair_skewed_elements",
    "repair_thin_chains",
    "repair_under_resolved_channels",
    "smooth_mesh_laplacian",
    "trim_dead_ends",
    "widen_thin_elements_at_centroid",
]
