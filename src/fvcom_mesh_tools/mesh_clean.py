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
UnderResolvedMode = Literal["widen", "delete", "medial", "none"]

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
    min_channel_elements: int = 1,
    bbox: tuple[float, float, float, float] | None = None,
    tol_deg: float | None = None,
    land_ibtype: int = 0,
    open_merge_coast_gap: int = 0,
    repair_flipped: bool = True,
    max_repair_passes: int = 5,
) -> tuple[Fort14Mesh, dict[str, Any]]:
    """Apply Phase E: widen, delete, or medial-axis-remesh
    under-resolved channel elements.

    The target flag is :func:`under_resolved_channels_flag` (detector 6)
    with the supplied ``min_w_h`` threshold and same medial-axis
    parameters used by ``fmesh-mesh-check``. ``mode='widen'`` (default)
    inserts a centroid in every flagged element; boundaries are
    preserved. ``mode='delete'`` removes the flagged elements and
    re-derives boundaries via bbox proximity (``bbox`` and ``tol_deg``
    must be provided). ``mode='medial'`` rebuilds each face-face-
    connected channel of at least ``min_channel_elements`` flagged
    members by sampling its diameter spine through the channel's
    face-face graph at ``h_local_median`` spacing and re-triangulating
    the patch (rim ∪ existing-interior ∪ new-spine) via Delaunay
    pruned by the rim polygon. ``mode='none'`` is a no-op.

    Centroid insertion (``widen``) shrinks each new sub-triangle's
    median edge length to ~0.577 × the parent's, while the geometric
    channel width (set by the bank polylines) is unchanged. So the
    post-widen w/h ratio is ~1.73 × the pre-widen ratio:
    borderline-flagged elements (ratio just below ``min_w_h``) cross
    the threshold, but very narrow channels (ratio well below) stay
    flagged. Read ``widen`` as "lift local resolution one step" rather
    than "guarantee 3 cells across every narrow channel"; the latter
    is what ``medial`` does.

    Empirically on the PoC #19 cleaned Tokyo-Bay mesh, ``widen`` cuts
    3,178 → 3,032 flagged (4.6 % reduction). PoC #37 quantifies the
    cost / quality story for ``medial``: with
    ``min_channel_elements=10`` the medial-axis cost is 0.66× the
    centroid-widen cost on the same input — and the surviving
    components have mean ``long_axis_m / h_local_median ≈ 6.5``
    (genuinely ribbon-like channels where centroid widen cannot
    deliver "3 cells across" no matter how many times it is run).

    Detector 6 typically flags thousands of elements on real meshes,
    so ``widen`` adds one interior node and two net elements per
    flagged element — the output mesh size grows proportionally.
    ``delete`` is useful for stripping cosmetic narrow inlets but is
    destructive on meshes where the inlets matter. ``medial`` only
    operates on components above the ``min_channel_elements`` filter,
    so it skips the small isolated clusters that dominate detector
    6's raw count and concentrates the cost where it actually
    improves cross-channel resolution.

    The ``repair_flipped=True`` (default) flag wraps the assembled
    output with the same flipped-element rollback safety net Phase G
    uses (``repair_flipped_elements``). It applies only to ``medial``
    where re-triangulation can in principle invert a sliver near the
    rim; for ``widen`` / ``delete`` no inversion is possible by
    construction.
    """
    if mode == "none":
        return mesh, {"mode": "none", "n_flagged": 0, "skipped": True}

    flag, _metric = under_resolved_channels_flag(
        mesh,
        min_w_h=min_w_h,
        sample_ds_m=sample_ds_m,
        arc_separation_factor=arc_separation_factor,
        opposite_bank_cos_max=opposite_bank_cos_max,
        min_channel_elements=int(min_channel_elements),
    )
    h_local = np.asarray(_metric["h_local_m"], dtype=float)
    n_flagged = int(flag.sum())
    info: dict[str, Any] = {
        "mode": mode,
        "min_w_h": float(min_w_h),
        "sample_ds_m": float(sample_ds_m),
        "min_channel_elements": int(min_channel_elements),
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

    if mode == "medial":
        new_mesh, m_info = _repair_under_resolved_channels_medial(
            mesh, flag, h_local,
            bbox=bbox, tol_deg=tol_deg,
            land_ibtype=land_ibtype,
            open_merge_coast_gap=open_merge_coast_gap,
            repair_flipped=repair_flipped,
            max_repair_passes=int(max_repair_passes),
        )
        info.update(m_info)
        return new_mesh, info

    raise ValueError(f"unknown under_resolved_mode: {mode!r}")


# ---------------------------------------------------------------------------
# Phase E Stage 1: medial-axis potential analysis (no re-meshing)
# ---------------------------------------------------------------------------


_EARTH_R_M_FOR_ANALYSIS: float = 6_371_000.0


def _analysis_xy_metric(nodes: np.ndarray) -> np.ndarray:
    """Convert ``(NP, 2)`` lon/lat to a local equirectangular metric
    frame (x in metres east, y in metres north) for in-plane distance
    work. Cheap; uncertainty in the projection is far below the
    accuracy this analysis cares about.
    """
    if nodes.size == 0:
        return nodes.astype(float).copy()
    lon = nodes[:, 0]
    lat = nodes[:, 1]
    lat0 = float(lat.mean())
    lon0 = float(lon.mean())
    dy = (lat - lat0) * (np.pi / 180.0) * _EARTH_R_M_FOR_ANALYSIS
    dx = (
        (lon - lon0) * (np.pi / 180.0) * _EARTH_R_M_FOR_ANALYSIS
        * float(np.cos(np.deg2rad(lat0)))
    )
    return np.column_stack([dx, dy])


def _component_long_axis_m(xy_m: np.ndarray) -> float:
    """Length of the longest principal axis of ``xy_m``'s point cloud
    in metres (PCA-based). Used as a cheap "channel length" proxy
    for a connected channel component.
    """
    if xy_m.shape[0] < 2:
        return 0.0
    centre = xy_m.mean(axis=0)
    centred = xy_m - centre
    # 2x2 covariance, closed-form principal direction.
    cov = centred.T @ centred / max(xy_m.shape[0] - 1, 1)
    eigvals, eigvecs = np.linalg.eigh(cov)
    # Largest eigvec → principal axis. Project points onto it and
    # take the range as the long-axis extent (a tighter envelope
    # than 2*sqrt(eigval), which is just a Gaussian σ).
    principal = eigvecs[:, -1]
    proj = centred @ principal
    return float(proj.max() - proj.min())


def analyze_under_resolved_channels(
    mesh: Fort14Mesh,
    *,
    min_w_h: float = DEFAULT_MIN_W_H,
    sample_ds_m: float = 50.0,
    arc_separation_factor: float = 4.0,
    opposite_bank_cos_max: float = -0.8,
    min_channel_elements: int = 1,
    target_cells_across: int = 3,
) -> dict[str, Any]:
    """Stage 1 of the "true medial-axis Phase E" project.

    Runs detector 6 (``under_resolved_channels_flag``), splits the
    flagged elements into face-face-adjacent connected components
    ("channels"), and reports per-channel statistics that quantify
    the gap between the existing centroid-widen Phase E and the
    medial-axis-based Stage 2 not yet implemented.

    Per-channel record:

    * ``n_elements`` — flagged elements in the component.
    * ``n_nodes`` — distinct nodes touched.
    * ``h_local_median_m`` — median cell-edge length in metres.
    * ``long_axis_m`` — PCA principal-axis extent of the component's
      node cloud (a "channel length" proxy).
    * ``current_phase_e_new_nodes`` — node count if the existing
      centroid-widen Phase E ran on this component (= ``n_elements``).
    * ``medial_axis_new_nodes_estimate`` — node count required to
      lift the channel to ``target_cells_across`` cells across,
      assuming the medial axis sits along ``long_axis`` and the new
      nodes are placed at ``h_local_median_m`` spacing along
      ``(target_cells_across - 1)`` interior rows. **Estimate only**;
      the real number depends on the local channel-width profile
      and is computed in Stage 2.

    ``min_channel_elements`` (default 1 = no filter) is forwarded to
    the underlying detector; small components are dropped before the
    per-component tally is computed. PoC #37 sweeps this to
    characterise where the medial-axis estimate breaks even with the
    centroid-widen baseline.

    The aggregate dict includes ``n_components``,
    ``total_flagged_elements``, the existing-Phase-E vs medial-axis
    new-node counts summed across channels, and the implied
    "extra-nodes-vs-current" delta — the *upper bound* on how many
    nodes Stage 2 would have to insert beyond what Phase E
    centroid-widen already adds.
    """
    if target_cells_across < 2:
        raise ValueError(
            f"target_cells_across must be >= 2, got {target_cells_across}"
        )

    flag, metric_info = under_resolved_channels_flag(
        mesh,
        min_w_h=float(min_w_h),
        sample_ds_m=float(sample_ds_m),
        arc_separation_factor=float(arc_separation_factor),
        opposite_bank_cos_max=float(opposite_bank_cos_max),
        min_channel_elements=int(min_channel_elements),
    )
    n_flagged = int(flag.sum())
    info: dict[str, Any] = {
        "min_w_h": float(min_w_h),
        "min_channel_elements": int(min_channel_elements),
        "target_cells_across": int(target_cells_across),
        "total_flagged_elements": n_flagged,
        "n_components": 0,
        "components": [],
        "current_phase_e_new_nodes": 0,
        "medial_axis_new_nodes_estimate": 0,
        "delta_nodes_vs_current": 0,
    }
    if n_flagged == 0:
        return info

    # Subgraph adjacency restricted to flagged elements.
    flagged_idx = np.where(flag)[0]
    full_adj = face_face_adjacency(mesh.elements)
    sub = full_adj[flagged_idx][:, flagged_idx]
    n_comp, labels = connected_components(sub, directed=False, return_labels=True)
    info["n_components"] = int(n_comp)

    nodes_m = _analysis_xy_metric(np.asarray(mesh.nodes, dtype=float))
    h_local = np.asarray(metric_info["h_local_m"], dtype=float)
    interior_rows = int(target_cells_across) - 1

    cur_total = 0
    medial_total = 0
    for c in range(int(n_comp)):
        member_local = np.where(labels == c)[0]
        member_global = flagged_idx[member_local]
        if member_global.size == 0:
            continue
        elem_block = mesh.elements[member_global]
        node_ids = np.unique(elem_block.ravel())
        comp_xy = nodes_m[node_ids]
        long_axis_m = _component_long_axis_m(comp_xy)
        comp_h = h_local[member_global]
        comp_h_pos = comp_h[comp_h > 0]
        h_med = float(np.median(comp_h_pos)) if comp_h_pos.size else 0.0

        cur_phase_e = int(member_global.size)  # one centroid each
        if h_med > 0 and long_axis_m > 0:
            nodes_per_row = max(int(np.ceil(long_axis_m / h_med)), 1)
            medial_estimate = int(interior_rows * nodes_per_row)
        else:
            medial_estimate = 0

        info["components"].append({
            "component_id": int(c),
            "n_elements": int(member_global.size),
            "n_nodes": int(node_ids.size),
            "h_local_median_m": h_med,
            "long_axis_m": float(long_axis_m),
            "current_phase_e_new_nodes": cur_phase_e,
            "medial_axis_new_nodes_estimate": medial_estimate,
        })
        cur_total += cur_phase_e
        medial_total += medial_estimate

    info["current_phase_e_new_nodes"] = int(cur_total)
    info["medial_axis_new_nodes_estimate"] = int(medial_total)
    info["delta_nodes_vs_current"] = int(medial_total - cur_total)
    return info


# ---------------------------------------------------------------------------
# Phase E Stage 2: medial-axis CDT re-meshing
# ---------------------------------------------------------------------------


def _patch_rim_polygon(
    elements: np.ndarray, comp_elem_idx: np.ndarray,
) -> np.ndarray | None:
    """Walk the boundary of the face-face-component formed by
    ``comp_elem_idx`` and return its rim node IDs in CCW order.

    Returns ``None`` if the patch is not simply connected (multiple
    boundary loops) or if a node has odd degree on the boundary
    (degenerate topology). Identifies rim edges as those used by
    exactly one element in the patch.
    """
    if comp_elem_idx.size == 0:
        return None
    patch = elements[comp_elem_idx]
    edges = np.vstack([
        patch[:, [0, 1]],
        patch[:, [1, 2]],
        patch[:, [2, 0]],
    ])
    edges_sorted = np.sort(edges, axis=1)
    # An edge appears once if rim, twice if interior.
    edge_keys = (
        edges_sorted[:, 0].astype(np.int64) << 32
    ) | edges_sorted[:, 1].astype(np.int64)
    uniq, counts = np.unique(edge_keys, return_counts=True)
    rim_keys = uniq[counts == 1]
    if rim_keys.size == 0:
        return None
    rim_a = (rim_keys >> 32).astype(np.int64)
    rim_b = (rim_keys & 0xFFFFFFFF).astype(np.int64)
    rim_pairs = np.column_stack([rim_a, rim_b])

    # Build node → adjacent rim-node list (each rim node has degree 2
    # for a simply connected patch).
    adj: dict[int, list[int]] = {}
    for a, b in rim_pairs:
        adj.setdefault(int(a), []).append(int(b))
        adj.setdefault(int(b), []).append(int(a))
    if any(len(v) != 2 for v in adj.values()):
        # Pinch points or branching rim — not simply connected.
        return None

    # Walk the loop. Pick the lex-smallest node as start, and pick the
    # next node so that the walk is uniquely determined.
    start = min(adj.keys())
    walk = [start]
    prev = -1
    cur = start
    while True:
        nbrs = adj[cur]
        nxt = nbrs[0] if nbrs[0] != prev else nbrs[1]
        if nxt == start:
            break
        walk.append(nxt)
        prev = cur
        cur = nxt
        if len(walk) > rim_pairs.shape[0]:
            return None  # safety: walk longer than rim → not simple
    if len(walk) != rim_pairs.shape[0]:
        return None  # multi-loop rim (more than one connected boundary)
    return np.asarray(walk, dtype=np.int64)


def _ccw_orient(rim_xy: np.ndarray, rim_ids: np.ndarray) -> np.ndarray:
    """Reverse ``rim_ids`` if the polygon ``rim_xy`` is CW."""
    x = rim_xy[:, 0]
    y = rim_xy[:, 1]
    signed_area = 0.5 * float(np.sum(
        x * np.roll(y, -1) - np.roll(x, -1) * y,
    ))
    if signed_area < 0:
        return rim_ids[::-1].copy()
    return rim_ids


def _diameter_path(
    comp_elem_idx: np.ndarray, full_adj: Any,
) -> np.ndarray:
    """Two-BFS diameter path through the face-face subgraph restricted
    to ``comp_elem_idx``. Returns global element indices in path order.
    """
    if comp_elem_idx.size == 1:
        return comp_elem_idx.copy()
    # Build a local adjacency dict.
    member_set = set(int(i) for i in comp_elem_idx)
    sub_adj: dict[int, list[int]] = {i: [] for i in member_set}
    sub = full_adj[comp_elem_idx][:, comp_elem_idx]
    coo = sub.tocoo()
    for r, c in zip(coo.row, coo.col):
        if r >= c:
            continue
        a = int(comp_elem_idx[r])
        b = int(comp_elem_idx[c])
        sub_adj[a].append(b)
        sub_adj[b].append(a)

    def _bfs_farthest(start: int) -> tuple[int, dict[int, int]]:
        visited = {start: -1}
        queue = [start]
        order = []
        while queue:
            nxt: list[int] = []
            for node in queue:
                order.append(node)
                for nb in sub_adj[node]:
                    if nb not in visited:
                        visited[nb] = node
                        nxt.append(nb)
            queue = nxt
        return order[-1], visited

    s = int(comp_elem_idx[0])
    far1, _ = _bfs_farthest(s)
    far2, parents = _bfs_farthest(far1)
    path: list[int] = []
    node = far2
    while node != -1:
        path.append(node)
        node = parents[node]
    path.reverse()
    return np.asarray(path, dtype=np.int64)


def _resample_polyline(
    xy: np.ndarray, depths: np.ndarray, target_spacing_m: float,
    *, lat_centre_deg: float,
) -> tuple[np.ndarray, np.ndarray]:
    """Place samples uniformly along the polyline ``xy`` at
    approximately ``target_spacing_m`` spacing in metres. ``depths``
    travels with each polyline vertex and is linearly interpolated to
    each sample. Endpoints are dropped (a small offset of one spacing
    is left at each end so the spine does not collide with the rim).
    """
    if xy.shape[0] < 2:
        return np.empty((0, 2)), np.empty((0,))
    deg_per_m_lat = 1.0 / (EARTH_R_M * np.pi / 180.0)
    deg_per_m_lon = deg_per_m_lat / max(np.cos(np.deg2rad(lat_centre_deg)), 1e-6)
    seg_dx_m = (xy[1:, 0] - xy[:-1, 0]) / deg_per_m_lon
    seg_dy_m = (xy[1:, 1] - xy[:-1, 1]) / deg_per_m_lat
    seg_len_m = np.hypot(seg_dx_m, seg_dy_m)
    arc_m = np.concatenate([[0.0], np.cumsum(seg_len_m)])
    total_m = float(arc_m[-1])
    if total_m < 2.0 * target_spacing_m:
        return np.empty((0, 2)), np.empty((0,))
    # Skip the first and last spacing (endpoints near rim).
    s_first = target_spacing_m
    s_last = total_m - target_spacing_m
    if s_last <= s_first:
        return np.empty((0, 2)), np.empty((0,))
    n_samples = int(np.floor((s_last - s_first) / target_spacing_m)) + 1
    if n_samples < 1:
        return np.empty((0, 2)), np.empty((0,))
    sample_s = np.linspace(s_first, s_last, n_samples)
    out_xy = np.empty((n_samples, 2), dtype=float)
    out_depth = np.empty(n_samples, dtype=float)
    for k, s in enumerate(sample_s):
        # Locate the segment.
        j = int(np.searchsorted(arc_m, s, side="right")) - 1
        j = min(max(j, 0), xy.shape[0] - 2)
        seg_start = arc_m[j]
        seg_end = arc_m[j + 1]
        t = (s - seg_start) / max(seg_end - seg_start, 1e-12)
        out_xy[k] = xy[j] + t * (xy[j + 1] - xy[j])
        out_depth[k] = depths[j] + t * (depths[j + 1] - depths[j])
    return out_xy, out_depth


def _retriangulate_patch(
    rim_xy: np.ndarray, spine_xy: np.ndarray, n_rim: int,
) -> tuple[np.ndarray | None, str]:
    """Delaunay-triangulate ``(rim ∪ spine)`` and prune triangles whose
    centroid falls outside the rim polygon. Returns ``(triangles,
    reason)``. ``triangles`` is the (T, 3) index array into the
    combined ``[rim_xy; spine_xy]`` ordering, or ``None`` if the patch
    is rejected; ``reason`` is the rejection bucket label or "ok".
    Triangles are oriented CCW (positive signed area).
    """
    from scipy.spatial import Delaunay
    from shapely import contains_xy
    from shapely.geometry import Polygon

    points = np.vstack([rim_xy, spine_xy])
    if points.shape[0] < 3:
        return None, "delaunay_failed"
    try:
        tri = Delaunay(points)
    except Exception:  # noqa: BLE001 — scipy raises QhullError variants
        return None, "delaunay_failed"
    triangles = tri.simplices
    if triangles.size == 0:
        return None, "delaunay_failed"

    poly = Polygon(rim_xy)
    if not poly.is_valid:
        return None, "rim_polygon_invalid"

    centroids = points[triangles].mean(axis=1)
    inside = contains_xy(poly, centroids[:, 0], centroids[:, 1])
    triangles = triangles[inside]
    if triangles.size == 0:
        return None, "all_triangles_pruned"

    # Re-orient any retained triangle so its signed area is positive
    # (Delaunay output may be CW after the prune mask reordering on
    # certain numpy versions; cheap to enforce explicitly).
    p0 = points[triangles[:, 0]]
    p1 = points[triangles[:, 1]]
    p2 = points[triangles[:, 2]]
    cross = (p1[:, 0] - p0[:, 0]) * (p2[:, 1] - p0[:, 1]) - (
        p1[:, 1] - p0[:, 1]
    ) * (p2[:, 0] - p0[:, 0])
    cw = cross < 0
    if cw.any():
        triangles = triangles.copy()
        triangles[cw] = triangles[cw][:, [0, 2, 1]]
        # Recompute and abort if any triangle is degenerate.
        cross = (
            (points[triangles[:, 1], 0] - points[triangles[:, 0], 0])
            * (points[triangles[:, 2], 1] - points[triangles[:, 0], 1])
            - (points[triangles[:, 1], 1] - points[triangles[:, 0], 1])
            * (points[triangles[:, 2], 0] - points[triangles[:, 0], 0])
        )
    if (cross <= 0).any():
        return None, "degenerate_triangle"

    # Verify rim edges all appear.
    edges = set()
    for tri_ in triangles:
        a, b, c = int(tri_[0]), int(tri_[1]), int(tri_[2])
        edges.add((min(a, b), max(a, b)))
        edges.add((min(b, c), max(b, c)))
        edges.add((min(c, a), max(c, a)))
    rim_edges = set()
    for i in range(n_rim):
        a, b = i, (i + 1) % n_rim
        rim_edges.add((min(a, b), max(a, b)))
    if not rim_edges.issubset(edges):
        return None, "non_convex_rim_unsolved"
    return triangles, "ok"


def _repair_under_resolved_channels_medial(
    mesh: Fort14Mesh,
    flag: np.ndarray,
    h_local: np.ndarray,
    *,
    bbox: tuple[float, float, float, float] | None,
    tol_deg: float | None,
    land_ibtype: int,
    open_merge_coast_gap: int,
    repair_flipped: bool,
    max_repair_passes: int,
) -> tuple[Fort14Mesh, dict[str, Any]]:
    """Stage 2 driver: medial-axis CDT re-meshing of every face-face
    component of ``flag``. Components for which CDT fails are left in
    place; their original triangulation survives. The CDT routine
    returns CCW-oriented triangles only — components whose retriangu-
    lation would yield a degenerate or inverted triangle are rejected
    rather than committed, so the final mesh is flip-free by
    construction. ``repair_flipped`` / ``max_repair_passes`` are
    accepted for symmetry with Phase G but are no-ops here.
    """
    del repair_flipped, max_repair_passes  # accepted for symmetry, unused.
    flag = np.asarray(flag, dtype=bool)
    flagged_idx = np.where(flag)[0]
    info: dict[str, Any] = {
        "n_components": 0,
        "n_components_replaced": 0,
        "n_components_skipped": 0,
        "n_nodes_inserted": 0,
        "n_elements_removed": 0,
        "n_elements_inserted": 0,
        "skip_reasons": {
            "rim_walk_failed": 0,
            "spine_too_short": 0,
            "delaunay_failed": 0,
            "non_convex_rim_unsolved": 0,
            "rim_polygon_invalid": 0,
            "all_triangles_pruned": 0,
            "degenerate_triangle": 0,
        },
    }
    if flagged_idx.size == 0:
        info["skipped"] = True
        return mesh, info

    full_adj = face_face_adjacency(mesh.elements)
    sub = full_adj[flagged_idx][:, flagged_idx]
    n_comp, labels = connected_components(sub, directed=False, return_labels=True)
    info["n_components"] = int(n_comp)

    lat_centre = float(np.mean(mesh.nodes[:, 1]))
    deg_per_m_lat = 1.0 / (EARTH_R_M * np.pi / 180.0)
    deg_per_m_lon = deg_per_m_lat / max(np.cos(np.deg2rad(lat_centre)), 1e-6)

    centroids_all = mesh.nodes[mesh.elements].mean(axis=1)
    centroid_depths_all = mesh.depths[mesh.elements].mean(axis=1)

    # Per-component records used in the assembly pass.
    replace_records: list[dict[str, Any]] = []

    for c in range(int(n_comp)):
        comp_local = np.where(labels == c)[0]
        comp_elems = flagged_idx[comp_local]

        rim_node_ids = _patch_rim_polygon(mesh.elements, comp_elems)
        if rim_node_ids is None:
            info["skip_reasons"]["rim_walk_failed"] += 1
            info["n_components_skipped"] += 1
            continue
        rim_xy = mesh.nodes[rim_node_ids]
        rim_node_ids = _ccw_orient(rim_xy, rim_node_ids)
        rim_xy = mesh.nodes[rim_node_ids]
        n_rim = rim_node_ids.size

        # Spine through the component's centroids.
        diam_elems = _diameter_path(comp_elems, full_adj)
        spine_xy = centroids_all[diam_elems]
        spine_depths = centroid_depths_all[diam_elems]
        h_local_comp = h_local[comp_elems]
        positive_h = h_local_comp[h_local_comp > 0]
        h_med = float(np.median(positive_h)) if positive_h.size else 0.0
        if h_med <= 0:
            info["skip_reasons"]["spine_too_short"] += 1
            info["n_components_skipped"] += 1
            continue
        # Need spine longer than 2*h_med to have any interior sample.
        new_spine_xy, new_spine_depths = _resample_polyline(
            spine_xy, spine_depths, h_med, lat_centre_deg=lat_centre,
        )
        if new_spine_xy.shape[0] == 0:
            info["skip_reasons"]["spine_too_short"] += 1
            info["n_components_skipped"] += 1
            continue

        # Build the (rim + spine) point set in metres-equivalent coords.
        rim_xy_m = np.column_stack([
            (rim_xy[:, 0] - rim_xy[0, 0]) / deg_per_m_lon,
            (rim_xy[:, 1] - rim_xy[0, 1]) / deg_per_m_lat,
        ])
        spine_xy_m = np.column_stack([
            (new_spine_xy[:, 0] - rim_xy[0, 0]) / deg_per_m_lon,
            (new_spine_xy[:, 1] - rim_xy[0, 1]) / deg_per_m_lat,
        ])
        triangles, reason = _retriangulate_patch(rim_xy_m, spine_xy_m, n_rim)
        if triangles is None:
            info["skip_reasons"][reason] = (
                info["skip_reasons"].get(reason, 0) + 1
            )
            info["n_components_skipped"] += 1
            continue

        replace_records.append({
            "comp_elems": comp_elems,
            "rim_node_ids": rim_node_ids,
            "new_spine_xy": new_spine_xy,
            "new_spine_depths": new_spine_depths,
            "local_triangles": triangles,  # indices into [0..n_rim) ∪
                                            # [n_rim..n_rim+S)
            "n_rim": int(n_rim),
        })
        info["n_components_replaced"] += 1

    if not replace_records:
        info["skipped_all"] = True
        return mesh, info

    # Assembly: drop replaced elements, append new spine nodes, append
    # new triangles using global node IDs.
    elements_to_drop = np.concatenate(
        [r["comp_elems"] for r in replace_records]
    )
    keep_elem_mask = np.ones(mesh.n_elements, dtype=bool)
    keep_elem_mask[elements_to_drop] = False
    info["n_elements_removed"] = int(elements_to_drop.size)

    new_nodes_chunks: list[np.ndarray] = []
    new_depths_chunks: list[np.ndarray] = []
    new_elements_chunks: list[np.ndarray] = []
    new_node_offset = mesh.n_nodes
    for rec in replace_records:
        n_rim = rec["n_rim"]
        spine_xy = rec["new_spine_xy"]
        spine_depths = rec["new_spine_depths"]
        new_nodes_chunks.append(spine_xy)
        new_depths_chunks.append(spine_depths)
        # Map local indices: 0..n_rim-1 → rim_node_ids[i],
        # n_rim..n_rim+S-1 → new_node_offset + (k - n_rim).
        loc2global = np.empty(n_rim + spine_xy.shape[0], dtype=np.int64)
        loc2global[:n_rim] = rec["rim_node_ids"]
        loc2global[n_rim:] = (
            new_node_offset + np.arange(spine_xy.shape[0], dtype=np.int64)
        )
        new_elements_chunks.append(loc2global[rec["local_triangles"]])
        new_node_offset += spine_xy.shape[0]

    n_new_nodes = sum(c.shape[0] for c in new_nodes_chunks)
    info["n_nodes_inserted"] = int(n_new_nodes)
    info["n_elements_inserted"] = int(
        sum(c.shape[0] for c in new_elements_chunks)
    )

    new_nodes = np.vstack([mesh.nodes, *new_nodes_chunks])
    new_depths = np.concatenate([mesh.depths, *new_depths_chunks])
    surviving_elements = mesh.elements[keep_elem_mask]
    inserted_elements = np.vstack(new_elements_chunks).astype(
        mesh.elements.dtype,
    )
    new_elements = np.vstack([surviving_elements, inserted_elements])

    # Rim nodes' IDs are unchanged, so previously-stored boundary
    # segments remain pointwise valid. If ``bbox`` is supplied,
    # re-derive open / land segments from scratch (the cleanest
    # behaviour when the patch touched a boundary, since the new
    # triangulation may have fewer / more boundary edges per node).
    if bbox is not None and tol_deg is not None:
        new_mesh = Fort14Mesh(
            title=mesh.title,
            nodes=new_nodes,
            depths=new_depths,
            elements=new_elements,
            open_boundaries=[],
            land_boundaries=[],
        )
        new_mesh = rebuild_boundaries(
            new_mesh, bbox=bbox, tol_deg=tol_deg,
            land_ibtype=land_ibtype,
            open_merge_coast_gap=open_merge_coast_gap,
        )
    else:
        new_mesh = Fort14Mesh(
            title=mesh.title,
            nodes=new_nodes,
            depths=new_depths,
            elements=new_elements,
            open_boundaries=[
                np.asarray(s).copy() for s in mesh.open_boundaries
            ],
            land_boundaries=[
                (int(ib), np.asarray(s).copy())
                for ib, s in mesh.land_boundaries
            ],
        )

    return new_mesh, info


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

# Cap on the iterative rollback loop in :func:`_repair_flipped_elements`.
# Five passes is enough in practice (a flip propagates ≤ K rings away
# after K passes); the full-rollback fallback handles the rare case
# where it doesn't converge.
DEFAULT_SMOOTH_REPAIR_PASSES: int = 5


def _signed_areas_raw(nodes: np.ndarray, elements: np.ndarray) -> np.ndarray:
    """Per-element signed area for raw ``(NP, 2)`` / ``(NE, 3)`` arrays.

    Same convention as
    :func:`fvcom_mesh_tools.algorithms.perpendicularity.signed_areas`
    but takes plain numpy arrays so it can be used inside helpers
    that have not yet wrapped a result into a :class:`Fort14Mesh`.
    """
    p0 = nodes[elements[:, 0]]
    p1 = nodes[elements[:, 1]]
    p2 = nodes[elements[:, 2]]
    return 0.5 * (
        (p1[:, 0] - p0[:, 0]) * (p2[:, 1] - p0[:, 1])
        - (p1[:, 1] - p0[:, 1]) * (p2[:, 0] - p0[:, 0])
    )


def _repair_flipped_elements(
    pre_nodes: np.ndarray,
    post_nodes: np.ndarray,
    elements: np.ndarray,
    *,
    max_passes: int = DEFAULT_SMOOTH_REPAIR_PASSES,
) -> tuple[np.ndarray, dict[str, Any]]:
    """Roll back nodes belonging to negative-signed-area triangles.

    A Laplacian smoothing pass converges on edge-length stability but
    does not check signed area, so a topologically-valid input can
    converge to a state with locally-flipped triangles. This helper
    iteratively reverts the three nodes of every flipped triangle to
    their pre-smoothing positions until no flips remain or
    ``max_passes`` is reached. If flips persist, the function falls
    back to a full rollback (``post_nodes := pre_nodes``) so the
    caller is guaranteed a flip-free output.

    Assumes ``pre_nodes`` is itself flip-free; if the input is already
    pathological the full-rollback safety net does not help.

    Returns ``(repaired_nodes, info)``. ``info`` keys:

    - ``n_flipped_post_smooth``: signed-area-negative count in
      ``post_nodes`` before any rollback.
    - ``n_flipped_after_repair``: ditto after the rollback (always 0
      on success).
    - ``n_nodes_rolled_back``: count of nodes that ended up reverted
      to their pre-smoothing position.
    - ``n_rollback_passes``: number of rollback passes actually run
      (0 if no flips were detected, ``max_passes`` if the fallback
      fired).
    - ``full_rollback``: True iff the safety net wiped every node
      back to its pre-smoothing position.
    """
    repaired = np.asarray(post_nodes, dtype=float).copy()
    pre = np.asarray(pre_nodes, dtype=float)
    n_flipped_initial = int((_signed_areas_raw(repaired, elements) < 0).sum())
    rolled_back = np.zeros(pre.shape[0], dtype=bool)

    if n_flipped_initial == 0:
        return repaired, {
            "n_flipped_post_smooth": 0,
            "n_flipped_after_repair": 0,
            "n_nodes_rolled_back": 0,
            "n_rollback_passes": 0,
            "full_rollback": False,
        }

    passes_run = 0
    for pass_idx in range(max_passes):
        passes_run = pass_idx + 1
        sa = _signed_areas_raw(repaired, elements)
        flipped_mask = sa < 0
        if not flipped_mask.any():
            return repaired, {
                "n_flipped_post_smooth": n_flipped_initial,
                "n_flipped_after_repair": 0,
                "n_nodes_rolled_back": int(rolled_back.sum()),
                "n_rollback_passes": pass_idx,
                "full_rollback": False,
            }
        bad_nodes = np.unique(elements[flipped_mask].ravel())
        repaired[bad_nodes] = pre[bad_nodes]
        rolled_back[bad_nodes] = True

    # Final check after the loop.
    sa = _signed_areas_raw(repaired, elements)
    if not (sa < 0).any():
        return repaired, {
            "n_flipped_post_smooth": n_flipped_initial,
            "n_flipped_after_repair": 0,
            "n_nodes_rolled_back": int(rolled_back.sum()),
            "n_rollback_passes": passes_run,
            "full_rollback": False,
        }

    # Safety net: full rollback. Guaranteed flip-free as long as the
    # caller's pre_nodes was itself flip-free.
    repaired = pre.copy()
    rolled_back[:] = True
    return repaired, {
        "n_flipped_post_smooth": n_flipped_initial,
        "n_flipped_after_repair": 0,
        "n_nodes_rolled_back": int(rolled_back.sum()),
        "n_rollback_passes": max_passes,
        "full_rollback": True,
    }


# Public alias for callers outside this module who need the same
# flip-rollback algorithm (e.g. ``mesh_engine/oceanmesh.py`` when it
# wraps ``om.laplacian2`` in the build-time cleanup chain).
repair_flipped_elements = _repair_flipped_elements


def smooth_mesh_laplacian(
    mesh: Fort14Mesh,
    *,
    max_iter: int = DEFAULT_SMOOTH_LAPLACIAN_ITERS,
    tol: float = DEFAULT_SMOOTH_LAPLACIAN_TOL,
    pfix: np.ndarray | None = None,
    repair_flipped: bool = True,
    max_repair_passes: int = DEFAULT_SMOOTH_REPAIR_PASSES,
) -> tuple[Fort14Mesh, dict[str, Any]]:
    """Phase G: move every non-boundary node to the mean position of its
    connected neighbours (Laplacian smoothing).

    Wraps :func:`oceanmesh.laplacian2`. oceanmesh derives the boundary
    set from the mesh topology and pins it automatically, so the
    output's open / land boundary node coordinates and the boundary
    lists themselves are unchanged. Connectivity and depth indices
    are preserved (``laplacian2`` only moves vertices). ``pfix`` adds
    extra fixed coordinates beyond the topological boundary.

    ``oceanmesh.laplacian2`` converges on edge-length stability but
    does not check signed area, so a topologically-valid input can
    converge to a locally-flipped state. With ``repair_flipped=True``
    (default) any flipped triangles are repaired by rolling back
    their three nodes' positions; see :func:`_repair_flipped_elements`
    for the algorithm and its termination guarantees. The repair
    counts surface in the returned ``info`` dict
    (``n_flipped_post_smooth``, ``n_nodes_rolled_back``,
    ``full_rollback``). Set ``repair_flipped=False`` to surface the
    raw oceanmesh output (useful for diagnosing whether a particular
    mesh causes flipping).

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
    if max_repair_passes < 0:
        raise ValueError(f"max_repair_passes must be >= 0, got {max_repair_passes}")
    if mesh.n_elements == 0:
        return mesh, {
            "max_iter": int(max_iter),
            "tol": float(tol),
            "n_pfix": 0,
            "n_nodes_moved": 0,
            "displacement_max": 0.0,
            "displacement_mean": 0.0,
            "n_flipped_post_smooth": 0,
            "n_flipped_after_repair": 0,
            "n_nodes_rolled_back": 0,
            "n_rollback_passes": 0,
            "full_rollback": False,
            "skipped": True,
        }

    from oceanmesh import laplacian2  # GPL-3.0-or-later, lazy import

    pre = np.asarray(mesh.nodes, dtype=float).copy()
    elements = np.asarray(mesh.elements, dtype=int)
    new_vertices, _ = laplacian2(
        pre.copy(), elements,
        max_iter=int(max_iter), tol=float(tol), pfix=pfix,
    )
    new_vertices = np.asarray(new_vertices, dtype=float)

    if repair_flipped:
        new_vertices, repair_info = _repair_flipped_elements(
            pre, new_vertices, elements, max_passes=max_repair_passes,
        )
    else:
        sa = _signed_areas_raw(new_vertices, elements)
        n_flipped_raw = int((sa < 0).sum())
        repair_info = {
            "n_flipped_post_smooth": n_flipped_raw,
            "n_flipped_after_repair": n_flipped_raw,
            "n_nodes_rolled_back": 0,
            "n_rollback_passes": 0,
            "full_rollback": False,
        }

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
        **repair_info,
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
    under_resolved_min_channel_elements: int = 1,
    repair_skewed: bool = False,
    repair_skewed_min_angle_deg: float = DEFAULT_SKEWED_MIN_ANGLE_DEG,
    repair_skewed_max_angle_deg: float = DEFAULT_SKEWED_MAX_ANGLE_DEG,
    smooth_laplacian: bool = False,
    smooth_laplacian_iters: int = DEFAULT_SMOOTH_LAPLACIAN_ITERS,
    smooth_laplacian_tol: float = DEFAULT_SMOOTH_LAPLACIAN_TOL,
    smooth_repair_flipped: bool = True,
    smooth_max_repair_passes: int = DEFAULT_SMOOTH_REPAIR_PASSES,
    phase_h: bool = False,
    phase_h_alpha_target: float = 0.95,
    phase_h_min_angle_target: float = 20.0,
    phase_h_max_outer_rounds: int = 10,
    phase_h_max_topology_per_round: int = 10_000,
    phase_h_max_smooth_sweeps: int = 200,
    phase_h_coastline_paths: list | None = None,
    phase_h_max_snap_distance_m: float = 500.0,
    phase_h_lookahead: bool = False,
    phase_h_max_lookahead_per_round: int = 10_000,
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
    under_resolved_min_channel_elements:
        Forwarded to :func:`under_resolved_channels_flag`. Default 1
        (no filter); raise to N to ignore detector-6 flags whose
        face-face-connected component has fewer than N flagged
        elements. PoC #35 motivated this filter — most flagged
        clusters on real meshes are tiny (mean ~3 elements / channel)
        and not the long ribbon-like inlets Phase E targets.
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
    smooth_repair_flipped:
        Phase G safety net (default True). After
        ``oceanmesh.laplacian2`` returns, detect any negative
        signed-area triangles and roll the affected nodes back to
        their pre-smoothing positions; iterate until no flips
        remain. Set False to surface the raw oceanmesh result.
    smooth_max_repair_passes:
        Cap on the iterative-rollback loop. Default
        ``DEFAULT_SMOOTH_REPAIR_PASSES = 5``; full rollback to the
        pre-smoothing positions fires if convergence is not reached
        within this many passes.
    """
    if thin_chain_mode not in ("widen", "delete", "none"):
        raise ValueError(
            f"thin_chain_mode must be one of widen / delete / none, got {thin_chain_mode!r}"
        )
    if under_resolved_mode not in ("widen", "delete", "medial", "none"):
        raise ValueError(
            "under_resolved_mode must be one of widen / delete / medial / "
            f"none, got {under_resolved_mode!r}"
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
            "under_resolved_min_channel_elements":
                int(under_resolved_min_channel_elements),
            "repair_skewed": bool(repair_skewed),
            "repair_skewed_min_angle_deg": float(repair_skewed_min_angle_deg),
            "repair_skewed_max_angle_deg": float(repair_skewed_max_angle_deg),
            "smooth_laplacian": bool(smooth_laplacian),
            "smooth_laplacian_iters": int(smooth_laplacian_iters),
            "smooth_laplacian_tol": float(smooth_laplacian_tol),
            "smooth_repair_flipped": bool(smooth_repair_flipped),
            "smooth_max_repair_passes": int(smooth_max_repair_passes),
            "phase_h": bool(phase_h),
            "phase_h_alpha_target": float(phase_h_alpha_target),
            "phase_h_min_angle_target": float(phase_h_min_angle_target),
            "phase_h_max_outer_rounds": int(phase_h_max_outer_rounds),
            "phase_h_max_topology_per_round":
                int(phase_h_max_topology_per_round),
            "phase_h_max_smooth_sweeps": int(phase_h_max_smooth_sweeps),
            "phase_h_coastline_paths":
                [str(p) for p in (phase_h_coastline_paths or [])],
            "phase_h_max_snap_distance_m":
                float(phase_h_max_snap_distance_m),
            "phase_h_lookahead": bool(phase_h_lookahead),
            "phase_h_max_lookahead_per_round":
                int(phase_h_max_lookahead_per_round),
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
            min_channel_elements=under_resolved_min_channel_elements,
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
            repair_flipped=smooth_repair_flipped,
            max_repair_passes=smooth_max_repair_passes,
        )
        info["phases"].append({"name": "smooth_mesh_laplacian", **g_info})

    if phase_h:
        # Lazy import to avoid a module-level cycle (mesh_clean_phase_h
        # imports a couple of helpers from mesh_clean).
        from fvcom_mesh_tools.mesh_clean_phase_h import (
            build_coastline_projector,
            phase_h_optimize,
        )
        projector = None
        if phase_h_coastline_paths:
            projector = build_coastline_projector(
                phase_h_coastline_paths,
                max_snap_distance_m=float(phase_h_max_snap_distance_m),
                mean_latitude_deg=(
                    float(cur.nodes[:, 1].mean()) if cur.n_nodes else None
                ),
            )
        cur, h_info = phase_h_optimize(
            cur,
            alpha_target=float(phase_h_alpha_target),
            min_angle_target=float(phase_h_min_angle_target),
            max_smooth_sweeps=int(phase_h_max_smooth_sweeps),
            max_topology_per_round=int(phase_h_max_topology_per_round),
            max_outer_rounds=int(phase_h_max_outer_rounds),
            coastline_projector=projector,
            lookahead_enabled=bool(phase_h_lookahead),
            max_lookahead_per_round=int(phase_h_max_lookahead_per_round),
        )
        info["phases"].append({"name": "phase_h_optimize", **h_info})

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
    "DEFAULT_SMOOTH_REPAIR_PASSES",
    "ThinChainMode",
    "UnderResolvedMode",
    "analyze_under_resolved_channels",
    "clean_mesh",
    "keep_components",
    "rebuild_boundaries",
    "remove_elements",
    "repair_flipped_elements",
    "repair_overconnected_nodes",
    "repair_skewed_elements",
    "repair_thin_chains",
    "repair_under_resolved_channels",
    "smooth_mesh_laplacian",
    "trim_dead_ends",
    "widen_thin_elements_at_centroid",
]
