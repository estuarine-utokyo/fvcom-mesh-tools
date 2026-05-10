"""Mesh-clean operations: prune disjoint pools, trim dead-end elements,
and repair 1-cell-wide channels.

Three phases, applied in order:

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

After every Phase A / B / C-delete deletion the boundaries are
re-derived via
:func:`fvcom_mesh_tools.algorithms.classify_boundaries_by_bbox`, so the
output mesh has a consistent open / land segmentation in the same
style as ``fmesh-buildmesh``. Phase C-widen preserves boundaries
unchanged because centroid insertion only introduces interior nodes /
edges. Repair of over-connected nodes is **not** implemented here.
"""

from __future__ import annotations

from typing import Any, Literal

import numpy as np
from scipy.sparse.csgraph import connected_components

from fvcom_mesh_tools.algorithms import classify_boundaries_by_bbox
from fvcom_mesh_tools.diagnostics import (
    dead_end_elements_flag,
    face_face_adjacency,
    open_boundary_node_mask,
    thin_chain_elements_flag,
    thin_elements_flag,
)
from fvcom_mesh_tools.io import Fort14Mesh

ThinChainMode = Literal["widen", "delete", "none"]

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
    """
    if thin_chain_mode not in ("widen", "delete", "none"):
        raise ValueError(
            f"thin_chain_mode must be one of widen / delete / none, got {thin_chain_mode!r}"
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
    "ThinChainMode",
    "clean_mesh",
    "keep_components",
    "rebuild_boundaries",
    "remove_elements",
    "repair_thin_chains",
    "trim_dead_ends",
    "widen_thin_elements_at_centroid",
]
