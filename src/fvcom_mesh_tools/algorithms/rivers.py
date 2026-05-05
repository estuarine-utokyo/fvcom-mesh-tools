"""Snap user-supplied river-mouth points to land-boundary nodes and
re-segment the parent boundary so each river occupies its own segment
with a custom ibtype (typically 21 for FVCOM river inflow).

The routine here is intentionally light: it does not move any node,
it does not change the mesh connectivity, and it only touches the
``land_boundaries`` list of the input mesh. Open boundaries and
existing already-river-classified segments are skipped.
"""

from __future__ import annotations

from dataclasses import replace

import numpy as np

from fvcom_mesh_tools.io import Fort14Mesh

EARTH_R_M = 6_371_000.0


def _haversine_m(p_lonlat: np.ndarray, q_lonlat: np.ndarray) -> np.ndarray:
    """Great-circle distances (m) between rows of ``p_lonlat`` and the
    single point ``q_lonlat``. Inputs in degrees."""
    lon1 = np.deg2rad(p_lonlat[:, 0])
    lat1 = np.deg2rad(p_lonlat[:, 1])
    lon2 = np.deg2rad(q_lonlat[0])
    lat2 = np.deg2rad(q_lonlat[1])
    dlat = lat2 - lat1
    dlon = lon2 - lon1
    a = np.sin(dlat / 2) ** 2 + np.cos(lat1) * np.cos(lat2) * np.sin(dlon / 2) ** 2
    return 2 * EARTH_R_M * np.arcsin(np.sqrt(np.clip(a, 0.0, 1.0)))


def _split_segment(
    ids: np.ndarray, center_idx: int, n_river: int
) -> tuple[np.ndarray, np.ndarray, np.ndarray]:
    """Split a polyline ``ids`` so the river occupies ``n_river`` nodes
    centred on ``center_idx``. Returns ``(prefix, river, suffix)`` slices.

    The river slice is clamped to the segment endpoints, so the
    prefix or suffix may be empty.
    """
    n_seg = ids.shape[0]
    half = n_river // 2
    start = center_idx - half
    end = start + n_river
    if start < 0:
        start = 0
        end = min(n_river, n_seg)
    if end > n_seg:
        end = n_seg
        start = max(0, n_seg - n_river)
    return ids[:start].copy(), ids[start:end].copy(), ids[end:].copy()


def add_river_inflow_segments(
    mesh: Fort14Mesh,
    points: np.ndarray,
    *,
    n_nodes_per_river: int = 5,
    river_ibtype: int = 21,
    snap_tol_m: float | None = None,
) -> tuple[Fort14Mesh, dict]:
    """Reclassify a stretch of the land boundary closest to each river
    point as a separate segment with ``river_ibtype``.

    Parameters
    ----------
    mesh:
        Input mesh with classified land boundaries (typically the
        output of :func:`classify_boundaries_by_bbox`).
    points:
        ``(N, 2)`` array of ``(lon, lat)`` river-mouth coordinates in
        EPSG:4326.
    n_nodes_per_river:
        Number of consecutive land-boundary nodes to assign to each
        river segment. Centred on the closest hit; clamped at parent
        endpoints.
    river_ibtype:
        ibtype to write for each river segment (default 21 = FVCOM
        discharge boundary). Must differ from the parent segments'
        ibtype, otherwise the split would be invisible in fort.14.
    snap_tol_m:
        If set, river points whose nearest land node is farther away
        are skipped with a warning entry in ``info["skipped"]``.

    Returns
    -------
    (mesh_out, info)
        ``info`` keys: ``"rivers"`` (list of dicts: ``point``,
        ``snapped_node``, ``dist_m``, ``parent_segment``,
        ``river_n_nodes``), ``"skipped"`` (list of dicts:
        ``point``, ``dist_m``).
    """
    if points.ndim != 2 or points.shape[1] != 2:
        raise ValueError("points must be (N, 2) lon/lat array")
    if n_nodes_per_river < 1:
        raise ValueError("n_nodes_per_river must be >= 1")

    new_land: list[tuple[int, np.ndarray]] = list(mesh.land_boundaries)
    rivers_info: list[dict] = []
    skipped: list[dict] = []

    for pt in points:
        # Find the closest land node, restricting to non-river segments.
        best_seg = -1
        best_idx = -1
        best_dist = np.inf
        for s_idx, (ibtype, ids) in enumerate(new_land):
            if ibtype == river_ibtype:
                continue
            xy = mesh.nodes[ids]
            d = _haversine_m(xy, pt)
            j = int(np.argmin(d))
            if d[j] < best_dist:
                best_dist = float(d[j])
                best_seg = s_idx
                best_idx = j

        if best_seg < 0:
            skipped.append({"point": pt.tolist(), "dist_m": float("inf")})
            continue
        if snap_tol_m is not None and best_dist > snap_tol_m:
            skipped.append({"point": pt.tolist(), "dist_m": best_dist})
            continue

        parent_ibtype, parent_ids = new_land[best_seg]
        prefix, river, suffix = _split_segment(
            parent_ids, best_idx, n_nodes_per_river,
        )
        if river.size == 0:
            skipped.append({"point": pt.tolist(), "dist_m": best_dist})
            continue

        replacement: list[tuple[int, np.ndarray]] = []
        if prefix.size > 0:
            replacement.append((parent_ibtype, prefix))
        replacement.append((int(river_ibtype), river))
        if suffix.size > 0:
            replacement.append((parent_ibtype, suffix))

        new_land = (
            new_land[:best_seg] + replacement + new_land[best_seg + 1 :]
        )

        snapped_node_id = int(parent_ids[best_idx])
        rivers_info.append({
            "point": pt.tolist(),
            "snapped_node": snapped_node_id,
            "dist_m": best_dist,
            "parent_segment_index": best_seg,
            "river_n_nodes": int(river.size),
        })

    out = replace(mesh, land_boundaries=new_land)
    return out, {"rivers": rivers_info, "skipped": skipped}
