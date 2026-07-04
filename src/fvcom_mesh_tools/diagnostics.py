"""Detection of inadequate FVCOM meshes.

Seven detectors flag elements / nodes that are likely to cause
problems for FVCOM, especially in narrow water bodies (rivers, canals,
harbours):

    1. ``disjoint_components_flag`` — elements not in the largest
       triangle-dual connected component.
    2. ``dead_end_elements_flag`` — degree-1 elements with no
       open-boundary edge (corner-at-OB excluded).
    3. ``thin_elements_flag`` — elements whose three vertices are all on
       a boundary (open or land); a chain of these is a 1-cell channel.
    4. ``thin_chain_elements_flag`` — thin elements that are part of a
       connected chain of at least ``min_thin_chain`` elements (the
       channel-aware variant of the previous detector).
    5. ``overconnected_nodes_flag`` — nodes incident to more than
       ``max_nbr_elem`` elements; FVCOM allocates a fixed-size element
       neighbour list per node and rejects meshes that exceed its
       compile-time cap.
    6. ``unreachable_elements_flag`` — elements in dual-graph components
       that contain no open-boundary node; such regions cannot be driven
       from the open boundary and are de facto dead pools in FVCOM.
    7. ``under_resolved_channels_flag`` — elements whose local channel
       width (sum of distances to the two nearest non-adjacent boundary
       samples) divided by the local mesh size is below ``min_w_h``;
       catches 2- to 3-cell channels that the 1-cell ``thin`` detector
       misses.

The module exposes per-detector functions plus a single high-level
:func:`run_diagnostics` that returns a :class:`DiagnosticReport` with
every flag array attached. Repair is **out of scope**; downstream
pipelines consume the flagged ids.
"""

from __future__ import annotations

from dataclasses import dataclass, field
from pathlib import Path
from typing import Any

import numpy as np
import scipy.sparse as sp
from scipy.sparse.csgraph import connected_components
from scipy.spatial import cKDTree

from fvcom_mesh_tools.io import Fort14Mesh

# FVCOM's compile-time per-node element-neighbour cap. 8 is the legacy
# conservative value; newer 4.x builds raise this. Always overridable
# via :func:`run_diagnostics(max_nbr_elem=...)`.
DEFAULT_MAX_NBR_ELEM: int = 8

# Minimum length of a connected thin-element run that we treat as a
# 1-cell-wide channel. A single isolated thin element can be a normal
# corner artefact; a chain is strongly indicative of an under-resolved
# channel.
DEFAULT_MIN_THIN_CHAIN: int = 3

# Minimum cells across a channel that FVCOM users typically require for
# flow to be representable. A ratio below this is flagged as
# under-resolved.
DEFAULT_MIN_W_H: float = 3.0

# Boundary-sample spacing in metres for the channel-width detector.
# Smaller = more accurate but more memory / time.
DEFAULT_CHANNEL_SAMPLE_DS_M: float = 50.0

# Along-polyline arc separation factor used to filter out adjacent
# samples on the same polyline. Two boundary samples on the same
# polyline are considered "different banks" only if their along-polyline
# arc separation exceeds ``factor * d1`` where ``d1`` is the distance
# from the query point to the nearest sample. 4.0 is empirically
# sufficient to distinguish a narrow inlet (along-polyline distance
# >> straight-line distance across the inlet) from adjacent samples on
# a smooth coast (along-polyline distance comparable to straight-line).
DEFAULT_ARC_SEPARATION_FACTOR: float = 4.0

# Maximum allowed cos(angle) between (centroid → nearest sample) and
# (centroid → far-arc sample) on the same polyline; the same-polyline
# narrow-inlet candidate is accepted only when this cosine is below
# the threshold (i.e., the two vectors point in roughly opposite
# directions). The default ``-0.8`` requires an angle > 143°, which
# rejects coastal-corner false positives where the same polyline
# wraps around a peninsula tip.
DEFAULT_OPPOSITE_BANK_COS_MAX: float = -0.8

# detector 6 default min-channel-elements filter. ``1`` keeps every
# flagged element (legacy behaviour); raise to N to drop flagged
# elements whose face-face-connected channel component has fewer
# than N members. Useful to suppress the "noise" of small isolated
# clusters (river-mouth corners, jetty tips) flagged as
# under-resolved when only the long ribbon-like channels matter.
DEFAULT_MIN_CHANNEL_ELEMENTS: int = 1

# Earth radius used for the lon/lat → metric projection in the
# channel-width detector.
_EARTH_R_M: float = 6_371_000.0


# ---------------------------------------------------------------------------
# Mesh primitives
# ---------------------------------------------------------------------------


def face_face_adjacency(elements: np.ndarray) -> sp.csr_matrix:
    """Triangle-dual adjacency: two triangles are adjacent iff they share
    an edge (two nodes).

    Returns a symmetric ``(NE, NE)`` ``int8`` CSR matrix with 1 entries on
    adjacent pairs (and 0 on the diagonal).
    """
    ne = len(elements)
    if ne == 0:
        return sp.csr_matrix((0, 0), dtype=np.int8)
    e = np.vstack([
        elements[:, [0, 1]],
        elements[:, [1, 2]],
        elements[:, [2, 0]],
    ])
    e.sort(axis=1)
    elem_of = np.tile(np.arange(ne), 3)
    order = np.lexsort(e.T[::-1])
    e_sorted = e[order]
    elem_sorted = elem_of[order]
    same = np.all(e_sorted[:-1] == e_sorted[1:], axis=1)
    rows = elem_sorted[:-1][same]
    cols = elem_sorted[1:][same]
    if len(rows) == 0:
        return sp.csr_matrix((ne, ne), dtype=np.int8)
    data = np.ones(len(rows), dtype=np.int8)
    a = sp.coo_matrix((data, (rows, cols)), shape=(ne, ne))
    return (a + a.T).tocsr()


def node_valence(elements: np.ndarray, n_nodes: int) -> np.ndarray:
    """Per-node count of incident triangles."""
    counts = np.zeros(n_nodes, dtype=np.int64)
    np.add.at(counts, elements.ravel(), 1)
    return counts


def boundary_node_mask(mesh: Fort14Mesh) -> np.ndarray:
    """Boolean ``(n_nodes,)`` mask that is True where the node belongs to
    any open or land boundary segment.
    """
    mask = np.zeros(mesh.n_nodes, dtype=bool)
    for seg in mesh.open_boundaries:
        mask[np.asarray(seg, dtype=np.int64)] = True
    for _, seg in mesh.land_boundaries:
        mask[np.asarray(seg, dtype=np.int64)] = True
    return mask


def open_boundary_node_mask(mesh: Fort14Mesh) -> np.ndarray:
    mask = np.zeros(mesh.n_nodes, dtype=bool)
    for seg in mesh.open_boundaries:
        mask[np.asarray(seg, dtype=np.int64)] = True
    return mask


def _open_boundary_edge_codes(mesh: Fort14Mesh) -> np.ndarray:
    """Sorted unique codes for OB edges (encoded as ``a * n_nodes + b``,
    ``a < b``). Empty array if there are no open boundaries.
    """
    pairs: list[tuple[int, int]] = []
    n = mesh.n_nodes
    for seg in mesh.open_boundaries:
        seg = np.asarray(seg, dtype=np.int64)
        for i in range(len(seg) - 1):
            a, b = int(seg[i]), int(seg[i + 1])
            if a == b:
                continue
            pairs.append((min(a, b), max(a, b)))
    if not pairs:
        return np.empty(0, dtype=np.int64)
    arr = np.asarray(sorted(set(pairs)), dtype=np.int64)
    codes = arr[:, 0] * n + arr[:, 1]
    return np.sort(codes)


def _element_edge_codes(elements: np.ndarray, n_nodes: int) -> np.ndarray:
    """Per-element edge codes ``(NE, 3)`` with ``a < b`` and code
    ``a * n_nodes + b``.
    """
    ee = np.stack([
        np.sort(elements[:, [0, 1]], axis=1),
        np.sort(elements[:, [1, 2]], axis=1),
        np.sort(elements[:, [2, 0]], axis=1),
    ], axis=1)
    return ee[..., 0] * n_nodes + ee[..., 1]


def element_centroids(mesh: Fort14Mesh) -> np.ndarray:
    return mesh.nodes[mesh.elements].mean(axis=1)


# ---------------------------------------------------------------------------
# Detectors. Each returns a bool flag array shaped (NE,) or (NP,) plus, if
# useful, an ancillary array (component labels, valence, ...).
# ---------------------------------------------------------------------------


def disjoint_components_flag(
    adj: sp.csr_matrix,
) -> tuple[np.ndarray, np.ndarray, np.ndarray]:
    """Flag elements not in the largest dual-graph connected component.

    Returns ``(flag, labels, sizes)`` where ``labels`` is the per-element
    component id and ``sizes`` is the histogram of component sizes
    indexed by component id.
    """
    if adj.shape[0] == 0:
        empty_i = np.zeros(0, dtype=np.int64)
        return np.zeros(0, dtype=bool), empty_i, empty_i
    _n_comp, labels = connected_components(adj, directed=False, return_labels=True)
    sizes = np.bincount(labels)
    largest = int(sizes.argmax())
    return labels != largest, labels, sizes


def dead_end_elements_flag(adj: sp.csr_matrix, mesh: Fort14Mesh) -> np.ndarray:
    """Degree-1 elements that do not have an open-boundary edge.

    A degree-1 element has 1 internal edge and 2 boundary edges. Filtering
    out those whose 2 boundary edges include an open-boundary edge skips
    legitimate OB corner elements; what remains is the spit / dangling
    triangle that terminates a 1-cell-wide channel.
    """
    if adj.shape[0] == 0:
        return np.zeros(0, dtype=bool)
    deg = np.asarray(adj.sum(axis=1)).ravel()
    deg1 = deg == 1
    if not deg1.any():
        return deg1
    ob_codes = _open_boundary_edge_codes(mesh)
    if ob_codes.size == 0:
        return deg1
    elem_codes = _element_edge_codes(mesh.elements, mesh.n_nodes)
    flat = elem_codes.ravel()
    idx = np.searchsorted(ob_codes, flat)
    idx_clip = np.minimum(idx, len(ob_codes) - 1)
    has_ob = (ob_codes[idx_clip] == flat).reshape(-1, 3).any(axis=1)
    return deg1 & ~has_ob


def thin_elements_flag(mesh: Fort14Mesh) -> np.ndarray:
    """Elements whose three vertices are all on a boundary (any type)."""
    bdy = boundary_node_mask(mesh)
    return bdy[mesh.elements].all(axis=1)


def thin_chain_elements_flag(
    adj: sp.csr_matrix, thin_flag: np.ndarray, *, min_chain_length: int,
) -> np.ndarray:
    """Thin elements that belong to a connected chain of at least
    ``min_chain_length`` thin elements.

    The plain :func:`thin_elements_flag` fires on isolated boundary
    triangles too (legitimate corner artefacts). Restricting the dual
    graph to thin elements and requiring a minimum component size isolates
    runs that line up along a 1-cell-wide channel.
    """
    if adj.shape[0] == 0 or not thin_flag.any() or min_chain_length <= 1:
        return thin_flag.copy()
    thin_idx = np.where(thin_flag)[0]
    sub = adj[thin_idx][:, thin_idx]
    _n_comp, labels = connected_components(sub, directed=False, return_labels=True)
    sizes = np.bincount(labels)
    keep_local = sizes[labels] >= min_chain_length
    out = np.zeros_like(thin_flag)
    out[thin_idx] = keep_local
    return out


def overconnected_nodes_flag(
    elements: np.ndarray, n_nodes: int, *, max_nbr: int = DEFAULT_MAX_NBR_ELEM,
) -> tuple[np.ndarray, np.ndarray]:
    """Flag nodes whose valence exceeds ``max_nbr``.

    Returns ``(flag, valence)``.
    """
    val = node_valence(elements, n_nodes)
    return val > max_nbr, val


def unreachable_elements_flag(
    adj: sp.csr_matrix,
    mesh: Fort14Mesh,
    component_labels: np.ndarray | None = None,
) -> np.ndarray:
    """Elements in dual-graph components that contain no open-boundary
    node. If the mesh has no open boundaries, returns all-False.
    """
    if adj.shape[0] == 0:
        return np.zeros(0, dtype=bool)
    if not mesh.open_boundaries:
        return np.zeros(mesh.n_elements, dtype=bool)
    ob_node_mask = open_boundary_node_mask(mesh)
    elem_has_ob_node = ob_node_mask[mesh.elements].any(axis=1)
    if not elem_has_ob_node.any():
        return np.ones(mesh.n_elements, dtype=bool)
    if component_labels is None:
        _n, component_labels = connected_components(
            adj, directed=False, return_labels=True,
        )
    open_components = np.unique(component_labels[elem_has_ob_node])
    reachable = np.isin(component_labels, open_components)
    return ~reachable


# ---------------------------------------------------------------------------
# Detector 7: under-resolved channels (medial-axis-style w/h ratio)
# ---------------------------------------------------------------------------


def _to_metric(
    nodes_lonlat: np.ndarray, *, lat0: float, lon0: float,
) -> np.ndarray:
    """Flat-earth lon/lat → metric (x_m, y_m) projection about
    ``(lon0, lat0)``. Accurate enough for distance metrics over a
    single coastal basin.
    """
    cos_lat0 = np.cos(np.deg2rad(lat0))
    x = (nodes_lonlat[:, 0] - lon0) * np.deg2rad(1.0) * _EARTH_R_M * cos_lat0
    y = (nodes_lonlat[:, 1] - lat0) * np.deg2rad(1.0) * _EARTH_R_M
    return np.column_stack([x, y])


def _sample_polyline_with_arc(
    nodes_m: np.ndarray, seg_ids: np.ndarray, ds: float,
) -> tuple[np.ndarray, np.ndarray]:
    """Sample a polyline at ~``ds`` spacing in metres.

    Returns ``(samples_xy, arc_pos)`` where ``arc_pos[i]`` is the
    cumulative arc length (in metres) from the polyline's first node
    to ``samples_xy[i]``. Always emits at least 2 samples per edge.
    """
    seg_ids = np.asarray(seg_ids, dtype=np.int64)
    if seg_ids.size < 2:
        return np.empty((0, 2), dtype=np.float64), np.empty((0,), dtype=np.float64)
    pts: list[np.ndarray] = []
    arcs: list[np.ndarray] = []
    cum_arc = 0.0
    for i in range(len(seg_ids) - 1):
        p0 = nodes_m[seg_ids[i]]
        p1 = nodes_m[seg_ids[i + 1]]
        L = float(np.linalg.norm(p1 - p0))
        n = max(2, int(np.ceil(L / ds)) + 1)
        # Emit n samples along this edge; drop the last to avoid
        # duplicating the next edge's first sample (we add the polyline
        # endpoint at the very end).
        t = np.linspace(0.0, 1.0, n, endpoint=False)
        pts.append(p0 + t[:, None] * (p1 - p0))
        arcs.append(cum_arc + t * L)
        cum_arc += L
    # Append the final endpoint.
    pts.append(nodes_m[seg_ids[-1:]].astype(np.float64))
    arcs.append(np.array([cum_arc], dtype=np.float64))
    return np.vstack(pts), np.concatenate(arcs)


def _median_edge_length_m(
    nodes_m: np.ndarray, elements: np.ndarray,
) -> np.ndarray:
    """Per-element median of the 3 edge lengths in metres."""
    p0 = nodes_m[elements[:, 0]]
    p1 = nodes_m[elements[:, 1]]
    p2 = nodes_m[elements[:, 2]]
    e01 = np.linalg.norm(p1 - p0, axis=1)
    e12 = np.linalg.norm(p2 - p1, axis=1)
    e20 = np.linalg.norm(p0 - p2, axis=1)
    return np.median(np.column_stack([e01, e12, e20]), axis=1)


def channel_width_metric(
    mesh: Fort14Mesh,
    *,
    sample_ds_m: float = DEFAULT_CHANNEL_SAMPLE_DS_M,
    arc_separation_factor: float = DEFAULT_ARC_SEPARATION_FACTOR,
    opposite_bank_cos_max: float = DEFAULT_OPPOSITE_BANK_COS_MAX,
    k_far: int = 50,
    coords: str = "lonlat",
) -> dict[str, np.ndarray]:
    """Per-element channel width (m), local h (m), and the w/h ratio.

    The channel width estimate at an element centroid is the smaller of
    two candidates:

        1. *Cross-polyline*: ``d(P, polyline_A) + d(P, polyline_B)``
           where ``A`` and ``B`` are the two distinct polylines whose
           nearest samples to ``P`` are smallest. Catches the "channel
           between mainland and an island" case.
        2. *Same-polyline narrow inlet*: ``min over p of
           (d(P, polyline_p) + d_far_arc(P, polyline_p))`` where
           ``d_far_arc`` is the distance from ``P`` to the nearest
           sample on ``polyline_p`` whose along-polyline arc separation
           from the absolute-nearest sample is at least
           ``arc_separation_factor * d(P, polyline_p)``. Catches a
           narrow inlet whose two banks lie on the same continuous
           coastline polyline.

    A separate ``cKDTree`` is built per polyline, so cross-polyline
    distance is exact and same-polyline arc-filtered probing is local
    to each polyline (the previous combined-tree implementation could
    miss the "other" polyline when the nearest one was densely
    sampled).

    Returns a dict with keys ``"channel_width_m"``, ``"h_local_m"``,
    ``"w_h_ratio"``, and ``"sample_count"`` (total boundary samples).
    The ratio is set to ``+inf`` when the mesh has no boundary at all.

    ``coords`` selects the node-coordinate interpretation: ``"lonlat"``
    (default, legacy behaviour — flat-earth projection to metres) or
    ``"metric"`` (coordinates are already metres, e.g. UTM; used as-is).
    """
    if coords not in ("lonlat", "metric"):
        raise ValueError(f"coords must be 'lonlat' or 'metric', got {coords!r}")
    if mesh.n_elements == 0:
        return {
            "channel_width_m": np.zeros(0),
            "h_local_m": np.zeros(0),
            "w_h_ratio": np.zeros(0),
            "sample_count": 0,
        }
    if not mesh.open_boundaries and not mesh.land_boundaries:
        return {
            "channel_width_m": np.full(mesh.n_elements, np.inf),
            "h_local_m": np.zeros(mesh.n_elements),
            "w_h_ratio": np.full(mesh.n_elements, np.inf),
            "sample_count": 0,
        }

    if coords == "metric":
        nodes_m = np.asarray(mesh.nodes, dtype=np.float64)
    else:
        lat0 = float(mesh.nodes[:, 1].mean())
        lon0 = float(mesh.nodes[:, 0].mean())
        nodes_m = _to_metric(mesh.nodes, lat0=lat0, lon0=lon0)

    polylines: list[tuple[np.ndarray, np.ndarray]] = []
    for seg in mesh.open_boundaries:
        pts, arc = _sample_polyline_with_arc(nodes_m, seg, sample_ds_m)
        if pts.size:
            polylines.append((pts, arc))
    for _ib, seg in mesh.land_boundaries:
        pts, arc = _sample_polyline_with_arc(nodes_m, seg, sample_ds_m)
        if pts.size:
            polylines.append((pts, arc))

    n_poly = len(polylines)
    sample_count = int(sum(len(pts) for pts, _ in polylines))
    if n_poly == 0:
        return {
            "channel_width_m": np.full(mesh.n_elements, np.inf),
            "h_local_m": np.zeros(mesh.n_elements),
            "w_h_ratio": np.full(mesh.n_elements, np.inf),
            "sample_count": 0,
        }

    centroids = nodes_m[mesh.elements].mean(axis=1)
    NE = mesh.n_elements
    rows = np.arange(NE)
    d_per_poly = np.empty((NE, n_poly), dtype=np.float64)
    d_far_per_poly = np.full((NE, n_poly), np.inf, dtype=np.float64)

    for p, (pts, arcs) in enumerate(polylines):
        tree = cKDTree(pts)
        # workers=-1: the per-polyline KD queries dominate QA wall
        # time on large meshes; use every available core.
        d1, idx1 = tree.query(centroids, k=1, workers=-1)
        d_per_poly[:, p] = d1
        arc_at_nearest = arcs[idx1]
        K = min(k_far, len(pts))
        if K < 2:
            continue
        dK, idxK = tree.query(centroids, k=K, workers=-1)
        arc_K = arcs[idxK]
        arc_diff = np.abs(arc_K - arc_at_nearest[:, None])
        threshold = arc_separation_factor * d1[:, None]
        far_enough = arc_diff >= threshold
        far_enough[:, 0] = False
        has_far = far_enough.any(axis=1)
        first_far_k = np.argmax(far_enough, axis=1)
        d_far_p = dK[rows, first_far_k]
        far_arc_idx = idxK[rows, first_far_k]

        # Direction filter: only accept the far-arc sample as the "other
        # bank" if it sits on the opposite side of the centroid from the
        # nearest sample. cos(angle) between (centroid -> nearest) and
        # (centroid -> far-arc) must be below -0.5 (angle > 120°). This
        # rejects coastal-corner false positives where the same polyline
        # wraps around a peninsula tip — both samples point in roughly
        # the same direction from the centroid, so cos(angle) is positive.
        s1 = pts[idx1]
        s_far = pts[far_arc_idx]
        v1 = s1 - centroids
        v2 = s_far - centroids
        dot = (v1 * v2).sum(axis=1)
        norm_prod = d1 * d_far_p
        cos_angle = dot / np.where(norm_prod > 1e-9, norm_prod, 1e-9)
        direction_ok = cos_angle < opposite_bank_cos_max

        accepted = has_far & direction_ok
        d_far_per_poly[:, p] = np.where(accepted, d_far_p, np.inf)

    if n_poly >= 2:
        sorted_d = np.sort(d_per_poly, axis=1)
        cross_w = sorted_d[:, 0] + sorted_d[:, 1]
    else:
        cross_w = np.full(NE, np.inf)
    same_w = np.min(d_per_poly + d_far_per_poly, axis=1)
    channel_width = np.minimum(cross_w, same_w)

    h_local = _median_edge_length_m(nodes_m, mesh.elements)
    ratio = channel_width / np.where(h_local > 0, h_local, 1.0)
    return {
        "channel_width_m": channel_width,
        "h_local_m": h_local,
        "w_h_ratio": ratio,
        "sample_count": sample_count,
    }


def under_resolved_channels_flag(
    mesh: Fort14Mesh,
    *,
    min_w_h: float = DEFAULT_MIN_W_H,
    sample_ds_m: float = DEFAULT_CHANNEL_SAMPLE_DS_M,
    arc_separation_factor: float = DEFAULT_ARC_SEPARATION_FACTOR,
    opposite_bank_cos_max: float = DEFAULT_OPPOSITE_BANK_COS_MAX,
    min_channel_elements: int = DEFAULT_MIN_CHANNEL_ELEMENTS,
    coords: str = "lonlat",
) -> tuple[np.ndarray, dict[str, np.ndarray]]:
    """Element-level flag: True if ``w_h_ratio < min_w_h``.

    Returns ``(flag, metric_info)`` where ``metric_info`` is the dict
    returned by :func:`channel_width_metric`.

    ``min_channel_elements`` (default 1 = no filter) drops any flagged
    element whose face-face-connected component has fewer than
    ``min_channel_elements`` flagged members. Use to suppress small
    isolated clusters (river-mouth corners, jetty tips) when only
    longer ribbon-like channels are of interest. PoC #35 found that
    on Tokyo Bay the 3,178 default-flag elements split into 1,010
    components with mean ~3 elements / component; a filter of, say,
    10 keeps the dozens of meaningful long channels while pruning
    the 970+ small-cluster noise.
    """
    info = channel_width_metric(
        mesh,
        sample_ds_m=sample_ds_m,
        arc_separation_factor=arc_separation_factor,
        opposite_bank_cos_max=opposite_bank_cos_max,
        coords=coords,
    )
    flag = info["w_h_ratio"] < min_w_h
    if min_channel_elements <= 1 or not flag.any():
        return flag, info

    # Drop flagged elements whose channel component is too small.
    flagged_idx = np.where(flag)[0]
    full_adj = face_face_adjacency(mesh.elements)
    sub = full_adj[flagged_idx][:, flagged_idx]
    n_comp, labels = connected_components(sub, directed=False, return_labels=True)
    sizes = np.bincount(labels, minlength=int(n_comp))
    small_component = sizes < int(min_channel_elements)
    drop_mask = small_component[labels]   # True for flagged elements to drop
    flag = flag.copy()
    flag[flagged_idx[drop_mask]] = False
    return flag, info


# ---------------------------------------------------------------------------
# High-level driver
# ---------------------------------------------------------------------------


@dataclass
class DiagnosticReport:
    """Aggregated detector results for one mesh.

    All flag arrays are bool; element-level arrays have shape ``(NE,)``
    and ``overconnected_flag`` has shape ``(NP,)``. ``valence`` and
    ``component_labels`` are kept so callers can query the underlying
    metrics directly.
    """

    mesh_name: str
    mesh_path: Path | None

    n_nodes: int
    n_elements: int
    n_open_boundaries: int
    n_land_boundaries: int

    # Element-level flags (shape (NE,))
    disjoint_flag: np.ndarray
    dead_end_flag: np.ndarray
    thin_flag: np.ndarray
    thin_chain_flag: np.ndarray
    unreachable_flag: np.ndarray
    under_resolved_channels_flag: np.ndarray

    # Node-level flag (shape (NP,))
    overconnected_flag: np.ndarray

    # Auxiliaries
    component_labels: np.ndarray
    component_sizes: np.ndarray
    valence: np.ndarray
    channel_width_m: np.ndarray
    h_local_m: np.ndarray
    w_h_ratio: np.ndarray

    # Configuration
    max_nbr_elem: int
    min_thin_chain: int
    min_w_h: float

    # Cached coordinates so JSON / plotting consumers do not need the
    # original mesh again.
    nodes: np.ndarray = field(repr=False)
    elements: np.ndarray = field(repr=False)
    centroids: np.ndarray = field(repr=False)

    def any_flagged(self) -> bool:
        return any(
            arr.any()
            for arr in (
                self.disjoint_flag, self.dead_end_flag, self.thin_flag,
                self.thin_chain_flag, self.unreachable_flag,
                self.overconnected_flag, self.under_resolved_channels_flag,
            )
        )


def run_diagnostics(
    mesh: Fort14Mesh,
    *,
    name: str | None = None,
    path: Path | None = None,
    max_nbr_elem: int = DEFAULT_MAX_NBR_ELEM,
    min_thin_chain: int = DEFAULT_MIN_THIN_CHAIN,
    min_w_h: float = DEFAULT_MIN_W_H,
    channel_sample_ds_m: float = DEFAULT_CHANNEL_SAMPLE_DS_M,
    channel_arc_separation_factor: float = DEFAULT_ARC_SEPARATION_FACTOR,
    channel_opposite_bank_cos_max: float = DEFAULT_OPPOSITE_BANK_COS_MAX,
    min_channel_elements: int = DEFAULT_MIN_CHANNEL_ELEMENTS,
) -> DiagnosticReport:
    """Apply all seven detectors to ``mesh`` and return a
    :class:`DiagnosticReport`.

    Repair is out of scope; the report carries the per-element /
    per-node flag arrays unchanged.
    """
    adj = face_face_adjacency(mesh.elements)
    disjoint_flag, labels, sizes = disjoint_components_flag(adj)
    dead_end_flag = dead_end_elements_flag(adj, mesh)
    thin_flag = thin_elements_flag(mesh)
    thin_chain_flag = thin_chain_elements_flag(
        adj, thin_flag, min_chain_length=min_thin_chain,
    )
    overconn_flag, valence = overconnected_nodes_flag(
        mesh.elements, mesh.n_nodes, max_nbr=max_nbr_elem,
    )
    unreach_flag = unreachable_elements_flag(adj, mesh, component_labels=labels)
    under_flag, ch_info = under_resolved_channels_flag(
        mesh,
        min_w_h=min_w_h,
        sample_ds_m=channel_sample_ds_m,
        arc_separation_factor=channel_arc_separation_factor,
        opposite_bank_cos_max=channel_opposite_bank_cos_max,
        min_channel_elements=min_channel_elements,
    )

    return DiagnosticReport(
        mesh_name=name or (path.name if path else "mesh"),
        mesh_path=path,
        n_nodes=mesh.n_nodes,
        n_elements=mesh.n_elements,
        n_open_boundaries=len(mesh.open_boundaries),
        n_land_boundaries=len(mesh.land_boundaries),
        disjoint_flag=disjoint_flag,
        dead_end_flag=dead_end_flag,
        thin_flag=thin_flag,
        thin_chain_flag=thin_chain_flag,
        unreachable_flag=unreach_flag,
        under_resolved_channels_flag=under_flag,
        overconnected_flag=overconn_flag,
        component_labels=labels,
        component_sizes=sizes,
        valence=valence,
        channel_width_m=ch_info["channel_width_m"],
        h_local_m=ch_info["h_local_m"],
        w_h_ratio=ch_info["w_h_ratio"],
        max_nbr_elem=max_nbr_elem,
        min_thin_chain=min_thin_chain,
        min_w_h=float(min_w_h),
        nodes=mesh.nodes,
        elements=mesh.elements,
        centroids=element_centroids(mesh),
    )


# ---------------------------------------------------------------------------
# Output helpers
# ---------------------------------------------------------------------------


def _elem_records(
    flag: np.ndarray, centroids: np.ndarray, *, labels: np.ndarray | None = None,
) -> list[dict[str, Any]]:
    idx = np.where(flag)[0]
    out: list[dict[str, Any]] = []
    for i in idx:
        rec: dict[str, Any] = {
            "id": int(i),
            "centroid_lon": float(centroids[i, 0]),
            "centroid_lat": float(centroids[i, 1]),
        }
        if labels is not None:
            rec["component"] = int(labels[i])
        out.append(rec)
    return out


def _node_records(
    flag: np.ndarray, nodes: np.ndarray, *, valence: np.ndarray,
) -> list[dict[str, Any]]:
    idx = np.where(flag)[0]
    return [
        {
            "id": int(i),
            "lon": float(nodes[i, 0]),
            "lat": float(nodes[i, 1]),
            "valence": int(valence[i]),
        }
        for i in idx
    ]


def _under_resolved_records(
    flag: np.ndarray,
    centroids: np.ndarray,
    *,
    channel_width_m: np.ndarray,
    h_local_m: np.ndarray,
    w_h_ratio: np.ndarray,
) -> list[dict[str, Any]]:
    idx = np.where(flag)[0]
    return [
        {
            "id": int(i),
            "centroid_lon": float(centroids[i, 0]),
            "centroid_lat": float(centroids[i, 1]),
            "channel_width_m": float(channel_width_m[i]),
            "h_local_m": float(h_local_m[i]),
            "w_h_ratio": float(w_h_ratio[i]),
        }
        for i in idx
    ]


def report_to_dict(report: DiagnosticReport) -> dict[str, Any]:
    """JSON-friendly representation of a :class:`DiagnosticReport`."""
    return {
        "mesh": {
            "name": report.mesh_name,
            "path": str(report.mesh_path) if report.mesh_path else None,
            "n_nodes": int(report.n_nodes),
            "n_elements": int(report.n_elements),
            "n_open_boundaries": int(report.n_open_boundaries),
            "n_land_boundaries": int(report.n_land_boundaries),
        },
        "config": {
            "max_nbr_elem": int(report.max_nbr_elem),
            "min_thin_chain": int(report.min_thin_chain),
        },
        "detectors": {
            "disjoint_components": {
                "description": (
                    "Elements not in the largest dual-graph connected "
                    "component."
                ),
                "n_components": int(len(report.component_sizes)),
                "component_sizes": [int(s) for s in report.component_sizes.tolist()],
                "n_flagged": int(report.disjoint_flag.sum()),
                "elements": _elem_records(
                    report.disjoint_flag, report.centroids,
                    labels=report.component_labels,
                ),
            },
            "dead_end_elements": {
                "description": (
                    "Degree-1 elements with no open-boundary edge "
                    "(corner-at-OB cases excluded)."
                ),
                "n_flagged": int(report.dead_end_flag.sum()),
                "elements": _elem_records(report.dead_end_flag, report.centroids),
            },
            "thin_elements": {
                "description": (
                    "Elements whose 3 vertices are all on a boundary "
                    "(open or land)."
                ),
                "n_flagged": int(report.thin_flag.sum()),
                "elements": _elem_records(report.thin_flag, report.centroids),
            },
            "thin_chain_elements": {
                "description": (
                    "Thin elements forming a chain of at least "
                    f"{report.min_thin_chain} adjacent thin elements; "
                    "indicates a 1-cell-wide channel run."
                ),
                "min_chain_length": int(report.min_thin_chain),
                "n_flagged": int(report.thin_chain_flag.sum()),
                "elements": _elem_records(report.thin_chain_flag, report.centroids),
            },
            "overconnected_nodes": {
                "description": (
                    f"Nodes incident to more than {report.max_nbr_elem} "
                    "elements (FVCOM MAX_NBR_ELEM cap)."
                ),
                "max_nbr_elem": int(report.max_nbr_elem),
                "max_valence_observed": int(report.valence.max())
                    if report.valence.size else 0,
                "n_flagged": int(report.overconnected_flag.sum()),
                "nodes": _node_records(
                    report.overconnected_flag, report.nodes,
                    valence=report.valence,
                ),
            },
            "open_boundary_unreachable": {
                "description": (
                    "Elements in components without any open-boundary node."
                ),
                "n_flagged": int(report.unreachable_flag.sum()),
                "elements": _elem_records(
                    report.unreachable_flag, report.centroids,
                ),
            },
            "under_resolved_channels": {
                "description": (
                    "Elements whose local channel width (sum of distances "
                    "to the two nearest non-adjacent boundary samples) "
                    "divided by the median edge length is below "
                    f"min_w_h={report.min_w_h:.1f}. Catches 2- and 3-cell "
                    "narrow channels that the 1-cell thin-chain detector "
                    "misses."
                ),
                "min_w_h": float(report.min_w_h),
                "n_flagged": int(report.under_resolved_channels_flag.sum()),
                "w_h_ratio_p50": float(np.percentile(report.w_h_ratio, 50))
                if report.w_h_ratio.size else 0.0,
                "elements": _under_resolved_records(
                    report.under_resolved_channels_flag,
                    report.centroids,
                    channel_width_m=report.channel_width_m,
                    h_local_m=report.h_local_m,
                    w_h_ratio=report.w_h_ratio,
                ),
            },
        },
    }


def report_to_summary_text(report: DiagnosticReport) -> str:
    """Human-readable single-mesh summary."""
    sizes = report.component_sizes
    sizes_str = (
        str([int(s) for s in sizes.tolist()])
        if len(sizes) <= 6
        else str([int(s) for s in sizes[:6].tolist()] + ["..."])
    )
    max_v = int(report.valence.max()) if report.valence.size else 0
    lines = [
        f"mesh:    {report.mesh_path}" if report.mesh_path
        else f"mesh:    {report.mesh_name}",
        f"NP={report.n_nodes:,}  NE={report.n_elements:,}  "
        f"open={report.n_open_boundaries}  land={report.n_land_boundaries}",
        f"config:  max_nbr_elem={report.max_nbr_elem}  "
        f"min_thin_chain={report.min_thin_chain}  "
        f"min_w_h={report.min_w_h:.1f}",
        "",
        f"  1. disjoint comp:    n_components={len(sizes)}  "
        f"sizes={sizes_str}  flagged_elems={int(report.disjoint_flag.sum()):,}",
        f"  2. dead-end elems:   {int(report.dead_end_flag.sum()):,}",
        f"  3. thin elems:       {int(report.thin_flag.sum()):,}",
        f"  3b. thin chains "
        f"(>= {report.min_thin_chain}): "
        f"{int(report.thin_chain_flag.sum()):,}",
        f"  4. over-conn nodes:  {int(report.overconnected_flag.sum()):,}  "
        f"max_valence={max_v}",
        f"  5. unreachable:      {int(report.unreachable_flag.sum()):,}",
        f"  6. under-resolved channels (w/h < {report.min_w_h:.1f}): "
        f"{int(report.under_resolved_channels_flag.sum()):,}  "
        f"w_h p50="
        f"{float(np.percentile(report.w_h_ratio, 50)):.1f}"
        if report.w_h_ratio.size else
        f"  6. under-resolved channels (w/h < {report.min_w_h:.1f}): 0",
        "",
        f"any flagged: {report.any_flagged()}",
    ]
    return "\n".join(lines)


_OVERLAY_STYLE: list[tuple[str, str]] = [
    ("disjoint", "tab:pink"),
    ("dead-end", "tab:orange"),
    ("thin chain", "tab:cyan"),
    ("unreachable", "tab:red"),
    ("under-resolved", "tab:olive"),
]


def plot_report(report: DiagnosticReport, png: Path, *, dpi: int | None = None) -> None:
    """Render the mesh with detector flags overlaid as coloured patches.

    matplotlib is imported lazily so callers that only need text/JSON
    output do not pull it in.
    """
    import matplotlib.pyplot as plt
    from matplotlib.collections import PolyCollection

    from fvcom_mesh_tools.plotting import MESH_PNG_DPI
    if dpi is None:
        dpi = MESH_PNG_DPI

    fig, ax = plt.subplots(figsize=(11, 9), dpi=120)
    ax.triplot(
        report.nodes[:, 0], report.nodes[:, 1], report.elements,
        color="0.85", lw=0.1,
    )

    # Boundaries — drawn from the report's nodes plus a re-derived edge
    # collection. We keep this lightweight: the plot helper does not
    # carry the original Fort14Mesh.boundaries lists, only the per-element
    # / per-node flags. Boundary geometry is therefore reconstructed
    # implicitly via the thin/over-connected overlays; callers that want
    # rich boundary colouring can use the JSON output to drive a custom
    # plot.

    overlays: dict[str, np.ndarray] = {
        "disjoint": report.disjoint_flag,
        "dead-end": report.dead_end_flag,
        "thin chain": report.thin_chain_flag,
        "unreachable": report.unreachable_flag,
        "under-resolved": report.under_resolved_channels_flag,
    }
    proxies = []
    for label, color in _OVERLAY_STYLE:
        flag = overlays[label]
        if not flag.any():
            continue
        tris = report.nodes[report.elements[flag]]
        pc = PolyCollection(
            tris, facecolor=color, edgecolor=color,
            alpha=0.55, linewidths=0.2,
        )
        ax.add_collection(pc)
        proxies.append(plt.Line2D(
            [], [], marker="s", linestyle="",
            markerfacecolor=color, markeredgecolor=color,
            markersize=8, alpha=0.7,
            label=f"{label} ({int(flag.sum())})",
        ))

    if report.overconnected_flag.any():
        ax.scatter(
            report.nodes[report.overconnected_flag, 0],
            report.nodes[report.overconnected_flag, 1],
            color="black", marker="x", s=35, linewidths=1.0, zorder=5,
        )
        proxies.append(plt.Line2D(
            [], [], marker="x", linestyle="", color="black", markersize=8,
            label=f"overconnected nodes ({int(report.overconnected_flag.sum())})",
        ))

    ax.set_aspect("equal")
    ax.set_xlabel("lon (deg)")
    ax.set_ylabel("lat (deg)")
    ax.set_title(
        f"fmesh-mesh-check  {report.mesh_name}  "
        f"NP={report.n_nodes:,}  NE={report.n_elements:,}"
    )
    if proxies:
        ax.legend(handles=proxies, loc="lower right", fontsize=8)
    fig.tight_layout()
    fig.savefig(png, bbox_inches="tight", dpi=dpi)
    plt.close(fig)


__all__ = [
    "DEFAULT_ARC_SEPARATION_FACTOR",
    "DEFAULT_CHANNEL_SAMPLE_DS_M",
    "DEFAULT_MAX_NBR_ELEM",
    "DEFAULT_MIN_CHANNEL_ELEMENTS",
    "DEFAULT_MIN_THIN_CHAIN",
    "DEFAULT_MIN_W_H",
    "DEFAULT_OPPOSITE_BANK_COS_MAX",
    "DiagnosticReport",
    "boundary_node_mask",
    "channel_width_metric",
    "dead_end_elements_flag",
    "disjoint_components_flag",
    "element_centroids",
    "face_face_adjacency",
    "node_valence",
    "open_boundary_node_mask",
    "overconnected_nodes_flag",
    "plot_report",
    "report_to_dict",
    "report_to_summary_text",
    "run_diagnostics",
    "thin_chain_elements_flag",
    "thin_elements_flag",
    "under_resolved_channels_flag",
    "unreachable_elements_flag",
]
