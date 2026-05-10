"""Detection of inadequate FVCOM meshes.

Six detectors flag elements / nodes that are likely to cause problems
for FVCOM, especially in narrow water bodies (rivers, canals,
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

    # Node-level flag (shape (NP,))
    overconnected_flag: np.ndarray

    # Auxiliaries
    component_labels: np.ndarray
    component_sizes: np.ndarray
    valence: np.ndarray

    # Configuration
    max_nbr_elem: int
    min_thin_chain: int

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
                self.overconnected_flag,
            )
        )


def run_diagnostics(
    mesh: Fort14Mesh,
    *,
    name: str | None = None,
    path: Path | None = None,
    max_nbr_elem: int = DEFAULT_MAX_NBR_ELEM,
    min_thin_chain: int = DEFAULT_MIN_THIN_CHAIN,
) -> DiagnosticReport:
    """Apply all six detectors to ``mesh`` and return a
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
        overconnected_flag=overconn_flag,
        component_labels=labels,
        component_sizes=sizes,
        valence=valence,
        max_nbr_elem=max_nbr_elem,
        min_thin_chain=min_thin_chain,
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
        f"min_thin_chain={report.min_thin_chain}",
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
        "",
        f"any flagged: {report.any_flagged()}",
    ]
    return "\n".join(lines)


_OVERLAY_STYLE: list[tuple[str, str]] = [
    ("disjoint", "tab:pink"),
    ("dead-end", "tab:orange"),
    ("thin chain", "tab:cyan"),
    ("unreachable", "tab:red"),
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
    "DEFAULT_MAX_NBR_ELEM",
    "DEFAULT_MIN_THIN_CHAIN",
    "DiagnosticReport",
    "boundary_node_mask",
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
    "unreachable_elements_flag",
]
