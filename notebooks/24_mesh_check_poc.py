"""PoC #24: detect inadequate FVCOM meshes (no repair).

Applies 5 detectors to existing fort.14 outputs to surface mesh defects
that are typical in narrow water bodies (rivers, canals, harbours):

    1. Disjoint wet-domain components
       Triangle-dual ``connected_components``; flag elements not in the
       largest component.
    2. Dead-end elements
       Triangles with only one dual-graph neighbour AND no edge on the
       open boundary (excludes legitimate corner elements at the open
       boundary, isolates "spit" / 1-cell-wide channel terminations).
    3. Thin elements
       Triangles whose three vertices are all on a boundary (open or
       land). A chain of such elements is a 1-element-wide channel.
    4. Over-connected nodes
       Nodes with more than ``MAX_NBR_ELEM`` (=8) incident triangles.
       FVCOM allocates a fixed-size element neighbour list per node; a
       value exceeding the build's ``MAX_NBR_ELEM`` is a hard runtime
       error. 8 is the conservative legacy default; newer builds permit
       more, so the threshold is configurable.
    5. Open-boundary unreachable elements
       Elements whose dual-graph component contains no open-boundary
       node. Such regions cannot be driven from the open boundary and
       are de facto dead pools in FVCOM.

Inputs (must already exist; run notebooks 16, 19, 20 first if missing):
    outputs/16_tokyo_bay_with_rivers.14
    outputs/19_tokyo_bay_oceanmesh.14
    outputs/20_osaka_bay_oceanmesh.14

Outputs (per input mesh ``<name>``):
    outputs/24_mesh_check_summary.txt
    outputs/24_mesh_check_<name>_map.png
    outputs/24_mesh_check_<name>_diag.json
"""

from __future__ import annotations

import json
from pathlib import Path

import matplotlib

matplotlib.use("Agg")

import matplotlib.pyplot as plt  # noqa: E402
import numpy as np  # noqa: E402
import scipy.sparse as sp  # noqa: E402
from matplotlib.collections import PolyCollection  # noqa: E402
from scipy.sparse.csgraph import connected_components  # noqa: E402

from fvcom_mesh_tools.io import Fort14Mesh, read_fort14  # noqa: E402
from fvcom_mesh_tools.plotting import MESH_PNG_DPI  # noqa: E402

REPO_ROOT = Path(__file__).resolve().parent.parent
OUT_DIR = REPO_ROOT / "outputs"
SUMMARY_TXT = OUT_DIR / "24_mesh_check_summary.txt"

INPUTS: list[tuple[str, Path]] = [
    ("19_tokyo_bay_oceanmesh", OUT_DIR / "19_tokyo_bay_oceanmesh.14"),
    ("16_tokyo_bay_with_rivers", OUT_DIR / "16_tokyo_bay_with_rivers.14"),
    ("20_osaka_bay_oceanmesh", OUT_DIR / "20_osaka_bay_oceanmesh.14"),
]

# FVCOM compile-time per-node element-neighbour cap. 8 is the legacy
# conservative value; newer builds (4.x) raise this. Keep configurable.
MAX_NBR_ELEM = 8


# ---------------------------------------------------------------------------
# Mesh primitives
# ---------------------------------------------------------------------------


def face_face_adjacency(elements: np.ndarray) -> sp.csr_matrix:
    """Triangle-dual adjacency: two triangles adjacent iff they share an edge."""
    ne = len(elements)
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


def open_boundary_edge_codes(mesh: Fort14Mesh) -> np.ndarray:
    """Sorted unique codes for OB edges, encoded as ``a * n_nodes + b``."""
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


def element_edge_codes(elements: np.ndarray, n_nodes: int) -> np.ndarray:
    """Per-element edge codes (NE, 3), each entry ``a * n_nodes + b`` with a < b."""
    ee = np.stack([
        np.sort(elements[:, [0, 1]], axis=1),
        np.sort(elements[:, [1, 2]], axis=1),
        np.sort(elements[:, [2, 0]], axis=1),
    ], axis=1)  # (NE, 3, 2)
    return ee[..., 0] * n_nodes + ee[..., 1]


def boundary_node_mask(mesh: Fort14Mesh) -> np.ndarray:
    """Boolean (n_nodes,): True if the node is on any open or land boundary."""
    mask = np.zeros(mesh.n_nodes, dtype=bool)
    for seg in mesh.open_boundaries:
        mask[np.asarray(seg, dtype=np.int64)] = True
    for _, seg in mesh.land_boundaries:
        mask[np.asarray(seg, dtype=np.int64)] = True
    return mask


def element_centroids(mesh: Fort14Mesh) -> np.ndarray:
    return mesh.nodes[mesh.elements].mean(axis=1)


def node_valence(mesh: Fort14Mesh) -> np.ndarray:
    counts = np.zeros(mesh.n_nodes, dtype=np.int64)
    np.add.at(counts, mesh.elements.ravel(), 1)
    return counts


# ---------------------------------------------------------------------------
# Detectors (return a bool flag array shaped (NE,) or (NP,) plus extras)
# ---------------------------------------------------------------------------


def detect_disjoint_components(
    adj: sp.csr_matrix,
) -> tuple[np.ndarray, np.ndarray]:
    """Return (flag, labels). ``flag[i] = True`` if element ``i`` is NOT in
    the largest connected component of the dual graph.
    """
    if adj.shape[0] == 0:
        return np.zeros(0, dtype=bool), np.zeros(0, dtype=np.int64)
    n_comp, labels = connected_components(adj, directed=False, return_labels=True)
    sizes = np.bincount(labels)
    largest = int(sizes.argmax())
    return labels != largest, labels


def detect_dead_end_elements(adj: sp.csr_matrix, mesh: Fort14Mesh) -> np.ndarray:
    """Degree-1 elements that do NOT have an open-boundary edge.

    A degree-1 element has 1 internal edge and 2 boundary edges. Filtering
    out those whose 2 boundary edges include an open-boundary edge skips
    legitimate OB corner elements; what remains is the "spit" / dangling
    triangle terminating a 1-cell channel.
    """
    deg = np.asarray(adj.sum(axis=1)).ravel()
    deg1 = deg == 1
    if not deg1.any():
        return deg1
    ob_codes = open_boundary_edge_codes(mesh)
    if len(ob_codes) == 0:
        return deg1
    elem_codes = element_edge_codes(mesh.elements, mesh.n_nodes)  # (NE, 3)
    flat = elem_codes.ravel()
    idx = np.searchsorted(ob_codes, flat)
    idx_clip = np.minimum(idx, len(ob_codes) - 1)
    has_ob = (ob_codes[idx_clip] == flat).reshape(-1, 3).any(axis=1)
    return deg1 & ~has_ob


def detect_thin_elements(mesh: Fort14Mesh, bdy_mask: np.ndarray) -> np.ndarray:
    """Triangles whose three vertices are all on a boundary (any type)."""
    return bdy_mask[mesh.elements].all(axis=1)


def detect_overconnected_nodes(
    mesh: Fort14Mesh, max_nbr: int = MAX_NBR_ELEM,
) -> tuple[np.ndarray, np.ndarray]:
    """Return (flag, valence). Flag is True where valence > max_nbr."""
    val = node_valence(mesh)
    return val > max_nbr, val


def detect_unreachable_elements(
    adj: sp.csr_matrix, mesh: Fort14Mesh, labels: np.ndarray,
) -> np.ndarray:
    """Elements in dual-graph components that contain no open-boundary node."""
    if not mesh.open_boundaries:
        return np.zeros(mesh.n_elements, dtype=bool)
    ob_node_mask = np.zeros(mesh.n_nodes, dtype=bool)
    for seg in mesh.open_boundaries:
        ob_node_mask[np.asarray(seg, dtype=np.int64)] = True
    elem_has_ob_node = ob_node_mask[mesh.elements].any(axis=1)
    if not elem_has_ob_node.any():
        return np.ones(mesh.n_elements, dtype=bool)
    open_components = np.unique(labels[elem_has_ob_node])
    reachable = np.isin(labels, open_components)
    return ~reachable


# ---------------------------------------------------------------------------
# Reporting
# ---------------------------------------------------------------------------


def _elem_records(
    flag: np.ndarray, centroids: np.ndarray, *, labels: np.ndarray | None = None,
) -> list[dict]:
    idx = np.where(flag)[0]
    out = []
    for i in idx:
        rec = {
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
) -> list[dict]:
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


def build_diag(
    name: str,
    path: Path,
    mesh: Fort14Mesh,
    *,
    disjoint_flag: np.ndarray,
    labels: np.ndarray,
    dead_end_flag: np.ndarray,
    thin_flag: np.ndarray,
    overconn_flag: np.ndarray,
    valence: np.ndarray,
    unreach_flag: np.ndarray,
) -> dict:
    centroids = element_centroids(mesh)
    sizes = np.bincount(labels) if labels.size else np.array([], dtype=np.int64)
    return {
        "mesh": {
            "name": name,
            "path": str(path),
            "n_nodes": int(mesh.n_nodes),
            "n_elements": int(mesh.n_elements),
            "n_open_boundaries": len(mesh.open_boundaries),
            "n_land_boundaries": len(mesh.land_boundaries),
        },
        "detectors": {
            "disjoint_components": {
                "description": (
                    "Elements not in the largest dual-graph connected component."
                ),
                "n_components": int(len(sizes)),
                "component_sizes": [int(s) for s in sizes.tolist()],
                "n_flagged": int(disjoint_flag.sum()),
                "elements": _elem_records(disjoint_flag, centroids, labels=labels),
            },
            "dead_end_elements": {
                "description": (
                    "Degree-1 elements with no open-boundary edge "
                    "(corner-at-OB cases excluded)."
                ),
                "n_flagged": int(dead_end_flag.sum()),
                "elements": _elem_records(dead_end_flag, centroids),
            },
            "thin_elements": {
                "description": (
                    "Elements whose 3 vertices are all on a boundary "
                    "(open or land); chain of these = 1-cell-wide channel."
                ),
                "n_flagged": int(thin_flag.sum()),
                "elements": _elem_records(thin_flag, centroids),
            },
            "overconnected_nodes": {
                "description": (
                    f"Nodes incident to more than {MAX_NBR_ELEM} elements "
                    f"(FVCOM legacy MAX_NBR_ELEM={MAX_NBR_ELEM})."
                ),
                "max_valence_observed": int(valence.max()) if valence.size else 0,
                "n_flagged": int(overconn_flag.sum()),
                "nodes": _node_records(overconn_flag, mesh.nodes, valence=valence),
            },
            "open_boundary_unreachable": {
                "description": (
                    "Elements in components without any open-boundary node."
                ),
                "n_flagged": int(unreach_flag.sum()),
                "elements": _elem_records(unreach_flag, centroids),
            },
        },
    }


# ---------------------------------------------------------------------------
# Plotting
# ---------------------------------------------------------------------------


_OVERLAY_STYLE = [
    ("disjoint", "tab:pink"),
    ("dead-end", "tab:orange"),
    ("thin", "tab:cyan"),
    ("unreachable", "tab:red"),
]


def plot_diagnostics(
    name: str,
    mesh: Fort14Mesh,
    *,
    disjoint_flag: np.ndarray,
    dead_end_flag: np.ndarray,
    thin_flag: np.ndarray,
    overconn_flag: np.ndarray,
    unreach_flag: np.ndarray,
    png: Path,
) -> None:
    fig, ax = plt.subplots(figsize=(11, 9), dpi=120)
    ax.triplot(
        mesh.nodes[:, 0], mesh.nodes[:, 1], mesh.elements,
        color="0.85", lw=0.1,
    )
    # Boundaries
    for seg in mesh.open_boundaries:
        ax.plot(
            mesh.nodes[seg, 0], mesh.nodes[seg, 1], "-",
            color="tab:red", lw=0.9,
        )
    for ib, seg in mesh.land_boundaries:
        col = "tab:blue" if ib == 21 else "0.4"
        lw = 1.2 if ib == 21 else 0.4
        ax.plot(
            mesh.nodes[seg, 0], mesh.nodes[seg, 1], "-",
            color=col, lw=lw,
        )

    flags = {
        "disjoint": disjoint_flag,
        "dead-end": dead_end_flag,
        "thin": thin_flag,
        "unreachable": unreach_flag,
    }
    proxies = []
    for label, color in _OVERLAY_STYLE:
        flag = flags[label]
        if not flag.any():
            continue
        tris = mesh.nodes[mesh.elements[flag]]
        pc = PolyCollection(
            tris, facecolor=color, edgecolor=color, alpha=0.55, linewidths=0.2,
        )
        ax.add_collection(pc)
        proxies.append(plt.Line2D(
            [], [], marker="s", linestyle="",
            markerfacecolor=color, markeredgecolor=color,
            markersize=8, alpha=0.7,
            label=f"{label} ({int(flag.sum())})",
        ))

    if overconn_flag.any():
        ax.scatter(
            mesh.nodes[overconn_flag, 0], mesh.nodes[overconn_flag, 1],
            color="black", marker="x", s=35, linewidths=1.0, zorder=5,
        )
        proxies.append(plt.Line2D(
            [], [], marker="x", linestyle="", color="black",
            markersize=8,
            label=f"overconnected nodes ({int(overconn_flag.sum())})",
        ))

    proxies += [
        plt.Line2D([], [], color="tab:red", lw=1.0, label="open boundary"),
        plt.Line2D([], [], color="0.4", lw=0.6, label="land (ibtype=20)"),
        plt.Line2D([], [], color="tab:blue", lw=1.4, label="river (ibtype=21)"),
    ]

    ax.set_aspect("equal")
    ax.set_xlabel("lon (deg)")
    ax.set_ylabel("lat (deg)")
    ax.set_title(
        f"PoC #24 mesh-check  {name}  "
        f"NP={mesh.n_nodes:,}  NE={mesh.n_elements:,}"
    )
    ax.legend(handles=proxies, loc="lower right", fontsize=8)
    fig.tight_layout()
    fig.savefig(png, bbox_inches="tight", dpi=MESH_PNG_DPI)
    plt.close(fig)


# ---------------------------------------------------------------------------
# Driver
# ---------------------------------------------------------------------------


def analyse_one(name: str, path: Path) -> tuple[str, dict | None]:
    if not path.exists():
        return f"=== {name}: SKIP (missing {path.name}) ===", None
    mesh = read_fort14(path)
    adj = face_face_adjacency(mesh.elements)
    bdy = boundary_node_mask(mesh)

    disjoint_flag, labels = detect_disjoint_components(adj)
    dead_end_flag = detect_dead_end_elements(adj, mesh)
    thin_flag = detect_thin_elements(mesh, bdy)
    overconn_flag, valence = detect_overconnected_nodes(mesh)
    unreach_flag = detect_unreachable_elements(adj, mesh, labels)

    diag = build_diag(
        name, path, mesh,
        disjoint_flag=disjoint_flag, labels=labels,
        dead_end_flag=dead_end_flag, thin_flag=thin_flag,
        overconn_flag=overconn_flag, valence=valence,
        unreach_flag=unreach_flag,
    )

    png = OUT_DIR / f"24_mesh_check_{name}_map.png"
    plot_diagnostics(
        name, mesh,
        disjoint_flag=disjoint_flag, dead_end_flag=dead_end_flag,
        thin_flag=thin_flag, overconn_flag=overconn_flag,
        unreach_flag=unreach_flag,
        png=png,
    )
    diag_path = OUT_DIR / f"24_mesh_check_{name}_diag.json"
    diag_path.write_text(json.dumps(diag, indent=2), encoding="utf-8")

    sizes = diag["detectors"]["disjoint_components"]["component_sizes"]
    text = "\n".join([
        f"=== {name} ===",
        f"file: {path}",
        f"NP={mesh.n_nodes:,}  NE={mesh.n_elements:,}  "
        f"open={len(mesh.open_boundaries)}  land={len(mesh.land_boundaries)}",
        f"  1. disjoint comp:    n_components={len(sizes)}  "
        f"sizes={sizes if len(sizes) <= 6 else sizes[:6] + ['...']}  "
        f"flagged_elems={diag['detectors']['disjoint_components']['n_flagged']:,}",
        f"  2. dead-end elems:   "
        f"{diag['detectors']['dead_end_elements']['n_flagged']:,}",
        f"  3. thin elems:       "
        f"{diag['detectors']['thin_elements']['n_flagged']:,}",
        f"  4. over-conn nodes:  "
        f"{diag['detectors']['overconnected_nodes']['n_flagged']:,}  "
        f"max_valence="
        f"{diag['detectors']['overconnected_nodes']['max_valence_observed']}",
        f"  5. unreachable:      "
        f"{diag['detectors']['open_boundary_unreachable']['n_flagged']:,}",
        f"  wrote: {png.name}",
        f"  wrote: {diag_path.name}",
    ])
    print(text)
    return text, diag


def main() -> None:
    OUT_DIR.mkdir(parents=True, exist_ok=True)
    sections: list[str] = []
    for name, path in INPUTS:
        text, _ = analyse_one(name, path)
        sections.append(text)
    summary = "\n\n".join(sections) + "\n"
    SUMMARY_TXT.write_text(summary, encoding="utf-8")
    print(f"\n[24] wrote {SUMMARY_TXT}")


if __name__ == "__main__":
    main()
