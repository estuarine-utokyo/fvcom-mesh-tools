"""PoC #2: orthogonality / open-boundary edge perpendicularity.

Two metrics on tb_futtsu20220311.14:

1. MeshKernel `mesh2d_get_orthogonality` — Deltares' general per-edge
   orthogonality (ratio between the edge and the segment connecting the
   adjacent face circumcenters; 0 == perfectly orthogonal).

2. FVCOM-relevant open-boundary perpendicularity — for every interior
   edge incident on a node of the open boundary, the angle between the
   edge and the local boundary tangent. Deviation from 90 degrees is
   reported (0 deg == perpendicular).

Outputs:
    outputs/02_orthogonality_summary.txt
    outputs/02_mk_orthogonality_hist.png
    outputs/02_open_bdy_perpendicularity_hist.png
    outputs/02_open_bdy_perpendicularity_map.png
"""

from __future__ import annotations

from pathlib import Path

import matplotlib

matplotlib.use("Agg")

import matplotlib.pyplot as plt  # noqa: E402
import numpy as np  # noqa: E402
from meshkernel import Mesh2d, MeshKernel, ProjectionType  # noqa: E402

from fvcom_mesh_tools.io import Fort14Mesh, read_fort14  # noqa: E402
from fvcom_mesh_tools.plotting import MESH_PNG_DPI  # noqa: E402

# MeshKernel returns this sentinel for edges where the orthogonality cannot
# be computed (e.g. at the mesh boundary, where the second face is absent).
MK_MISSING = -999.0

REPO_ROOT = Path(__file__).resolve().parent.parent
MESH_PATH = REPO_ROOT / "data" / "mesh" / "reference" / "tokyo_bay" / "tb_futtsu20220311.14"
OUT_DIR = REPO_ROOT / "outputs"
SUMMARY_TXT = OUT_DIR / "02_orthogonality_summary.txt"
MK_HIST = OUT_DIR / "02_mk_orthogonality_hist.png"
PERP_HIST = OUT_DIR / "02_open_bdy_perpendicularity_hist.png"
PERP_MAP = OUT_DIR / "02_open_bdy_perpendicularity_map.png"


def unique_edges(elements: np.ndarray) -> np.ndarray:
    """Return unique mesh edges as a ``(M, 2)`` array of sorted node indices."""
    e = np.vstack([
        elements[:, [0, 1]],
        elements[:, [1, 2]],
        elements[:, [2, 0]],
    ])
    e.sort(axis=1)
    return np.unique(e, axis=0)


def mk_orthogonality(mesh: Fort14Mesh) -> tuple[np.ndarray, np.ndarray]:
    """Run MeshKernel's per-edge orthogonality on the whole mesh.

    Returns ``(values, valid_mask)`` where ``valid_mask`` is False at edges
    where MeshKernel could not compute the metric (sentinel ``MK_MISSING``).

    ``edge_nodes`` is computed and passed explicitly: MeshKernel 8.2.2's
    ``mesh2d_set`` segfaults if it has to derive the edge list itself.
    """
    edges = unique_edges(mesh.elements)

    m2d = Mesh2d(
        node_x=np.ascontiguousarray(mesh.nodes[:, 0], dtype=np.float64),
        node_y=np.ascontiguousarray(mesh.nodes[:, 1], dtype=np.float64),
        edge_nodes=edges.flatten().astype(np.int32),
        face_nodes=mesh.elements.flatten().astype(np.int32),
        nodes_per_face=np.full(mesh.n_elements, 3, dtype=np.int32),
    )
    mk = MeshKernel(projection=ProjectionType.SPHERICAL)
    mk.mesh2d_set(m2d)
    geom = mk.mesh2d_get_orthogonality()
    values = np.asarray(geom.values, dtype=np.float64)
    valid = values != MK_MISSING
    return values, valid


def open_bdy_perpendicularity(mesh: Fort14Mesh) -> dict:
    """For every interior edge incident on an open-boundary node, compute
    the deviation (deg) from a 90-deg angle relative to the boundary tangent.
    """
    if not mesh.open_boundaries:
        return {"perp_dev": np.array([]), "node_dev_mean": np.array([])}

    bdy = np.asarray(mesh.open_boundaries[0], dtype=np.int64)
    coords = mesh.nodes

    # Per-node tangent: central difference along the boundary, one-sided at ends.
    bdy_xy = coords[bdy]
    tangents = np.empty_like(bdy_xy)
    tangents[1:-1] = bdy_xy[2:] - bdy_xy[:-2]
    tangents[0] = bdy_xy[1] - bdy_xy[0]
    tangents[-1] = bdy_xy[-1] - bdy_xy[-2]
    norms = np.linalg.norm(tangents, axis=1, keepdims=True)
    tangents = tangents / np.where(norms == 0, 1.0, norms)

    # Inverse map node-id -> position along boundary; -1 if not on boundary.
    inv_map = np.full(mesh.n_nodes, -1, dtype=np.int64)
    inv_map[bdy] = np.arange(len(bdy))

    edges = unique_edges(mesh.elements)
    a_in = inv_map[edges[:, 0]] >= 0
    b_in = inv_map[edges[:, 1]] >= 0
    incident_mask = a_in ^ b_in
    inc = edges[incident_mask]
    inc_a_in = a_in[incident_mask]

    bdy_node = np.where(inc_a_in, inc[:, 0], inc[:, 1])
    int_node = np.where(inc_a_in, inc[:, 1], inc[:, 0])

    edge_vec = coords[int_node] - coords[bdy_node]
    edge_norms = np.linalg.norm(edge_vec, axis=1, keepdims=True)
    edge_vec = edge_vec / np.where(edge_norms == 0, 1.0, edge_norms)

    edge_tangent = tangents[inv_map[bdy_node]]
    cos_angle = np.clip((edge_vec * edge_tangent).sum(axis=1), -1.0, 1.0)
    angle_deg = np.degrees(np.arccos(np.abs(cos_angle)))   # 0..90
    perp_dev = 90.0 - angle_deg                              # 0 = perfectly perpendicular

    # Average deviation per boundary node (only nodes with incident interior edges).
    n_bdy = len(bdy)
    sum_dev = np.zeros(n_bdy, dtype=np.float64)
    cnt = np.zeros(n_bdy, dtype=np.int64)
    bdy_pos = inv_map[bdy_node]
    np.add.at(sum_dev, bdy_pos, perp_dev)
    np.add.at(cnt, bdy_pos, 1)
    node_dev_mean = np.where(cnt > 0, sum_dev / np.where(cnt == 0, 1, cnt), np.nan)

    return {
        "n_open_bdy_nodes": n_bdy,
        "n_incident_edges": int(incident_mask.sum()),
        "perp_dev": perp_dev,
        "node_dev_mean": node_dev_mean,
        "bdy_xy": bdy_xy,
    }


def _stats(name: str, x: np.ndarray) -> str:
    if x.size == 0:
        return f"{name}: empty"
    return (
        f"{name}: n={x.size:,}  "
        f"min={np.nanmin(x):.4f}  "
        f"p50={np.nanpercentile(x, 50):.4f}  "
        f"p95={np.nanpercentile(x, 95):.4f}  "
        f"max={np.nanmax(x):.4f}  "
        f"mean={np.nanmean(x):.4f}"
    )


def write_summary(
    mesh: Fort14Mesh,
    mk_ortho: np.ndarray,
    mk_valid: np.ndarray,
    perp: dict,
) -> str:
    lines = [
        f"file:                       {MESH_PATH}",
        f"NP={mesh.n_nodes:,}   NE={mesh.n_elements:,}",
        "",
        "[MeshKernel mesh2d_get_orthogonality (per-edge, 0 = orthogonal)]",
        f"  total edges: {len(mk_ortho):,}",
        f"  valid:       {int(mk_valid.sum()):,}   "
        f"missing (sentinel -999): {int((~mk_valid).sum()):,}",
        _stats("  values", mk_ortho[mk_valid]),
        "",
        f"[FVCOM open-boundary perpendicularity (deg from 90 deg, "
        f"{perp.get('n_incident_edges', 0):,} edges across "
        f"{perp.get('n_open_bdy_nodes', 0):,} bdy nodes)]",
        _stats("  per-edge dev (deg)", perp["perp_dev"]),
        _stats("  per-node mean dev", perp["node_dev_mean"]),
    ]
    return "\n".join(lines)


def plot_hist(values: np.ndarray, title: str, xlabel: str, png: Path) -> None:
    fig, ax = plt.subplots(figsize=(7, 4), dpi=120)
    finite = values[np.isfinite(values)]
    ax.hist(finite, bins=80, color="tab:blue", alpha=0.85)
    ax.set_title(title)
    ax.set_xlabel(xlabel)
    ax.set_ylabel("count")
    ax.grid(True, lw=0.3, color="0.85")
    fig.tight_layout()
    fig.savefig(png, bbox_inches="tight")
    plt.close(fig)


def plot_perp_map(mesh: Fort14Mesh, perp: dict, png: Path) -> None:
    bdy_xy = perp["bdy_xy"]
    node_dev = perp["node_dev_mean"]

    fig, ax = plt.subplots(figsize=(8, 7), dpi=120)
    ax.scatter(
        mesh.nodes[:, 0], mesh.nodes[:, 1],
        s=0.2, color="0.85", linewidths=0,
    )
    sc = ax.scatter(
        bdy_xy[:, 0], bdy_xy[:, 1],
        c=node_dev, cmap="magma", s=12, edgecolors="none",
        vmin=0, vmax=max(30.0, float(np.nanpercentile(node_dev, 95))),
    )
    cb = fig.colorbar(sc, ax=ax, shrink=0.7)
    cb.set_label("mean perpendicularity dev (deg, 0 = perpendicular)")

    ax.set_aspect("equal")
    ax.set_xlabel("lon (deg)")
    ax.set_ylabel("lat (deg)")
    ax.set_title("open boundary nodes colored by mean perp deviation")
    ax.grid(True, lw=0.3, color="0.9")
    fig.tight_layout()
    fig.savefig(png, bbox_inches="tight", dpi=MESH_PNG_DPI)
    plt.close(fig)


def main() -> None:
    if not MESH_PATH.exists():
        raise SystemExit(f"mesh not found: {MESH_PATH}")
    OUT_DIR.mkdir(parents=True, exist_ok=True)

    print(f"[02] reading {MESH_PATH}")
    mesh = read_fort14(MESH_PATH)

    print("[02] running MeshKernel orthogonality...")
    mk_ortho, mk_valid = mk_orthogonality(mesh)

    print("[02] computing open-boundary perpendicularity...")
    perp = open_bdy_perpendicularity(mesh)

    text = write_summary(mesh, mk_ortho, mk_valid, perp)
    print(text)
    SUMMARY_TXT.write_text(text + "\n", encoding="utf-8")
    print(f"[02] wrote {SUMMARY_TXT}")

    plot_hist(
        mk_ortho[mk_valid],
        "MeshKernel orthogonality (valid edges, 0 = orthogonal)",
        "orthogonality value",
        MK_HIST,
    )
    print(f"[02] wrote {MK_HIST}")

    plot_hist(
        perp["perp_dev"],
        "Open-boundary edge perpendicularity (per incident interior edge)",
        "deviation from 90 deg (deg)",
        PERP_HIST,
    )
    print(f"[02] wrote {PERP_HIST}")

    plot_perp_map(mesh, perp, PERP_MAP)
    print(f"[02] wrote {PERP_MAP}")


if __name__ == "__main__":
    main()
