"""PoC #4: custom open-boundary edge-perpendicularity fix.

Moves the first-ring interior neighbours of every open-boundary node so
that each incident interior edge becomes perpendicular to the local
boundary tangent. Boundary nodes (both open and land) are kept fixed,
which preserves the boundary shape and the coastline. Edge length is
preserved when projecting an interior node, so triangle area shrinks/
grows only as required by the rotation.

Outputs:
    outputs/04_perp_fix_summary.txt
    outputs/04_perp_before_after_hist.png
    outputs/04_open_bdy_first_ring.png
    outputs/04_tb_futtsu20220311_perp_fixed.14
"""

from __future__ import annotations

from pathlib import Path

import matplotlib

matplotlib.use("Agg")

import matplotlib.pyplot as plt  # noqa: E402
import numpy as np  # noqa: E402

from fvcom_mesh_tools.io import Fort14Mesh, read_fort14, write_fort14  # noqa: E402

REPO_ROOT = Path(__file__).resolve().parent.parent
MESH_PATH = REPO_ROOT / "data" / "mesh" / "reference" / "tokyo_bay" / "tb_futtsu20220311.14"
OUT_DIR = REPO_ROOT / "outputs"
SUMMARY_TXT = OUT_DIR / "04_perp_fix_summary.txt"
PERP_HIST_BA = OUT_DIR / "04_perp_before_after_hist.png"
FIRST_RING_PNG = OUT_DIR / "04_open_bdy_first_ring.png"
FIXED_F14 = OUT_DIR / "04_tb_futtsu20220311_perp_fixed.14"

ALPHA = 1.0          # target blend (1.0 = full projection)
N_ITERS = 1          # iterations (set >1 with alpha<1 for damped relaxation)
SEGMENT_INDEX = 0    # which open-boundary segment to operate on


def unique_edges(elements: np.ndarray) -> np.ndarray:
    e = np.vstack([
        elements[:, [0, 1]], elements[:, [1, 2]], elements[:, [2, 0]],
    ])
    e.sort(axis=1)
    return np.unique(e, axis=0)


def boundary_tangents(bdy_xy: np.ndarray) -> np.ndarray:
    """Per-node unit tangent: central difference; one-sided at the ends."""
    tangents = np.empty_like(bdy_xy)
    tangents[1:-1] = bdy_xy[2:] - bdy_xy[:-2]
    tangents[0] = bdy_xy[1] - bdy_xy[0]
    tangents[-1] = bdy_xy[-1] - bdy_xy[-2]
    norms = np.linalg.norm(tangents, axis=1, keepdims=True)
    return tangents / np.where(norms == 0, 1.0, norms)


def fixed_node_mask(mesh: Fort14Mesh) -> np.ndarray:
    """Boolean mask of nodes that must stay put: every open and land boundary node."""
    fixed = np.zeros(mesh.n_nodes, dtype=bool)
    for ids in mesh.open_boundaries:
        fixed[ids] = True
    for _ibtype, ids in mesh.land_boundaries:
        fixed[ids] = True
    return fixed


def open_bdy_perp(mesh: Fort14Mesh, segment_index: int = 0) -> np.ndarray:
    """Per-edge deviation from 90 deg between interior edges incident on the
    selected open-boundary segment and the local boundary tangent."""
    if not mesh.open_boundaries:
        return np.array([])
    bdy = np.asarray(mesh.open_boundaries[segment_index], dtype=np.int64)
    bdy_xy = mesh.nodes[bdy]
    tangents = boundary_tangents(bdy_xy)

    inv_map = np.full(mesh.n_nodes, -1, dtype=np.int64)
    inv_map[bdy] = np.arange(len(bdy))

    edges = unique_edges(mesh.elements)
    a_in = inv_map[edges[:, 0]] >= 0
    b_in = inv_map[edges[:, 1]] >= 0
    incident = a_in ^ b_in
    inc = edges[incident]
    inc_a_in = a_in[incident]
    bdy_node = np.where(inc_a_in, inc[:, 0], inc[:, 1])
    int_node = np.where(inc_a_in, inc[:, 1], inc[:, 0])

    edge_vec = mesh.nodes[int_node] - mesh.nodes[bdy_node]
    edge_norms = np.linalg.norm(edge_vec, axis=1, keepdims=True)
    edge_vec = edge_vec / np.where(edge_norms == 0, 1.0, edge_norms)
    edge_tangent = tangents[inv_map[bdy_node]]
    cos_angle = np.clip((edge_vec * edge_tangent).sum(axis=1), -1.0, 1.0)
    angle = np.degrees(np.arccos(np.abs(cos_angle)))
    return 90.0 - angle


def signed_areas(mesh: Fort14Mesh) -> np.ndarray:
    p0 = mesh.nodes[mesh.elements[:, 0]]
    p1 = mesh.nodes[mesh.elements[:, 1]]
    p2 = mesh.nodes[mesh.elements[:, 2]]
    return 0.5 * (
        (p1[:, 0] - p0[:, 0]) * (p2[:, 1] - p0[:, 1])
        - (p1[:, 1] - p0[:, 1]) * (p2[:, 0] - p0[:, 0])
    )


def perp_fix_step(
    nodes: np.ndarray,
    elements: np.ndarray,
    bdy: np.ndarray,
    fixed: np.ndarray,
    alpha: float,
) -> np.ndarray:
    """One projection step. Returns a new node array."""
    bdy_xy = nodes[bdy]
    tangents = boundary_tangents(bdy_xy)
    perp = np.column_stack([-tangents[:, 1], tangents[:, 0]])  # 90 deg CCW

    inv_map = np.full(len(nodes), -1, dtype=np.int64)
    inv_map[bdy] = np.arange(len(bdy))

    edges = unique_edges(elements)
    a_in = inv_map[edges[:, 0]] >= 0
    b_in = inv_map[edges[:, 1]] >= 0
    incident = a_in ^ b_in
    inc = edges[incident]
    inc_a_in = a_in[incident]
    bdy_node = np.where(inc_a_in, inc[:, 0], inc[:, 1])
    int_node = np.where(inc_a_in, inc[:, 1], inc[:, 0])

    bdy_pos = nodes[bdy_node]
    int_pos = nodes[int_node]
    edge_vec = int_pos - bdy_pos
    edge_len = np.linalg.norm(edge_vec, axis=1, keepdims=True)

    perp_at_bdy = perp[inv_map[bdy_node]]
    sign = np.sign((edge_vec * perp_at_bdy).sum(axis=1, keepdims=True))
    sign = np.where(sign == 0, 1.0, sign)

    target_pos = bdy_pos + edge_len * (perp_at_bdy * sign)

    accum = np.zeros_like(nodes)
    cnt = np.zeros(len(nodes), dtype=np.int64)
    np.add.at(accum, int_node, target_pos)
    np.add.at(cnt, int_node, 1)

    new_nodes = nodes.copy()
    moved_idx = np.where(cnt > 0)[0]
    if moved_idx.size:
        avg_target = accum[moved_idx] / cnt[moved_idx, None]
        new_nodes[moved_idx] = (1.0 - alpha) * nodes[moved_idx] + alpha * avg_target

    new_nodes[fixed] = nodes[fixed]   # boundary nodes never move
    return new_nodes


def fix_open_boundary_perpendicular(
    mesh: Fort14Mesh,
    alpha: float = 1.0,
    n_iters: int = 1,
    segment_index: int = 0,
) -> tuple[Fort14Mesh, dict]:
    bdy = np.asarray(mesh.open_boundaries[segment_index], dtype=np.int64)
    fixed = fixed_node_mask(mesh)

    nodes = mesh.nodes.copy()
    for _ in range(n_iters):
        nodes = perp_fix_step(nodes, mesh.elements, bdy, fixed, alpha)

    out = Fort14Mesh(
        title=mesh.title + " (boundary-perp fixed)",
        nodes=nodes,
        depths=mesh.depths.copy(),
        elements=mesh.elements.copy(),
        open_boundaries=[a.copy() for a in mesh.open_boundaries],
        land_boundaries=[(ib, a.copy()) for (ib, a) in mesh.land_boundaries],
    )
    info = {
        "alpha": alpha,
        "n_iters": n_iters,
        "n_first_ring": int(((mesh.nodes != nodes).any(axis=1)).sum()),
    }
    return out, info


def _stats(name: str, x: np.ndarray) -> str:
    if x.size == 0:
        return f"{name}: empty"
    return (
        f"{name}: n={x.size:,}  min={np.nanmin(x):.4f}  "
        f"p50={np.nanpercentile(x, 50):.4f}  p95={np.nanpercentile(x, 95):.4f}  "
        f"max={np.nanmax(x):.4f}  mean={np.nanmean(x):.4f}"
    )


def plot_hist(perp_b: np.ndarray, perp_a: np.ndarray, png: Path) -> None:
    fig, ax = plt.subplots(figsize=(7, 4), dpi=120)
    bins = np.linspace(0, 90, 60)
    ax.hist(perp_b, bins=bins, alpha=0.6, label="before", color="tab:red")
    ax.hist(perp_a, bins=bins, alpha=0.6, label="after", color="tab:blue")
    ax.set_xlabel("perp deviation from 90 deg (deg)")
    ax.set_ylabel("count")
    ax.set_title(
        "FVCOM open-boundary edge perpendicularity: before vs custom fix"
    )
    ax.legend()
    ax.grid(True, lw=0.3, color="0.85")
    fig.tight_layout()
    fig.savefig(png, bbox_inches="tight")
    plt.close(fig)


def plot_first_ring(
    before: Fort14Mesh,
    after: Fort14Mesh,
    segment_index: int,
    png: Path,
) -> None:
    """Zoomed plot of the open-boundary segment with first-ring nodes
    before/after the projection."""
    bdy = before.open_boundaries[segment_index]

    inv_map = np.full(before.n_nodes, -1, dtype=np.int64)
    inv_map[bdy] = np.arange(len(bdy))
    edges = unique_edges(before.elements)
    a_in = inv_map[edges[:, 0]] >= 0
    b_in = inv_map[edges[:, 1]] >= 0
    incident = a_in ^ b_in
    inc = edges[incident]
    inc_a_in = a_in[incident]
    bdy_node = np.where(inc_a_in, inc[:, 0], inc[:, 1])
    int_node = np.where(inc_a_in, inc[:, 1], inc[:, 0])

    fig, ax = plt.subplots(figsize=(8, 7), dpi=120)
    ax.scatter(
        before.nodes[:, 0], before.nodes[:, 1],
        s=0.2, color="0.92", linewidths=0,
    )
    ax.plot(
        before.nodes[bdy, 0], before.nodes[bdy, 1], "-",
        color="tab:gray", lw=1.2, label="open boundary",
    )
    # Edges before/after
    for b, e in zip(bdy_node, int_node):
        ax.plot(
            [before.nodes[b, 0], before.nodes[e, 0]],
            [before.nodes[b, 1], before.nodes[e, 1]],
            "-", color="tab:red", lw=0.4, alpha=0.7,
        )
        ax.plot(
            [after.nodes[b, 0], after.nodes[e, 0]],
            [after.nodes[b, 1], after.nodes[e, 1]],
            "-", color="tab:blue", lw=0.4, alpha=0.7,
        )

    ax.set_aspect("equal")
    ax.set_xlabel("lon (deg)")
    ax.set_ylabel("lat (deg)")
    ax.set_title("First-ring incident edges: before (red) vs after (blue)")
    proxy_before = plt.Line2D([], [], color="tab:red", lw=1.2, label="before")
    proxy_after = plt.Line2D([], [], color="tab:blue", lw=1.2, label="after")
    proxy_bdy = plt.Line2D([], [], color="tab:gray", lw=1.2, label="open boundary")
    ax.legend(handles=[proxy_before, proxy_after, proxy_bdy], loc="best", fontsize=8)
    ax.grid(True, lw=0.3, color="0.9")
    fig.tight_layout()
    fig.savefig(png, bbox_inches="tight")
    plt.close(fig)


def main() -> None:
    if not MESH_PATH.exists():
        raise SystemExit(f"mesh not found: {MESH_PATH}")
    OUT_DIR.mkdir(parents=True, exist_ok=True)

    print(f"[04] reading {MESH_PATH}")
    before = read_fort14(MESH_PATH)
    print(f"[04] NP={before.n_nodes:,}  NE={before.n_elements:,}")

    perp_b = open_bdy_perp(before, segment_index=SEGMENT_INDEX)
    areas_b = signed_areas(before)
    n_flipped_b = int((areas_b <= 0).sum())

    print(f"[04] running perpendicularity fix: alpha={ALPHA}, n_iters={N_ITERS}")
    after, info = fix_open_boundary_perpendicular(
        before, alpha=ALPHA, n_iters=N_ITERS, segment_index=SEGMENT_INDEX,
    )
    print(f"[04] moved {info['n_first_ring']:,} interior nodes")

    perp_a = open_bdy_perp(after, segment_index=SEGMENT_INDEX)
    areas_a = signed_areas(after)
    n_flipped_a = int((areas_a <= 0).sum())

    disp = np.linalg.norm(after.nodes - before.nodes, axis=1)

    lines = [
        f"file:  {MESH_PATH}",
        f"NP={before.n_nodes:,}  NE={before.n_elements:,}",
        f"alpha={ALPHA}  n_iters={N_ITERS}  segment_index={SEGMENT_INDEX}",
        f"interior nodes moved: {info['n_first_ring']:,}",
        "",
        "[Per-node displacement (deg)]",
        _stats("  ||after - before||", disp),
        "",
        "[Element validity (signed area <= 0 indicates a flip)]",
        f"  before flipped triangles: {n_flipped_b:,} / {before.n_elements:,}",
        f"  after  flipped triangles: {n_flipped_a:,} / {after.n_elements:,}",
        _stats("  signed area before", areas_b),
        _stats("  signed area after ", areas_a),
        "",
        "[FVCOM open-boundary perpendicularity (deg from 90)]",
        _stats("  before per-edge", perp_b),
        _stats("  after  per-edge", perp_a),
    ]
    summary = "\n".join(lines)
    print(summary)
    SUMMARY_TXT.write_text(summary + "\n", encoding="utf-8")
    print(f"[04] wrote {SUMMARY_TXT}")

    write_fort14(after, FIXED_F14)
    print(f"[04] wrote {FIXED_F14}")

    plot_hist(perp_b, perp_a, PERP_HIST_BA)
    print(f"[04] wrote {PERP_HIST_BA}")

    plot_first_ring(before, after, SEGMENT_INDEX, FIRST_RING_PNG)
    print(f"[04] wrote {FIRST_RING_PNG}")


if __name__ == "__main__":
    main()
