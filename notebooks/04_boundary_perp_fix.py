"""PoC #4: custom open-boundary edge-perpendicularity fix.

Thin notebook driver around
:func:`fvcom_mesh_tools.algorithms.align_open_boundary_first_ring`. The
algorithm moves the first-ring interior neighbours of every open-boundary
node so each incident edge becomes as perpendicular to the local boundary
tangent as possible while preserving original edge length. Boundary nodes
(open + land) are kept fixed.

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

from fvcom_mesh_tools.algorithms import (  # noqa: E402
    align_open_boundary_first_ring,
    open_bdy_perpendicularity,
    signed_areas,
    unique_edges,
)
from fvcom_mesh_tools.io import Fort14Mesh, read_fort14, write_fort14  # noqa: E402
from fvcom_mesh_tools.plotting import MESH_PNG_DPI  # noqa: E402

REPO_ROOT = Path(__file__).resolve().parent.parent
MESH_PATH = REPO_ROOT / "data" / "mesh" / "reference" / "tokyo_bay" / "tb_futtsu20220311.14"
OUT_DIR = REPO_ROOT / "outputs"
SUMMARY_TXT = OUT_DIR / "04_perp_fix_summary.txt"
PERP_HIST_BA = OUT_DIR / "04_perp_before_after_hist.png"
FIRST_RING_PNG = OUT_DIR / "04_open_bdy_first_ring.png"
FIXED_F14 = OUT_DIR / "04_tb_futtsu20220311_perp_fixed.14"

ALPHA = 1.0
N_ITERS = 1
SMOOTH_ITERS = 0
SEGMENT_INDEX = 0


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
        color="tab:gray", lw=1.2,
    )
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
    fig.savefig(png, bbox_inches="tight", dpi=MESH_PNG_DPI)
    plt.close(fig)


def main() -> None:
    if not MESH_PATH.exists():
        raise SystemExit(f"mesh not found: {MESH_PATH}")
    OUT_DIR.mkdir(parents=True, exist_ok=True)

    print(f"[04] reading {MESH_PATH}")
    before = read_fort14(MESH_PATH)
    print(f"[04] NP={before.n_nodes:,}  NE={before.n_elements:,}")

    perp_b = open_bdy_perpendicularity(before, segment_index=SEGMENT_INDEX)
    n_flipped_b = int((signed_areas(before) <= 0).sum())

    print(
        f"[04] running align_open_boundary_first_ring: "
        f"alpha={ALPHA}, n_iters={N_ITERS}, smooth_iters={SMOOTH_ITERS}"
    )
    after, info = align_open_boundary_first_ring(
        before,
        alpha=ALPHA,
        n_iters=N_ITERS,
        smooth_iters=SMOOTH_ITERS,
        segment_index=SEGMENT_INDEX,
    )
    print(f"[04] moved {info['moved']:,} interior nodes")
    print(f"[04] movable first-ring by parent count: {info['movable_first_ring_by_parent_count']}")

    perp_a = open_bdy_perpendicularity(after, segment_index=SEGMENT_INDEX)
    n_flipped_a = int((signed_areas(after) <= 0).sum())
    disp = np.linalg.norm(after.nodes - before.nodes, axis=1)

    lines = [
        f"file:  {MESH_PATH}",
        f"NP={before.n_nodes:,}  NE={before.n_elements:,}",
        f"alpha={ALPHA}  n_iters={N_ITERS}  smooth_iters={SMOOTH_ITERS}  "
        f"segment_index={SEGMENT_INDEX}",
        f"interior nodes moved: {info['moved']:,}",
        f"movable first-ring by parent count: {info['movable_first_ring_by_parent_count']}",
        "",
        "[Per-node displacement (deg)]",
        _stats("  ||after - before||", disp),
        "",
        "[Element validity]",
        f"  flipped triangles before/after: {n_flipped_b:,} / {n_flipped_a:,}",
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
