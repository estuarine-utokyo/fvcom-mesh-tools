"""PoC #9: edge-swap quality pass on the fmesh-buildmesh Tokyo Bay output.

PoC #8 showed that damped Laplacian smoothing alone barely shifts mean
quality on this mesh - slivers are *topological*, not metric. This PoC
runs edge swapping (Lawson / min-angle flip) until convergence and
re-measures. Per-pass swap counts are printed so the convergence
behaviour is auditable from the log.

Outputs:
    outputs/09_tokyo_bay_swapped.14
    outputs/09_edge_swap_summary.txt
    outputs/09_edge_swap_alpha_hist.png
    outputs/09_edge_swap_min_angle_hist.png
    outputs/09_edge_swap_mesh.png
"""

from __future__ import annotations

import time
from pathlib import Path

import matplotlib

matplotlib.use("Agg")

import matplotlib.pyplot as plt  # noqa: E402
import numpy as np  # noqa: E402

from fvcom_mesh_tools.algorithms import (  # noqa: E402
    alpha_quality,
    min_interior_angle,
    open_bdy_perpendicularity,
    signed_areas,
    swap_edges_for_quality,
)
from fvcom_mesh_tools.io import read_fort14, write_fort14  # noqa: E402
from fvcom_mesh_tools.plotting import MESH_PNG_DPI  # noqa: E402

REPO_ROOT = Path(__file__).resolve().parent.parent
INPUT_F14 = REPO_ROOT / "outputs" / "07_tokyo_bay_buildmesh.14"
OUT_DIR = REPO_ROOT / "outputs"
SWAPPED_F14 = OUT_DIR / "09_tokyo_bay_swapped.14"
SUMMARY_TXT = OUT_DIR / "09_edge_swap_summary.txt"
ALPHA_PNG = OUT_DIR / "09_edge_swap_alpha_hist.png"
MIN_ANGLE_PNG = OUT_DIR / "09_edge_swap_min_angle_hist.png"
MESH_PNG = OUT_DIR / "09_edge_swap_mesh.png"

MAX_ITERS = 30


def _stats(name: str, x: np.ndarray, fmt: str = ".4f") -> str:
    if x.size == 0:
        return f"{name}: empty"
    return (
        f"{name}: n={x.size:,}  "
        f"min={np.nanmin(x):{fmt}}  "
        f"p05={np.nanpercentile(x, 5):{fmt}}  "
        f"p50={np.nanpercentile(x, 50):{fmt}}  "
        f"p95={np.nanpercentile(x, 95):{fmt}}  "
        f"max={np.nanmax(x):{fmt}}  "
        f"mean={np.nanmean(x):{fmt}}"
    )


def plot_alpha_hist(b: np.ndarray, a: np.ndarray, png: Path) -> None:
    fig, ax = plt.subplots(figsize=(8, 4.5), dpi=120)
    bins = np.linspace(0.0, 1.0, 41)
    ax.hist(b, bins=bins, alpha=0.6, label="before swap", color="tab:red")
    ax.hist(a, bins=bins, alpha=0.6, label="after swap", color="tab:blue")
    ax.set_xlabel("alpha-quality (1 = equilateral)")
    ax.set_ylabel("count")
    ax.set_title(f"Edge swap: max_iters={MAX_ITERS}")
    ax.legend()
    ax.grid(True, lw=0.3, color="0.9")
    fig.tight_layout()
    fig.savefig(png, bbox_inches="tight")
    plt.close(fig)


def plot_min_angle_hist(b: np.ndarray, a: np.ndarray, png: Path) -> None:
    fig, ax = plt.subplots(figsize=(8, 4.5), dpi=120)
    bins = np.linspace(0.0, 60.0, 61)
    ax.hist(b, bins=bins, alpha=0.6, label="before swap", color="tab:red")
    ax.hist(a, bins=bins, alpha=0.6, label="after swap", color="tab:blue")
    ax.axvline(20.0, color="0.4", ls="--", lw=0.7, label="FVCOM 20-deg threshold")
    ax.set_xlabel("triangle minimum interior angle (deg)")
    ax.set_ylabel("count")
    ax.set_title("Min-interior-angle distribution after edge swap")
    ax.legend()
    ax.grid(True, lw=0.3, color="0.9")
    fig.tight_layout()
    fig.savefig(png, bbox_inches="tight")
    plt.close(fig)


def plot_mesh(nodes: np.ndarray, elements: np.ndarray, depths: np.ndarray, png: Path) -> None:
    fig, ax = plt.subplots(figsize=(8, 7), dpi=120)
    tri = ax.tripcolor(
        nodes[:, 0], nodes[:, 1], elements,
        facecolors=depths[elements].mean(axis=1),
        cmap="viridis_r", edgecolors="0.85", linewidth=0.05,
    )
    fig.colorbar(tri, ax=ax, label="depth (m, +down)")
    ax.set_aspect("equal")
    ax.set_xlabel("lon (deg)")
    ax.set_ylabel("lat (deg)")
    ax.set_title(f"Tokyo Bay after edge-swap pass (max_iters={MAX_ITERS})")
    fig.tight_layout()
    fig.savefig(png, bbox_inches="tight", dpi=MESH_PNG_DPI)
    plt.close(fig)


def main() -> None:
    if not INPUT_F14.exists():
        raise SystemExit(
            f"input not found: {INPUT_F14}\n"
            f"Run notebooks/07_buildmesh_e2e.py first."
        )
    OUT_DIR.mkdir(parents=True, exist_ok=True)

    print(f"[09] reading {INPUT_F14}")
    before = read_fort14(INPUT_F14)
    print(f"[09] NP={before.n_nodes:,}  NE={before.n_elements:,}")

    q_b = alpha_quality(before)
    a_b = min_interior_angle(before)
    perp_b = open_bdy_perpendicularity(before, segment_index=0)

    t0 = time.perf_counter()
    after, info = swap_edges_for_quality(before, max_iters=MAX_ITERS)
    t_swap = time.perf_counter() - t0
    print(f"[09] swap_edges_for_quality: {t_swap:.2f} s")
    print(f"[09] swaps per iteration: {info['swaps_per_iter']}")
    print(f"[09] total swaps: {info['total_swaps']:,}")

    q_a = alpha_quality(after)
    a_a = min_interior_angle(after)
    perp_a = open_bdy_perpendicularity(after, segment_index=0)
    flipped_a = int((signed_areas(after) <= 0).sum())

    write_fort14(after, SWAPPED_F14)

    lines = [
        f"input:  {INPUT_F14}",
        f"output: {SWAPPED_F14}",
        f"max_iters={MAX_ITERS}",
        f"edge-swap wall time: {t_swap:.2f} s",
        f"swaps per iter: {info['swaps_per_iter']}",
        f"total swaps:    {info['total_swaps']:,}",
        "",
        "[Element validity]",
        f"  flipped triangles after: {flipped_a:,}",
        "",
        "[Triangle alpha-quality (1 = equilateral)]",
        _stats("  before", q_b),
        _stats("  after ", q_a),
        f"  before frac alpha < 0.3: {(q_b < 0.3).mean() * 100:.2f} %",
        f"  after  frac alpha < 0.3: {(q_a < 0.3).mean() * 100:.2f} %",
        "",
        "[Triangle minimum interior angle (deg)]",
        _stats("  before", a_b, ".2f"),
        _stats("  after ", a_a, ".2f"),
        f"  before frac min-angle < 20 deg: {(a_b < 20).mean() * 100:.2f} %",
        f"  after  frac min-angle < 20 deg: {(a_a < 20).mean() * 100:.2f} %",
        "",
        "[Open-boundary perpendicularity (deg from 90)]",
        _stats("  before", perp_b, ".4f"),
        _stats("  after ", perp_a, ".4f"),
    ]
    summary = "\n".join(lines)
    print(summary)
    SUMMARY_TXT.write_text(summary + "\n", encoding="utf-8")
    print(f"[09] wrote {SUMMARY_TXT}")
    print(f"[09] wrote {SWAPPED_F14}")

    plot_alpha_hist(q_b, q_a, ALPHA_PNG)
    print(f"[09] wrote {ALPHA_PNG}")
    plot_min_angle_hist(a_b, a_a, MIN_ANGLE_PNG)
    print(f"[09] wrote {MIN_ANGLE_PNG}")
    plot_mesh(after.nodes, after.elements, after.depths, MESH_PNG)
    print(f"[09] wrote {MESH_PNG}")

    assert flipped_a == 0, f"{flipped_a} flipped triangles after swap"
    assert q_a.mean() >= q_b.mean() - 1e-6, "edge swap regressed mean alpha"
    delta_alpha = float(q_a.mean() - q_b.mean())
    delta_bad = float((a_a < 20).mean() - (a_b < 20).mean())
    print(
        f"[09] delta mean alpha = {delta_alpha:+.4f}  "
        f"delta frac<20deg = {delta_bad * 100:+.2f} %"
    )


if __name__ == "__main__":
    main()
