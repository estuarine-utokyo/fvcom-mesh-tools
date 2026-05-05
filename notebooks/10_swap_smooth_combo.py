"""PoC #10: alternating edge-swap + Laplacian smoothing on PoC #7 output.

PoC #8 found Laplacian-only smoothing did nothing. PoC #9 found edge
swapping alone monotonically improved quality but stalled with ~18 %
of triangles still below the FVCOM 20-degree min-angle threshold. The
classic recipe is to *alternate*: swap re-triangulates around bad
positions, smooth then drifts the new node neighbourhood toward better
positions, swap again unlocks new candidates, etc.

Outputs:
    outputs/10_tokyo_bay_quality.14
    outputs/10_swap_smooth_summary.txt
    outputs/10_swap_smooth_min_angle_hist.png
    outputs/10_swap_smooth_alpha_hist.png
    outputs/10_swap_smooth_mesh.png
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
    laplacian_smooth,
    min_interior_angle,
    open_bdy_perpendicularity,
    signed_areas,
    swap_edges_for_quality,
)
from fvcom_mesh_tools.io import read_fort14, write_fort14  # noqa: E402

REPO_ROOT = Path(__file__).resolve().parent.parent
INPUT_F14 = REPO_ROOT / "outputs" / "07_tokyo_bay_buildmesh.14"
OUT_DIR = REPO_ROOT / "outputs"
OUTPUT_F14 = OUT_DIR / "10_tokyo_bay_quality.14"
SUMMARY_TXT = OUT_DIR / "10_swap_smooth_summary.txt"
ALPHA_PNG = OUT_DIR / "10_swap_smooth_alpha_hist.png"
MIN_ANGLE_PNG = OUT_DIR / "10_swap_smooth_min_angle_hist.png"
MESH_PNG = OUT_DIR / "10_swap_smooth_mesh.png"

OUTER_ROUNDS = 6
SMOOTH_ITERS_PER_ROUND = 5
SMOOTH_ALPHA = 0.5
MAX_SWAP_ITERS_PER_ROUND = 10


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
    ax.hist(b, bins=bins, alpha=0.6, label="before", color="tab:red")
    ax.hist(a, bins=bins, alpha=0.6, label="after swap+smooth", color="tab:blue")
    ax.set_xlabel("alpha-quality (1 = equilateral)")
    ax.set_ylabel("count")
    ax.set_title(
        f"swap+smooth: outer={OUTER_ROUNDS}, smooth_iters={SMOOTH_ITERS_PER_ROUND}, "
        f"alpha={SMOOTH_ALPHA}"
    )
    ax.legend()
    ax.grid(True, lw=0.3, color="0.9")
    fig.tight_layout()
    fig.savefig(png, bbox_inches="tight")
    plt.close(fig)


def plot_min_angle_hist(b: np.ndarray, a: np.ndarray, png: Path) -> None:
    fig, ax = plt.subplots(figsize=(8, 4.5), dpi=120)
    bins = np.linspace(0.0, 60.0, 61)
    ax.hist(b, bins=bins, alpha=0.6, label="before", color="tab:red")
    ax.hist(a, bins=bins, alpha=0.6, label="after swap+smooth", color="tab:blue")
    ax.axvline(20.0, color="0.4", ls="--", lw=0.7, label="FVCOM 20-deg threshold")
    ax.set_xlabel("triangle minimum interior angle (deg)")
    ax.set_ylabel("count")
    ax.set_title("Min-interior-angle distribution")
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
    ax.set_title(
        f"Tokyo Bay after swap+smooth quality pass (rounds={OUTER_ROUNDS})"
    )
    fig.tight_layout()
    fig.savefig(png, bbox_inches="tight")
    plt.close(fig)


def main() -> None:
    if not INPUT_F14.exists():
        raise SystemExit(
            f"input not found: {INPUT_F14}\n"
            f"Run notebooks/07_buildmesh_e2e.py first."
        )
    OUT_DIR.mkdir(parents=True, exist_ok=True)

    print(f"[10] reading {INPUT_F14}")
    before = read_fort14(INPUT_F14)
    print(f"[10] NP={before.n_nodes:,}  NE={before.n_elements:,}")

    q_b = alpha_quality(before)
    a_b = min_interior_angle(before)
    perp_b = open_bdy_perpendicularity(before, segment_index=0)

    mesh = before
    history: list[dict] = []
    t0 = time.perf_counter()
    for r in range(OUTER_ROUNDS):
        # 1) Edge-swap pass (re-triangulate)
        mesh, swap_info = swap_edges_for_quality(
            mesh, max_iters=MAX_SWAP_ITERS_PER_ROUND
        )
        n_swaps = swap_info["total_swaps"]

        # 2) Laplacian smoothing pass (drift nodes)
        mesh, smooth_info = laplacian_smooth(
            mesh,
            n_iters=SMOOTH_ITERS_PER_ROUND,
            alpha=SMOOTH_ALPHA,
            prevent_flips=True,
        )

        ma_now = min_interior_angle(mesh)
        aq_now = alpha_quality(mesh)
        history.append({
            "round": r + 1,
            "swaps": n_swaps,
            "smooth_reverts": int(sum(smooth_info["reverts"])),
            "frac_below_20": float((ma_now < 20).mean()),
            "alpha_mean": float(aq_now.mean()),
        })
        print(
            f"[10] round {r + 1}: swaps={n_swaps:,}  "
            f"smooth_reverts={int(sum(smooth_info['reverts'])):,}  "
            f"frac<20deg={(ma_now < 20).mean() * 100:.2f} %  "
            f"alpha_mean={aq_now.mean():.4f}"
        )
    t_total = time.perf_counter() - t0
    print(f"[10] total wall: {t_total:.2f} s")

    after = mesh
    q_a = alpha_quality(after)
    a_a = min_interior_angle(after)
    perp_a = open_bdy_perpendicularity(after, segment_index=0)
    flipped_a = int((signed_areas(after) <= 0).sum())

    write_fort14(after, OUTPUT_F14)

    lines = [
        f"input:  {INPUT_F14}",
        f"output: {OUTPUT_F14}",
        f"rounds={OUTER_ROUNDS}  smooth_iters={SMOOTH_ITERS_PER_ROUND}  "
        f"smooth_alpha={SMOOTH_ALPHA}  swap_max={MAX_SWAP_ITERS_PER_ROUND}",
        f"wall time: {t_total:.2f} s",
        "",
        "[Convergence]",
    ]
    for h in history:
        lines.append(
            f"  round {h['round']}: swaps={h['swaps']:,}  "
            f"reverts={h['smooth_reverts']:,}  "
            f"frac<20deg={h['frac_below_20'] * 100:.2f} %  "
            f"alpha_mean={h['alpha_mean']:.4f}"
        )
    lines += [
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
    print(f"[10] wrote {SUMMARY_TXT}")
    print(f"[10] wrote {OUTPUT_F14}")

    plot_alpha_hist(q_b, q_a, ALPHA_PNG)
    print(f"[10] wrote {ALPHA_PNG}")
    plot_min_angle_hist(a_b, a_a, MIN_ANGLE_PNG)
    print(f"[10] wrote {MIN_ANGLE_PNG}")
    plot_mesh(after.nodes, after.elements, after.depths, MESH_PNG)
    print(f"[10] wrote {MESH_PNG}")

    assert flipped_a == 0, f"{flipped_a} flipped triangles after combo"


if __name__ == "__main__":
    main()
