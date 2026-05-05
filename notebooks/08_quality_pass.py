"""PoC #8: Laplacian quality pass on the fmesh-buildmesh Tokyo Bay output.

The minimal pipeline produces a valid fort.14 (PoC #7) but with
20%-of-triangles below the 20-degree min-interior-angle threshold,
which blocks FVCOM stability checks. Run a damped Laplacian smoother
on every interior node, with all open and land boundary nodes pinned
and a per-iteration no-flip rollback, and re-measure.

Outputs:
    outputs/08_tokyo_bay_smoothed.14
    outputs/08_quality_pass_summary.txt
    outputs/08_quality_pass_alpha_hist.png
    outputs/08_quality_pass_min_angle_hist.png
    outputs/08_quality_pass_mesh.png
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
)
from fvcom_mesh_tools.io import read_fort14, write_fort14  # noqa: E402

REPO_ROOT = Path(__file__).resolve().parent.parent
INPUT_F14 = REPO_ROOT / "outputs" / "07_tokyo_bay_buildmesh.14"
OUT_DIR = REPO_ROOT / "outputs"
SMOOTHED_F14 = OUT_DIR / "08_tokyo_bay_smoothed.14"
SUMMARY_TXT = OUT_DIR / "08_quality_pass_summary.txt"
ALPHA_PNG = OUT_DIR / "08_quality_pass_alpha_hist.png"
MIN_ANGLE_PNG = OUT_DIR / "08_quality_pass_min_angle_hist.png"
MESH_PNG = OUT_DIR / "08_quality_pass_mesh.png"

N_ITERS = 30
ALPHA = 0.5


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
    ax.hist(b, bins=bins, alpha=0.6, label="before smooth", color="tab:red")
    ax.hist(a, bins=bins, alpha=0.6, label="after smooth", color="tab:blue")
    ax.set_xlabel("alpha-quality (1 = equilateral)")
    ax.set_ylabel("count")
    ax.set_title(f"Laplacian smoothing: n_iters={N_ITERS}, alpha={ALPHA}")
    ax.legend()
    ax.grid(True, lw=0.3, color="0.9")
    fig.tight_layout()
    fig.savefig(png, bbox_inches="tight")
    plt.close(fig)


def plot_min_angle_hist(b: np.ndarray, a: np.ndarray, png: Path) -> None:
    fig, ax = plt.subplots(figsize=(8, 4.5), dpi=120)
    bins = np.linspace(0.0, 60.0, 61)
    ax.hist(b, bins=bins, alpha=0.6, label="before smooth", color="tab:red")
    ax.hist(a, bins=bins, alpha=0.6, label="after smooth", color="tab:blue")
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
    ax.set_title(f"Tokyo Bay after Laplacian quality pass (n_iters={N_ITERS}, alpha={ALPHA})")
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

    print(f"[08] reading {INPUT_F14}")
    before = read_fort14(INPUT_F14)
    print(f"[08] NP={before.n_nodes:,}  NE={before.n_elements:,}")

    q_b = alpha_quality(before)
    a_b = min_interior_angle(before)
    perp_b = open_bdy_perpendicularity(before, segment_index=0)

    t0 = time.perf_counter()
    after, info = laplacian_smooth(before, n_iters=N_ITERS, alpha=ALPHA, prevent_flips=True)
    t_smooth = time.perf_counter() - t0
    print(f"[08] laplacian_smooth: {t_smooth:.2f} s")
    print(f"[08] reverts per iteration: {info['reverts']}")
    print(f"[08] degenerate after: {info['degenerate_remaining']}")

    q_a = alpha_quality(after)
    a_a = min_interior_angle(after)
    perp_a = open_bdy_perpendicularity(after, segment_index=0)
    flipped_a = int((signed_areas(after) <= 0).sum())

    movable_count = int(info["moved"].sum())
    disp = np.linalg.norm(after.nodes - before.nodes, axis=1)

    write_fort14(after, SMOOTHED_F14)

    lines = [
        f"input:  {INPUT_F14}",
        f"output: {SMOOTHED_F14}",
        f"n_iters={N_ITERS}  alpha={ALPHA}",
        f"laplacian_smooth wall time: {t_smooth:.2f} s",
        "",
        "[Move counts]",
        f"  movable nodes touched: {movable_count:,}",
        _stats("  ||after - before||", disp, ".6f"),
        f"  reverts per iter: {info['reverts']}",
        "",
        "[Element validity]",
        f"  flipped triangles after: {flipped_a:,}",
        f"  degenerate remaining (info): {info['degenerate_remaining']}",
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
    print(f"[08] wrote {SUMMARY_TXT}")
    print(f"[08] wrote {SMOOTHED_F14}")

    plot_alpha_hist(q_b, q_a, ALPHA_PNG)
    print(f"[08] wrote {ALPHA_PNG}")
    plot_min_angle_hist(a_b, a_a, MIN_ANGLE_PNG)
    print(f"[08] wrote {MIN_ANGLE_PNG}")
    plot_mesh(after.nodes, after.elements, after.depths, MESH_PNG)
    print(f"[08] wrote {MESH_PNG}")

    # Sanity: the no-flip rollback must work; we should never produce
    # an invalid mesh even when the move is otherwise unhelpful.
    assert flipped_a == 0, f"{flipped_a} flipped triangles after smooth"

    # Empirical finding: damped Laplacian alone barely shifts the mean
    # alpha-quality on this Tokyo Bay mesh. Slivers are topologically
    # locked, not metrically locked - they need an edge-swap pass
    # (PoC #9). Print the verdict but do not fail.
    delta = float(q_a.mean() - q_b.mean())
    if delta > 1e-3:
        print(f"[08] OUTCOME: improved (mean alpha {delta:+.4f}).")
    elif abs(delta) <= 1e-3:
        print(f"[08] OUTCOME: neutral (mean alpha {delta:+.4f}).")
    else:
        print(
            f"[08] OUTCOME: regressed (mean alpha {delta:+.4f}). "
            f"Topological fixes (edge swapping) are required next."
        )


if __name__ == "__main__":
    main()
