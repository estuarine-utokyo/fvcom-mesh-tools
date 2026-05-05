"""PoC #15: longest-edge bisection refinement on PoC #14 mesh.

PoC #14 produced a Tokyo Bay mesh with frac<20deg = 2.84 % - the
plateau set by initial sizing once swap+smooth has run. This PoC
applies longest-edge bisection (Rivara-style) to bad triangles to
push that fraction down further. Each pass stops early if it would
*increase* the bad-triangle count, since unrestricted bisection can
oscillate once the easy wins are taken.

Outputs:
    outputs/15_tokyo_bay_refined.14
    outputs/15_refine_summary.txt
    outputs/15_refine_min_angle_hist.png
    outputs/15_refine_alpha_hist.png
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
    refine_bad_triangles,
    signed_areas,
)
from fvcom_mesh_tools.io import read_fort14, write_fort14  # noqa: E402

REPO_ROOT = Path(__file__).resolve().parent.parent
INPUT_F14 = REPO_ROOT / "outputs" / "14_boundary_cleanup_islands_open50.14"
OUT_DIR = REPO_ROOT / "outputs"
REFINED_F14 = OUT_DIR / "15_tokyo_bay_refined.14"
SUMMARY_TXT = OUT_DIR / "15_refine_summary.txt"
MIN_ANGLE_PNG = OUT_DIR / "15_refine_min_angle_hist.png"
ALPHA_PNG = OUT_DIR / "15_refine_alpha_hist.png"

THRESHOLDS = [15.0, 18.0, 20.0]
MAX_PASSES = 5


def _stats(name: str, x: np.ndarray, fmt: str = ".4f") -> str:
    if x.size == 0:
        return f"{name}: empty"
    return (
        f"{name}: n={x.size:,}  "
        f"min={np.nanmin(x):{fmt}}  "
        f"p05={np.nanpercentile(x, 5):{fmt}}  "
        f"p50={np.nanpercentile(x, 50):{fmt}}  "
        f"p95={np.nanpercentile(x, 95):{fmt}}  "
        f"mean={np.nanmean(x):{fmt}}"
    )


def plot_min_angle(b: np.ndarray, after_runs: dict[float, np.ndarray], png: Path) -> None:
    fig, ax = plt.subplots(figsize=(9, 4.5), dpi=120)
    bins = np.linspace(0.0, 60.0, 61)
    ax.hist(b, bins=bins, alpha=0.55, label="before refine", color="tab:red")
    cmap = plt.cm.viridis(np.linspace(0.2, 0.85, len(after_runs)))
    for c, (thr, a) in zip(cmap, after_runs.items()):
        ax.hist(a, bins=bins, alpha=0.45, label=f"refine thr={thr:.0f}", color=c)
    ax.axvline(20.0, color="0.4", ls="--", lw=0.7, label="20-deg threshold")
    ax.set_xlabel("triangle minimum interior angle (deg)")
    ax.set_ylabel("count")
    ax.set_title("Refine threshold sweep: min-angle distribution")
    ax.legend()
    ax.grid(True, lw=0.3, color="0.9")
    fig.tight_layout()
    fig.savefig(png, bbox_inches="tight")
    plt.close(fig)


def plot_alpha(b: np.ndarray, after_runs: dict[float, np.ndarray], png: Path) -> None:
    fig, ax = plt.subplots(figsize=(9, 4.5), dpi=120)
    bins = np.linspace(0.0, 1.0, 41)
    ax.hist(b, bins=bins, alpha=0.55, label="before refine", color="tab:red")
    cmap = plt.cm.viridis(np.linspace(0.2, 0.85, len(after_runs)))
    for c, (thr, a) in zip(cmap, after_runs.items()):
        ax.hist(a, bins=bins, alpha=0.45, label=f"refine thr={thr:.0f}", color=c)
    ax.set_xlabel("alpha-quality (1 = equilateral)")
    ax.set_ylabel("count")
    ax.set_title("Refine threshold sweep: alpha-quality distribution")
    ax.legend()
    ax.grid(True, lw=0.3, color="0.9")
    fig.tight_layout()
    fig.savefig(png, bbox_inches="tight")
    plt.close(fig)


def main() -> None:
    if not INPUT_F14.exists():
        raise SystemExit(
            f"input not found: {INPUT_F14}\n"
            "Run notebooks/14_boundary_cleanup.py first."
        )
    OUT_DIR.mkdir(parents=True, exist_ok=True)

    print(f"[15] reading {INPUT_F14}")
    before = read_fort14(INPUT_F14)
    print(f"[15] NP={before.n_nodes:,}  NE={before.n_elements:,}")

    a_b = min_interior_angle(before)
    q_b = alpha_quality(before)

    after_min_angles: dict[float, np.ndarray] = {}
    after_alphas: dict[float, np.ndarray] = {}
    sections = [
        f"input:  {INPUT_F14}",
        f"max_passes: {MAX_PASSES}",
        "",
        "[Reference] NP=95,551 NE=182,603 alpha=0.979 frac<20deg=0%",
        "",
        "=== before refine ===",
        f"  NP={before.n_nodes:,}  NE={before.n_elements:,}",
        _stats("  alpha    ", q_b),
        f"  alpha < 0.3   : {(q_b < 0.3).mean() * 100:.2f} %",
        _stats("  min-angle", a_b, ".2f"),
        f"  min-angle <20 : {(a_b < 20).mean() * 100:.2f} %",
        "",
    ]

    for thr in THRESHOLDS:
        print(f"[15] refine threshold={thr} deg ...")
        t0 = time.perf_counter()
        out, info = refine_bad_triangles(
            before, min_angle_threshold=thr, max_passes=MAX_PASSES,
        )
        wall = time.perf_counter() - t0
        a = min_interior_angle(out)
        q = alpha_quality(out)
        after_min_angles[thr] = a
        after_alphas[thr] = q
        sections += [
            f"=== refine threshold={thr:g} deg (wall {wall:.2f} s) ===",
            f"  passes={info['passes']}  stop_reason={info['stop_reason']}",
            f"  nodes_inserted={info['total_nodes_inserted']:,}  "
            f"swaps={info['total_swaps']:,}",
            f"  NP={out.n_nodes:,}  NE={out.n_elements:,}  "
            f"flipped={int((signed_areas(out) <= 0).sum())}",
            _stats("  alpha    ", q),
            f"  alpha < 0.3   : {(q < 0.3).mean() * 100:.2f} %",
            _stats("  min-angle", a, ".2f"),
            f"  min-angle <20 : {(a < 20).mean() * 100:.2f} %",
            "",
        ]
        if thr == 20.0:
            write_fort14(out, REFINED_F14)
            print(f"[15] wrote {REFINED_F14} (threshold=20)")

    summary = "\n".join(sections)
    print(summary)
    SUMMARY_TXT.write_text(summary + "\n", encoding="utf-8")
    print(f"[15] wrote {SUMMARY_TXT}")

    plot_min_angle(a_b, after_min_angles, MIN_ANGLE_PNG)
    print(f"[15] wrote {MIN_ANGLE_PNG}")
    plot_alpha(q_b, after_alphas, ALPHA_PNG)
    print(f"[15] wrote {ALPHA_PNG}")


if __name__ == "__main__":
    main()
