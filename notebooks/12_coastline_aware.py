"""PoC #12: coastline-aware sizing via ``Hfun.add_feature``.

PoC #11 confirmed the swap+smooth quality pass plateaus at frac<20deg
~ 17 % because the underlying size function is uniform (hmin = 200 m,
hmax = 5000 m). The remaining slivers are concentrated along the
coast where gmsh has to stitch a uniform 200 m mesh against an
irregular shoreline. This PoC drives ``Hfun.add_feature`` with the
MLIT C23 Tokyo Bay coastline shapefile so the size function refines
adaptively near land.

Three configurations are run for comparison:

1. ``baseline`` - no coastline, no quality pass (PoC #7 mesh).
2. ``coast``    - coastline-aware sizing, no quality pass.
3. ``coast+qp`` - coastline-aware sizing + 6-round quality pass.

Outputs:
    outputs/12_baseline.14
    outputs/12_coast.14
    outputs/12_coast_qp.14
    outputs/12_coastline_summary.txt
    outputs/12_coastline_min_angle_hist.png
    outputs/12_coastline_alpha_hist.png
    outputs/12_coastline_meshes.png
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
)
from fvcom_mesh_tools.cli.buildmesh import main as buildmesh_main  # noqa: E402
from fvcom_mesh_tools.io import read_fort14  # noqa: E402
from fvcom_mesh_tools.plotting import MESH_PNG_DPI  # noqa: E402

REPO_ROOT = Path(__file__).resolve().parent.parent
DEM = REPO_ROOT / "data" / "bathymetry" / "tokyo_bay" / "dem_00_01_change.nc"
COASTLINE = (
    REPO_ROOT
    / "data"
    / "coastline"
    / "tokyo_bay"
    / "MLIT_C23"
    / "C23-06_TOKYOBAY.shp"
)
OUT_DIR = REPO_ROOT / "outputs"
F14_BASELINE = OUT_DIR / "12_baseline.14"
F14_COAST = OUT_DIR / "12_coast.14"
F14_COAST_QP = OUT_DIR / "12_coast_qp.14"
SUMMARY_TXT = OUT_DIR / "12_coastline_summary.txt"
MIN_ANGLE_PNG = OUT_DIR / "12_coastline_min_angle_hist.png"
ALPHA_PNG = OUT_DIR / "12_coastline_alpha_hist.png"
MESHES_PNG = OUT_DIR / "12_coastline_meshes.png"

HMIN_M = 200.0
HMAX_M = 5000.0
COAST_TARGET = 200.0
COAST_EXP_RATE = 0.005
QP_ROUNDS = 6


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


def _summarise(label: str, mesh, t_wall: float) -> tuple[str, dict]:
    q = alpha_quality(mesh)
    a = min_interior_angle(mesh)
    perp = open_bdy_perpendicularity(mesh, segment_index=0)
    flipped = int((signed_areas(mesh) <= 0).sum())
    out = [
        f"=== {label} ===",
        f"  wall: {t_wall:.2f} s   NP={mesh.n_nodes:,}  NE={mesh.n_elements:,}  flipped={flipped}",
        _stats("  alpha    ", q),
        f"  alpha < 0.3 frac    : {(q < 0.3).mean() * 100:.2f} %",
        _stats("  min-angle", a, ".2f"),
        f"  min-angle < 20 frac : {(a < 20).mean() * 100:.2f} %",
        _stats("  open-perp", perp, ".4f"),
    ]
    return "\n".join(out), {
        "alpha_mean": float(q.mean()),
        "alpha_lt_03": float((q < 0.3).mean()),
        "ma_p50": float(np.percentile(a, 50)),
        "ma_lt_20": float((a < 20).mean()),
        "flipped": flipped,
        "np": mesh.n_nodes,
        "ne": mesh.n_elements,
    }


def _run_buildmesh(
    label: str,
    out_path: Path,
    *,
    coastline: bool,
    quality_pass: int,
) -> tuple[float, ...]:
    args = [
        str(DEM), str(out_path),
        "--hmin", str(HMIN_M),
        "--hmax", str(HMAX_M),
        "--zmax", "0.0",
        "--interp-method", "linear",
        "--land-ibtype", "20",
        "--quality-pass", str(quality_pass),
        "--smooth-iters", "5",
        "--smooth-alpha", "0.5",
        "--perpfix-iters", "1",
        "--quiet",
    ]
    if coastline:
        args += [
            "--coastline", str(COASTLINE),
            "--coast-target-size", str(COAST_TARGET),
            "--coast-expansion-rate", str(COAST_EXP_RATE),
        ]
    print(f"[12] {label}: running fmesh-buildmesh ...")
    t0 = time.perf_counter()
    rc = buildmesh_main(args)
    t = time.perf_counter() - t0
    if rc != 0:
        raise SystemExit(f"buildmesh exited {rc}")
    print(f"[12] {label}: {t:.2f} s")
    return (t,)


def plot_min_angle_hists(meshes: dict, png: Path) -> None:
    fig, ax = plt.subplots(figsize=(9, 4.5), dpi=120)
    bins = np.linspace(0.0, 60.0, 61)
    colours = {"baseline": "tab:red", "coast": "tab:orange", "coast+qp": "tab:blue"}
    for label, mesh in meshes.items():
        ax.hist(
            min_interior_angle(mesh), bins=bins, alpha=0.55,
            label=label, color=colours[label],
        )
    ax.axvline(20.0, color="0.4", ls="--", lw=0.7, label="FVCOM 20-deg threshold")
    ax.set_xlabel("triangle minimum interior angle (deg)")
    ax.set_ylabel("count")
    ax.set_title("Coastline-aware sizing: min-angle distribution")
    ax.legend()
    ax.grid(True, lw=0.3, color="0.9")
    fig.tight_layout()
    fig.savefig(png, bbox_inches="tight")
    plt.close(fig)


def plot_alpha_hists(meshes: dict, png: Path) -> None:
    fig, ax = plt.subplots(figsize=(9, 4.5), dpi=120)
    bins = np.linspace(0.0, 1.0, 41)
    colours = {"baseline": "tab:red", "coast": "tab:orange", "coast+qp": "tab:blue"}
    for label, mesh in meshes.items():
        ax.hist(
            alpha_quality(mesh), bins=bins, alpha=0.55,
            label=label, color=colours[label],
        )
    ax.set_xlabel("alpha-quality (1 = equilateral)")
    ax.set_ylabel("count")
    ax.set_title("Coastline-aware sizing: alpha-quality distribution")
    ax.legend()
    ax.grid(True, lw=0.3, color="0.9")
    fig.tight_layout()
    fig.savefig(png, bbox_inches="tight")
    plt.close(fig)


def plot_meshes(meshes: dict, png: Path) -> None:
    fig, axes = plt.subplots(1, 3, figsize=(15, 6), dpi=120, sharex=True, sharey=True)
    titles = {
        "baseline": f"baseline (NP={meshes['baseline'].n_nodes:,})",
        "coast": f"coast only (NP={meshes['coast'].n_nodes:,})",
        "coast+qp": f"coast + qp (NP={meshes['coast+qp'].n_nodes:,})",
    }
    for ax, label in zip(axes, ["baseline", "coast", "coast+qp"]):
        m = meshes[label]
        ax.triplot(
            m.nodes[:, 0], m.nodes[:, 1], m.elements,
            color="0.4", lw=0.15,
        )
        ax.set_aspect("equal")
        ax.set_title(titles[label])
        ax.set_xlabel("lon (deg)")
        ax.grid(True, lw=0.3, color="0.9")
    axes[0].set_ylabel("lat (deg)")
    fig.tight_layout()
    fig.savefig(png, bbox_inches="tight", dpi=MESH_PNG_DPI)
    plt.close(fig)


def main() -> None:
    if not DEM.exists():
        raise SystemExit(f"DEM not found: {DEM}")
    if not COASTLINE.exists():
        raise SystemExit(f"coastline not found: {COASTLINE}")
    OUT_DIR.mkdir(parents=True, exist_ok=True)

    walls: dict[str, float] = {}
    walls["baseline"], = _run_buildmesh("baseline", F14_BASELINE, coastline=False, quality_pass=0)
    walls["coast"], = _run_buildmesh("coast", F14_COAST, coastline=True, quality_pass=0)
    walls["coast+qp"], = _run_buildmesh(
        "coast+qp", F14_COAST_QP, coastline=True, quality_pass=QP_ROUNDS,
    )

    meshes = {
        "baseline": read_fort14(F14_BASELINE),
        "coast": read_fort14(F14_COAST),
        "coast+qp": read_fort14(F14_COAST_QP),
    }

    sections = []
    rows = {}
    for label, mesh in meshes.items():
        text, row = _summarise(label, mesh, walls[label])
        sections.append(text)
        rows[label] = row

    summary = (
        f"DEM:       {DEM}\n"
        f"coastline: {COASTLINE}\n"
        f"hmin={HMIN_M:g} m  hmax={HMAX_M:g} m\n"
        f"coast_target={COAST_TARGET:g} m  coast_expansion_rate={COAST_EXP_RATE:g}\n"
        f"quality_pass rounds = {QP_ROUNDS}\n\n"
        + "\n\n".join(sections)
        + "\n\n[Comparison vs reference (NP=95,551, alpha=0.979, frac<20deg=0)]\n"
        + "\n".join(
            f"  {label:9s}: alpha_mean={r['alpha_mean']:.4f}  "
            f"frac<20deg={r['ma_lt_20'] * 100:.2f} %  "
            f"NP={r['np']:,}  NE={r['ne']:,}"
            for label, r in rows.items()
        )
    )
    print(summary)
    SUMMARY_TXT.write_text(summary + "\n", encoding="utf-8")
    print(f"[12] wrote {SUMMARY_TXT}")

    plot_min_angle_hists(meshes, MIN_ANGLE_PNG)
    print(f"[12] wrote {MIN_ANGLE_PNG}")
    plot_alpha_hists(meshes, ALPHA_PNG)
    print(f"[12] wrote {ALPHA_PNG}")
    plot_meshes(meshes, MESHES_PNG)
    print(f"[12] wrote {MESHES_PNG}")

    # Sanity: every mesh must be valid and non-empty.
    for label, mesh in meshes.items():
        assert (signed_areas(mesh) > 0).all(), f"{label}: flipped triangles in output"
        assert mesh.n_elements > 0


if __name__ == "__main__":
    main()
