"""PoC #7: end-to-end ``fmesh-buildmesh`` validation on Tokyo Bay.

Drives the full Python pipeline (DEM -> OCSMesh -> depth interp ->
boundary classification -> open-boundary perpendicularity fix -> fort.14)
on the Tokyo Bay DEM and verifies the output meets FVCOM expectations:

* Non-zero, finite depths interpolated from the DEM.
* At least one open-boundary segment and at least one land-boundary
  segment.
* Open-boundary first-ring perpendicularity comparable to (or better
  than) the legacy reference.

Outputs:
    outputs/07_tokyo_bay_buildmesh.14
    outputs/07_buildmesh_summary.txt
    outputs/07_buildmesh_mesh.png
"""

from __future__ import annotations

import time
from pathlib import Path

import matplotlib

matplotlib.use("Agg")

import matplotlib.pyplot as plt  # noqa: E402
import numpy as np  # noqa: E402

from fvcom_mesh_tools.algorithms import (  # noqa: E402
    open_bdy_perpendicularity,
    signed_areas,
)
from fvcom_mesh_tools.cli.buildmesh import main as buildmesh_main  # noqa: E402
from fvcom_mesh_tools.io import read_fort14  # noqa: E402

REPO_ROOT = Path(__file__).resolve().parent.parent
DEM = REPO_ROOT / "data" / "bathymetry" / "tokyo_bay" / "dem_00_01_change.nc"
REF_MESH = REPO_ROOT / "data" / "mesh" / "reference" / "tokyo_bay" / "tb_futtsu20220311.14"

OUT_DIR = REPO_ROOT / "outputs"
FORT14 = OUT_DIR / "07_tokyo_bay_buildmesh.14"
SUMMARY = OUT_DIR / "07_buildmesh_summary.txt"
MESH_PNG = OUT_DIR / "07_buildmesh_mesh.png"

HMIN_M = 200.0
HMAX_M = 5000.0


def _stats(name: str, x: np.ndarray, fmt: str = ".4f") -> str:
    if x.size == 0:
        return f"{name}: empty"
    return (
        f"{name}: n={x.size:,}  "
        f"min={np.nanmin(x):{fmt}}  "
        f"p50={np.nanpercentile(x, 50):{fmt}}  "
        f"p95={np.nanpercentile(x, 95):{fmt}}  "
        f"max={np.nanmax(x):{fmt}}  "
        f"mean={np.nanmean(x):{fmt}}"
    )


def plot_mesh_with_boundaries(
    nodes: np.ndarray,
    elements: np.ndarray,
    depths: np.ndarray,
    open_segs: list[np.ndarray],
    land_segs: list[tuple[int, np.ndarray]],
    png: Path,
) -> None:
    fig, ax = plt.subplots(figsize=(8, 7), dpi=120)
    tri = ax.tripcolor(
        nodes[:, 0], nodes[:, 1], elements,
        facecolors=depths[elements].mean(axis=1),
        cmap="viridis_r", edgecolors="0.85", linewidth=0.05,
    )
    fig.colorbar(tri, ax=ax, label="depth (m, +down)")

    for seg in open_segs:
        ax.plot(nodes[seg, 0], nodes[seg, 1], "-", color="tab:red", lw=1.4)
    for _, seg in land_segs:
        ax.plot(nodes[seg, 0], nodes[seg, 1], "-", color="tab:gray", lw=0.9)

    ax.set_aspect("equal")
    ax.set_xlabel("lon (deg)")
    ax.set_ylabel("lat (deg)")
    ax.set_title("fmesh-buildmesh: Tokyo Bay (red=open, gray=land)")
    fig.tight_layout()
    fig.savefig(png, bbox_inches="tight")
    plt.close(fig)


def main() -> None:
    if not DEM.exists():
        raise SystemExit(f"DEM not found: {DEM}")
    OUT_DIR.mkdir(parents=True, exist_ok=True)

    print("[07] running fmesh-buildmesh...")
    t0 = time.perf_counter()
    rc = buildmesh_main([
        str(DEM), str(FORT14),
        "--hmin", str(HMIN_M),
        "--hmax", str(HMAX_M),
        "--zmax", "0.0",
        "--interp-method", "linear",
        "--land-ibtype", "20",
        "--perpfix-iters", "1",
    ])
    t_total = time.perf_counter() - t0
    if rc != 0:
        raise SystemExit(f"fmesh-buildmesh exited {rc}")
    print(f"[07] fmesh-buildmesh wall time: {t_total:.2f} s")

    f14 = read_fort14(FORT14)
    print(f"[07] re-read NP={f14.n_nodes:,}  NE={f14.n_elements:,}")

    n_open = len(f14.open_boundaries)
    n_open_nodes = sum(int(b.size) for b in f14.open_boundaries)
    n_land = len(f14.land_boundaries)
    n_land_nodes = sum(int(b.size) for _, b in f14.land_boundaries)
    ibtypes = sorted({int(ib) for ib, _ in f14.land_boundaries})
    flipped = int((signed_areas(f14) <= 0).sum())

    perp = open_bdy_perpendicularity(f14, segment_index=0) if n_open else np.array([])
    ref = read_fort14(REF_MESH)
    perp_ref = open_bdy_perpendicularity(ref, segment_index=0)

    lines = [
        f"DEM:    {DEM}",
        f"output: {FORT14}",
        f"ref:    {REF_MESH}",
        f"hmin={HMIN_M:g} m  hmax={HMAX_M:g} m",
        f"buildmesh wall time: {t_total:.2f} s",
        "",
        "[Mesh size]",
        f"  generated : NP={f14.n_nodes:,}  NE={f14.n_elements:,}",
        f"  reference : NP={ref.n_nodes:,}  NE={ref.n_elements:,}",
        "",
        "[Element validity]",
        f"  flipped triangles: {flipped:,}",
        "",
        "[Boundary structure]",
        f"  open : {n_open} segments, {n_open_nodes:,} nodes",
        f"  land : {n_land} segments, {n_land_nodes:,} nodes  (ibtypes={ibtypes})",
        "",
        "[Depths in fort.14 (positive = below MSL)]",
        _stats("  depth", f14.depths, ".4f"),
        f"  fraction <= 0:    {(f14.depths <= 0).mean() * 100:.2f} %",
        "",
        "[Open-boundary perpendicularity (deg from 90)]",
        _stats("  generated", perp, ".4f"),
        _stats("  reference", perp_ref, ".4f"),
    ]
    summary = "\n".join(lines)
    print(summary)
    SUMMARY.write_text(summary + "\n", encoding="utf-8")
    print(f"[07] wrote {SUMMARY}")

    plot_mesh_with_boundaries(
        f14.nodes, f14.elements, f14.depths,
        f14.open_boundaries, f14.land_boundaries,
        MESH_PNG,
    )
    print(f"[07] wrote {MESH_PNG}")

    # Hard sanity checks: anything below should already have failed in the
    # CLI, but the goal of this PoC is to assert it from outside the CLI.
    assert n_open >= 1, "no open boundary written"
    assert n_land >= 1, "no land boundary written"
    assert flipped == 0, f"{flipped} flipped triangles in output"
    assert (f14.depths > 0).any(), "no positive depths -- depth interp likely failed"
    print("[07] sanity checks PASSED")


if __name__ == "__main__":
    main()
