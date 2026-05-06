"""PoC #19: full fmesh-buildmesh pipeline through the oceanmesh engine.

PoC #18 confirmed oceanmesh produces near-equilateral meshes (alpha
0.96, frac<20deg 0.03 %) on Tokyo Bay's raw output. This PoC drives
the same Tokyo Bay inputs as PoC #16 through ``fmesh-buildmesh``
with ``--engine oceanmesh`` and the full post-processing chain
(boundary classify, river inflow, perpfix, fort.14 write) to confirm
end-to-end success and report final numbers comparable to PoC #16.

Outputs:
    outputs/19_tokyo_bay_oceanmesh.14
    outputs/19_tokyo_bay_oceanmesh_summary.txt
    outputs/19_tokyo_bay_oceanmesh_mesh.png
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
    signed_areas,
)
from fvcom_mesh_tools.cli.buildmesh import main as buildmesh_main  # noqa: E402
from fvcom_mesh_tools.io import load_river_points, read_fort14  # noqa: E402

REPO_ROOT = Path(__file__).resolve().parent.parent
DEM = REPO_ROOT / "data" / "bathymetry" / "tokyo_bay" / "dem_00_01_change.nc"
COASTLINE = (
    REPO_ROOT / "data" / "coastline" / "tokyo_bay"
    / "MLIT_C23" / "C23-06_TOKYOBAY.shp"
)
RIVERS = REPO_ROOT / "data" / "rivers" / "tokyo_bay" / "tokyo_bay_rivers.csv"
OUT_DIR = REPO_ROOT / "outputs"
F14_OUT = OUT_DIR / "19_tokyo_bay_oceanmesh.14"
SUMMARY_TXT = OUT_DIR / "19_tokyo_bay_oceanmesh_summary.txt"
MESH_PNG = OUT_DIR / "19_tokyo_bay_oceanmesh_mesh.png"

HMIN_M = 200.0
HMAX_M = 5000.0


def plot_mesh(mesh, river_pts: np.ndarray, png: Path) -> None:
    fig, ax = plt.subplots(figsize=(9, 7), dpi=120)
    ax.triplot(
        mesh.nodes[:, 0], mesh.nodes[:, 1], mesh.elements,
        color="0.7", lw=0.1,
    )
    plotted_open = False
    for seg in mesh.open_boundaries:
        ax.plot(
            mesh.nodes[seg, 0], mesh.nodes[seg, 1], "-",
            color="tab:red", lw=1.2,
            label="open" if not plotted_open else None,
        )
        plotted_open = True
    plotted_land = False
    plotted_river = False
    for ib, seg in mesh.land_boundaries:
        if ib == 21:
            ax.plot(
                mesh.nodes[seg, 0], mesh.nodes[seg, 1], "-",
                color="tab:blue", lw=2.0,
                label="river (ibtype=21)" if not plotted_river else None,
            )
            plotted_river = True
        else:
            ax.plot(
                mesh.nodes[seg, 0], mesh.nodes[seg, 1], "-",
                color="0.4", lw=0.6,
                label=f"land (ibtype={ib})" if not plotted_land else None,
            )
            plotted_land = True
    ax.scatter(
        river_pts[:, 0], river_pts[:, 1], marker="x",
        color="tab:purple", s=60, label="river input pt",
    )
    ax.set_aspect("equal")
    ax.set_xlabel("lon (deg)")
    ax.set_ylabel("lat (deg)")
    ax.set_title("PoC #19  fmesh-buildmesh --engine oceanmesh")
    ax.legend(loc="lower right", fontsize=8)
    fig.tight_layout()
    fig.savefig(png, bbox_inches="tight")
    plt.close(fig)


def main() -> None:
    for p in (DEM, COASTLINE, RIVERS):
        if not p.exists():
            raise SystemExit(f"required input missing: {p}")
    OUT_DIR.mkdir(parents=True, exist_ok=True)

    args = [
        str(DEM), str(F14_OUT),
        "--engine", "oceanmesh",
        "--hmin", str(HMIN_M), "--hmax", str(HMAX_M), "--zmax", "0.0",
        "--interp-method", "linear",
        "--coastline", str(COASTLINE),
        "--om-slope-parameter", "20",
        "--om-gradation", "0.15",
        "--om-max-iter", "50",
        "--om-seed", "0",
        "--open-merge-coast-gap", "50",
        "--river-inflow-points", str(RIVERS),
        "--river-segment-nodes", "5",
        "--river-ibtype", "21",
        "--river-snap-tol-m", "5000",
        "--quality-pass", "0",
        "--refine-min-angle", "0",
        "--land-ibtype", "20",
        "--perpfix-iters", "1",
    ]
    print("[19] running fmesh-buildmesh --engine oceanmesh ...")
    t0 = time.perf_counter()
    rc = buildmesh_main(args)
    wall = time.perf_counter() - t0
    if rc != 0:
        raise SystemExit(f"buildmesh exited {rc}")
    print(f"[19] wall: {wall:.2f} s")

    river_pts = load_river_points([RIVERS])
    mesh = read_fort14(F14_OUT)

    n_open = len(mesh.open_boundaries)
    n_land = len(mesh.land_boundaries)
    n_river_segs = sum(1 for ib, _ in mesh.land_boundaries if ib == 21)
    n_river_nodes = sum(b.size for ib, b in mesh.land_boundaries if ib == 21)
    n_land20 = sum(1 for ib, _ in mesh.land_boundaries if ib == 20)
    flipped = int((signed_areas(mesh) <= 0).sum())
    q = alpha_quality(mesh)
    a = min_interior_angle(mesh)

    summary_lines = [
        f"DEM:       {DEM}",
        f"coastline: {COASTLINE}",
        f"rivers:    {RIVERS}",
        "engine:    oceanmesh",
        f"wall:      {wall:.2f} s",
        "",
        f"NP={mesh.n_nodes:,}  NE={mesh.n_elements:,}  flipped={flipped}",
        "",
        "[Boundary structure]",
        f"  open boundaries     : {n_open}",
        f"  land segments total : {n_land}",
        f"    ibtype=20 (coast) : {n_land20}",
        f"    ibtype=21 (river) : {n_river_segs} ({n_river_nodes} nodes)",
        "",
        "[Quality]",
        f"  alpha mean          : {q.mean():.4f}",
        f"  alpha < 0.3         : {(q < 0.3).mean() * 100:.2f} %",
        f"  min-angle p50       : {np.percentile(a, 50):.2f}",
        f"  min-angle < 20 deg  : {(a < 20).mean() * 100:.2f} %",
        "",
        "[River input points]",
    ]
    for px, py in river_pts.tolist():
        summary_lines.append(f"  {px:.4f}, {py:.4f}")
    summary = "\n".join(summary_lines)
    print(summary)
    SUMMARY_TXT.write_text(summary + "\n", encoding="utf-8")
    print(f"[19] wrote {SUMMARY_TXT}")

    plot_mesh(mesh, river_pts, MESH_PNG)
    print(f"[19] wrote {MESH_PNG}")

    assert flipped == 0
    assert n_river_segs > 0, "no river ibtype=21 segments produced"
    assert n_open > 0, "no open boundary segments produced"


if __name__ == "__main__":
    main()
