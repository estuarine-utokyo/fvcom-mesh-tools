"""PoC #20: full fmesh-buildmesh pipeline through oceanmesh on Osaka Bay.

Cross-basin de-risk of the new default engine (PoC #19 made
``--engine oceanmesh`` the default in fmesh-buildmesh on Tokyo Bay).
Inputs mirror PoC #17 (Osaka Bay × OCSMesh): SRTM15+ subset DEM,
GSHHS-f L1 global coastline, 4 rivers (Yodo / Yamato / Mukou /
Kanzaki).

If ``--engine oceanmesh`` ports cleanly here, we have direct evidence
that the new default works on geometries unrelated to its training
basin (Osaka has two straits + a large internal island, no shared
recipe with Tokyo Bay).

Outputs:
    outputs/20_osaka_bay_oceanmesh.14
    outputs/20_osaka_bay_oceanmesh_summary.txt
    outputs/20_osaka_bay_oceanmesh_mesh.png
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
from fvcom_mesh_tools.io import (  # noqa: E402
    load_river_points,
    read_fort14,
    subset_dem_to_geotiff,
)

REPO_ROOT = Path(__file__).resolve().parent.parent
SRTM = REPO_ROOT / "data" / "bathymetry" / "SRTM15plus" / "SRTM15+.nc"
COASTLINE = (
    REPO_ROOT / "data" / "coastline" / "GSHHS" / "f" / "GSHHS_f_L1.shp"
)
RIVERS = REPO_ROOT / "data" / "rivers" / "osaka_bay" / "osaka_bay_rivers.csv"
DEM_SUBSET = (
    REPO_ROOT / "data" / "bathymetry" / "osaka_bay" / "srtm15_osaka_bay.tif"
)
OUT_DIR = REPO_ROOT / "outputs"
F14_OUT = OUT_DIR / "20_osaka_bay_oceanmesh.14"
SUMMARY_TXT = OUT_DIR / "20_osaka_bay_oceanmesh_summary.txt"
MESH_PNG = OUT_DIR / "20_osaka_bay_oceanmesh_mesh.png"

BBOX = (134.90, 34.20, 135.55, 34.85)
HMIN_M = 200.0
HMAX_M = 5000.0


def plot_mesh(mesh, river_pts, png: Path) -> None:
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
    ax.set_title("PoC #20  Osaka Bay  --engine oceanmesh")
    ax.legend(loc="lower right", fontsize=8)
    fig.tight_layout()
    fig.savefig(png, bbox_inches="tight")
    plt.close(fig)


def main() -> None:
    for p in (SRTM, COASTLINE, RIVERS):
        if not p.exists():
            raise SystemExit(f"required input missing: {p}")
    OUT_DIR.mkdir(parents=True, exist_ok=True)
    DEM_SUBSET.parent.mkdir(parents=True, exist_ok=True)

    if not DEM_SUBSET.exists():
        print(f"[20] subsetting SRTM15+ to {BBOX} -> {DEM_SUBSET}")
        t0 = time.perf_counter()
        subset_dem_to_geotiff(SRTM, DEM_SUBSET, BBOX, src_var="z")
        print(f"[20] subset: {time.perf_counter() - t0:.2f} s")

    args = [
        str(DEM_SUBSET), str(F14_OUT),
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
    print("[20] running fmesh-buildmesh --engine oceanmesh on Osaka Bay ...")
    t0 = time.perf_counter()
    rc = buildmesh_main(args)
    wall = time.perf_counter() - t0
    if rc != 0:
        raise SystemExit(f"buildmesh exited {rc}")
    print(f"[20] wall: {wall:.2f} s")

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
        f"DEM subset: {DEM_SUBSET}  (from {SRTM.name})",
        f"bbox:       {BBOX}",
        f"coastline:  {COASTLINE}",
        f"rivers:     {RIVERS}",
        f"engine:     oceanmesh",
        f"wall:       {wall:.2f} s",
        "",
        f"NP={mesh.n_nodes:,}  NE={mesh.n_elements:,}  flipped={flipped}",
        "",
        "[Boundary structure]",
        f"  open boundaries     : {n_open}",
        "  open segment sizes  : "
        + ", ".join(str(s.size) for s in mesh.open_boundaries),
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
    print(f"[20] wrote {SUMMARY_TXT}")

    plot_mesh(mesh, river_pts, MESH_PNG)
    print(f"[20] wrote {MESH_PNG}")

    assert flipped == 0
    assert n_river_segs > 0, "no river ibtype=21 segments produced"
    assert n_open > 0, "no open boundary segments produced"


if __name__ == "__main__":
    main()
