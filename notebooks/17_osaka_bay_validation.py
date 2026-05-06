"""PoC #17: second-basin validation of fmesh-buildmesh on Osaka Bay.

Cross-checks that the pipeline tuned on Tokyo Bay (PoC #16) ports to a
basin with very different topology:

* two straits (Akashi N-W, Tomogashima S-W) instead of one,
* a large internal island (Awaji) on the west side,
* and a different river roster (Yodo, Yamato, Mukou, Kanzaki).

The DEM source is the global SRTM15+ raster (no embedded CRS, so we
subset and re-emit a CF-compliant GeoTIFF first via
``fmesh-subset-dem``); coastlines come from GSHHS-f L1 (continental +
large-island full resolution); river-mouth points are in
``data/rivers/osaka_bay/osaka_bay_rivers.csv``. The same
``fmesh-buildmesh`` flag set as PoC #16 is used; if Osaka Bay needs
different parameters, the gap should surface here.

Outputs:
    outputs/17_osaka_bay.14
    outputs/17_osaka_bay_summary.txt
    outputs/17_osaka_bay_mesh.png
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
from fvcom_mesh_tools.dem.subset import to_geotiff  # noqa: E402
from fvcom_mesh_tools.io import load_river_points, read_fort14  # noqa: E402

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
F14_OUT = OUT_DIR / "17_osaka_bay.14"
SUMMARY_TXT = OUT_DIR / "17_osaka_bay_summary.txt"
MESH_PNG = OUT_DIR / "17_osaka_bay_mesh.png"

# Osaka Bay bbox: includes the bay proper, east coast of Awaji,
# Akashi Strait at the NW corner, Tomogashima Strait at the SW corner.
# The DEM bbox edges that cut through Awaji land contain no water and
# therefore produce no open-boundary mesh edges - only the water-filled
# strait sections will be classified as open.
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
    ax.set_title("PoC #17 Osaka Bay")
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
        print(f"[17] subsetting SRTM15+ to {BBOX} -> {DEM_SUBSET}")
        t0 = time.perf_counter()
        to_geotiff(SRTM, DEM_SUBSET, BBOX, src_var="z")
        print(f"[17] subset: {time.perf_counter() - t0:.2f} s")

    args = [
        str(DEM_SUBSET), str(F14_OUT),
        "--hmin", str(HMIN_M), "--hmax", str(HMAX_M), "--zmax", "0.0",
        "--interp-method", "linear",
        "--coastline", str(COASTLINE),
        "--coast-target-size", "200", "--coast-expansion-rate", "0.005",
        "--min-polygon-area-m2", "1000000",
        "--min-island-area-m2", "100000",
        "--open-merge-coast-gap", "50",
        "--river-inflow-points", str(RIVERS),
        "--river-segment-nodes", "5",
        "--river-ibtype", "21",
        "--river-snap-tol-m", "5000",
        "--quality-pass", "6",
        "--refine-min-angle", "20",
        "--land-ibtype", "20",
        "--perpfix-iters", "1",
    ]
    print("[17] running fmesh-buildmesh on Osaka Bay ...")
    t0 = time.perf_counter()
    rc = buildmesh_main(args)
    wall = time.perf_counter() - t0
    if rc != 0:
        raise SystemExit(f"buildmesh exited {rc}")
    print(f"[17] wall: {wall:.2f} s")

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
    print(f"[17] wrote {SUMMARY_TXT}")

    plot_mesh(mesh, river_pts, MESH_PNG)
    print(f"[17] wrote {MESH_PNG}")

    assert flipped == 0
    assert n_river_segs > 0, "no river ibtype=21 segments produced"
    assert n_open > 0, "no open boundary segments produced"


if __name__ == "__main__":
    main()
