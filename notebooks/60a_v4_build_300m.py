"""PoC #59j — v3 coastline-fidelity build at 100 m nearshore.

The user requires the mesh boundary to reproduce the coastline
faithfully; depth adjustment (2 m floor incl. tidal flats) comes
afterwards and must not move the shoreline. Investigation of the v1
lineage showed the boundary IS the C23 coastline, but resampled by
``om.Shoreline`` at h0 = hmin = 200 m (fidelity p50 48 m / p90 130 m
against C23) and further nibbled by repair deletions.

This PoC rebuilds the raw mesh at **hmin = 100 m** with two more
fidelity levers pulled:

* DEM = the real 30 m survey-derived grid
  ``depth_0030-11+12+13+14+15.nc`` (covers 139.565-140.172E,
  35.10-35.86N) instead of the v1 scenario DEM — and the domain now
  extends EAST to the DEM edge (~140.172E), so the Chiba-side
  140.10E artificial clip of v1 disappears (default decision,
  flagged for review);
* island-pruning threshold scales with h0^2, so 100 m retains 4x
  smaller islands automatically.

Deliverables: raw mesh, v1-vs-v3 fidelity comparison stats (distance
of real-coast boundary nodes to C23, artificial-cut bands excluded),
and paired zoom figures. The full finishing chain (structural ->
UTM -> metric finish -> OBC arc -> QA) runs as a follow-up once the
fidelity numbers justify the cost.
"""

from __future__ import annotations

import json
import os
import time
from pathlib import Path

import numpy as np

from fvcom_mesh_tools.cli.buildmesh import main as buildmesh_main
from fvcom_mesh_tools.io import Fort14Mesh, read_fort14

REPO = Path(__file__).resolve().parents[1]
DATA_DIR = Path(os.environ["DATA_DIR"])  # fail loudly if missing
DEM = DATA_DIR / "geodata" / "bathymetry" / "tokyo_bay" / "depth_0030-11+12+13+14+15.nc"
COASTLINE = REPO / "data" / "coastline" / "tokyo_bay" / "MLIT_C23" / "C23-06_TOKYOBAY.shp"
RIVERS = REPO / "data" / "rivers" / "tokyo_bay" / "tokyo_bay_rivers.csv"

RAW_OUT = REPO / "outputs" / "60a_v4_raw_300m.14"
STATS_JSON = REPO / "outputs" / "60a_fidelity_stats.json"
FIG_DIR = REPO / "outputs" / "figures"
V1_MESH = REPO / "outputs" / "59e_gate_passed.14"

HMIN_M = 300.0
HMAX_M = 5000.0
GRADATION = 0.10
SLOPE_PARAMETER = 20.0
MAX_ITER = 50
# The build region comes from the DEM extent (139.565-140.172E,
# 35.10-35.86N), so the domain automatically reaches the real Chiba
# coast — the v1 Chiba clip at 140.10E was the old pipeline bbox and
# disappears here.

# Zoom windows (lon/lat) for visual comparison: Banzu flat off
# Kisarazu, and the Keihin port area.
ZOOMS = {
    "banzu": (139.87, 35.36, 139.98, 35.44),
    "keihin": (139.75, 35.48, 139.88, 35.56),
}


def _coast_fidelity(mesh: Fort14Mesh, coords: str) -> dict:
    """Distance stats from real-coast land-boundary nodes to C23.

    Artificial-cut bands (southern/western/eastern domain edges) are
    excluded so the numbers measure coastline fidelity only.
    """
    import geopandas as gpd
    import shapely
    from pyproj import Transformer
    from shapely.strtree import STRtree

    if coords == "metric":
        tr = Transformer.from_crs("EPSG:32654", "EPSG:4326", always_xy=True)
        lon, lat = tr.transform(mesh.nodes[:, 0], mesh.nodes[:, 1])
        tr_fwd = None
    else:
        lon, lat = mesh.nodes[:, 0], mesh.nodes[:, 1]
        tr_fwd = Transformer.from_crs("EPSG:4326", "EPSG:32654", always_xy=True)

    land_nodes = np.unique(np.concatenate(
        [np.asarray(s) for _ib, s in mesh.land_boundaries]
    )) if mesh.land_boundaries else np.empty(0, np.int64)
    band = 0.01  # ~1.1 km: artificial-cut exclusion band at domain edges
    keep = (
        (lat[land_nodes] > lat[land_nodes].min() + band)
        & (lon[land_nodes] > lon[land_nodes].min() + band)
        & (lon[land_nodes] < lon[land_nodes].max() - band)
    )
    nodes = land_nodes[keep]

    gdf = gpd.read_file(COASTLINE).set_crs(4326, allow_override=True).to_crs(32654)
    tree = STRtree(list(gdf.geometry))
    if coords == "metric":
        px, py = mesh.nodes[nodes, 0], mesh.nodes[nodes, 1]
    else:
        px, py = tr_fwd.transform(lon[nodes], lat[nodes])
    pts = shapely.points(px, py)
    nearest = tree.nearest(pts)
    d = shapely.distance(pts, np.array(list(gdf.geometry))[nearest])
    return {
        "n_nodes": int(nodes.size),
        "p50_m": float(np.percentile(d, 50)),
        "p90_m": float(np.percentile(d, 90)),
        "p99_m": float(np.percentile(d, 99)),
        "max_m": float(d.max()),
        "n_gt_100m": int((d > 100).sum()),
        "n_gt_200m": int((d > 200).sum()),
    }


def main() -> int:
    t0 = time.perf_counter()
    for p in (DEM, COASTLINE, RIVERS):
        if not p.exists():
            raise SystemExit(f"required input missing: {p}")

    print(f"[60a] building v4 raw mesh: hmin={HMIN_M} m, DEM={DEM.name}", flush=True)
    rc = buildmesh_main([
        str(DEM), str(RAW_OUT),
        "--engine", "oceanmesh",
        "--hmin", str(HMIN_M), "--hmax", str(HMAX_M), "--zmax", "0.0",
        "--interp-method", "linear",
        "--coastline", str(COASTLINE),
        "--om-slope-parameter", str(SLOPE_PARAMETER),
        "--om-gradation", str(GRADATION),
        "--om-max-iter", str(MAX_ITER),
        "--om-seed", "0",
        "--bbox-tol-m", "150",
        "--open-merge-coast-gap", "50",
        "--river-inflow-points", str(RIVERS),
        "--river-segment-nodes", "5",
        "--river-ibtype", "21",
        "--river-snap-tol-m", "5000",
        "--quality-pass", "0",
        "--refine-min-angle", "0",
        "--land-ibtype", "20",
        "--perpfix-iters", "0",
    ])
    build_wall = time.perf_counter() - t0
    print(f"[60a] build wall: {build_wall:.0f} s (exit {rc})", flush=True)
    if rc != 0:
        raise SystemExit(f"buildmesh failed with exit {rc}")

    v3 = read_fort14(RAW_OUT)
    print(f"[60a] v4 raw: NP={v3.n_nodes:,} NE={v3.n_elements:,}", flush=True)
    stats = {
        "v3_raw_100m": {
            "n_nodes": v3.n_nodes, "n_elements": v3.n_elements,
            "build_wall_s": build_wall,
            "fidelity": _coast_fidelity(v3, "lonlat"),
        },
    }
    if V1_MESH.exists():
        v1 = read_fort14(V1_MESH)
        stats["v1_59e_200m"] = {
            "n_nodes": v1.n_nodes, "n_elements": v1.n_elements,
            "fidelity": _coast_fidelity(v1, "metric"),
        }
    STATS_JSON.write_text(json.dumps(stats, indent=2), encoding="utf-8")
    print(json.dumps(stats, indent=2), flush=True)

    # Paired zoom figures (v3 lonlat; v1 already has UTM figures).
    os.environ.setdefault("MPLBACKEND", "Agg")
    from fvcom_mesh_tools.plotting import plot_mesh_overview

    for name, zb in ZOOMS.items():
        plot_mesh_overview(
            v3, FIG_DIR / f"60a_v4_{name}.png", crs="EPSG:4326",
            cell_m=None, coast=(139.4, 35.05, 140.4, 35.85),
            zoom=zb, dpi=300,
            title=f"v4 raw 100 m — {name}",
        )
        if V1_MESH.exists():
            from pyproj import Transformer

            tr = Transformer.from_crs("EPSG:4326", "EPSG:32654", always_xy=True)
            x0, y0 = tr.transform(zb[0], zb[1])
            x1, y1 = tr.transform(zb[2], zb[3])
            plot_mesh_overview(
                v1, FIG_DIR / f"60a_v1_{name}.png", crs="EPSG:32654",
                cell_m=None, coast=(139.4, 35.05, 140.4, 35.85),
                zoom=(x0, y0, x1, y1), dpi=300,
                title=f"v1 (200 m lineage) — {name}",
            )
    print(f"[60a] wall total: {time.perf_counter() - t0:.0f} s", flush=True)
    return 0


if __name__ == "__main__":
    raise SystemExit(main())
