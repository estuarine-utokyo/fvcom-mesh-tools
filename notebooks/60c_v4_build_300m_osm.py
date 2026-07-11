"""PoC #60c — v4 build at hmin 300 m with the OSM coastline (new default).

User decision 2026-07-04: OSM (osmdata land polygons minus inland
water, via the xcoast cache; data date 2026-03) replaces MLIT C23 as
the DEFAULT coastline source — newer, and worldwide-applicable. The
shoreline fed to ``om.Shoreline`` is
``outputs/osm_shoreline/osm_true_land_tokyo_bay.shp`` (60b clip).

Resolution policy per the same message: do NOT chase resolution —
build a rational coarse mesh (hmin 300 m) and fix the two observed
defects downstream (60d): fragmented sub-grid river channels
(excluded by the w/h < 2 criterion, rivers become mouth inflows) and
boundary/OBC conformity (snap-to-polyline constraints).
"""

from __future__ import annotations

import json
import os
import time
from pathlib import Path

import numpy as np

from fvcom_mesh_tools.cli.buildmesh import main as buildmesh_main
from fvcom_mesh_tools.io import read_fort14

REPO = Path(__file__).resolve().parents[1]
DATA_DIR = Path(os.environ["DATA_DIR"])
DEM = DATA_DIR / "geodata" / "bathymetry" / "tokyo_bay" / "depth_0030-11+12+13+14+15.nc"
COASTLINE = REPO / "outputs" / "osm_shoreline" / "osm_true_land_tokyo_bay.shp"
RIVERS = REPO / "data" / "rivers" / "tokyo_bay" / "tokyo_bay_rivers.csv"

RAW_OUT = REPO / "outputs" / "60c_v4_raw_300m_osm.14"
STATS_JSON = REPO / "outputs" / "60c_fidelity_stats.json"

HMIN_M = 300.0
HMAX_M = 5000.0
GRADATION = 0.10
SLOPE_PARAMETER = 20.0
MAX_ITER = 50


def _fidelity_vs_osm(mesh) -> dict:
    """Distance from land-boundary nodes to the OSM land-polygon
    BOUNDARY lines (polygon distance is 0 inside — boundaries only),
    artificial-cut bands excluded."""
    import geopandas as gpd
    import shapely
    from pyproj import Transformer
    from shapely.strtree import STRtree

    lon, lat = mesh.nodes[:, 0], mesh.nodes[:, 1]
    land_nodes = np.unique(np.concatenate(
        [np.asarray(s) for _ib, s in mesh.land_boundaries]
    ))
    band = 0.01
    keep = (
        (lat[land_nodes] > lat[land_nodes].min() + band)
        & (lon[land_nodes] > lon[land_nodes].min() + band)
        & (lon[land_nodes] < lon[land_nodes].max() - band)
    )
    nodes = land_nodes[keep]
    gdf = gpd.read_file(COASTLINE).to_crs(32654)
    lines = [g.boundary for g in gdf.geometry if g is not None and not g.is_empty]
    tree = STRtree(lines)
    tr = Transformer.from_crs("EPSG:4326", "EPSG:32654", always_xy=True)
    px, py = tr.transform(lon[nodes], lat[nodes])
    pts = shapely.points(px, py)
    d = shapely.distance(pts, np.array(lines, dtype=object)[tree.nearest(pts)])
    return {
        "n_nodes": int(nodes.size),
        "p50_m": float(np.percentile(d, 50)),
        "p90_m": float(np.percentile(d, 90)),
        "n_gt_150m": int((d > 150).sum()),
    }


def main() -> int:
    t0 = time.perf_counter()
    for p in (DEM, COASTLINE, RIVERS):
        if not p.exists():
            raise SystemExit(f"required input missing: {p}")
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
        "--bbox-tol-m", "450",
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
    print(f"[60c] build wall: {time.perf_counter() - t0:.0f} s (exit {rc})",
          flush=True)
    if rc != 0:
        raise SystemExit(rc)
    mesh = read_fort14(RAW_OUT)
    stats = {
        "v4_raw_300m_osm": {
            "n_nodes": mesh.n_nodes,
            "n_elements": mesh.n_elements,
            "fidelity_vs_osm": _fidelity_vs_osm(mesh),
        },
    }
    STATS_JSON.write_text(json.dumps(stats, indent=2), encoding="utf-8")
    print(json.dumps(stats, indent=2), flush=True)
    return 0


if __name__ == "__main__":
    raise SystemExit(main())
