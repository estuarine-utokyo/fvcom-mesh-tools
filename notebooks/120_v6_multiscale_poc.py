# PoC #120: Tokyo Bay via OM2D-parity MULTISCALE oceanmesh (port
# acceptance test). Two nests, per-nest coastline detail (same
# engineered shapefile read at different h0 — the OM2D geodata
# pattern), varres_3r-parity sizing, finalize_sizing with the
# two-sided CFL limiter. NO simplify_outside_region anywhere: the
# outer nest's coarseness comes from h0=1 km + coarsen_outside.
import json
import os
from pathlib import Path

import numpy as np

os.environ.setdefault("MPLBACKEND", "Agg")

import geopandas as gpd  # noqa: E402
import oceanmesh as om  # noqa: E402
from oceanmesh import Region  # noqa: E402

REPO = Path(__file__).resolve().parents[1]
OUT = REPO / "outputs" / "pipeline_v6"
OUT.mkdir(parents=True, exist_ok=True)
PREP = OUT / "prep"
PREP.mkdir(exist_ok=True)

import yaml  # noqa: E402

recipe = yaml.safe_load((REPO / "recipes" / "tokyo_bay_v5u.yaml").read_text())
BBOX = tuple(recipe["prep"]["bbox"])  # lon_min, lat_min, lon_max, lat_max
OBC_LINE = [tuple(q) for q in recipe["prep"]["obc_line"]]
# Interest polygon override (PoC #121 diagnosis): the v5u polygon's
# EAST side dips to 35.02-35.18, putting the Futtsu bank / east
# mouth waters inside the fine nest — the sample's fine/coarse
# transition sits at ~35.20-35.25 (band p50: 541 m at 35.20-35.30
# vs 1146 m at 35.10-35.20). South edge raised to the narrows on
# both sides; single-variable change vs PoC #120.
INTEREST = [
    (139.58, 35.32),
    (139.58, 35.74),
    (140.16, 35.74),
    (140.16, 35.32),
    (139.88, 35.25),
    (139.80, 35.10),
]
DEM_PATH = os.path.expandvars(recipe["build"]["dem"])
assert "$" not in DEM_PATH

# --- prep: engineered shoreline WITHOUT outside simplification -------
from fvcom_mesh_tools.prep import fetch_true_land, open_land  # noqa: E402
from fvcom_mesh_tools.prep.shoreline import (  # noqa: E402
    cut_domain_at_obc_line,
    extend_obc_ends_perpendicular,
)

shp_path = PREP / "land_opened.shp"
if not shp_path.exists():
    print("[prep] fetching true land ...", flush=True)
    land = fetch_true_land(
        BBOX, min_water_area_deg2=float(
            recipe["prep"].get("min_water_area", 1e-6))
    )
    opened = open_land(
        land,
        r_open_m=float(recipe["prep"].get("r_open_m", 150.0)),
        min_island_area_m2=float(recipe["prep"].get("min_island", 3.6e5)),
        clip_bbox=BBOX,
    )
    eff = extend_obc_ends_perpendicular(OBC_LINE, opened)
    opened = cut_domain_at_obc_line(opened, eff, BBOX)
    opened.to_file(shp_path)
    (PREP / "obc_line_effective.json").write_text(json.dumps(eff))
    print(f"[prep] wrote {shp_path} ({len(opened)} polygons)", flush=True)
else:
    eff = json.loads((PREP / "obc_line_effective.json").read_text())
    print("[prep] reusing existing shoreline", flush=True)

DEG = 1.0 / 111194.9266

# --- nest definitions (varres_3r parity, adapted to walled domain) ---
# Outer nest: whole domain, coarse geodata (h0 = 1 km), sizing
# min 1000 / max 2000, grade 0.1, wl 30.
# Inner nest: interest polygon, detailed geodata (h0 = 100 m),
# sizing min 400 / max 500 (sample interior p50 ~540), grade 0.1,
# feature r=3.
region_outer = Region((BBOX[0], BBOX[2], BBOX[1], BBOX[3]), 4326)
poly_inner = np.asarray(INTEREST + [INTEREST[0]], dtype=float)

print("[nest1] outer Shoreline h0=1000 m ...", flush=True)
shore_o = om.Shoreline(str(shp_path), region_outer.bbox, 1000.0 * DEG)
sdf_o = om.signed_distance_function(shore_o)
dem = om.DEM(DEM_PATH, bbox=region_outer)

edge_o_feat = om.feature_sizing_function(
    shore_o, sdf_o, r=3, max_edge_length=2000.0 * DEG)
edge_o_wl = om.wavelength_sizing_function(
    dem, wl=30, min_edgelength=1000.0 * DEG,
    max_edge_length=2000.0 * DEG)
edge_o, dt_o = om.finalize_sizing(
    [edge_o_feat, edge_o_wl], dem=dem,
    hmin=1000.0 * DEG, max_edge_length=2000.0,
    gradation=0.10,
    courant={"timestep": 16.0, "max": 0.5},
)
print(f"[nest1] finalize done (dt={dt_o})", flush=True)

print("[nest2] inner Shoreline h0=100 m, polygon boubox ...", flush=True)
shore_i = om.Shoreline(str(shp_path), poly_inner, 100.0 * DEG)
sdf_i = om.signed_distance_function(shore_i)
edge_i_feat = om.feature_sizing_function(
    shore_i, sdf_i, r=3, max_edge_length=500.0 * DEG)
edge_i, dt_i = om.finalize_sizing(
    [edge_i_feat], dem=dem,
    hmin=400.0, max_edge_length=500.0,
    gradation=0.10,
    courant={"timestep": 16.0, "max": 0.5},
)
print(f"[nest2] finalize done (dt={dt_i})", flush=True)

print("[multiscale] generating ...", flush=True)
points, cells = om.generate_multiscale_mesh(
    [sdf_o, sdf_i], [edge_o, edge_i], seed=0,
)
print(f"[multiscale] raw NP={len(points):,} NE={len(cells):,}", flush=True)

# cleanup (OM2D msh.clean-equivalent, packaged)
points, cells = om.make_mesh_boundaries_traversable(points, cells)
points, cells = om.delete_faces_connected_to_one_face(points, cells)
points, cells = om.laplacian2(points, cells)
from oceanmesh.mesh_improve import area_length_quality  # noqa: E402

q = area_length_quality(points, cells)
print(f"[clean] NP={len(points):,} NE={len(cells):,} "
      f"AL-qual min/mean = {q.min():.3f}/{q.mean():.3f}", flush=True)

# depth + fort.14 with auto boundaries
b = om.interp_bathymetry(points, cells, dem, method="cell-averaging",
                         min_depth=1.0)
bc = om.make_bc_auto(points, cells, depth=b, classifier="depth",
                     depth_lim=5.0, cut_lim=6)
f14 = OUT / "tokyo_bay_v6_raw.14"
om.write_fort14(str(f14), points, cells, depth=b, boundaries=bc)
print(f"[out] {f14}", flush=True)

# --- A/B metrics vs the goto2023 sample -------------------------------
def edge_stats(pts, tris, label):
    e = np.vstack([tris[:, [0, 1]], tris[:, [1, 2]], tris[:, [2, 0]]])
    e = np.unique(np.sort(e, axis=1), axis=0)
    L = np.linalg.norm(pts[e[:, 0]] - pts[e[:, 1]], axis=1) / DEG
    lat = 0.5 * (pts[e[:, 0], 1] + pts[e[:, 1], 1])
    print(f"=== {label}: edges by lat band (m) ===", flush=True)
    for lo, hi in [(34.95, 35.10), (35.10, 35.20), (35.20, 35.30),
                   (35.30, 35.40), (35.40, 35.55), (35.55, 35.70)]:
        sel = (lat >= lo) & (lat < hi)
        if sel.sum():
            print(f"  {lo:.2f}-{hi:.2f}: n={int(sel.sum()):5d} "
                  f"p10/p50/p90 = "
                  f"{np.percentile(L[sel], [10, 50, 90]).round(0)}",
                  flush=True)

edge_stats(points, cells, "v6 multiscale")

# figures
import matplotlib.pyplot as plt  # noqa: E402

for name, zoom in {
    "full": None,
    "mouth": (139.60, 34.95, 140.05, 35.25),
    "yokosuka": (139.60, 35.15, 139.78, 35.32),
    "tokyoport": (139.72, 35.53, 139.90, 35.68),
}.items():
    fig, ax = plt.subplots(figsize=(11, 11))
    ax.triplot(points[:, 0], points[:, 1], cells, lw=0.3,
               color="steelblue")
    if zoom:
        ax.set_xlim(zoom[0], zoom[2])
        ax.set_ylim(zoom[1], zoom[3])
    ax.set_aspect(1 / np.cos(np.deg2rad(35.4)))
    ax.set_title(f"tokyo_bay_v6 multiscale raw - {name}")
    fig.savefig(OUT / f"v6_raw_{name}.png", dpi=180,
                bbox_inches="tight")
    plt.close(fig)
print("[out] figures done", flush=True)
