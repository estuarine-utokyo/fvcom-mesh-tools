# PoC #111: reconstruct the v5u sizing field exactly as the engine
# composes it and locate the fine-size source that defeats the
# 1.2 km interest/OBC floors at the Miura junction (observed mesh
# p50 366 m there instead of ~1.2 km).
import json
import os
from pathlib import Path

import numpy as np

os.environ.setdefault("MPLBACKEND", "Agg")

import oceanmesh as om  # noqa: E402

REPO = Path(__file__).resolve().parents[1]
OUT = REPO / "outputs" / "pipeline_v5u"
PREP = OUT / "prep"

import yaml  # noqa: E402

recipe = yaml.safe_load((REPO / "recipes" / "tokyo_bay_v5u.yaml").read_text())
bcfg = recipe["build"]
dem_path = os.path.expandvars(bcfg["dem"])
assert "$" not in dem_path, f"unresolved env var in {bcfg['dem']}"

bbox = recipe["prep"]["bbox"]
region_bbox = (bbox[0], bbox[2], bbox[1], bbox[3])  # om order lon lon lat lat
region = om.Region(extent=region_bbox, crs=4326)

hmin_m = float(bcfg["hmin_m"])
DEG = 1.0 / 111194.9266
hmin_deg = hmin_m * DEG
hmax_deg = float(bcfg["hmax_m"]) * DEG

shp = PREP / "land_opened.shp"
h0_deg = float(bcfg.get("shoreline_h0_m", hmin_m)) * DEG
shore = om.Shoreline(str(shp), region.bbox, h0_deg, crs=4326)
sdf = om.signed_distance_function(shore)
edge_feat = om.feature_sizing_function(
    shore, sdf, max_edge_length=hmax_deg, crs=4326,
)
n_below = int((edge_feat.values < hmin_deg).sum())
edge_feat.values = np.maximum(edge_feat.values, hmin_deg)
print(f"feature floor: raised {n_below}")

combined = edge_feat

# OBC coarsening ramp (recipe prep.obc_line -- same as build stage)
import shapely  # noqa: E402
from shapely.geometry import LineString, Polygon  # noqa: E402

obc_line = recipe["prep"]["obc_line"]
arc = LineString([(q[0], q[1]) for q in obc_line])
xg, yg = combined.create_grid()
print("grid shapes:", xg.shape, yg.shape, combined.values.shape)
pts = shapely.points(xg.ravel(), yg.ravel())
d_deg = shapely.distance(pts, arc).reshape(xg.shape)
d_m = d_deg / DEG
size_deg = float(bcfg.get("obc_coarsen_size_m", 1600.0)) * DEG
R = float(bcfg.get("obc_coarsen_radius_m", 10000.0))
ramp = size_deg * np.clip(1.0 - d_m / R, 0.0, 1.0)
print("obc ramp raised:", int((combined.values < ramp).sum()))
combined.values = np.maximum(combined.values, ramp)

poly = Polygon([(q[0], q[1]) for q in recipe["build"]["interest_region"]])
d_out = shapely.distance(pts, poly).reshape(xg.shape) / DEG
floor2 = (float(bcfg["outside_min_m"]) * DEG) * np.clip(
    d_out / float(bcfg["outside_blend_m"]), 0.0, 1.0)
print("interest floor raised:", int((combined.values < floor2).sum()))
combined.values = np.maximum(combined.values, floor2)

pre_grad = combined.values.copy()
edge = om.enforce_mesh_gradation(combined, gradation=float(bcfg.get("gradation", 0.15)))

# --- evaluate at probes ---------------------------------------------------
probes = {
    "junction_zigzag": (139.677, 35.117),
    "junction_coast": (139.670, 35.130),
    "arc_mid": (139.706, 35.055),
    "narrows_uraga": (139.72, 35.24),
    "banzu": (139.92, 35.39),
}
for name, (px, py) in probes.items():
    h = float(np.asarray(edge.eval(np.array([[px, py]]))).ravel()[0]) / DEG
    print(f"h({name}) = {h:.0f} m")

# min sizing within 3 km of the junction: WHERE is the fine source?
jx, jy = probes["junction_zigzag"]
dj = np.hypot(xg - jx, yg - jy) / DEG
mask = dj < 3000.0
vals = edge.values.copy()
vv = np.where(mask, vals, np.inf)
k = np.unravel_index(np.argmin(vv), vv.shape)
print(f"min h within 3 km of junction: {vals[k]/DEG:.0f} m at "
      f"({xg[k]:.4f},{yg[k]:.4f}); pre-gradation there: {pre_grad[k]/DEG:.0f} m")

# heatmap around the mouth
import matplotlib.pyplot as plt  # noqa: E402

fig, axes = plt.subplots(1, 2, figsize=(16, 8))
for ax, (title, V) in zip(
    axes,
    [("pre-gradation", pre_grad), ("final edge", edge.values)],
):
    sel = (xg > 139.60) & (xg < 139.90) & (yg > 34.98) & (yg < 35.30)
    im = ax.scatter(xg[sel], yg[sel], c=V[sel] / DEG, s=1,
                    cmap="viridis", vmin=200, vmax=1800)
    ax.plot(*arc.xy, "r-", lw=1.5)
    ax.set_title(title)
    ax.set_aspect(1 / np.cos(np.deg2rad(35.1)))
    plt.colorbar(im, ax=ax, label="h (m)")
fig.savefig(REPO / "outputs" / "figures" / "111_sizing_mouth.png",
            dpi=150, bbox_inches="tight")
print("saved outputs/figures/111_sizing_mouth.png")
