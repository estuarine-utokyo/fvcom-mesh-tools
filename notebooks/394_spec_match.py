"""How closely does the mesh follow the sizing SPEC it was given?

The spec is the input condition, not a preference:

    h_spec(x) = min( MAXEL, max( H0 + GRADE * distance_to_coast,
                                 sqrt(g*H) * DT / CR ) )

where the second term is the CFL floor. This script measures, per
element, achieved / h_spec, so that

    ratio > 1   coarser than the spec asked for  (under-resolved)
    ratio < 1   finer than the spec asked for    (nodes spent for nothing)

Changing H0 is changing the spec; making the ratio 1 is the algorithm's
job (owner 2026-09-20).

Outputs a map, a histogram, band statistics and a JSON summary.

Usage: python notebooks/394_spec_match.py [mesh.14] [label]
"""

from __future__ import annotations

import json
import os
import sys
from pathlib import Path

import geopandas as gpd
import matplotlib

matplotlib.use("Agg")
import matplotlib.pyplot as plt  # noqa: E402
import numpy as np  # noqa: E402
import shapely  # noqa: E402
from pyproj import Transformer  # noqa: E402
from shapely.ops import unary_union  # noqa: E402

sys.path.insert(0, str(Path(__file__).resolve().parents[1] / "src"))
from fvcom_mesh_tools.plotting import add_atlas_grid, use_readable_style  # noqa: E402

use_readable_style()
ROOT = Path(__file__).resolve().parents[1]
FIG = ROOT / "outputs/figures"
OUT = ROOT / "outputs/sample_repro"
COSW = float(np.cos(np.deg2rad(35.35)))
tr = Transformer.from_crs(32654, 4326, always_xy=True)

MESH = Path(sys.argv[1]) if len(sys.argv) > 1 else OUT / "sample_repro_final.14"
LABEL = sys.argv[2] if len(sys.argv) > 2 else "certified"

H0 = float(os.environ.get("SR_H0", 290.0))
MAXEL = float(os.environ.get("SR_MAXEL", 1400.0))
GRADE = float(os.environ.get("SR_GRADE", 0.165))
DT = float(os.environ.get("SR_DT", 15.0))
CR = float(os.environ.get("SR_CRMIN", 0.45))

lines = MESH.read_text().splitlines()
NE, NN = map(int, lines[1].split()[:2])
nod = np.array([lines[2 + i].split()[1:4] for i in range(NN)], float)
tri = np.array([lines[2 + NN + i].split()[2:5] for i in range(NE)], int) - 1
xy_m, depth = nod[:, :2], nod[:, 2]
lon, lat = tr.transform(xy_m[:, 0], xy_m[:, 1])
po = np.column_stack([lon, lat])
cen_m = xy_m[tri].mean(1)
cen = po[tri].mean(1)

# achieved size: mean edge length of the element (metres, UTM)
p = xy_m[tri]
edges = np.stack([np.hypot(*(p[:, 1] - p[:, 0]).T),
                  np.hypot(*(p[:, 2] - p[:, 1]).T),
                  np.hypot(*(p[:, 0] - p[:, 2]).T)], axis=1)
achieved = edges.mean(1)

# spec: distance to the ORIGINAL coast, plus the CFL floor from the depth
land = unary_union(list(gpd.read_file(ROOT / "outputs/tb_varres_3r/land_osm_wide.shp").geometry))
to_utm = Transformer.from_crs(4326, 32654, always_xy=True)
coast = shapely.ops.transform(lambda x, y: to_utm.transform(x, y), land.boundary)
dist = shapely.distance(shapely.points(cen_m[:, 0], cen_m[:, 1]), coast)
h_dist = H0 + GRADE * dist
h_cfl = np.sqrt(9.81 * np.maximum(depth[tri].max(1), 0.1)) * DT / CR
h_spec = np.minimum(MAXEL, np.maximum(h_dist, h_cfl))
ratio = achieved / h_spec
binding = np.where(h_dist >= h_cfl, "distance", "cfl")
binding = np.where(np.minimum(h_dist, h_cfl) >= MAXEL, "maxel", binding)

print(f"[394] {LABEL}: {NN} nodes, {NE} elements; spec H0={H0:.0f} grade={GRADE} "
      f"maxel={MAXEL:.0f} dt={DT:.0f} Cr={CR}", flush=True)
print(f"[394] achieved/spec: median {np.median(ratio):.2f}  p05 {np.percentile(ratio, 5):.2f} "
      f" p95 {np.percentile(ratio, 95):.2f}", flush=True)
summary = {"label": LABEL, "nodes": NN, "elements": NE,
           "ratio_median": float(np.median(ratio)),
           "ratio_p05": float(np.percentile(ratio, 5)),
           "ratio_p95": float(np.percentile(ratio, 95)),
           "spec": {"H0": H0, "grade": GRADE, "maxel": MAXEL, "dt": DT, "cr": CR}}
for name, sel in (("coastal band (<=1 h0 from shore)", dist <= H0),
                  ("shallow (<=10 m)", depth[tri].max(1) <= 10),
                  ("deep (>20 m)", depth[tri].max(1) > 20)):
    if sel.sum():
        print(f"[394]   {name:34s} n={sel.sum():5d}  achieved {np.median(achieved[sel]):6.0f} m"
              f"  spec {np.median(h_spec[sel]):6.0f} m  ratio {np.median(ratio[sel]):.2f}",
              flush=True)
        summary[name] = dict(n=int(sel.sum()), achieved_median=float(np.median(achieved[sel])),
                             spec_median=float(np.median(h_spec[sel])),
                             ratio_median=float(np.median(ratio[sel])))
for b in ("distance", "cfl", "maxel"):
    sel = binding == b
    if sel.sum():
        print(f"[394]   binding={b:9s} n={sel.sum():5d}  ratio median "
              f"{np.median(ratio[sel]):.2f}", flush=True)
        summary[f"binding_{b}"] = dict(n=int(sel.sum()),
                                       ratio_median=float(np.median(ratio[sel])))
(OUT / f"spec_match_{LABEL}.json").write_text(json.dumps(summary, indent=1) + "\n")

gl = gpd.GeoSeries([land], crs="EPSG:4326")
fig, axes = plt.subplots(1, 2, figsize=(21, 15), width_ratios=[1.35, 1])
ax = axes[0]
gl.plot(ax=ax, color="0.9", edgecolor="0.7", linewidth=0.4)
tpc = ax.tripcolor(po[:, 0], po[:, 1], tri, facecolors=ratio, cmap="RdYlBu_r",
                   vmin=0.6, vmax=1.6, edgecolors="none")
ax.triplot(po[:, 0], po[:, 1], tri, color="0.35", linewidth=0.15, alpha=0.5)
ax.set_xlim(po[:, 0].min() - 0.01, po[:, 0].max() + 0.01)
ax.set_ylim(po[:, 1].min() - 0.01, po[:, 1].max() + 0.01)
ax.set_aspect(1.0 / COSW)
add_atlas_grid(ax, crs="EPSG:4326")
fig.colorbar(tpc, ax=ax, shrink=0.6, label="achieved edge / spec size")
ax.set_title(f"{LABEL}: how closely the mesh follows its sizing spec\n"
             f"red = coarser than specified, blue = finer than specified",
             fontsize=14)
ax = axes[1]
for b, c in (("distance", "#2c7bb6"), ("cfl", "#d7191c"), ("maxel", "0.4")):
    sel = binding == b
    if sel.sum():
        ax.hist(ratio[sel], bins=np.linspace(0.4, 2.2, 60), histtype="step", linewidth=2,
                color=c, label=f"{b} binds (n={sel.sum()})")
ax.axvline(1.0, color="0.2", linestyle="--")
ax.set_xlabel("achieved / spec")
ax.set_ylabel("elements")
ax.legend()
ax.set_title("Distribution by binding constraint", fontsize=14)
fig.tight_layout()
out = FIG / f"394_spec_match_{LABEL}.png"
fig.savefig(out, dpi=110, bbox_inches="tight")
print(f"[394] saved {out}", flush=True)
