"""OW05 visual-judgment sheet (owner 2026-07-15): a large zoom
with fort.14 NODE numbers, the pipeline shoreline AND the
original OSM coastline, so the human/AI can see which nodes sit
off the terrain line and where unused water lies."""

from pathlib import Path

import geopandas as gpd
import matplotlib
matplotlib.use("Agg")
import matplotlib.pyplot as plt
import numpy as np
from pyproj import Transformer
from shapely.ops import unary_union

import sys
sys.path.insert(0, str(Path(__file__).resolve().parents[1] / "src"))
from fvcom_mesh_tools.io import read_fort14
from fvcom_mesh_tools.plotting import use_readable_style

use_readable_style()
ROOT = Path(__file__).resolve().parents[1]
mesh = read_fort14(ROOT / "outputs/sample_repro/sample_repro_final.14")
tr = Transformer.from_crs(32654, 4326, always_xy=True)
lon, lat = tr.transform(mesh.nodes[:, 0], mesh.nodes[:, 1])
po = np.column_stack([lon, lat])
tri = mesh.elements
COSW = float(np.cos(np.deg2rad(35.35)))
CX, CY = 139.8395, 35.6310
HW = 0.0062

pipe = gpd.read_file(ROOT / "outputs/sample_repro/land_channel_adj.shp")
orig = gpd.read_file(ROOT / "outputs/tb_varres_3r/land_osm_wide.shp")

fig, ax = plt.subplots(figsize=(16, 15))
gpd.GeoSeries([unary_union(list(orig.geometry))], crs="EPSG:4326").plot(
    ax=ax, color="0.90", edgecolor="saddlebrown", linewidth=1.8,
    zorder=1, label="_")
gpd.GeoSeries([unary_union(list(pipe.geometry))], crs="EPSG:4326").boundary.plot(
    ax=ax, color="green", linewidth=1.6, linestyle="--", zorder=3)
ax.triplot(po[:, 0], po[:, 1], tri, color="steelblue",
           linewidth=1.0, zorder=4)

inw = ((np.abs(po[:, 0] - CX) < HW / COSW * 1.05)
       & (np.abs(po[:, 1] - CY) < HW * 1.05))
for i in np.nonzero(inw)[0]:
    ax.annotate(str(i + 1), po[i], fontsize=11, color="black",
                ha="left", va="bottom", zorder=6,
                xytext=(2, 2), textcoords="offset points",
                bbox=dict(boxstyle="round,pad=0.12",
                          fc="white", ec="none", alpha=0.75))
cent = po[tri].mean(axis=1)
inc = ((np.abs(cent[:, 0] - CX) < HW / COSW)
       & (np.abs(cent[:, 1] - CY) < HW))
for j in np.nonzero(inc)[0]:
    ax.annotate(f"[{j+1}]", cent[j], fontsize=9, color="purple",
                ha="center", va="center", zorder=5, alpha=0.85)
for cid in (735, 736):
    q = po[tri[cid - 1]]
    ax.fill(q[:, 0], q[:, 1], facecolor="crimson", alpha=0.12,
            edgecolor="crimson", linewidth=2.0, zorder=4.5)
ax.set_xlim(CX - HW / COSW, CX + HW / COSW)
ax.set_ylim(CY - HW, CY + HW)
ax.set_aspect(1.0 / COSW)
ax.set_title("OW05 detail: mesh + node ids (black), cell ids "
             "[purple], ORIGINAL OSM coast (brown/grey), "
             "PIPELINE shoreline (green dashed)")
fig.tight_layout()
out = ROOT / "outputs/figures/ow05_detail.png"
fig.savefig(out, dpi=150, bbox_inches="tight")
print(f"[367] saved {out}", flush=True)
