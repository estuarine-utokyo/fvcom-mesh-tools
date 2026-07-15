"""Whole-domain 3-way mesh comparison (owner 2026-07-16): the
full-bay view of goto2023 sample | manual-edit baseline (3dfa621)
| certified normalization config (run 6210307). Companion to the
per-site sheet 376."""

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
from fvcom_mesh_tools.plotting import add_atlas_grid, use_readable_style

use_readable_style()
ROOT = Path(__file__).resolve().parents[1]
OUT = ROOT / "outputs/sample_repro"
COSW = float(np.cos(np.deg2rad(35.35)))
tr = Transformer.from_crs(32654, 4326, always_xy=True)


def load14_ll(p):
    ls = Path(p).read_text().splitlines()
    ne, nn = map(int, ls[1].split()[:2])
    nod = np.array([ls[2 + i].split()[1:3] for i in range(nn)],
                   float)
    tri = np.array([ls[2 + nn + i].split()[2:5]
                    for i in range(ne)], int) - 1
    lon, lat = tr.transform(nod[:, 0], nod[:, 1])
    return np.column_stack([lon, lat]), tri


MESHES = [
    ("goto2023 sample", *load14_ll(OUT / "sample_original.14")),
    ("BEFORE automation (3dfa621, manual edits)",
     *load14_ll(OUT / "sample_repro_final_A_3dfa621.14")),
    ("AFTER automation (normalize cert, run 6210307)",
     *load14_ll(OUT / "sample_repro_final.14")),
]
land = unary_union(list(gpd.read_file(
    ROOT / "outputs/tb_varres_3r/land_osm_wide.shp").geometry))
gl = gpd.GeoSeries([land], crs="EPSG:4326")

# common frame: union of all three meshes plus a small margin
x0 = min(po[:, 0].min() for _, po, _ in MESHES) - 0.01
x1 = max(po[:, 0].max() for _, po, _ in MESHES) + 0.01
y0 = min(po[:, 1].min() for _, po, _ in MESHES) - 0.01
y1 = max(po[:, 1].max() for _, po, _ in MESHES) + 0.01

h_in = 13.5
w_in = h_in * ((x1 - x0) * COSW / (y1 - y0)) + 1.2
fig, axes = plt.subplots(1, 3, figsize=(w_in * 3, h_in + 0.8))
for c, (title, po, tri) in enumerate(MESHES):
    ax = axes[c]
    gl.plot(ax=ax, color="0.90", edgecolor="0.65", linewidth=0.4)
    ax.triplot(po[:, 0], po[:, 1], tri, color="steelblue",
               linewidth=0.35, alpha=0.95)
    ax.set_xlim(x0, x1)
    ax.set_ylim(y0, y1)
    ax.set_aspect(1.0 / COSW)
    add_atlas_grid(ax, crs="EPSG:4326")
    ax.set_title(f"{title}\nNP={len(po):,}  NE={len(tri):,}",
                 fontsize=14)
fig.tight_layout()
out = ROOT / "outputs/figures/full_domain_3way.png"
fig.savefig(out, dpi=150, bbox_inches="tight")
print(f"[378] saved {out}", flush=True)
