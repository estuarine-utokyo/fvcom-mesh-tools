"""Sample / before / after 3-way zooms at the RESOLVED one-wide
sites (owner 2026-07-15). Same 2x window as 366; each row is one
fixed site with the fixing mechanism in the title."""

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
from fvcom_mesh_tools.gridref import TOKYO_BAY_GRID
from fvcom_mesh_tools.plotting import add_atlas_grid, use_readable_style

use_readable_style()
ROOT = Path(__file__).resolve().parents[1]
OUT = ROOT / "outputs/sample_repro"
COSW = float(np.cos(np.deg2rad(35.35)))
HALF = 0.024
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
    ("BEFORE (2026-07-12)",
     *load14_ll(OUT / "sample_repro_final_baseline_pre_edit001.14")),
    ("AFTER (current)",
     *load14_ll(OUT / "sample_repro_final.14")),
]
land = unary_union(list(gpd.read_file(
    ROOT / "outputs/tb_varres_3r/land_osm_wide.shp").geometry))
gl = gpd.GeoSeries([land], crs="EPSG:4326")

SITES = [
    ("OW05", "I6-d2", 139.8395, 35.6310,
     "data carve 700 m x 2 km (edit_005)"),
    ("OW13", "I6-c3", 139.8241, 35.6292,
     "water-fringe widen + split"),
    ("OW10", "F9-c5", 139.6736, 35.4555,
     "water-fringe widen + split"),
    ("OW23", "H8-c5", 139.7763, 35.5077,
     "gate-first widen + split"),
    ("OW04", "H8-a5", 139.7506, 35.5082,
     "off-line midpoint widen + split"),
    ("OW16", "G8-e4", 139.7475, 35.5140,
     "thin-short-stub closure (do-not-mesh)"),
]

fig, axes = plt.subplots(len(SITES), 3,
                         figsize=(6.4 * 3, 6.9 * len(SITES)))
for r, (sid, gref, cx, cy, how) in enumerate(SITES):
    for c, (title, po, tri) in enumerate(MESHES):
        ax = axes[r, c]
        gl.plot(ax=ax, color="0.88", edgecolor="0.6",
                linewidth=0.5)
        ax.triplot(po[:, 0], po[:, 1], tri, color="steelblue",
                   linewidth=0.55, alpha=0.9)
        ax.plot(cx, cy, marker="o", ms=26, mfc="none",
                mec="darkorange", mew=2.0, zorder=7)
        ax.set_xlim(cx - HALF / COSW, cx + HALF / COSW)
        ax.set_ylim(cy - HALF, cy + HALF)
        ax.set_aspect(1.0 / COSW)
        add_atlas_grid(ax, crs="EPSG:4326")
        ttl = f"{sid} [{gref}]  {title}"
        if c == 2:
            ttl += f"\nfixed by: {how}"
        ax.set_title(ttl, fontsize=13)
fig.tight_layout()
out = ROOT / "outputs/figures/fixed_sites_3way.png"
fig.savefig(out, dpi=130, bbox_inches="tight")
print(f"[368] saved {out}", flush=True)
