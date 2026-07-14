"""Sample / pre-campaign / current 3-way zoom at the accepted
residual sites (owner 2026-07-15), at 2x the linear window of the
351 panels (half = 0.024 deg ~ 2.7 km) so the surroundings stay
visible.

Panels per site:
  left   goto2023 sample
  middle before the widen/prune/CFL campaign
         (sample_repro_final_baseline_pre_edit001.14, 2026-07-12)
  right  current final mesh, flagged one-wide cells hatched
"""

import json
import os
from pathlib import Path

import geopandas as gpd
import matplotlib
matplotlib.use("Agg")
import matplotlib.pyplot as plt
import numpy as np
import shapely
from pyproj import Transformer
from shapely.ops import unary_union

import sys
sys.path.insert(0, str(Path(__file__).resolve().parents[1] / "src"))
from fvcom_mesh_tools.gridref import TOKYO_BAY_GRID
from fvcom_mesh_tools.plotting import add_atlas_grid, use_readable_style

use_readable_style()
ROOT = Path(__file__).resolve().parents[1]
OUT = ROOT / "outputs/sample_repro"
FIG = ROOT / "outputs/figures"
COSW = float(np.cos(np.deg2rad(35.35)))
HALF = 0.024          # deg lat; 2x the 351 window (4x area)

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
    ("BEFORE campaign (2026-07-12)",
     *load14_ll(OUT / "sample_repro_final_baseline_pre_edit001.14")),
    ("AFTER (current final)",
     *load14_ll(OUT / "sample_repro_final.14")),
]
land = unary_union(list(gpd.read_file(
    ROOT / "outputs/tb_varres_3r/land_osm_wide.shp").geometry))
gl = gpd.GeoSeries([land], crs="EPSG:4326")

ow = json.loads((OUT / "one_wide_cells.json").read_text())
sites = ow["confirmed_sites"]

fig, axes = plt.subplots(len(sites), 3,
                         figsize=(6.4 * 3, 6.8 * len(sites)))
axes = np.atleast_2d(axes)
for r, s in enumerate(sites):
    cx, cy = s["center_lonlat"]
    for c, (title, po, tri) in enumerate(MESHES):
        ax = axes[r, c]
        gl.plot(ax=ax, color="0.88", edgecolor="0.6",
                linewidth=0.5)
        ax.triplot(po[:, 0], po[:, 1], tri, color="steelblue",
                   linewidth=0.55, alpha=0.9)
        if c == 2:
            cells = np.array(s["cell_ids_fort14"]) - 1
            for cid in cells:
                q = po[tri[cid]]
                ax.fill(q[:, 0], q[:, 1], facecolor="none",
                        edgecolor="crimson", hatch="///",
                        linewidth=1.6, zorder=5)
                cc = q.mean(axis=0)
                ax.annotate(str(cid + 1), cc, color="crimson",
                            fontsize=10, ha="center",
                            va="center", zorder=6,
                            path_effects=None)
        ax.plot(cx, cy, marker="o", ms=26, mfc="none",
                mec="darkorange", mew=2.0, zorder=7)
        ax.set_xlim(cx - HALF / COSW, cx + HALF / COSW)
        ax.set_ylim(cy - HALF, cy + HALF)
        ax.set_aspect(1.0 / COSW)
        add_atlas_grid(ax, crs="EPSG:4326")
        ax.set_title(f"{s['site']} [{s['gridref']}]  {title}",
                     fontsize=14)
fig.tight_layout()
FIG.mkdir(parents=True, exist_ok=True)
fig.savefig(FIG / "residual_before_after.png", dpi=140,
            bbox_inches="tight")
print(f"[366] saved {FIG / 'residual_before_after.png'}",
      flush=True)
