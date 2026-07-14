"""Urayasu wall/basin: sample / before / after 3-way (owner
2026-07-15). BEFORE highlights, in red, the offending cells --
centroids inside the enclosed-basin+wall region (the edit_006
fill polygon) that the mesh wrongly treated as open water."""

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
from fvcom_mesh_tools.plotting import add_atlas_grid, use_readable_style

use_readable_style()
ROOT = Path(__file__).resolve().parents[1]
OUT = ROOT / "outputs/sample_repro"
COSW = float(np.cos(np.deg2rad(35.35)))
CX, CY = 139.8300, 35.6320
HW = 0.0100
tr = Transformer.from_crs(32654, 4326, always_xy=True)

BASIN = shapely.Polygon([
    (139.8238, 35.6326), (139.8348, 35.6326),
    (139.8348, 35.6395), (139.8238, 35.6395)])


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
    ("BEFORE = edit rev2 (wall pierced, basin opened)",
     *load14_ll(OUT / "wall_before_demo.14")),
    ("AFTER (rev4: wall kept, basin closed)",
     *load14_ll(OUT / "sample_repro_final.14")),
]
land = unary_union(list(gpd.read_file(
    ROOT / "outputs/tb_varres_3r/land_osm_wide.shp").geometry))
gl = gpd.GeoSeries([land], crs="EPSG:4326")

fig, axes = plt.subplots(1, 3, figsize=(6.8 * 3, 7.4))
for c, (title, po, tri) in enumerate(MESHES):
    ax = axes[c]
    gl.plot(ax=ax, color="0.88", edgecolor="saddlebrown",
            linewidth=0.9)
    ax.triplot(po[:, 0], po[:, 1], tri, color="steelblue",
               linewidth=0.6, alpha=0.9)
    cent = po[tri].mean(axis=1)
    bad = [j for j in range(len(tri))
           if abs(cent[j, 0] - CX) < 2 * HW
           and abs(cent[j, 1] - CY) < 2 * HW
           and BASIN.covers(shapely.Point(cent[j]))]
    for j in bad:
        q = po[tri[j]]
        ax.fill(q[:, 0], q[:, 1], facecolor="crimson",
                alpha=0.25, edgecolor="crimson",
                linewidth=1.6, zorder=5)
    xs, ys = BASIN.exterior.xy
    ax.plot(xs, ys, color="darkorange", linewidth=1.6,
            linestyle=":", zorder=6)
    ax.set_xlim(CX - HW / COSW, CX + HW / COSW)
    ax.set_ylim(CY - HW, CY + HW)
    ax.set_aspect(1.0 / COSW)
    add_atlas_grid(ax, crs="EPSG:4326")
    n_bad = len(bad)
    ax.set_title(f"{title}\ncells inside basin+wall region: "
                 f"{n_bad}", fontsize=13)
fig.tight_layout()
out = ROOT / "outputs/figures/wall_3way.png"
fig.savefig(out, dpi=150, bbox_inches="tight")
print(f"[370] saved {out}", flush=True)
