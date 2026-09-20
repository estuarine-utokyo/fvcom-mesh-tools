"""Show the meshes of the coastal-resolution sweep side by side.

Columns: certified | A (grade 0.12) | B (H0 250) | C (both).
Rows: whole domain, then zooms at representative sites.

Mesh boundary edges (the mesh's own shoreline) are heavier than interior
edges; the atlas grid is dark grey dashed; land is light grey.

Usage: python notebooks/393_sizing_compare.py [outputs/sizing_392]
"""

from __future__ import annotations

import sys
from pathlib import Path

import geopandas as gpd
import matplotlib

matplotlib.use("Agg")
import matplotlib.pyplot as plt  # noqa: E402
import numpy as np  # noqa: E402
from pyproj import Transformer  # noqa: E402
from shapely.ops import unary_union  # noqa: E402

sys.path.insert(0, str(Path(__file__).resolve().parents[1] / "src"))
from fvcom_mesh_tools.plotting import add_atlas_grid, use_readable_style  # noqa: E402

use_readable_style()
ROOT = Path(__file__).resolve().parents[1]
SWEEP = Path(sys.argv[1] if len(sys.argv) > 1 else ROOT / "outputs/sizing_392")
FIG = ROOT / "outputs/figures"
COSW = float(np.cos(np.deg2rad(35.35)))
tr = Transformer.from_crs(32654, 4326, always_xy=True)


def load(path: Path):
    lines = Path(path).read_text().splitlines()
    ne, nn = map(int, lines[1].split()[:2])
    nod = np.array([lines[2 + i].split()[1:3] for i in range(nn)], float)
    tri = np.array([lines[2 + nn + i].split()[2:5] for i in range(ne)], int) - 1
    lon, lat = tr.transform(nod[:, 0], nod[:, 1])
    e = np.vstack([tri[:, [0, 1]], tri[:, [1, 2]], tri[:, [2, 0]]])
    e.sort(axis=1)
    uniq, cnt = np.unique(e, axis=0, return_counts=True)
    return np.column_stack([lon, lat]), tri, uniq[cnt == 1], nn, ne


CASES = [("certified", ROOT / "outputs/sample_repro/sample_repro_final.14")]
CASES += [(d.name, d / "sample_repro_final.14") for d in sorted(SWEEP.iterdir())
          if d.is_dir() and (d / "sample_repro_final.14").exists()]
CASES = [(t, p) for t, p in CASES if Path(p).exists()]
meshes = [(t, *load(p)) for t, p in CASES]

land = unary_union(list(gpd.read_file(ROOT / "outputs/tb_varres_3r/land_osm_wide.shp").geometry))
gl = gpd.GeoSeries([land], crs="EPSG:4326")

SITES = [("whole domain", 139.87, 35.35, 0.42),
         ("Keihin canals", 139.6800, 35.4550, 0.020),
         ("Urayasu / OW05", 139.8300, 35.6300, 0.020),
         ("Kisarazu shallows", 139.8600, 35.4000, 0.030)]

fig, axes = plt.subplots(len(SITES), len(meshes),
                         figsize=(6.0 * len(meshes), 6.6 * len(SITES)))
for r, (sid, cx, cy, half) in enumerate(SITES):
    for c, (title, po, tri, bnd, nn, ne) in enumerate(meshes):
        ax = axes[r, c]
        gl.plot(ax=ax, color="0.88", edgecolor="0.7", linewidth=0.4)
        lw = 0.25 if r == 0 else 0.9
        ax.triplot(po[:, 0], po[:, 1], tri, color="0.5", linewidth=lw, alpha=0.9)
        seg = po[bnd]
        ax.plot(np.c_[seg[:, 0, 0], seg[:, 1, 0], np.full(len(seg), np.nan)].ravel(),
                np.c_[seg[:, 0, 1], seg[:, 1, 1], np.full(len(seg), np.nan)].ravel(),
                color="0.2", linewidth=0.6 if r == 0 else 1.8, zorder=3)
        ax.set_xlim(cx - half / COSW, cx + half / COSW)
        ax.set_ylim(cy - half, cy + half)
        ax.set_aspect(1.0 / COSW)
        add_atlas_grid(ax, crs="EPSG:4326")
        ax.set_title(f"{sid}\n{title}" + (f"   NP={nn:,} NE={ne:,}" if r == 0 else ""),
                     fontsize=12)
fig.tight_layout()
out = FIG / f"393_mesh_compare_{SWEEP.name}.png"
fig.savefig(out, dpi=110, bbox_inches="tight")
print(f"[393] saved {out}", flush=True)
