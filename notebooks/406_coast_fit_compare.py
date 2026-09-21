"""Show what the coastline fit actually moved, site by site.

Draws the same windows twice -- the mesh before the fit and after it -- over
the coastline the generator was given, so the gain is visible rather than
only tabulated.

Usage:
  python notebooks/406_coast_fit_compare.py before.14 after.14 <dir with land_channel_adj.shp> [tag]
"""

from __future__ import annotations

import sys
from pathlib import Path

import geopandas as gpd
import matplotlib

matplotlib.use("Agg")
import matplotlib.pyplot as plt  # noqa: E402
import numpy as np  # noqa: E402
import shapely  # noqa: E402
from matplotlib.patches import Rectangle  # noqa: E402
from pyproj import Transformer  # noqa: E402
from shapely.ops import unary_union  # noqa: E402

sys.path.insert(0, str(Path(__file__).resolve().parents[1] / "src"))
from fvcom_mesh_tools.coast_fit import boundary_edges  # noqa: E402
from fvcom_mesh_tools.plotting import add_atlas_grid, use_readable_style  # noqa: E402

use_readable_style()
ROOT = Path(__file__).resolve().parents[1]
FIG = ROOT / "outputs/figures"
FIG.mkdir(parents=True, exist_ok=True)
COSW = float(np.cos(np.deg2rad(35.35)))
tr = Transformer.from_crs(32654, 4326, always_xy=True)

BEFORE, AFTER, SHPDIR = Path(sys.argv[1]), Path(sys.argv[2]), Path(sys.argv[3])
TAG = sys.argv[4] if len(sys.argv) > 4 else "coast_fit"


def load(path: Path):
    lines = path.read_text().splitlines()
    ne, nn = map(int, lines[1].split()[:2])
    xy = np.array([lines[2 + i].split()[1:3] for i in range(nn)], float)
    tri = np.array([lines[2 + nn + i].split()[2:5] for i in range(ne)], int) - 1
    lo, la = tr.transform(xy[:, 0], xy[:, 1])
    return xy, np.column_stack([lo, la]), tri


xy_b, po_b, tri_b = load(BEFORE)
xy_a, po_a, tri_a = load(AFTER)
if tri_b.shape != tri_a.shape or not (tri_b == tri_a).all():
    print("[406] WARNING: the two meshes do not share a triangulation; "
          "the fit is supposed to leave it untouched", flush=True)

land = unary_union(list(gpd.read_file(SHPDIR / "land_channel_adj.shp")
                        .to_crs(32654).geometry))
gl = gpd.GeoSeries([land], crs=32654).to_crs(4326)

be = boundary_edges(tri_b)
nod = np.unique(be.ravel())
d_b = shapely.distance(shapely.points(xy_b[nod, 0], xy_b[nod, 1]), land.boundary)
moved = np.linalg.norm(xy_a[nod] - xy_b[nod], axis=1)
print(f"[406] {len(nod)} boundary nodes, {int((moved > 1.0).sum())} moved more than 1 m, "
      f"largest move {moved.max():.1f} m", flush=True)

# Windows around the nodes the fit moved the furthest.
HALF = 0.008
order = np.argsort(-moved)
windows = []
used = np.zeros(len(nod), bool)
for k in order:
    if used[k] or moved[k] < 1.0:
        continue
    c = po_b[nod[k]]
    members = np.where((~used)
                       & (np.abs(po_b[nod][:, 1] - c[1]) < HALF)
                       & (np.abs(po_b[nod][:, 0] - c[0]) * COSW < HALF))[0]
    used[members] = True
    windows.append(dict(centre=c, best=float(moved[k]),
                        n=int((moved[members] > 1.0).sum()),
                        gain=float(d_b[members].mean())))
    if len(windows) >= 9:
        break
print(f"[406] {len(windows)} windows", flush=True)

BEF, AFT = "#d7191c", "#2c7bb6"


def draw(ax, po, tri, colour, xlim, ylim, label):
    gl.plot(ax=ax, color="0.88", edgecolor="0.45", linewidth=1.2, zorder=1)
    ax.triplot(po[:, 0], po[:, 1], tri, color="0.55", linewidth=0.8, alpha=0.9, zorder=3)
    b = boundary_edges(tri)
    seg = po[b]
    ax.plot(np.c_[seg[:, 0, 0], seg[:, 1, 0], np.full(len(seg), np.nan)].ravel(),
            np.c_[seg[:, 0, 1], seg[:, 1, 1], np.full(len(seg), np.nan)].ravel(),
            color=colour, linewidth=2.4, zorder=5, label=label)
    ax.set_xlim(*xlim)
    ax.set_ylim(*ylim)
    ax.set_aspect(1.0 / COSW)
    add_atlas_grid(ax, crs="EPSG:4326")


nrow = len(windows)
fig, axes = plt.subplots(nrow, 2, figsize=(13.0, 6.0 * nrow), squeeze=False)
for n, (w, row) in enumerate(zip(windows, axes), 1):
    cx, cy = w["centre"]
    xlim = (cx - HALF / COSW, cx + HALF / COSW)
    ylim = (cy - HALF, cy + HALF)
    draw(row[0], po_b, tri_b, BEF, xlim, ylim, "before")
    draw(row[1], po_a, tri_a, AFT, xlim, ylim, "after")
    row[0].set_title(f"[{n}] before the fit  ({cx:.4f}, {cy:.4f})\n"
                     f"mean offset here {w['gain']:.0f} m", fontsize=12, color=BEF)
    row[1].set_title(f"[{n}] after the fit  |  {w['n']} nodes moved, "
                     f"largest {w['best']:.0f} m", fontsize=12, color=AFT)
fig.tight_layout()
p = FIG / f"406_coast_fit_compare_{TAG}.png"
fig.savefig(p, dpi=120, bbox_inches="tight")
print(f"[406] saved {p}", flush=True)

# Overview of where the fit did its work.
fig, ax = plt.subplots(figsize=(12.5, 15.5))
gl.plot(ax=ax, color="0.9", edgecolor="0.6", linewidth=0.5)
ax.triplot(po_a[:, 0], po_a[:, 1], tri_a, color="0.55", linewidth=0.3, alpha=0.9)
seg = po_a[boundary_edges(tri_a)]
ax.plot(np.c_[seg[:, 0, 0], seg[:, 1, 0], np.full(len(seg), np.nan)].ravel(),
        np.c_[seg[:, 0, 1], seg[:, 1, 1], np.full(len(seg), np.nan)].ravel(),
        color="0.2", linewidth=0.7, zorder=4)
for n, w in enumerate(windows, 1):
    cx, cy = w["centre"]
    ax.add_patch(Rectangle((cx - HALF / COSW, cy - HALF), 2 * HALF / COSW, 2 * HALF,
                           fill=False, edgecolor=AFT, linewidth=2.2, zorder=7))
    ax.text(cx, cy + HALF + 0.004, str(n), color=AFT, fontsize=16, fontweight="bold",
            ha="center", va="bottom", zorder=8)
ax.set_xlim(po_a[:, 0].min() - 0.01, po_a[:, 0].max() + 0.01)
ax.set_ylim(po_a[:, 1].min() - 0.01, po_a[:, 1].max() + 0.01)
ax.set_aspect(1.0 / COSW)
add_atlas_grid(ax, crs="EPSG:4326")
ax.set_title(f"Mesh after the coastline fit ({AFTER.name})\n"
             f"{int((moved > 1.0).sum())} of {len(nod)} boundary nodes moved, "
             f"largest {moved.max():.0f} m  |  boxes = the windows of "
             f"406_coast_fit_compare_{TAG}.png", fontsize=14)
fig.tight_layout()
p = FIG / f"406_coast_fit_overview_{TAG}.png"
fig.savefig(p, dpi=130, bbox_inches="tight")
print(f"[406] saved {p}", flush=True)
