"""Measure, and map, how far the mesh land boundary sits from the coastline.

Owner observation 2026-09-21: several stretches of the mesh shoreline are
visibly off the OSM coastline even though a human could place them by eye.
This script turns that observation into a number and a picture.

Two references are measured, because they answer different questions:

  normalized land   the polygon oceanmesh was actually given
                    (land_channel_adj.shp: OSM land after the recipe edits,
                    the normalize pass and the waterway carving).  The gap
                    against THIS reference is pure mesh-generation error and
                    is what a fitting pass can remove.
  raw OSM land      land_osm_wide.shp, before any of our edits.  The extra
                    gap against this reference is what our own preprocessing
                    deliberately changed and must NOT be "fixed" blindly.

The sign convention is `+ = the mesh node is short of the coast (in water)`,
`- = the node has crossed onto land`.

Usage:
  python notebooks/404_coast_fit.py [mesh.14] [dir holding land_channel_adj.shp]
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
from fvcom_mesh_tools.plotting import add_atlas_grid, use_readable_style  # noqa: E402

use_readable_style()
ROOT = Path(__file__).resolve().parents[1]
OUT = ROOT / "outputs/sample_repro"
FIG = ROOT / "outputs/figures"
FIG.mkdir(parents=True, exist_ok=True)
COSW = float(np.cos(np.deg2rad(35.35)))

MESH = Path(sys.argv[1]).resolve() if len(sys.argv) > 1 else OUT / "sample_repro_final.14"
SHPDIR = Path(sys.argv[2]).resolve() if len(sys.argv) > 2 else OUT
TAG = MESH.stem

lines = MESH.read_text().splitlines()
NE, NN = map(int, lines[1].split()[:2])
xy = np.array([lines[2 + i].split()[1:3] for i in range(NN)], float)
tri = np.array([lines[2 + NN + i].split()[2:5] for i in range(NE)], int) - 1
tr = Transformer.from_crs(32654, 4326, always_xy=True)
lon, lat = tr.transform(xy[:, 0], xy[:, 1])
po = np.column_stack([lon, lat])

# Boundary edges (used by exactly one element).  The open boundary is the
# straight OBC arc and is an INPUT, so it is excluded from the coast fit.
e = np.sort(np.vstack([tri[:, [0, 1]], tri[:, [1, 2]], tri[:, [2, 0]]]), axis=1)
uniq, cnt = np.unique(e, axis=0, return_counts=True)
bnd = uniq[cnt == 1]

nopen = int(lines[2 + NN + NE].split()[0])
obc: set[int] = set()
row = 2 + NN + NE + 2
for _ in range(nopen):
    n = int(lines[row].split()[0])
    obc.update(int(lines[row + 1 + k].split()[0]) - 1 for k in range(n))
    row += 1 + n
land_edges = np.array([b for b in bnd if not (b[0] in obc and b[1] in obc)])
nod = np.unique(land_edges.ravel())
print(f"[404] {MESH.name}: {len(bnd)} boundary edges, "
      f"{len(land_edges)} land boundary edges, {len(nod)} land boundary nodes",
      flush=True)

# Local target size per boundary node: the mean length of its boundary edges.
elen = np.linalg.norm(xy[land_edges[:, 0]] - xy[land_edges[:, 1]], axis=1)
hsum = np.zeros(NN)
hcnt = np.zeros(NN)
for (a, b), L in zip(land_edges, elen):
    hsum[[a, b]] += L
    hcnt[[a, b]] += 1
hloc = hsum[nod] / np.maximum(hcnt[nod], 1)

REFS = [("normalized land (given to oceanmesh)", SHPDIR / "land_channel_adj.shp"),
        ("raw OSM land", ROOT / "outputs/tb_varres_3r/land_osm_wide.shp")]
signed = {}
for tag, shp in REFS:
    if not shp.exists():
        print(f"[404] missing {shp}", flush=True)
        continue
    g = gpd.read_file(shp).to_crs(32654)
    poly = unary_union(list(g.geometry))
    P = shapely.points(xy[nod, 0], xy[nod, 1])
    d = shapely.distance(P, poly.boundary)
    s = np.where(shapely.contains(poly, P), -1.0, 1.0) * d
    signed[tag] = s
    r = np.abs(s) / np.maximum(hloc, 1.0)
    print(f"[404] vs {tag}", flush=True)
    print(f"[404]   |offset|  median {np.median(np.abs(s)):6.1f}  "
          f"p90 {np.percentile(np.abs(s), 90):6.1f}  "
          f"p99 {np.percentile(np.abs(s), 99):6.1f}  max {np.abs(s).max():7.1f} m",
          flush=True)
    print(f"[404]   signed    p10 {np.percentile(s, 10):7.1f}  median {np.median(s):6.1f}  "
          f"p90 {np.percentile(s, 90):6.1f} m   on land {int((s < 0).sum())}"
          f" ({100 * (s < 0).mean():.0f}%)", flush=True)
    for f in (0.25, 0.5, 1.0):
        print(f"[404]   nodes |offset| > {f:4.2f} * local h : "
              f"{int((r > f).sum()):5d} / {len(nod)}", flush=True)

REF = REFS[0][0]
if REF not in signed:
    sys.exit("[404] the normalized land shapefile is required for the map")
off = signed[REF]
ratio = np.abs(off) / np.maximum(hloc, 1.0)

# Sites worth a look: a node more than half a local edge length off the
# coastline it was meant to follow.
THRESH = 0.5
flag = np.where(ratio > THRESH)[0]
print(f"[404] {len(flag)} nodes beyond {THRESH} * local h", flush=True)

land_g = gpd.read_file(REFS[0][1])
osm_g = gpd.read_file(REFS[1][1])
gl = gpd.GeoSeries([unary_union(list(land_g.geometry))], crs=land_g.crs).to_crs(4326)
go = gpd.GeoSeries([unary_union(list(osm_g.geometry))], crs=osm_g.crs).to_crs(4326)

HALF = 0.008
windows = []
used = np.zeros(len(flag), bool)
order = np.argsort(-ratio[flag])
for k in order:
    if used[k]:
        continue
    c = po[nod[flag[k]]]
    members = np.where((~used)
                       & (np.abs(po[nod[flag]][:, 1] - c[1]) < HALF)
                       & (np.abs(po[nod[flag]][:, 0] - c[0]) * COSW < HALF))[0]
    used[members] = True
    windows.append(dict(centre=np.array([po[nod[flag[members]]][:, 0].mean(),
                                         po[nod[flag[members]]][:, 1].mean()]),
                        members=flag[members],
                        worst=float(np.abs(off[flag[members]]).max())))
windows.sort(key=lambda w: -w["worst"])
windows = windows[:12]
print(f"[404] {len(windows)} zoom windows", flush=True)

BND = bnd
RED = "#d7191c"


def base(ax, xlim, ylim, lw, show_osm):
    gl.plot(ax=ax, color="0.88", edgecolor="0.55", linewidth=0.6, zorder=1)
    if show_osm:
        go.boundary.plot(ax=ax, color="#2c7bb6", linewidth=1.4, zorder=2,
                         linestyle=(0, (4, 2)))
    ax.triplot(po[:, 0], po[:, 1], tri, color="0.5", linewidth=lw, alpha=0.9, zorder=3)
    seg = po[BND]
    ax.plot(np.c_[seg[:, 0, 0], seg[:, 1, 0], np.full(len(seg), np.nan)].ravel(),
            np.c_[seg[:, 0, 1], seg[:, 1, 1], np.full(len(seg), np.nan)].ravel(),
            color="0.2", linewidth=0.7 if lw <= 0.5 else 2.0, zorder=4)
    ax.set_xlim(*xlim)
    ax.set_ylim(*ylim)
    ax.set_aspect(1.0 / COSW)
    add_atlas_grid(ax, crs="EPSG:4326")


fig, ax = plt.subplots(figsize=(12.5, 15.5))
base(ax, (po[:, 0].min() - 0.01, po[:, 0].max() + 0.01),
     (po[:, 1].min() - 0.01, po[:, 1].max() + 0.01), 0.3, False)
for n, w in enumerate(windows, 1):
    cx, cy = w["centre"]
    ax.add_patch(Rectangle((cx - HALF / COSW, cy - HALF), 2 * HALF / COSW, 2 * HALF,
                           fill=False, edgecolor=RED, linewidth=2.2, zorder=7))
    ax.text(cx, cy + HALF + 0.004, str(n), color=RED, fontsize=16, fontweight="bold",
            ha="center", va="bottom", zorder=8)
ax.set_title(
    f"Mesh shoreline vs the coastline it was given ({MESH.name})\n"
    f"|offset| median {np.median(np.abs(off)):.0f} m, p90 {np.percentile(np.abs(off), 90):.0f} m, "
    f"max {np.abs(off).max():.0f} m  |  {len(flag)} of {len(nod)} boundary nodes "
    f"beyond {THRESH:g} x local edge length", fontsize=14)
fig.tight_layout()
p = FIG / f"404_coast_fit_map_{TAG}.png"
fig.savefig(p, dpi=130, bbox_inches="tight")
print(f"[404] saved {p}", flush=True)

ncol = 3
nrow = int(np.ceil(max(len(windows), 1) / ncol))
fig, axes = plt.subplots(nrow, ncol, figsize=(6.4 * ncol, 6.8 * nrow), squeeze=False)
for n, (ax, w) in enumerate(zip(np.ravel(axes), windows), 1):
    cx, cy = w["centre"]
    base(ax, (cx - HALF / COSW, cx + HALF / COSW), (cy - HALF, cy + HALF), 1.0, True)
    for j in w["members"]:
        i = nod[j]
        ax.plot([po[i, 0]], [po[i, 1]], marker="o", markersize=7,
                color=RED, zorder=9)
    # Every land boundary edge touching a flagged node, in red.
    sel = np.array([(a in set(nod[w["members"]])) or (b in set(nod[w["members"]]))
                    for a, b in land_edges])
    for a, b in land_edges[sel]:
        ax.plot(po[[a, b], 0], po[[a, b], 1], color=RED, linewidth=3.0, zorder=8)
    ax.set_title(f"[{n}] worst offset {w['worst']:.0f} m, {len(w['members'])} nodes\n"
                 f"({cx:.4f}, {cy:.4f})  blue dashed = raw OSM coastline",
                 fontsize=12, color=RED)
for ax in np.ravel(axes)[len(windows):]:
    ax.set_axis_off()
fig.tight_layout()
p = FIG / f"404_coast_fit_zooms_{TAG}.png"
fig.savefig(p, dpi=120, bbox_inches="tight")
print(f"[404] saved {p}", flush=True)
