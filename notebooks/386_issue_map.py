"""Map every flagged site of the current mesh, prominently.

Owner request 2026-09-20: the overview must show the zoom WINDOWS as
boxes (not markers), and each zoom must colour the offending mesh --
red edges over a light fill -- so the reader sees which elements are
meant.

Inputs, all written by the comparators on the same mesh:
  outputs/sample_repro/over_resolution.json   364: stray elements, wall edges
  outputs/sample_repro/land_breaches.json     342: elements on original land
  outputs/sample_repro/one_wide_cells.json    346: one-wide + advisory cells

Usage: python notebooks/386_issue_map.py
"""

from __future__ import annotations

import json
import sys
from pathlib import Path

import geopandas as gpd
import matplotlib

matplotlib.use("Agg")
import matplotlib.pyplot as plt  # noqa: E402
import numpy as np  # noqa: E402
from matplotlib.patches import Rectangle  # noqa: E402
from pyproj import Transformer  # noqa: E402
from shapely.ops import unary_union  # noqa: E402

sys.path.insert(0, str(Path(__file__).resolve().parents[1] / "src"))
from fvcom_mesh_tools.plotting import add_atlas_grid, use_readable_style  # noqa: E402

use_readable_style()
ROOT = Path(__file__).resolve().parents[1]
OUT = ROOT / "outputs/sample_repro"
FIG = ROOT / "outputs/figures"
COSW = float(np.cos(np.deg2rad(35.35)))
tr = Transformer.from_crs(32654, 4326, always_xy=True)

MESH = OUT / "sample_repro_final.14"
lines = MESH.read_text().splitlines()
NE, NN = map(int, lines[1].split()[:2])
xy = np.array([lines[2 + i].split()[1:3] for i in range(NN)], float)
tri = np.array([lines[2 + NN + i].split()[2:5] for i in range(NE)], int) - 1
lon, lat = tr.transform(xy[:, 0], xy[:, 1])
po = np.column_stack([lon, lat])
cen = po[tri].mean(1)

ovr = json.loads((OUT / "over_resolution.json").read_text())
brc = json.loads((OUT / "land_breaches.json").read_text())
ow = json.loads((OUT / "one_wide_cells.json").read_text())

# (class, element ids, edge node pairs, position)
items = []
for rec in ow["confirmed_sites"]:
    ids = [i - 1 for i in rec["cell_ids_fort14"]]
    items.append(("1WIDE", ids, [], cen[ids].mean(0), f"{rec['site']} [{rec['gridref']}]"))
for e, (i, u) in enumerate(zip(brc["elements"], brc["unintended"])):
    items.append(("LAND-UNINTENDED" if u else "LAND-INTENDED", [i], [], cen[i],
                  "on original land"))
for i in ovr.get("severe_elements", []):
    items.append(("NARROW-SEVERE", [i], [], cen[i], "water narrower than 0.5 h"))
for i in ovr.get("narrow_elements", ovr.get("stray_elements", [])):
    if i in set(ovr.get("severe_elements", [])):
        continue
    items.append(("NARROW", [i], [], cen[i], "water narrower than one row"))
for a, b in ovr["wall_edges"]:
    items.append(("WALL", [], [(a, b)], po[[a, b]].mean(0), "crosses a thin land wall"))

STYLE = {
    "WALL": dict(color="#d7191c", label="WALL crossing (364 gate fail)"),
    "NARROW-SEVERE": dict(color="#e66101", label="water meshed at w/h<0.5 (364 gate fail)"),
    "NARROW": dict(color="#fdb863", label="water meshed at w/h<1 (364 advisory)"),
    "1WIDE": dict(color="#7b3294", label="confirmed one-element-wide (346)"),
    "LAND-UNINTENDED": dict(color="#2c7bb6", label="element on land, UNINTENDED (342)"),
    "LAND-INTENDED": dict(color="#92c5de", label="element on land, intended widening (342)"),
}
counts = {k: sum(1 for it in items if it[0] == k) for k in STYLE}
print("[386] " + "  ".join(f"{k}={v}" for k, v in counts.items()), flush=True)

# Cluster the items that matter into zoom windows, worst class first.
PRIORITY = ["1WIDE", "WALL", "NARROW-SEVERE", "LAND-UNINTENDED", "NARROW"]
HALF = 0.010  # degrees of latitude; the window is 2*HALF tall
windows = []
used = [False] * len(items)
for cls in PRIORITY:
    idx = [i for i, it in enumerate(items) if it[0] == cls]
    for i in idx:
        if used[i]:
            continue
        c = items[i][3]
        members = [j for j, it in enumerate(items)
                   if not used[j] and abs(it[3][1] - c[1]) < HALF
                   and abs(it[3][0] - c[0]) * COSW < HALF]
        for j in members:
            used[j] = True
        centre = np.array([np.mean([items[j][3][0] for j in members]),
                           np.mean([items[j][3][1] for j in members])])
        classes = sorted({items[j][0] for j in members},
                         key=lambda k: PRIORITY.index(k) if k in PRIORITY else 9)
        windows.append(dict(centre=centre, members=members, lead=cls, classes=classes))
windows = windows[:12]
print(f"[386] {len(windows)} zoom windows", flush=True)

land = unary_union(list(gpd.read_file(ROOT / "outputs/tb_varres_3r/land_osm_wide.shp").geometry))
gl = gpd.GeoSeries([land], crs="EPSG:4326")


# Mesh boundary edges (an edge used by exactly one element): the mesh's own
# shoreline. Drawn heavier than the interior so the reader can see which mesh
# edges follow the coast (owner 2026-09-20); modest in the overview so the
# mesh itself stays visible.
_e = np.vstack([tri[:, [0, 1]], tri[:, [1, 2]], tri[:, [2, 0]]])
_e.sort(axis=1)
_uniq, _cnt = np.unique(_e, axis=0, return_counts=True)
BND = _uniq[_cnt == 1]


def base(ax, xlim, ylim, lw):
    gl.plot(ax=ax, color="0.88", edgecolor="0.7", linewidth=0.4)
    ax.triplot(po[:, 0], po[:, 1], tri, color="0.5", linewidth=lw, alpha=0.9)
    seg = po[BND]
    ax.plot(np.c_[seg[:, 0, 0], seg[:, 1, 0], np.full(len(seg), np.nan)].ravel(),
            np.c_[seg[:, 0, 1], seg[:, 1, 1], np.full(len(seg), np.nan)].ravel(),
            color="0.2", linewidth=0.7 if lw <= 0.5 else 2.0, zorder=3)
    ax.set_xlim(*xlim)
    ax.set_ylim(*ylim)
    ax.set_aspect(1.0 / COSW)
    add_atlas_grid(ax, crs="EPSG:4326")


def highlight(ax, member_ids):
    """Fill the flagged elements and draw their edges in the class colour."""
    for j in member_ids:
        cls, ids, edges, _pos, _lbl = items[j]
        col = STYLE[cls]["color"]
        for i in ids:
            poly = po[tri[i]]
            ax.fill(poly[:, 0], poly[:, 1], color=col, alpha=0.30, zorder=5)
            ax.plot(np.append(poly[:, 0], poly[0, 0]), np.append(poly[:, 1], poly[0, 1]),
                    color=col, linewidth=2.6, zorder=6)
        for a, b in edges:
            ax.plot(po[[a, b], 0], po[[a, b], 1], color=col, linewidth=3.2, zorder=6)


# --- overview: boxes, numbered, no markers -------------------------------
fig, ax = plt.subplots(figsize=(12.5, 15.5))
base(ax, (po[:, 0].min() - 0.01, po[:, 0].max() + 0.01),
     (po[:, 1].min() - 0.01, po[:, 1].max() + 0.01), 0.3)
for n, w in enumerate(windows, 1):
    cx, cy = w["centre"]
    col = STYLE[w["lead"]]["color"]
    ax.add_patch(Rectangle((cx - HALF / COSW, cy - HALF), 2 * HALF / COSW, 2 * HALF,
                           fill=False, edgecolor=col, linewidth=2.2, zorder=7))
    # Plain text: a label box would hide the mesh underneath (owner 2026-09-20).
    ax.text(cx, cy + HALF + 0.004, str(n), color=col, fontsize=16, fontweight="bold",
            ha="center", va="bottom", zorder=8)
handles = [plt.Line2D([], [], color=st["color"], lw=3, label=st["label"])
           for st in STYLE.values()]
ax.legend(handles=handles, loc="lower right", fontsize=11, framealpha=0.95,
          title="box colour = worst class inside")
ax.set_title(
    "Flagged sites of the certified mesh (sample_repro_final.14)\n"
    f"boxes 1-{len(windows)} are the zoom windows of 386_issue_zooms.png  |  "
    + "  ".join(f"{k} {v}" for k, v in counts.items()),
    fontsize=14,
)
fig.tight_layout()
p = FIG / "386_issue_map.png"
fig.savefig(p, dpi=130, bbox_inches="tight")
print(f"[386] saved {p}", flush=True)

# --- zooms: same windows, flagged elements coloured ----------------------
ncol = 3
nrow = int(np.ceil(len(windows) / ncol))
fig, axes = plt.subplots(nrow, ncol, figsize=(6.4 * ncol, 6.8 * nrow))
for n, (ax, w) in enumerate(zip(np.ravel(axes), windows), 1):
    cx, cy = w["centre"]
    base(ax, (cx - HALF / COSW, cx + HALF / COSW), (cy - HALF, cy + HALF), 1.0)
    highlight(ax, w["members"])
    per = {c: sum(1 for j in w["members"] if items[j][0] == c) for c in w["classes"]}
    ax.set_title(f"[{n}] " + ", ".join(f"{c} x{v}" for c, v in per.items())
                 + f"\n({cx:.4f}, {cy:.4f})",
                 fontsize=12, color=STYLE[w["lead"]]["color"])
for ax in np.ravel(axes)[len(windows):]:
    ax.set_axis_off()
fig.tight_layout()
p = FIG / "386_issue_zooms.png"
fig.savefig(p, dpi=120, bbox_inches="tight")
print(f"[386] saved {p}", flush=True)
