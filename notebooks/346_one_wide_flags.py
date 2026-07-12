# One-wide channel cell reporter (owner spec 2026-07-12): where a
# waterway is carried by cells spanning bank to bank (all three
# nodes on the land boundary), mark each such cell IN COLOUR WITH
# ITS CELL NUMBER so the next manual-editing round can address
# them by ID. This is a REPORT, not a fixer -- sub-cell-width
# channels meshed 1-wide are accepted by design (minimum mesh size
# is inviolable; banks are only widened by explicit arc edits).
import json
import os
from pathlib import Path

import numpy as np
os.environ.setdefault("MPLBACKEND", "Agg")
import geopandas as gpd
import matplotlib.pyplot as plt
from pyproj import Transformer
from scipy.sparse import coo_matrix
from scipy.sparse.csgraph import connected_components
from shapely.ops import unary_union

from fvcom_mesh_tools.io import read_fort14
from fvcom_mesh_tools.plotting import use_readable_style

use_readable_style()
OUT = Path("outputs/sample_repro")
COSW = float(np.cos(np.deg2rad(35.35)))

m = read_fort14(str(OUT / "sample_repro_final.14"))
tr = Transformer.from_crs("EPSG:32654", "EPSG:4326", always_xy=True)
lon, lat = tr.transform(m.nodes[:, 0], m.nodes[:, 1])
P = np.column_stack([lon, lat])
T = m.elements

# boundary nodes = nodes of edges used by exactly one element
keys, eids = [], np.tile(np.arange(len(T)), 3)
for a, b in ((0, 1), (1, 2), (2, 0)):
    lo = np.minimum(T[:, a], T[:, b]).astype(np.int64)
    hi = np.maximum(T[:, a], T[:, b]).astype(np.int64)
    keys.append(lo * m.n_nodes + hi)
keys = np.concatenate(keys)
o = np.argsort(keys, kind="stable")
ks, es = keys[o], eids[o]
same = np.zeros(len(ks), bool)
same[1:] |= ks[1:] == ks[:-1]
same[:-1] |= ks[1:] == ks[:-1]
bkeys = ks[~same]
bnode = np.zeros(m.n_nodes, bool)
bnode[(bkeys // m.n_nodes)] = True
bnode[(bkeys % m.n_nodes)] = True

obc = np.zeros(m.n_nodes, bool)
for ob in m.open_boundaries:
    obc[np.asarray(ob) - 1] = True

allb = bnode[T].all(axis=1)
touches_obc = obc[T].any(axis=1)
flag = np.where(allb & ~touches_obc)[0]
print(f"[1wide] bank-to-bank cells (all 3 nodes on land "
      f"boundary, no OBC): {len(flag)}", flush=True)

# cluster flagged cells by node-sharing so nearby cells report as
# one site
lab = -np.ones(len(T), int)
if len(flag):
    nid = {}
    rows, cols = [], []
    for fi, e in enumerate(flag):
        for v in T[e]:
            nid.setdefault(v, []).append(fi)
    for lst in nid.values():
        for a2 in lst[1:]:
            rows.append(lst[0])
            cols.append(a2)
    g = coo_matrix((np.ones(len(rows)), (rows, cols)),
                   shape=(len(flag), len(flag)))
    _, cl = connected_components(g + g.T, directed=False)
    lab[flag] = cl

cent = P[T].mean(axis=1)
sites = []
for c in sorted(set(lab[flag])) if len(flag) else []:
    els = flag[lab[flag] == c]
    cc = cent[els].mean(axis=0)
    sites.append({
        "site": f"OW{len(sites) + 1:02d}",
        "cell_ids_fort14": (els + 1).tolist(),
        "center_lonlat": [round(float(cc[0]), 4),
                          round(float(cc[1]), 4)],
    })
(OUT / "one_wide_cells.json").write_text(json.dumps(sites,
                                                    indent=1))
for s in sites:
    print(f"[1wide] {s['site']}: cells {s['cell_ids_fort14']} at "
          f"({s['center_lonlat'][0]}, {s['center_lonlat'][1]})",
          flush=True)

land = unary_union(list(gpd.read_file(
    "outputs/tb_varres_3r/land_osm_wide.shp").geometry))
gser = gpd.GeoSeries([land], crs="EPSG:4326")
show = sites[:12]
if show:
    ncol = min(3, len(show))
    nrow = int(np.ceil(len(show) / ncol))
    fig, axes = plt.subplots(nrow, ncol,
                             figsize=(5.6 * ncol, 5.6 * nrow),
                             squeeze=False)
    axes = axes.ravel()
    for ax, s in zip(axes, show):
        cx, cy = s["center_lonlat"]
        half = 0.012
        gser.plot(ax=ax, color="0.9", edgecolor="0.6", lw=0.5)
        ax.triplot(lon, lat, T, lw=0.5, color="steelblue")
        els = np.asarray(s["cell_ids_fort14"]) - 1
        ax.tripcolor(lon, lat, T,
                     facecolors=np.isin(np.arange(len(T)), els)
                     .astype(float), cmap="Reds", alpha=0.55,
                     vmin=0.0, vmax=1.3)
        for e in els:
            ax.annotate(str(e + 1), cent[e],
                        ha="center", va="center",
                        fontsize=12, fontweight="bold",
                        color="darkred")
        ax.set_xlim(cx - half / COSW, cx + half / COSW)
        ax.set_ylim(cy - half, cy + half)
        ax.set_aspect(1 / COSW)
        ax.set_xticks([])
        ax.set_yticks([])
        ax.set_title(f"{s['site']}  {len(els)} cell(s)  "
                     f"({cx:.3f}, {cy:.3f})")
    for ax in axes[len(show):]:
        ax.axis("off")
    fig.suptitle("one-wide channel cells (red + fort.14 cell "
                 "number): candidates for the next manual edit")
    fig.tight_layout(rect=(0, 0, 1, 0.95))
    fig.savefig("outputs/figures/one_wide_cells.png", dpi=170,
                bbox_inches="tight")
    print("[1wide] saved outputs/figures/one_wide_cells.png",
          flush=True)
else:
    print("[1wide] no one-wide cells -- no figure", flush=True)
