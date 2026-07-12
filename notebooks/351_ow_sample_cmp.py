# Sample comparison for confirmed one-wide sites (owner request
# 2026-07-12): for every confirmed bank-to-bank cell site, show
# the goto2023 sample and our mesh side by side at the same
# window, and MEASURE the sample's structure at each flagged cell:
# is the site inside the sample's meshed water at all, and how
# many sample cells does the local cross-section pass through?
import json
import os
from pathlib import Path

import numpy as np
os.environ.setdefault("MPLBACKEND", "Agg")
import geopandas as gpd
import matplotlib.pyplot as plt
import shapely
from pyproj import Transformer
from shapely.geometry import LineString
from shapely.ops import unary_union
from shapely.strtree import STRtree

from fvcom_mesh_tools.io import read_fort14
from fvcom_mesh_tools.plotting import (
    add_atlas_grid,
    use_readable_style,
)

use_readable_style()
OUT = Path("outputs/sample_repro")
COSW = float(np.cos(np.deg2rad(35.35)))
SCALE = 0.5 * (111e3 * COSW + 111e3)

tr = Transformer.from_crs("EPSG:32654", "EPSG:4326",
                          always_xy=True)
G = os.path.expanduser('~/Github/TB-FVCOM/input/goto2023/grid/')
gd = open(G + 'TokyoBay_grd.dat').read().split('\n')
nn = int(gd[0].split('=')[1]); ne = int(gd[1].split('=')[1])
Ts = np.array([[int(w) for w in gd[2 + i].split()[1:4]]
               for i in range(ne)]) - 1
Ps = np.array([[float(w) for w in gd[2 + ne + i].split()[1:3]]
               for i in range(nn)])
lon_s, lat_s = tr.transform(Ps[:, 0], Ps[:, 1])
P_s = np.column_stack([lon_s, lat_s])
spolys = [shapely.Polygon(P_s[t]) for t in Ts]
stree = STRtree(spolys)

m = read_fort14(str(OUT / "sample_repro_final.14"))
lon_o, lat_o = tr.transform(m.nodes[:, 0], m.nodes[:, 1])
P_o = np.column_stack([lon_o, lat_o])
T_o = m.elements
cent_o = P_o[T_o].mean(axis=1)

# boundary edges of OUR mesh (for the local channel direction at
# a flagged cell = along its boundary edge)
keys, eids = [], np.tile(np.arange(len(T_o)), 3)
edges = []
for a, b in ((0, 1), (1, 2), (2, 0)):
    lo = np.minimum(T_o[:, a], T_o[:, b]).astype(np.int64)
    hi = np.maximum(T_o[:, a], T_o[:, b]).astype(np.int64)
    keys.append(lo * m.n_nodes + hi)
keys = np.concatenate(keys)
o = np.argsort(keys, kind="stable")
ks = keys[o]
same = np.zeros(len(ks), bool)
same[1:] |= ks[1:] == ks[:-1]
same[:-1] |= ks[1:] == ks[:-1]
bedges = set(ks[~same].tolist())

data = json.loads((OUT / "one_wide_cells.json").read_text())
sites = data["confirmed_sites"]

land = unary_union(list(gpd.read_file(
    "outputs/tb_varres_3r/land_osm_wide.shp").geometry))
gser = gpd.GeoSeries([land], crs="EPSG:4326")


def sample_section_count(e):
    """Cross-section through cell e's centroid, perpendicular to
    its boundary edge (the bank): how many SAMPLE cells does it
    substantially cross, and is the centroid sample-covered?"""
    c = cent_o[e]
    tang = None
    for a, b in ((0, 1), (1, 2), (2, 0)):
        vlo, vhi = sorted((T_o[e][a], T_o[e][b]))
        if vlo * m.n_nodes + vhi in bedges:
            tang = P_o[vhi] - P_o[vlo]
            break
    if tang is None:
        tang = P_o[T_o[e][1]] - P_o[T_o[e][0]]
    tang = tang / (np.hypot(*tang) + 1e-15)
    nv = np.array([-tang[1], tang[0]])
    half = 500.0 / SCALE
    sec = LineString([c - nv * half, c + nv * half])
    fr, tot = [], 0.0
    for j in stree.query(sec):
        seg = spolys[j].intersection(sec)
        if not seg.is_empty and seg.length > 0:
            fr.append(seg.length)
            tot += seg.length
    nsub = sum(1 for ln in fr if tot > 0 and ln / tot > 0.3)
    pt = shapely.Point(c)
    covered = any(spolys[j].covers(pt) for j in stree.query(pt))
    return covered, nsub


rows = [s for s in sites]
fig, axes = plt.subplots(len(rows), 2,
                         figsize=(11.5, 5.6 * len(rows)),
                         squeeze=False)
for r, s in enumerate(rows):
    cx, cy = s["center_lonlat"]
    half = 0.012
    els = np.asarray(s["cell_ids_fort14"]) - 1
    meas = [sample_section_count(e) for e in els]
    verdicts = "; ".join(
        f"{e + 1}: sample {'covers' if cov else 'NOT meshed'}"
        f", {n} cells across" for e, (cov, n) in zip(els, meas))
    print(f"[cmp] {s['site']} [{s['gridref']}]  {verdicts}",
          flush=True)
    for k, (title, lo, la, T, col) in enumerate([
            ("goto2023 sample", lon_s, lat_s, Ts, "0.35"),
            ("ours (red = one-wide)", lon_o, lat_o, T_o,
             "steelblue")]):
        ax = axes[r, k]
        gser.plot(ax=ax, color="0.9", edgecolor="0.6", lw=0.5,
                  zorder=1)
        ax.triplot(lo, la, T, lw=0.6, color=col, zorder=3)
        if k == 1:
            for e in els:
                ring = P_o[np.append(T_o[e], T_o[e][0])]
                ax.fill(ring[:, 0], ring[:, 1], color="crimson",
                        alpha=0.18, zorder=2)
                ax.plot(ring[:, 0], ring[:, 1], color="crimson",
                        lw=1.6, zorder=4)
                ax.annotate(str(e + 1), cent_o[e], ha="center",
                            va="center", fontsize=12,
                            fontweight="bold", color="darkred",
                            zorder=5)
        ax.set_xlim(cx - half / COSW, cx + half / COSW)
        ax.set_ylim(cy - half, cy + half)
        ax.set_aspect(1 / COSW)
        ax.set_xticks([])
        ax.set_yticks([])
        add_atlas_grid(ax, crs="EPSG:4326")
        ax.set_title(f"{s['site']} [{s['gridref']}]  {title}")
fig.suptitle("confirmed one-wide sites vs the goto2023 sample")
fig.tight_layout(rect=(0, 0, 1, 0.97))
fig.savefig("outputs/figures/one_wide_sample_cmp.png", dpi=170,
            bbox_inches="tight")
print("[cmp] saved outputs/figures/one_wide_sample_cmp.png",
      flush=True)
