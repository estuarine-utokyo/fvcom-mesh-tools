# Site-exact verification for edit_001 (Haneda D-runway channel,
# W10): sample vs baseline (pre-edit) vs edited mesh, zoomed at
# the owner-reported coordinates. Complements the 342 comparator
# gate -- a fix claim needs BOTH the gate PASS and this zoom.
import os
from pathlib import Path

import numpy as np
os.environ.setdefault("MPLBACKEND", "Agg")
import geopandas as gpd
import matplotlib.pyplot as plt
from pyproj import Transformer
from shapely.ops import unary_union

from fvcom_mesh_tools.io import read_fort14

OUT = Path("outputs/sample_repro")
CX, CY = 139.8043, 35.5259          # W10 center (owner site)
HALF = 0.022
COSW = float(np.cos(np.deg2rad(35.35)))

tr = Transformer.from_crs("EPSG:32654", "EPSG:4326", always_xy=True)
G = os.path.expanduser('~/Github/TB-FVCOM/input/goto2023/grid/')
gd = open(G + 'TokyoBay_grd.dat').read().split('\n')
nn = int(gd[0].split('=')[1]); ne = int(gd[1].split('=')[1])
Ts = np.array([[int(w) for w in gd[2 + i].split()[1:4]]
               for i in range(ne)]) - 1
Ps = np.array([[float(w) for w in gd[2 + ne + i].split()[1:3]]
               for i in range(nn)])
lon_s, lat_s = tr.transform(Ps[:, 0], Ps[:, 1])

land = unary_union(list(gpd.read_file(
    "outputs/tb_varres_3r/land_osm_wide.shp").geometry))
gser = gpd.GeoSeries([land], crs="EPSG:4326")

panels = [("goto2023 sample", None, lon_s, lat_s, Ts, "0.35")]
for title, path in [
        ("baseline (pre-edit)",
         OUT / "sample_repro_final_baseline_pre_edit001.14"),
        ("with edit_001 arc band",
         OUT / "sample_repro_final.14")]:
    m = read_fort14(str(path))
    lo, la = tr.transform(m.nodes[:, 0], m.nodes[:, 1])
    panels.append((f"{title}  NP={m.n_nodes}", path, lo, la,
                   m.elements, "steelblue"))

fig, axes = plt.subplots(1, 3, figsize=(15, 5.6))
for ax, (title, _, lo, la, T, col) in zip(axes, panels):
    gser.plot(ax=ax, color="0.9", edgecolor="0.6", lw=0.5)
    ax.triplot(lo, la, T, lw=0.5, color=col)
    ax.set_xlim(CX - HALF / COSW, CX + HALF / COSW)
    ax.set_ylim(CY - HALF, CY + HALF)
    ax.set_aspect(1 / COSW)
    ax.plot([CX], [CY], marker="+", ms=14, mew=2, color="red")
    ax.set_title(title, fontsize=10)
fig.suptitle("W10 Haneda D-runway channel: site-exact "
             f"verification at ({CX}, {CY})", fontsize=12)
fig.tight_layout(rect=(0, 0, 1, 0.95))
fig.savefig("outputs/figures/haneda_edit001_verify.png", dpi=190,
            bbox_inches="tight")
print("[verify] saved outputs/figures/haneda_edit001_verify.png",
      flush=True)

# MEASURED through-passage check (connectivity-verification rule:
# never claim a fix from figures). Subgraph = our elements whose
# centroids lie within the edit corridor (+width margin); the
# passage is open iff the elements nearest the two arc ends are
# connected inside that subgraph (a detour around the bay cannot
# satisfy this).
import json

import shapely
from scipy.sparse import coo_matrix
from scipy.sparse.csgraph import connected_components
from shapely.geometry import LineString

ed = json.loads(Path(
    "recipes/edits/sample_repro/"
    "edit_001_haneda_d_runway.json").read_text())
arc = np.asarray(ed["arc"], float)
scale = 0.5 * (111e3 * COSW + 111e3)
tube = LineString(arc).buffer(ed["width_m"] / scale)

m = read_fort14(str(OUT / "sample_repro_final.14"))
lo, la = tr.transform(m.nodes[:, 0], m.nodes[:, 1])
T = m.elements
cent = np.column_stack([lo, la])[T].mean(axis=1)
inside = shapely.covers(
    tube, shapely.points(cent[:, 0], cent[:, 1]))
sub = np.where(inside)[0]
if len(sub) == 0:
    raise SystemExit("[verify] FAIL: no elements in corridor")
keys, eids = [], np.tile(np.arange(len(T)), 3)
for a, b in ((0, 1), (1, 2), (2, 0)):
    lonn = np.minimum(T[:, a], T[:, b]).astype(np.int64)
    hinn = np.maximum(T[:, a], T[:, b]).astype(np.int64)
    keys.append(lonn * m.n_nodes + hinn)
keys = np.concatenate(keys)
o = np.argsort(keys, kind="stable")
keys, eids = keys[o], eids[o]
same = keys[1:] == keys[:-1]
pairs = np.column_stack([eids[:-1][same], eids[1:][same]])
insub = inside[pairs[:, 0]] & inside[pairs[:, 1]]
pr = pairs[insub]
g = coo_matrix((np.ones(len(pr)), (pr[:, 0], pr[:, 1])),
               shape=(len(T), len(T)))
_, lab = connected_components(g + g.T, directed=False)
iw = sub[np.argmin(np.hypot((cent[sub, 0] - arc[0, 0]) * COSW,
                            cent[sub, 1] - arc[0, 1]))]
ie = sub[np.argmin(np.hypot((cent[sub, 0] - arc[-1, 0]) * COSW,
                            cent[sub, 1] - arc[-1, 1]))]
ok = lab[iw] == lab[ie]
print(f"[verify] corridor elements={len(sub)}, west elem {iw} "
      f"<-> east elem {ie}: "
      f"{'CONNECTED (passage open)' if ok else 'DISCONNECTED'}",
      flush=True)
raise SystemExit(0 if ok else 1)
