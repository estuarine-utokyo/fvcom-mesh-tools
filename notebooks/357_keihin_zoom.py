# H8 lower-left (Keihin canal, Kawasaki): sample vs ours zoom +
# measured through-connectivity along the canal chain.
import os
import numpy as np
os.environ.setdefault("MPLBACKEND", "Agg")
import geopandas as gpd
import matplotlib.pyplot as plt
import shapely
from pyproj import Transformer
from shapely.ops import unary_union
from fvcom_mesh_tools.io import read_fort14
from fvcom_mesh_tools.plotting import add_atlas_grid, use_readable_style

use_readable_style()
COSW = float(np.cos(np.deg2rad(35.35)))
CX, CY, HALF = 139.757, 35.505, 0.024

tr = Transformer.from_crs("EPSG:32654", "EPSG:4326", always_xy=True)
G = os.path.expanduser('~/Github/TB-FVCOM/input/goto2023/grid/')
gd = open(G + 'TokyoBay_grd.dat').read().split('\n')
nn = int(gd[0].split('=')[1]); ne = int(gd[1].split('=')[1])
Ts = np.array([[int(w) for w in gd[2 + i].split()[1:4]]
               for i in range(ne)]) - 1
Ps = np.array([[float(w) for w in gd[2 + ne + i].split()[1:3]]
               for i in range(nn)])
lon_s, lat_s = tr.transform(Ps[:, 0], Ps[:, 1])

m = read_fort14("outputs/sample_repro/sample_repro_final.14")
lon_o, lat_o = tr.transform(m.nodes[:, 0], m.nodes[:, 1])
T = m.elements
land = unary_union(list(gpd.read_file(
    "outputs/tb_varres_3r/land_osm_wide.shp").geometry))
gser = gpd.GeoSeries([land], crs="EPSG:4326")

fig, axes = plt.subplots(1, 2, figsize=(15, 8))
for ax, (ttl, lo, la, tt, col) in zip(axes, [
        ("goto2023 sample", lon_s, lat_s, Ts, "0.35"),
        ("ours (anchor-rule fix)", lon_o, lat_o, T,
         "steelblue")]):
    gser.plot(ax=ax, color="0.9", edgecolor="0.6", lw=0.5,
              zorder=1)
    ax.triplot(lo, la, tt, lw=0.6, color=col, zorder=3)
    ax.set_xlim(CX - HALF / COSW, CX + HALF / COSW)
    ax.set_ylim(CY - HALF, CY + HALF)
    ax.set_aspect(1 / COSW)
    add_atlas_grid(ax, crs="EPSG:4326")
    ax.set_title(ttl)
fig.suptitle("Keihin canal, Kawasaki section (H8 lower-left)")
fig.tight_layout(rect=(0, 0, 1, 0.94))
fig.savefig("outputs/figures/keihin_h8_zoom.png", dpi=180,
            bbox_inches="tight")
print("[keihin] saved", flush=True)

# measured through-connectivity: dual-graph path constrained to
# the canal corridor between west (G9) and east (Haneda) ends
P = np.column_stack([lon_o, lat_o])
cent = P[T].mean(axis=1)
# waypoints FOLLOW the real canal (the first probe cut across
# the island at 35.51-35.52 and split its own tube -> false
# DISCONNECTED): Keihin canal east along ~35.505, then Daishi
# canal NE to the Tama mouth
band = shapely.LineString([(139.712, 35.489), (139.735, 35.502),
                           (139.757, 35.5055), (139.775, 35.504),
                           (139.7815, 35.508), (139.786, 35.514),
                           (139.792, 35.521),
                           (139.796, 35.526)]).buffer(0.006)
inside = shapely.covers(band, shapely.points(cent[:, 0],
                                             cent[:, 1]))
sub = np.where(inside)[0]
keys = []
for a, b in ((0, 1), (1, 2), (2, 0)):
    lo2 = np.minimum(T[:, a], T[:, b]).astype(np.int64)
    hi2 = np.maximum(T[:, a], T[:, b]).astype(np.int64)
    keys.append(lo2 * m.n_nodes + hi2)
keys = np.concatenate(keys)
eids = np.tile(np.arange(len(T)), 3)
o = np.argsort(keys, kind="stable")
ks, es = keys[o], eids[o]
same = ks[1:] == ks[:-1]
pairs = np.column_stack([es[:-1][same], es[1:][same]])
pr = pairs[inside[pairs[:, 0]] & inside[pairs[:, 1]]]
from scipy.sparse import coo_matrix
from scipy.sparse.csgraph import connected_components
g = coo_matrix((np.ones(len(pr)), (pr[:, 0], pr[:, 1])),
               shape=(len(T), len(T)))
_, lab = connected_components(g + g.T, directed=False)
iw = sub[np.argmin(np.hypot((cent[sub, 0] - 139.712) * COSW,
                            cent[sub, 1] - 35.489))]
ie = sub[np.argmin(np.hypot((cent[sub, 0] - 139.796) * COSW,
                            cent[sub, 1] - 35.526))]
ok = lab[iw] == lab[ie]
print(f"[keihin] corridor cells={len(sub)}, west {iw} <-> east "
      f"{ie}: {'CONNECTED' if ok else 'DISCONNECTED'}",
      flush=True)
# where is the break? report each component's extent along the
# canal and the gap between consecutive components
import collections
comp = collections.defaultdict(list)
for e in sub:
    comp[int(lab[e])].append(e)
line = shapely.LineString(np.asarray(band.exterior.coords)[:0])
axis = shapely.LineString([(139.712, 35.489), (139.735, 35.502),
                           (139.757, 35.5055), (139.775, 35.504),
                           (139.7815, 35.508), (139.786, 35.514),
                           (139.792, 35.521), (139.796, 35.526)])
for cid, els in sorted(comp.items(),
                       key=lambda kv: -len(kv[1])):
    ss = [float(axis.project(shapely.Point(cent[e])))
          for e in els]
    print(f"[keihin]   comp {cid}: n={len(els)} span "
          f"s=[{min(ss)*111e3:.0f}, {max(ss)*111e3:.0f}] m "
          f"center=({np.mean(cent[els,0]):.4f},"
          f"{np.mean(cent[els,1]):.4f})", flush=True)
raise SystemExit(0 if ok else 1)
