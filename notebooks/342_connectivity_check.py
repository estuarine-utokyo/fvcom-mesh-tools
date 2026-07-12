# Connectivity comparator (owner mandate 2026-07-12): the sample's
# water connectivity is the reference truth for Tokyo Bay.
#   SEVERANCE: sample-meshed water we failed to mesh, clustered;
#     clusters whose removal disconnects the SAMPLE dual graph are
#     severed passages -> CRITICAL.
#   BREACH: our elements sitting on ORIGINAL land (widen artifacts),
#     clustered; clusters whose removal disconnects OUR dual graph
#     are fabricated passages through land barriers -> CRITICAL.
# Exit code 1 on any critical finding.
import os
import sys
import numpy as np
os.environ.setdefault("MPLBACKEND", "Agg")
import matplotlib.pyplot as plt
import geopandas as gpd
import shapely
from pyproj import Transformer
from scipy.sparse import coo_matrix
from scipy.sparse.csgraph import connected_components
from shapely.ops import unary_union
from shapely.strtree import STRtree
from fvcom_mesh_tools.io import read_fort14
from fvcom_mesh_tools.plotting import _add_coast, add_atlas_grid

tr = Transformer.from_crs("EPSG:32654", "EPSG:4326", always_xy=True)
G = os.path.expanduser('~/Github/TB-FVCOM/input/goto2023/grid/')
gd = open(G + 'TokyoBay_grd.dat').read().split('\n')
nn = int(gd[0].split('=')[1]); ne = int(gd[1].split('=')[1])
Ts = np.array([[int(w) for w in gd[2 + i].split()[1:4]]
               for i in range(ne)]) - 1
Ps = np.array([[float(w) for w in gd[2 + ne + i].split()[1:3]]
               for i in range(nn)])
lon_s, lat_s = tr.transform(Ps[:, 0], Ps[:, 1])
Pll_s = np.column_stack([lon_s, lat_s])

m = read_fort14("outputs/sample_repro/sample_repro_final.14")
lon_o, lat_o = tr.transform(m.nodes[:, 0], m.nodes[:, 1])
Pll_o = np.column_stack([lon_o, lat_o])
T_o = m.elements


def dual_pairs(T, n_nodes):
    keys = []
    for a, b in ((0, 1), (1, 2), (2, 0)):
        lo = np.minimum(T[:, a], T[:, b]).astype(np.int64)
        hi = np.maximum(T[:, a], T[:, b]).astype(np.int64)
        keys.append(lo * n_nodes + hi)
    keys = np.concatenate(keys)
    eids = np.tile(np.arange(len(T)), 3)
    o = np.argsort(keys, kind="stable")
    keys, eids = keys[o], eids[o]
    same = keys[1:] == keys[:-1]
    return np.column_stack([eids[:-1][same], eids[1:][same]])


def components_without(pairs, n, drop):
    keep = np.ones(n, bool)
    keep[list(drop)] = False
    pr = pairs[keep[pairs[:, 0]] & keep[pairs[:, 1]]]
    g = coo_matrix((np.ones(len(pr)), (pr[:, 0], pr[:, 1])),
                   shape=(n, n))
    ncomp, lab = connected_components(g + g.T, directed=False)
    return ncomp - int((~keep).sum()), lab


# ---------- SEVERANCE: sample water not covered by our mesh ------
our_polys = [shapely.Polygon(Pll_o[t]) for t in T_o]
tree = STRtree(our_polys)
cent_s = Pll_s[Ts].mean(axis=1)
pts = shapely.points(cent_s[:, 0], cent_s[:, 1])
covered = np.zeros(len(Ts), bool)
for i, pt in enumerate(pts):
    for j in tree.query(pt):
        if our_polys[j].covers(pt):
            covered[i] = True
            break
missing = np.where(~covered)[0]
print(f"[conn] sample elements NOT covered by our mesh: "
      f"{len(missing)} / {len(Ts)}", flush=True)

pairs_s = dual_pairs(Ts, nn)
n0, _ = components_without(pairs_s, len(Ts), [])
severed = []
if len(missing):
    # cluster missing elements (face adjacency)
    sub = set(missing.tolist())
    prs = pairs_s[np.isin(pairs_s, missing).all(axis=1)]
    g = coo_matrix((np.ones(len(prs)),
                    (prs[:, 0], prs[:, 1])),
                   shape=(len(Ts), len(Ts)))
    _, lab = connected_components(g + g.T, directed=False)
    for c in set(int(lab[i]) for i in missing):
        cl = [i for i in missing if lab[i] == c]
        nafter, _ = components_without(pairs_s, len(Ts), cl)
        cc = cent_s[cl].mean(axis=0)
        rec = (len(cl), float(cc[0]), float(cc[1]),
               nafter > n0)
        if nafter > n0:
            severed.append(rec)
        print(f"[conn]   missing x{rec[0]} at ({rec[1]:.4f}, "
              f"{rec[2]:.4f}) severs_sample={rec[3]}", flush=True)

# ---------- BREACH: our elements on ORIGINAL land ----------------
land = unary_union(list(gpd.read_file(
    "outputs/tb_varres_3r/land_osm_wide.shp").geometry))
cent_o = Pll_o[T_o].mean(axis=1)
opts = shapely.points(cent_o[:, 0], cent_o[:, 1])
# tolerance: centroid deeper than ~100 m into original land
land_shrunk = land.buffer(-100.0 / 111e3)
onland = np.array([land_shrunk.covers(p) for p in opts])
breach_idx = np.where(onland)[0]
print(f"[conn] our elements on ORIGINAL land: {len(breach_idx)}",
      flush=True)
pairs_o = dual_pairs(T_o, m.n_nodes)
m0, _ = components_without(pairs_o, len(T_o), [])
breaches = []
if len(breach_idx):
    prs = pairs_o[np.isin(pairs_o, breach_idx).all(axis=1)]
    g = coo_matrix((np.ones(len(prs)),
                    (prs[:, 0], prs[:, 1])),
                   shape=(len(T_o), len(T_o)))
    _, lab = connected_components(g + g.T, directed=False)
    for c in set(int(lab[i]) for i in breach_idx):
        cl = [i for i in breach_idx if lab[i] == c]
        nafter, _ = components_without(pairs_o, len(T_o), cl)
        cc = cent_o[cl].mean(axis=0)
        rec = (len(cl), float(cc[0]), float(cc[1]), nafter > m0)
        if nafter > m0:
            breaches.append(rec)
        print(f"[conn]   on-land x{rec[0]} at ({rec[1]:.4f}, "
              f"{rec[2]:.4f}) bridges={rec[3]}", flush=True)

print(f"[conn] CRITICAL severed passages: {len(severed)}",
      flush=True)
print(f"[conn] CRITICAL land breaches:    {len(breaches)}",
      flush=True)

# map figure
fig, ax = plt.subplots(figsize=(9, 12))
_add_coast(ax, (139.0, 34.5, 141.3, 36.2), "EPSG:4326")
ax.triplot(lon_o, lat_o, T_o, lw=0.25, color="steelblue")
if len(missing):
    ax.scatter(cent_s[missing, 0], cent_s[missing, 1], s=12,
               color="orange", zorder=5,
               label=f"sample water we closed ({len(missing)})")
for nrec, x, y, crit in severed:
    ax.plot([x], [y], marker="x", ms=16, mew=3, color="red",
            zorder=6)
if len(breach_idx):
    ax.scatter(cent_o[breach_idx, 0], cent_o[breach_idx, 1], s=20,
               marker="s", color="purple", zorder=5,
               label=f"our mesh on land ({len(breach_idx)})")
ax.set_xlim(139.57, 140.15); ax.set_ylim(34.93, 35.78)
add_atlas_grid(ax, crs="EPSG:4326")
ax.set_aspect(1 / np.cos(np.deg2rad(35.35)))
ax.legend(loc="lower right")
ax.set_title("connectivity comparison vs sample + original land\n"
             "red X = severed sample passage (critical)")
fig.savefig("outputs/figures/connectivity_check.png", dpi=190,
            bbox_inches="tight")
print("[conn] saved map", flush=True)
sys.exit(1 if (severed or breaches) else 0)
