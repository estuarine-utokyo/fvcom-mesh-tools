# Issue REGISTRY builder (owner-approved manual-editing loop,
# 2026-07-12): auto-extract remaining issues vs the goto2023
# sample, assign stable IDs, and for each missing-waterway cluster
# propose a channel-ARC edit (arc polyline + width) ready to be
# reviewed and copied into outputs/sample_repro/edits/.
#
# Outputs:
#   outputs/sample_repro/issue_registry.json
#   outputs/sample_repro/edits_proposed/W##_*.json  (arc proposals)
#   outputs/figures/issue_registry_atlas.png        (ID'd zooms)
import json
import os
from pathlib import Path

import numpy as np
os.environ.setdefault("MPLBACKEND", "Agg")
import geopandas as gpd
import matplotlib.pyplot as plt
import shapely
from pyproj import Transformer
from scipy.sparse import coo_matrix
from scipy.sparse.csgraph import connected_components
from shapely.ops import unary_union
from shapely.strtree import STRtree

from fvcom_mesh_tools.channel_arcs import arc_from_points
from fvcom_mesh_tools.io import read_fort14
from fvcom_mesh_tools.plotting import use_readable_style

use_readable_style()

OUT = Path("outputs/sample_repro")
PROP = OUT / "edits_proposed"
PROP.mkdir(parents=True, exist_ok=True)
COSW = float(np.cos(np.deg2rad(35.35)))
H0 = 290.0                       # field hmin (m); mesh ~1.2x
WIDTH_2ROW = round(1.75 * 1.2 * H0)   # 2 standard rows ~ 610 m

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

m = read_fort14(str(OUT / "sample_repro_final.14"))
lon_o, lat_o = tr.transform(m.nodes[:, 0], m.nodes[:, 1])
T_o = m.elements

land = unary_union(list(gpd.read_file(
    "outputs/tb_varres_3r/land_osm_wide.shp").geometry))


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


# ---- coverage (same criterion as the 342 comparator) ------------
our_polys = [shapely.Polygon(np.column_stack([lon_o, lat_o])[t])
             for t in T_o]
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
print(f"[reg] uncovered sample elements: {len(missing)}", flush=True)

pairs_s = dual_pairs(Ts, nn)
g_all = coo_matrix((np.ones(len(pairs_s)),
                    (pairs_s[:, 0], pairs_s[:, 1])),
                   shape=(len(Ts), len(Ts)))
n0, _ = connected_components(g_all + g_all.T, directed=False)

# adjacency restricted to missing elements -> clusters
prs = pairs_s[np.isin(pairs_s, missing).all(axis=1)]
g = coo_matrix((np.ones(len(prs)), (prs[:, 0], prs[:, 1])),
               shape=(len(Ts), len(Ts)))
_, lab = connected_components(g + g.T, directed=False)

# neighbour lists for ring expansion
nbr = [[] for _ in range(len(Ts))]
for a, b in pairs_s:
    nbr[a].append(b)
    nbr[b].append(a)

# sample local edge length per element (metres, rough)
P_m = np.column_stack([Pll_s[:, 0] * COSW * 111e3,
                       Pll_s[:, 1] * 111e3])
el_edge = np.zeros(len(Ts))
for k, t in enumerate(Ts):
    e = P_m[t]
    el_edge[k] = np.mean([np.linalg.norm(e[0] - e[1]),
                          np.linalg.norm(e[1] - e[2]),
                          np.linalg.norm(e[2] - e[0])])

clusters = []
for c in sorted(set(int(lab[i]) for i in missing)):
    cl = np.array([i for i in missing if lab[i] == c])
    keep = np.ones(len(Ts), bool)
    keep[cl] = False
    pr = pairs_s[keep[pairs_s[:, 0]] & keep[pairs_s[:, 1]]]
    g2 = coo_matrix((np.ones(len(pr)), (pr[:, 0], pr[:, 1])),
                    shape=(len(Ts), len(Ts)))
    ncomp, _ = connected_components(g2 + g2.T, directed=False)
    severs = bool(ncomp - len(cl) > n0)
    clusters.append((cl, severs))
clusters.sort(key=lambda t: -len(t[0]))

# STABLE IDs: match clusters to the previous registry by center
# proximity (< ~600 m) so an ID keeps meaning the same waterway
# across runs (the owner refers to issues by ID). Unmatched
# clusters get fresh sequential numbers.
prev_path = OUT / "issue_registry.json"
prev = (json.loads(prev_path.read_text())
        if prev_path.exists() else [])
_used: set[str] = set()
_next = 1 + max(
    [int(p["id"][1:]) for p in prev if p["id"][1:].isdigit()],
    default=0)


def stable_id(cc):
    global _next
    best, bestd = None, 0.0055
    for p in prev:
        if p["id"] in _used:
            continue
        d = float(np.hypot((p["center_lonlat"][0] - cc[0]) * COSW,
                           p["center_lonlat"][1] - cc[1]))
        if d < bestd:
            bestd, best = d, p["id"]
    if best is not None:
        _used.add(best)
        return best
    sid = f"W{_next:02d}"
    _next += 1
    return sid


registry = []
for rank, (cl, severs) in enumerate(clusters, start=1):
    cc = cent_s[cl].mean(axis=0)
    # 2-ring expansion along the sample dual graph: reaches the
    # open water the channel connects to, so the proposed arc's
    # ENDS land in meshable water on both sides
    ring = set(cl.tolist())
    for _ in range(2):
        ring |= {j for i in ring for j in nbr[i]}
    ring = np.array(sorted(ring))
    pm = np.column_stack([cent_s[ring, 0] * COSW, cent_s[ring, 1]])
    arc_sc = arc_from_points(pm, smooth_passes=2)
    arc_ll = np.column_stack([arc_sc[:, 0] / COSW, arc_sc[:, 1]])
    ext_m = np.ptp(P_m[Ts[cl]].reshape(-1, 2), axis=0)
    rec = {
        "id": stable_id(cc),
        "kind": "missing-waterway",
        "n_sample_elements": int(len(cl)),
        "center_lonlat": [round(float(cc[0]), 4),
                          round(float(cc[1]), 4)],
        "extent_m": [round(float(ext_m[0])), round(float(ext_m[1]))],
        "sample_local_edge_m": round(float(el_edge[cl].mean())),
        "severs_sample": severs,
        "proposed_action": "arc-band widen"
        if len(cl) >= 3 or severs else "review (tiny)",
    }
    registry.append(rec)
    prop = {
        "id": rec["id"],
        "note": f"auto-proposed arc for missing waterway at "
                f"({cc[0]:.4f}, {cc[1]:.4f}), "
                f"n={len(cl)}, severs_sample={severs}",
        "arc": [[round(float(x), 5), round(float(y), 5)]
                for x, y in arc_ll],
        "width_m": WIDTH_2ROW,
        "min_gap_m": 150.0,
    }
    (PROP / f"{rec['id']}.json").write_text(
        json.dumps(prop, indent=1))
    print(f"[reg] {rec['id']}: n={len(cl):3d} at ({cc[0]:.4f}, "
          f"{cc[1]:.4f}) severs={severs} "
          f"edge~{rec['sample_local_edge_m']} m "
          f"arc_pts={len(arc_ll)}", flush=True)

(OUT / "issue_registry.json").write_text(
    json.dumps(registry, indent=1))
print(f"[reg] registry: {len(registry)} issues -> "
      f"{OUT / 'issue_registry.json'}", flush=True)

# ---- ID'd zoom atlas -------------------------------------------
show = registry[:12]
ncol = 3
nrow = int(np.ceil(len(show) / ncol))
fig, axes = plt.subplots(nrow, ncol,
                         figsize=(5.6 * ncol, 5.6 * nrow))
axes = np.atleast_1d(axes).ravel()
gxy = gpd.GeoSeries([land], crs="EPSG:4326")
for ax, rec in zip(axes, show):
    cx, cy = rec["center_lonlat"]
    half = max(rec["extent_m"][0], rec["extent_m"][1], 1200) \
        * 0.75 / 111e3
    gxy.plot(ax=ax, color="0.88", edgecolor="0.55", lw=0.4)
    ax.triplot(Pll_s[:, 0], Pll_s[:, 1], Ts, lw=0.35,
               color="0.45")
    ax.triplot(lon_o, lat_o, T_o, lw=0.55, color="steelblue")
    cl_pts = [i for i in missing
              if abs(cent_s[i, 0] - cx) < 2 * half
              and abs(cent_s[i, 1] - cy) < 2 * half]
    ax.scatter(cent_s[cl_pts, 0], cent_s[cl_pts, 1], s=14,
               color="orange", zorder=6)
    prop = json.loads((PROP / f"{rec['id']}.json").read_text())
    a = np.array(prop["arc"])
    ax.plot(a[:, 0], a[:, 1], "-", color="crimson", lw=1.6,
            zorder=7)
    ax.set_xlim(cx - half / COSW, cx + half / COSW)
    ax.set_ylim(cy - half, cy + half)
    ax.set_aspect(1 / COSW)
    ax.set_xticks([]); ax.set_yticks([])
    ax.set_title(f"{rec['id']}  n={rec['n_sample_elements']}  "
                 f"({cx:.3f}, {cy:.3f})\n"
                 f"severs={rec['severs_sample']}  "
                 f"edge~{rec['sample_local_edge_m']} m")
for ax in axes[len(show):]:
    ax.axis("off")
fig.suptitle("issue registry: sample waterways our mesh misses\n"
             "(gray=sample, blue=ours, orange=uncovered, "
             "red=proposed arc)")
fig.tight_layout(rect=(0, 0, 1, 0.97))
Path("outputs/figures").mkdir(parents=True, exist_ok=True)
fig.savefig("outputs/figures/issue_registry_atlas.png", dpi=170,
            bbox_inches="tight")
print("[reg] saved atlas", flush=True)
