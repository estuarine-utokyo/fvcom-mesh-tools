"""How well does the extended mesh follow the coastline? (read-only analysis)

    python notebooks/455_coast_fit_analysis.py BUILD_DIR OUT_DIR

BUILD_DIR holds the 445 build (``<case>_grd.dat``, ``report.json``, ``generate/
land_with_base.shp``); OUT_DIR gets ``coast_fit_analysis.json`` and ``worst_nodes.csv``.

Two views of the same question, both against the land polygon the build itself used
(``land_with_base.shp``, UTM 54N):

* node view -- the signed offset of every NEW land-boundary node (+ in water, - on
  land), against its local edge length ``h``;
* coast view -- for points sampled along the coastline inside the new domain, the
  distance to the nearest boundary edge of the mesh: what the mesh leaves out or adds.
"""
from __future__ import annotations

import json
import sys
from pathlib import Path

import geopandas as gpd
import numpy as np
import shapely
from pyproj import Transformer
from scipy.spatial import cKDTree

REPO = Path(__file__).resolve().parents[1]
sys.path.insert(0, str(REPO / "src"))

from fvcom_mesh_tools.coast_fit import boundary_edges, boundary_offsets  # noqa: E402
from fvcom_mesh_tools.io.fvcom_native import read_fvcom_case  # noqa: E402

BUILD, OUT = (Path(a).resolve() for a in sys.argv[1:3])
OUT.mkdir(parents=True, exist_ok=True)
report = json.loads((BUILD / "report.json").read_text())
prefix = sorted(BUILD.glob("*_grd.dat"))[0].name[:-len("_grd.dat")]
mesh = read_fvcom_case(BUILD / f"{prefix}_grd.dat", BUILD / f"{prefix}_dep.dat",
                       BUILD / f"{prefix}_obc.dat")
xy = np.asarray(mesh.nodes)[:, :2]
NB = int(report["merge"]["frozen_base"]["n_base_nodes"])
land = shapely.union_all(list(gpd.read_file(BUILD / "generate" / "land_with_base.shp")
                              .to_crs(32654).geometry))
to_ll = Transformer.from_crs(32654, 4326, always_xy=True)

# ----------------------------------------------------------------- node view
obc = np.concatenate([np.asarray(c, int) for c in mesh.open_boundaries])
ids, off = boundary_offsets(xy, mesh.elements, land, exclude=obc)
be = boundary_edges(mesh.elements)
elen = np.linalg.norm(xy[be[:, 0]] - xy[be[:, 1]], axis=1)
hsum = np.zeros(len(xy))
hcnt = np.zeros(len(xy))
np.add.at(hsum, be[:, 0], elen)
np.add.at(hsum, be[:, 1], elen)
np.add.at(hcnt, be[:, 0], 1)
np.add.at(hcnt, be[:, 1], 1)
h = hsum[ids] / np.maximum(hcnt[ids], 1)
new = ids >= NB


def stats(sel):
    o, hh = off[sel], h[sel]
    a = np.abs(o)
    return {
        "n": int(sel.sum()),
        "in_water": int((o > 0).sum()), "on_land": int((o < 0).sum()),
        "abs_p50_m": float(np.median(a)), "abs_p90_m": float(np.percentile(a, 90)),
        "abs_p99_m": float(np.percentile(a, 99)), "abs_max_m": float(a.max()),
        "on_land_more_than_1m": int((o < -1).sum()),
        "on_land_more_than_10m": int((o < -10).sum()),
        "on_land_more_than_tenth_h": int((o < -0.1 * hh).sum()),
        "in_water_more_than_1m": int((o > 1).sum()),
        "in_water_more_than_tenth_h": int((o > 0.1 * hh).sum()),
        "abs_over_h_p90": float(np.percentile(a / hh, 90)),
        "abs_over_h_max": float((a / hh).max()),
        "h_median_m": float(np.median(hh)),
    }


res = {"case": prefix, "n_nodes": int(len(xy)), "n_base_nodes": NB,
       "new_boundary_nodes": stats(new)}
# by size class of the local edge: where is the misfit large?
bins = [0, 300, 600, 1200, 2400, 1e9]
res["new_by_h_class"] = []
for lo, hi in zip(bins[:-1], bins[1:]):
    s = new & (h >= lo) & (h < hi)
    if s.sum():
        res["new_by_h_class"].append({"h_from_m": lo, "h_to_m": hi, **stats(s)})

# the worst new nodes, with enough to look them up on the atlas
worst = np.argsort(off[new])[:40]
nid = ids[new][worst]
lon, lat = to_ll.transform(xy[nid, 0], xy[nid, 1])
rows = ["node,lon,lat,offset_m,h_m,depth_m"]
for k, n in enumerate(nid):
    j = int(np.flatnonzero(ids == n)[0])
    rows.append(f"{n + 1},{lon[k]:.5f},{lat[k]:.5f},{off[j]:.1f},{h[j]:.0f},"
                f"{float(mesh.depths[n]):.1f}")
(OUT / "worst_nodes.csv").write_text("\n".join(rows) + "\n")

# ----------------------------------------------------------------- coast view
# the new mesh's boundary edges (not the base's, not the open boundary)
nb_edges = be[(be[:, 0] >= NB) | (be[:, 1] >= NB)]
obc_set = set(obc.tolist())
nb_edges = np.array([e for e in nb_edges if not (e[0] in obc_set and e[1] in obc_set)])
tree_pts = []
for a, b in nb_edges:   # sample each boundary edge so the nearest-distance is to the segment
    n = max(2, int(np.ceil(np.linalg.norm(xy[a] - xy[b]) / 25.0)) + 1)
    t = np.linspace(0, 1, n)[:, None]
    tree_pts.append(xy[a] * (1 - t) + xy[b] * t)
mesh_pts = np.vstack(tree_pts)
tree = cKDTree(mesh_pts)
# the coast the NEW mesh is answerable for: samples within NEAR m of a new node (so the
# base's own coast, far coasts and inland water are left out)
NEAR = 2000.0
new_tree = cKDTree(xy[NB:])
coast = land.boundary
parts = [g for g in getattr(coast, "geoms", [coast]) if g.length > 0]
samples = []
for g in parts:
    L = g.length
    if L < 100:
        continue
    samples.append(np.array([g.interpolate(s).coords[0] for s in np.arange(0, L, 100.0)]))
samples = np.vstack(samples)
all_tree = cKDTree(xy)
dn, nearest = all_tree.query(samples)
# the base's own coast belongs to the base: a sample whose nearest node is a base node is out
near = (new_tree.query(samples)[0] <= NEAR) & (nearest >= NB)
samples = samples[near]
d, _ = tree.query(samples)
omitted = d > 600
res["coast_view"] = {
    "near_m": NEAR, "n_samples": int(len(samples)), "spacing_m": 100,
    "coastline_km_near_new_mesh": float(len(samples) * 0.1),
    "dist_p50_m": float(np.median(d)), "dist_p90_m": float(np.percentile(d, 90)),
    "dist_p99_m": float(np.percentile(d, 99)), "dist_max_m": float(d.max()),
    "share_within_50m": float((d <= 50).mean()), "share_within_150m": float((d <= 150).mean()),
    "share_over_300m": float((d > 300).mean()), "share_over_600m": float(omitted.mean()),
}
# where the omitted coast is: coarse 5 km cells, most omitted samples first
if omitted.any():
    ll = np.array(to_ll.transform(samples[omitted, 0], samples[omitted, 1])).T
    cell = np.floor(samples[omitted] / 5000).astype(int)
    keys, inv = np.unique(cell, axis=0, return_inverse=True)
    cnt = np.bincount(inv.ravel())
    rows = ["cell_x5km,cell_y5km,omitted_samples,lon,lat,max_dist_m"]
    for k in np.argsort(-cnt)[:30]:
        sel = inv.ravel() == k
        rows.append(f"{keys[k][0]},{keys[k][1]},{cnt[k]},{ll[sel, 0].mean():.4f},"
                    f"{ll[sel, 1].mean():.4f},{d[omitted][sel].max():.0f}")
    (OUT / "omitted_coast_cells.csv").write_text("\n".join(rows) + "\n")
(OUT / "coast_fit_analysis.json").write_text(json.dumps(res, indent=1))
print(json.dumps(res, indent=1))


# ------------------------------------------------------------------- figure
import matplotlib  # noqa: E402

matplotlib.use("Agg")
import matplotlib.pyplot as plt  # noqa: E402

from fvcom_mesh_tools.plotting import add_atlas_grid, use_readable_style  # noqa: E402

use_readable_style()
fig, ax = plt.subplots(figsize=(11, 10))
segs = np.array([[to_ll.transform(*xy[a]), to_ll.transform(*xy[b])] for a, b in nb_edges])
ax.add_collection(matplotlib.collections.LineCollection(segs, colors="0.55", linewidths=0.5,
                                                        label="new mesh boundary"))
if omitted.any():
    o = np.array(to_ll.transform(samples[omitted, 0], samples[omitted, 1])).T
    ax.plot(o[:, 0], o[:, 1], ".", ms=2.5, color="tab:red",
            label=f"coast >600 m from the mesh ({int(omitted.sum())} samples of 100 m)")
bad = ids[new][off[new] < -10]
bl = np.array(to_ll.transform(xy[bad, 0], xy[bad, 1])).T
ax.plot(bl[:, 0], bl[:, 1], "o", ms=4, mfc="none", color="tab:blue",
        label=f"new boundary node >10 m on land ({len(bad)})")
ax.set_xlim(137.8, 141.3)
ax.set_ylim(33.2, 35.8)
ax.set_aspect(1 / np.cos(np.deg2rad(34.8)))
add_atlas_grid(ax)
ax.legend(loc="lower right")
ax.set_xlabel("lon (deg E)")
ax.set_ylabel("lat (deg N)")
ax.set_title("Where the extended mesh departs from the coastline")
fig.savefig(OUT / "coast_fit_map.png", dpi=140, bbox_inches="tight")


# ------------------------------------------------------------------- zooms
ZOOMS = {"suruga_head": (138.40, 34.88, 138.62, 35.10), "oshima": (139.30, 34.65, 139.50, 34.85)}
els = np.asarray(mesh.elements)
for name, (x0, y0, x1, y1) in ZOOMS.items():
    px0, py0 = Transformer.from_crs(4326, 32654, always_xy=True).transform(x0, y0)
    px1, py1 = Transformer.from_crs(4326, 32654, always_xy=True).transform(x1, y1)
    fig, ax = plt.subplots(figsize=(10, 9))
    sel = els[(xy[els, 0].max(axis=1) > px0) & (xy[els, 0].min(axis=1) < px1)
              & (xy[els, 1].max(axis=1) > py0) & (xy[els, 1].min(axis=1) < py1)]
    ax.triplot(xy[:, 0] / 1e3, xy[:, 1] / 1e3, sel, color="0.6", lw=0.35)
    for g in getattr(land, "geoms", [land]):
        if g.intersects(shapely.box(px0, py0, px1, py1)):
            ax.fill(*(np.asarray(g.exterior.coords).T / 1e3), color="tan", alpha=0.45, zorder=0)
            for hole in g.interiors:
                ax.fill(*(np.asarray(hole.coords).T / 1e3), color="white", zorder=1)
            ax.plot(*(np.asarray(g.exterior.coords).T / 1e3), color="saddlebrown", lw=1.0)
    badn = ids[new][off[new] < -10]
    ax.plot(xy[badn, 0] / 1e3, xy[badn, 1] / 1e3, "o", ms=6, mfc="none", color="tab:blue",
            label="boundary node >10 m on land")
    ax.set_xlim(px0 / 1e3, px1 / 1e3)
    ax.set_ylim(py0 / 1e3, py1 / 1e3)
    ax.set_aspect(1)
    ax.set_title(f"{name}: mesh (grey), OSM land (tan), nodes on land (blue)")
    ax.set_xlabel("x UTM54N (km)")
    ax.set_ylabel("y UTM54N (km)")
    ax.legend(loc="lower right")
    fig.savefig(OUT / f"zoom_{name}.png", dpi=130, bbox_inches="tight")
