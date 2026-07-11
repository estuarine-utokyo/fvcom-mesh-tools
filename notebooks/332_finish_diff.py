# BEFORE (oceanmesh output: sample_repro_utm.14) vs AFTER (FVCOM
# finishing: sample_repro_final.14): what did the finishing chain
# change, and what issues remain.
import os
import numpy as np
os.environ.setdefault("MPLBACKEND", "Agg")
import matplotlib.pyplot as plt
from pyproj import Transformer
from scipy.spatial import cKDTree
from fvcom_mesh_tools.io import read_fort14
from fvcom_mesh_tools.plotting import _add_coast, add_atlas_grid

OUT = "outputs/figures"
mb = read_fort14("outputs/sample_repro/sample_repro_utm.14")
ma = read_fort14("outputs/sample_repro/sample_repro_final.14")
tr = Transformer.from_crs("EPSG:32654", "EPSG:4326", always_xy=True)
lon_b, lat_b = tr.transform(mb.nodes[:, 0], mb.nodes[:, 1])
lon_a, lat_a = tr.transform(ma.nodes[:, 0], ma.nodes[:, 1])
Pb_ll = np.column_stack([lon_b, lat_b])
Pa_ll = np.column_stack([lon_a, lat_a])

# node correspondence: after -> nearest before
d_ab, idx_ab = cKDTree(mb.nodes).query(ma.nodes)
moved = d_ab > 0.01           # nodes displaced by the finishing (m)
# removed BEFORE-nodes: those not matched by any after-node
matched_b = np.zeros(mb.n_nodes, bool)
matched_b[idx_ab[d_ab < 0.01]] = True
removed_b = np.where(~matched_b)[0]

# element connectivity change: canonical triples in BEFORE index space
map_a2b = idx_ab.copy()        # nearest-before id for every after node
eb = {tuple(sorted(t)) for t in mb.elements.tolist()}
ea = {tuple(sorted(map_a2b[list(t)])) for t in ma.elements.tolist()}
gone = eb - ea                 # elements destroyed by finishing
new = ea - eb                  # elements created by finishing
print(f"[diff] nodes: before={mb.n_nodes:,} after={ma.n_nodes:,} "
      f"removed={len(removed_b)} moved(>1cm)={int(moved.sum())} "
      f"disp p50/max = "
      f"{np.percentile(d_ab[moved], 50) if moved.any() else 0:.1f}/"
      f"{d_ab.max():.1f} m", flush=True)
print(f"[diff] elements: before={len(mb.elements):,} "
      f"after={len(ma.elements):,} rewired: -{len(gone)} +{len(new)}",
      flush=True)

# change sites for zooms: cluster moved nodes + rewired elements
site_pts = [Pa_ll[i] for i in np.where(moved)[0]]
cent_b = {t: Pb_ll[list(t)].mean(axis=0) for t in gone}
site_pts += list(cent_b.values())
site_pts = np.asarray(site_pts)
# greedy clustering at 2.5 km
sites = []
used = np.zeros(len(site_pts), bool)
order = np.argsort(-np.r_[d_ab[moved],
                          np.full(len(cent_b), 1e9)])  # flips first
for k in order:
    if used[k]:
        continue
    c = site_pts[k]
    dd = np.hypot((site_pts[:, 0] - c[0]) * 91e3,
                  (site_pts[:, 1] - c[1]) * 111e3)
    used |= dd < 2500
    sites.append(c)
print(f"[diff] change sites (clusters): {len(sites)}", flush=True)

COAST = (139.0, 34.5, 141.3, 36.2)
NZ = min(6, len(sites))
fig, axes = plt.subplots(2, 3, figsize=(19, 12))
for ax, c in zip(axes.ravel(), sites[:NZ]):
    R = 0.022
    x0, x1, y0, y1 = c[0] - R, c[0] + R, c[1] - R * 0.82, c[1] + R * 0.82
    _add_coast(ax, COAST, "EPSG:4326")
    ax.triplot(lon_b, lat_b, mb.elements, lw=1.6, color="0.72",
               zorder=2)
    ax.triplot(lon_a, lat_a, ma.elements, lw=0.9, color="crimson",
               zorder=3)
    mv = moved & (np.abs(Pa_ll[:, 0] - c[0]) < R) \
        & (np.abs(Pa_ll[:, 1] - c[1]) < R)
    for i in np.where(mv)[0]:
        j = idx_ab[i]
        ax.annotate("", xy=Pa_ll[i], xytext=Pb_ll[j],
                    arrowprops=dict(arrowstyle="->", color="blue",
                                    lw=1.8), zorder=5)
    ob = np.asarray(ma.open_boundaries[0], int)
    ax.plot(Pa_ll[ob, 0], Pa_ll[ob, 1], color="red", lw=2.0,
            marker="o", ms=3, zorder=4)
    ax.set_xlim(x0, x1); ax.set_ylim(y0, y1)
    add_atlas_grid(ax, crs="EPSG:4326")
    ax.set_aspect(1 / np.cos(np.deg2rad(c[1])))
    ax.set_title(f"({c[0]:.3f}, {c[1]:.3f})", fontsize=10)
fig.suptitle("finishing diff: BEFORE oceanmesh (gray) -> AFTER "
             "FVCOM finishing (red); arrows = node moves", y=0.995)
fig.savefig(f"{OUT}/finish_diff_zooms.png", dpi=190,
            bbox_inches="tight")
print("saved zooms", flush=True)

# overview map of all change locations
fig2, ax = plt.subplots(figsize=(9, 12))
_add_coast(ax, COAST, "EPSG:4326")
ax.triplot(lon_a, lat_a, ma.elements, lw=0.25, color="steelblue",
           zorder=2)
mvi = np.where(moved)[0]
ax.scatter(Pa_ll[mvi, 0], Pa_ll[mvi, 1], s=60, facecolors="none",
           edgecolors="blue", lw=1.6, zorder=5,
           label=f"moved nodes ({len(mvi)})")
if len(cent_b):
    cb = np.asarray(list(cent_b.values()))
    ax.scatter(cb[:, 0], cb[:, 1], s=90, marker="s",
               facecolors="none", edgecolors="crimson", lw=1.6,
               zorder=5, label=f"rewired elements ({len(gone)})")
if len(removed_b):
    ax.scatter(Pb_ll[removed_b, 0], Pb_ll[removed_b, 1], s=90,
               marker="x", color="black", zorder=5,
               label=f"removed nodes ({len(removed_b)})")
ob = np.asarray(ma.open_boundaries[0], int)
ax.plot(Pa_ll[ob, 0], Pa_ll[ob, 1], color="red", lw=2, zorder=4)
ax.set_xlim(139.57, 140.15); ax.set_ylim(34.93, 35.78)
add_atlas_grid(ax, crs="EPSG:4326")
ax.set_aspect(1 / np.cos(np.deg2rad(35.35)))
ax.legend(loc="lower right")
ax.set_title("FVCOM finishing: all change sites")
fig2.savefig(f"{OUT}/finish_diff_map.png", dpi=190,
             bbox_inches="tight")
print("saved map", flush=True)
