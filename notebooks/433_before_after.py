# Before and after: two refinements of the same port drawn over the raw OSM.
#
# Each column is one refinement output directory.  Grey is the raw OSM land;
# the mesh is drawn by plotting.draw_mesh: thin grey edges, every solid
# boundary (coast, quay, wall) black, the open boundary red.  Where a run
# saved its wall stages (walls_stages.npz), the centrelines extracted from
# OSM are dashed orange and the pieces kept inside the hole dotted purple,
# so a wall that is missing can be traced to the stage that lost it.
#
#   FMESH_LAND=<land.shp> python notebooks/433_before_after.py <before dir> <after dir>
import os
import sys
from pathlib import Path

import numpy as np

os.environ.setdefault("MPLBACKEND", "Agg")
import geopandas as gpd  # noqa: E402
import matplotlib.pyplot as plt  # noqa: E402
import shapely  # noqa: E402
from matplotlib.collections import LineCollection  # noqa: E402

ROOT = Path(__file__).resolve().parents[1]
sys.path.insert(0, str(ROOT / "src"))
from fvcom_mesh_tools.io.fort14 import read_fort14  # noqa: E402
from fvcom_mesh_tools.plotting import (  # noqa: E402
    OPEN_BOUNDARY_COLOR,
    SOLID_BOUNDARY_COLOR,
    draw_mesh,
)

MESH_EPSG = 32654
VIEWS = [(393010.0, 3909480.0, 1100.0), (392500.0, 3909000.0, 350.0),
         (392460.0, 3908790.0, 150.0)]


def boundary(m):
    xy, tri = m.nodes[:, :2], m.elements
    e = np.sort(np.vstack([tri[:, [0, 1]], tri[:, [1, 2]], tri[:, [2, 0]]]), axis=1)
    u, c = np.unique(e, axis=0, return_counts=True)
    b = u[c == 1]
    seg = np.stack([xy[b[:, 0]], xy[b[:, 1]]], axis=1)
    key = np.round(np.sort(seg.view([("x", float), ("y", float)]), axis=1)
                   .view(float).reshape(len(b), 4), 3)
    _, inv, cnt = np.unique(key, axis=0, return_inverse=True, return_counts=True)
    return seg, cnt[inv.ravel()] == 2


def tilt_table(d, land_raw, seg, wall):
    """Each structure the filter took from the land, and how the wall follows it.

    The footprint is the raw OSM land the width filter removed; its long axis
    is the minimum rotated rectangle's.  The wall's direction is the principal
    axis of the wall edges inside the footprint.
    """
    filt = shapely.union_all(gpd.read_file(d / "shoreline_filtered.shp").geometry.values)
    resid = shapely.difference(land_raw, filt)
    wl = shapely.MultiLineString(list(seg[wall])) if wall.any() else None
    rows = []
    for g in getattr(resid, "geoms", [resid]):
        if g.area < 50 or not g.intersects(box):
            continue
        cc = np.asarray(g.minimum_rotated_rectangle.exterior.coords)[:4]
        s1, s2 = cc[1] - cc[0], cc[2] - cc[1]
        length = max(np.linalg.norm(s1), np.linalg.norm(s2))
        # a structure, not a sliver the filter shaved off a quay
        if g.area / length < 2.0 or length < 3.0 * g.area / length:
            continue
        ax_ = s1 if np.linalg.norm(s1) >= np.linalg.norm(s2) else s2
        axis = np.degrees(np.arctan2(ax_[1], ax_[0])) % 180
        tilt = np.nan
        if wl is not None:
            w_in = shapely.intersection(wl, g.buffer(3.0))
            if not w_in.is_empty and w_in.length > 5:
                pts = np.vstack([np.asarray(q.coords) for q in getattr(w_in, "geoms", [w_in])])
                v = np.linalg.svd(pts - pts.mean(0))[2][0]
                wd = np.degrees(np.arctan2(v[1], v[0])) % 180
                tilt = min(abs(wd - axis), 180 - abs(wd - axis))
        p = g.representative_point()
        rows.append((p.x, p.y, length, g.area / length, tilt))
    print(f"[433] {d.name}: structures the filter took from the land")
    print(f"      {'x':>9} {'y':>10} {'length':>7} {'width':>6} {'wall tilt':>9}")
    for x, y, L, w, t in rows:
        print(f"      {x:9.0f} {y:10.0f} {L:7.0f} {w:6.1f} "
              + ("   (none)" if np.isnan(t) else f"{t:9.1f}"))


dirs = [Path(a).resolve() for a in sys.argv[1:3]]
cx, cy, half = VIEWS[0]
box = shapely.box(cx - half, cy - half, cx + half, cy + half)
land = gpd.read_file(os.environ["FMESH_LAND"]).to_crs(MESH_EPSG)
land = land[land.intersects(box)]
land_raw = shapely.union_all(land.geometry.values).intersection(box)

fig, axes = plt.subplots(3, 2, figsize=(17, 25.5), constrained_layout=True)
for col, d in enumerate(dirs):
    m = read_fort14(next(d.glob("*.14")))
    seg, wall = boundary(m)
    tilt_table(d, land_raw, seg, wall)
    stages = np.load(d / "walls_stages.npz") if (d / "walls_stages.npz").exists() else None
    for row, (x0, y0, h) in enumerate(VIEWS):
        ax = axes[row, col]
        land.plot(ax=ax, color="0.8", edgecolor="0.5", lw=0.5, zorder=0)
        if (d / "shoreline_filtered.shp").exists():
            gpd.read_file(d / "shoreline_filtered.shp").boundary.plot(
                ax=ax, color="tab:green", lw=0.8, ls="--", zorder=2)
        draw_mesh(ax, m.nodes, m.elements, m.open_boundaries, zorder=3)
        if stages is not None:
            src = [stages[k] for k in stages.files if k.startswith("src")]
            pcs = [stages[k] for k in stages.files if k.startswith("piece")]
            ax.add_collection(LineCollection(src, colors="tab:orange", lw=1.0,
                                             linestyles="--", zorder=5))
            ax.add_collection(LineCollection(pcs, colors="tab:purple", lw=1.4,
                                             linestyles=":", zorder=6))
        ax.set_xlim(x0 - h, x0 + h)
        ax.set_ylim(y0 - h, y0 + h)
        ax.set_aspect("equal")
        ax.set_title(f"{'before' if col == 0 else 'after'}: {d.name}  "
                     f"(NP={m.n_nodes:,}, NE={m.n_elements:,})", fontsize=10)
ax = axes[0, 0]
ax.plot([], [], color="0.8", lw=6, label="OSM land (raw)")
ax.plot([], [], color="tab:green", ls="--", label="land after the width filter")
ax.plot([], [], color=SOLID_BOUNDARY_COLOR, lw=1.4,
        label="solid boundary (coast, quay, wall)")
ax.plot([], [], color=OPEN_BOUNDARY_COLOR, lw=1.8, label="open boundary")
ax.plot([], [], color="tab:orange", ls="--", label="wall centreline from OSM")
ax.plot([], [], color="tab:purple", ls=":", label="wall piece kept in the hole")
ax.legend(loc="upper right", fontsize=9)
png = dirs[1] / "before_after.png"
fig.savefig(png, dpi=130)
print(f"wrote {png}")
