# Before and after: two refinements of the same port drawn over the raw OSM.
#
# Each column is one refinement output directory.  Grey is the raw OSM land,
# blue the mesh coast, red the mesh walls, thin grey the mesh.  Where a run
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

MESH_EPSG = 32654
VIEWS = [(393010.0, 3909480.0, 1100.0), (392500.0, 3909000.0, 350.0)]


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


dirs = [Path(a).resolve() for a in sys.argv[1:3]]
cx, cy, half = VIEWS[0]
box = shapely.box(cx - half, cy - half, cx + half, cy + half)
land = gpd.read_file(os.environ["FMESH_LAND"]).to_crs(MESH_EPSG)
land = land[land.intersects(box)]

fig, axes = plt.subplots(2, 2, figsize=(17, 17), constrained_layout=True)
for col, d in enumerate(dirs):
    m = read_fort14(next(d.glob("*.14")))
    seg, wall = boundary(m)
    stages = np.load(d / "walls_stages.npz") if (d / "walls_stages.npz").exists() else None
    for row, (x0, y0, h) in enumerate(VIEWS):
        ax = axes[row, col]
        land.plot(ax=ax, color="0.8", edgecolor="0.5", lw=0.5, zorder=0)
        ax.triplot(m.nodes[:, 0], m.nodes[:, 1], m.elements, lw=0.2, color="0.5", zorder=1)
        ax.add_collection(LineCollection(seg[~wall], colors="tab:blue", lw=1.2, zorder=3))
        ax.add_collection(LineCollection(seg[wall], colors="tab:red", lw=2.2, zorder=4))
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
ax.plot([], [], color="tab:blue", label="mesh coast")
ax.plot([], [], color="tab:red", lw=2, label="mesh wall")
ax.plot([], [], color="tab:orange", ls="--", label="wall centreline from OSM")
ax.plot([], [], color="tab:purple", ls=":", label="wall piece kept in the hole")
ax.legend(loc="upper right", fontsize=9)
png = dirs[1] / "before_after.png"
fig.savefig(png, dpi=130)
print(f"wrote {png}")
