# Are the gaps in the walls in OSM, or did the pipeline make them?
#
# The raw OSM land (FMESH_LAND) is drawn under the delivered mesh, whose
# solid boundaries -- coast, quay, wall -- are black.  Where a black line
# stops but the grey OSM land runs on, the gap is the pipeline's; where the
# grey stops too, it is OSM's.
#
#   FMESH_LAND=<land.shp> python notebooks/432_walls_vs_osm.py <refinement output dir>
import os
import sys
from pathlib import Path

os.environ.setdefault("MPLBACKEND", "Agg")
import geopandas as gpd  # noqa: E402
import matplotlib.pyplot as plt  # noqa: E402
import shapely  # noqa: E402

ROOT = Path(__file__).resolve().parents[1]
sys.path.insert(0, str(ROOT / "src"))
from fvcom_mesh_tools.io.fort14 import read_fort14  # noqa: E402
from fvcom_mesh_tools.plotting import SOLID_BOUNDARY_COLOR, draw_mesh  # noqa: E402

MESH_EPSG = 32654
VIEWS = [(393010.0, 3909480.0, 1100.0), (392500.0, 3909000.0, 350.0)]

out = Path(sys.argv[1]).resolve()
m = read_fort14(next(out.glob("*.14")))
xy, tri = m.nodes[:, :2], m.elements

cx, cy, half = VIEWS[0]
box = shapely.box(cx - half, cy - half, cx + half, cy + half)
land = gpd.read_file(os.environ["FMESH_LAND"]).to_crs(MESH_EPSG)
land = land[land.intersects(box)]
filt = gpd.read_file(out / "shoreline_filtered.shp")
filt = filt[filt.intersects(box)]

fig, axes = plt.subplots(1, 2, figsize=(17, 8.5), constrained_layout=True)
for ax, (x0, y0, h) in zip(axes, VIEWS):
    land.plot(ax=ax, color="0.75", edgecolor="0.45", lw=0.6, zorder=0)
    filt.boundary.plot(ax=ax, color="tab:green", lw=0.8, ls="--", zorder=1)
    draw_mesh(ax, xy, tri, m.open_boundaries, zorder=2)
    ax.set_xlim(x0 - h, x0 + h)
    ax.set_ylim(y0 - h, y0 + h)
    ax.set_aspect("equal")
axes[0].plot([], [], color="0.75", lw=6, label="OSM land (raw)")
axes[0].plot([], [], color="tab:green", ls="--", label="land after the width filter")
axes[0].plot([], [], color=SOLID_BOUNDARY_COLOR, lw=1.4,
             label="solid boundary (coast, quay, wall)")
axes[0].legend(loc="upper right", fontsize=9)
png = out / "walls_vs_osm.png"
fig.savefig(png, dpi=150)
print(f"wrote {png}")
