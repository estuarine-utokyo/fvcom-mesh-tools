# The walls extract_walls would build from a harbour's OSM land, drawn over
# the source and the current mesh, before any of it goes near a mesher.
#
#   python notebooks/428_walls_preview.py <lon> <lat> <radius_m> <h0> <mesh dir>
import json
import os
import sys
from pathlib import Path

import geopandas as gpd
import matplotlib
import numpy as np
import shapely

matplotlib.use("Agg")
import matplotlib.pyplot as plt  # noqa: E402
from pyproj import Transformer  # noqa: E402

from fvcom_mesh_tools.io.fort14 import read_fort14  # noqa: E402
from fvcom_mesh_tools.patch import filter_shoreline  # noqa: E402
from fvcom_mesh_tools.walls import extract_walls  # noqa: E402

lon, lat, R, h0 = (float(a) for a in sys.argv[1:5])
mdir = Path(sys.argv[5])
x, y = Transformer.from_crs("EPSG:4326", "EPSG:32654", always_xy=True).transform(lon, lat)
disc = shapely.Point(x, y).buffer(R, quad_segs=64)
land = gpd.read_file(os.environ["FMESH_LAND"]).to_crs(32654)
keep = [g for g in land.geometry if shapely.intersects(disc.buffer(2000.0), g)]
area, frep = filter_shoreline(keep, h0, elements_per_feature=2)
walls, wrep = extract_walls(keep, area, h0)
inside = [w for w in walls if shapely.intersects(w, disc)]
print(json.dumps({k: v for k, v in frep.items()}, indent=1))
print(json.dumps({k: v for k, v in wrep.items() if k != "dropped"}, indent=1))
print(f"walls touching the region: {len(inside)}")
for w in sorted(inside, key=lambda w: -w.length):
    c = np.asarray(w.coords)
    print(f"  {w.length:7.1f} m, {len(c)} vertices, from ({c[0, 0]:.0f}, {c[0, 1]:.0f}) "
          f"to ({c[-1, 0]:.0f}, {c[-1, 1]:.0f})")
mesh = read_fort14(next(mdir.glob("*.14")))
fig, axes = plt.subplots(1, 2, figsize=(17, 8.5))
for ax, (cx, cy, half) in zip(axes, [(x, y, 1.2 * R), (x - 0.55 * R, y - 0.55 * R, 0.45 * R)]):
    ax.triplot(mesh.nodes[:, 0], mesh.nodes[:, 1], mesh.elements, lw=0.25, color="0.75")
    for g in keep:
        for part in getattr(g, "geoms", [g]):
            for ring in [part.exterior, *part.interiors]:
                q = np.asarray(ring.coords)
                ax.plot(q[:, 0], q[:, 1], color="tab:blue", lw=0.9)
    for part in getattr(area, "geoms", [area]):
        q = np.asarray(part.exterior.coords)
        ax.fill(q[:, 0], q[:, 1], color="tab:blue", alpha=0.15, lw=0)
    for w in walls:
        q = np.asarray(w.coords)
        ax.plot(q[:, 0], q[:, 1], color="tab:red", lw=2.2)
    q = np.asarray(disc.exterior.coords)
    ax.plot(q[:, 0], q[:, 1], "--", color="tab:green", lw=1.4)
    ax.set_xlim(cx - half, cx + half)
    ax.set_ylim(cy - half, cy + half)
    ax.set_aspect("equal")
axes[0].plot([], [], color="tab:blue", label="OSM land")
axes[0].fill([], [], color="tab:blue", alpha=0.15, label=f"kept as land (>= {2 * h0:g} m)")
axes[0].plot([], [], color="tab:red", lw=2.2, label="walls")
axes[0].plot([], [], color="0.75", label="current mesh")
axes[0].legend(loc="upper right", fontsize=9)
axes[0].set_title("the region")
axes[1].set_title("lower-left harbour")
fig.tight_layout()
out = mdir / "walls_preview.png"
fig.savefig(out, dpi=150)
print(f"wrote {out}")
