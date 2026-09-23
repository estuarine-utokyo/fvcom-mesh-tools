# What did the width filter take from a harbour, and was it structures?
#
# The Kimitsu port run filtered OSM at 3 x h0 = 90 m.  The owner saw that the
# small harbour at the lower left of the region -- breakwaters, piers, a basin
# -- was not resolved though a 30 m mesh looks able to.  This measures what the
# filter removed there, piece by piece: its length, its mean width and whether
# it is LINEAR (a breakwater) or an AREA (a spit, a basin).
#
#   python notebooks/425_structure_survey.py <lon> <lat> <radius_m> <h0>
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

lon, lat, R, h0 = (float(a) for a in sys.argv[1:5])
to_m = Transformer.from_crs("EPSG:4326", "EPSG:32654", always_xy=True)
x, y = to_m.transform(lon, lat)
disc = shapely.Point(x, y).buffer(R, quad_segs=64)
land = gpd.read_file(os.environ["FMESH_LAND"]).to_crs(32654)
keep = [g for g in land.geometry if shapely.intersects(disc.buffer(500.0), g)]
before = shapely.union_all(keep)
print(f"region {lon}, {lat}, r {R:g} m, h0 {h0:g} m; {len(keep)} land polygon(s)")


def pieces(g):
    return [q for q in getattr(g, "geoms", [g]) if not q.is_empty and q.area > 1.0]


def describe(q):
    rr = shapely.minimum_rotated_rectangle(q)
    c = np.asarray(rr.exterior.coords)[:4]
    s = sorted([np.linalg.norm(c[1] - c[0]), np.linalg.norm(c[2] - c[1])])
    length = s[1]
    width = q.area / length if length > 0 else 0.0
    return length, width


for mult in (1.0, 2.0, 3.0):
    r = 0.5 * mult * h0
    opened = before.buffer(-r).buffer(r)
    closed = opened.buffer(r).buffer(-r)
    land_lost = [q for q in pieces(shapely.difference(before, opened))
                 if shapely.intersects(disc, q)]
    water_lost = [q for q in pieces(shapely.difference(closed, opened))
                  if shapely.intersects(disc, q)]
    lin = [(describe(q), q) for q in land_lost]
    linear = [t for t in lin if t[0][0] > 5 * max(t[0][1], 1e-6) and t[0][0] > h0]
    print(f"\n--- filter at {mult:g} x h0 (removes < {2 * r:.0f} m)")
    print(f"  land removed in region: {len(land_lost)} piece(s), "
          f"{sum(q.area for q in land_lost):.0f} m2; LINEAR (length > 5 x width "
          f"and > h0): {len(linear)}")
    for (L, W), q in sorted(linear, key=lambda t: -t[0][0])[:12]:
        cx, cy = q.centroid.x, q.centroid.y
        print(f"    linear piece: length {L:6.0f} m, mean width {W:5.1f} m, "
              f"at {cx:.0f}, {cy:.0f}")
    print(f"  water closed in region: {len(water_lost)} piece(s), "
          f"{sum(q.area for q in water_lost):.0f} m2")
    for q in sorted(water_lost, key=lambda q: -q.area)[:6]:
        L, W = describe(q)
        print(f"    closed water: area {q.area:8.0f} m2, length {L:5.0f} m, "
              f"mean width {W:5.1f} m, at {q.centroid.x:.0f}, {q.centroid.y:.0f}")

r = 1.5 * h0
opened = before.buffer(-r).buffer(r)
lost = [q for q in pieces(shapely.difference(before, opened)) if shapely.intersects(disc, q)]
mesh = read_fort14(next(Path(os.environ["FMESH_MESH_DIR"]).glob("*.14")))
fig, axes = plt.subplots(1, 2, figsize=(17, 8.5))
for ax, half in zip(axes, (1.2 * R, 0.45 * R)):
    cx, cy = (x, y) if half > 0.5 * R else (x - 0.55 * R, y - 0.55 * R)
    ax.triplot(mesh.nodes[:, 0], mesh.nodes[:, 1], mesh.elements, lw=0.25, color="0.7")
    for g in pieces(before):
        for ring in [g.exterior, *g.interiors]:
            q = np.asarray(ring.coords)
            ax.plot(q[:, 0], q[:, 1], color="tab:blue", lw=1.0)
    for g in lost:
        for part in getattr(g, "geoms", [g]):
            q = np.asarray(part.exterior.coords)
            ax.fill(q[:, 0], q[:, 1], color="tab:red", alpha=0.6, lw=0)
    q = np.asarray(disc.exterior.coords)
    ax.plot(q[:, 0], q[:, 1], "--", color="tab:green", lw=1.5)
    ax.set_xlim(cx - half, cx + half)
    ax.set_ylim(cy - half, cy + half)
    ax.set_aspect("equal")
ax = axes[0]
ax.plot([], [], color="tab:blue", label="OSM land")
ax.fill([], [], color="tab:red", alpha=0.6, label=f"removed by the {2 * r:.0f} m filter")
ax.plot([], [], color="0.7", label="delivered mesh")
ax.legend(loc="upper right", fontsize=9)
axes[0].set_title("the region")
axes[1].set_title("lower-left harbour")
fig.tight_layout()
out = Path(os.environ["FMESH_MESH_DIR"]) / "structures_removed.png"
fig.savefig(out, dpi=150)
print(f"\nwrote {out}")
