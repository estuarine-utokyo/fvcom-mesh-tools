# Size a refinement site before declaring it, the way banzu_nori.yaml was.
#
#   python notebooks/423_site_survey.py <x> <y> <radius_m> [crs]
#
# `crs` is `utm` (EPSG:32654, the default) or `ll`.  Reports what a recipe
# writer needs to know and cannot guess: how complex the shoreline is, how
# coarse the base mesh is there, how deep the water is, and which rung of the
# depth ladder will answer.
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

from fvcom_mesh_tools.dem import tokyo_bay as tb  # noqa: E402
from fvcom_mesh_tools.io.fvcom_native import read_fvcom_case  # noqa: E402

MESH_EPSG = 32654
x, y, R = float(sys.argv[1]), float(sys.argv[2]), float(sys.argv[3])
crs = sys.argv[4] if len(sys.argv) > 4 else "utm"
to_ll = Transformer.from_crs(f"EPSG:{MESH_EPSG}", "EPSG:4326", always_xy=True)
to_m = Transformer.from_crs("EPSG:4326", f"EPSG:{MESH_EPSG}", always_xy=True)
if crs == "ll":
    lon, lat = x, y
    x, y = to_m.transform(lon, lat)
else:
    lon, lat = to_ll.transform(x, y)
print(f"centre {lon:.5f}, {lat:.5f}  ({x:.0f}, {y:.0f} in UTM 54N), radius {R:g} m")

G = Path(os.environ.get("FMESH_BASE_GRD", os.path.expanduser(
    "~/Github/TB-FVCOM/input/goto2023/grid/TokyoBay_grd.dat")))
D = Path(os.environ.get("FMESH_BASE_DEP", os.path.expanduser(
    "~/Github/TB-FVCOM/input/goto2023/grid/TokyoBay_dep_m7001tp_rfac0p2_cap300.dat")))
base = read_fvcom_case(G, D, None)
disc = shapely.Point(x, y).buffer(R, quad_segs=96)

LAND = Path(os.environ["FMESH_LAND"]).resolve()
land = gpd.read_file(LAND).to_crs(MESH_EPSG)
rings = []
for g in land.geometry:
    for r in [g.exterior, *g.interiors]:
        ln = shapely.LineString(np.asarray(r.coords))
        if shapely.intersects(disc, ln):
            rings.append(ln)
inside = shapely.intersection(shapely.MultiLineString(rings), disc) if rings \
    else shapely.MultiLineString([])
verts = shapely.get_coordinates(inside)
print(f"\\nOSM shoreline inside: {len(rings)} ring(s) touched, "
      f"{inside.length / 1000:.2f} km of line, {len(verts):,} vertices "
      f"({inside.length / max(len(verts) - 1, 1):.1f} m per vertex)")
# How tight the tightest feature is: the largest erosion that keeps water.
water = shapely.difference(disc, shapely.union_all(
    [g for g in land.geometry if shapely.intersects(disc, g)]))
print(f"water inside the disc: {water.area / 1e6:.4f} km2 of "
      f"{disc.area / 1e6:.4f} km2 ({100 * water.area / disc.area:.0f} %)")
for w in (5.0, 10.0, 20.0, 40.0, 80.0):
    er = water.buffer(-w)
    print(f"  water surviving a {w:5.1f} m erosion: "
          f"{er.area / 1e6:.4f} km2, {len(getattr(er, 'geoms', [er]))} part(s)")

e = np.sort(np.vstack([base.elements[:, [0, 1]], base.elements[:, [1, 2]],
                       base.elements[:, [2, 0]]]), axis=1)
u, c = np.unique(e, axis=0, return_counts=True)
bnd = u[c == 1]
bxy = base.nodes[:, :2]
mid = 0.5 * (bxy[bnd[:, 0]] + bxy[bnd[:, 1]])
near = shapely.contains(disc, shapely.points(mid[:, 0], mid[:, 1]))
lens = np.linalg.norm(bxy[bnd[near, 0]] - bxy[bnd[near, 1]], axis=1)
print(f"\\nbase mesh coastline inside: {int(near.sum())} edge(s)"
      + (f", {lens.min():.0f}-{lens.max():.0f} m, median {np.median(lens):.0f} m"
         if near.any() else " -- the base has no boundary here"))
cen = base.nodes[base.elements][:, :, :2].mean(axis=1)
sel = shapely.contains(disc, shapely.points(cen[:, 0], cen[:, 1]))
if sel.any():
    s = base.elements[sel]
    el = np.concatenate([np.linalg.norm(bxy[s[:, (i + 1) % 3]] - bxy[s[:, i]], axis=1)
                         for i in range(3)])
    print(f"base elements inside: {int(sel.sum())}, edges median {np.median(el):.0f} m, "
          f"depths {base.depths[np.unique(s)].min():.2f}-"
          f"{base.depths[np.unique(s)].max():.2f} m")

k = int(np.ceil(np.sqrt(4000)))
gx, gy = np.meshgrid(np.linspace(x - R, x + R, k), np.linspace(y - R, y + R, k))
m = shapely.contains(water, shapely.points(gx.ravel(), gy.ravel()))
if m.any():
    slo, sla = to_ll.transform(gx.ravel()[m], gy.ravel()[m])
    d, rung, dist = tb.sample(slo, sla)
    print(f"\\nladder over the water: {int(m.sum())} sample(s); "
          + ", ".join(f"{n} {int((rung == i).sum())}"
                      for i, n in enumerate(tb.LADDER))
          + f", extrapolated {int((rung < 0).sum())}")
    print(f"  depth {d.min():.2f} to {d.max():.2f} m, median {np.median(d):.2f}; "
          f"{100 * (d <= 0).mean():.0f} % at or above the datum")
    sd = tb.sounding_distance(slo, sla)
    print(f"  nearest real M7001 sounding: median {np.median(sd):.0f} m, "
          f"max {sd.max():.0f} m")

fig, ax = plt.subplots(figsize=(10, 9))
ax.triplot(bxy[:, 0], bxy[:, 1], base.elements, lw=0.3, color="0.75")
for ln in rings:
    q = np.asarray(ln.coords)
    ax.plot(q[:, 0], q[:, 1], color="tab:blue", lw=1.0)
ax.plot([], [], color="tab:blue", lw=1.0, label="OSM shoreline")
ax.plot([], [], color="0.75", lw=0.3, label="base mesh")
q = np.asarray(disc.exterior.coords)
ax.plot(q[:, 0], q[:, 1], color="tab:green", lw=1.8, ls="--",
        label=f"candidate region (r = {R:g} m)")
ax.set_xlim(x - 1.6 * R, x + 1.6 * R)
ax.set_ylim(y - 1.6 * R, y + 1.6 * R)
ax.set_aspect("equal")
ax.legend(loc="upper right", fontsize=9)
ax.set_title(f"site survey {lon:.4f}, {lat:.4f}  r = {R:g} m")
fig.tight_layout()
out = Path(os.environ.get("LR_OUT", "outputs")) / f"site_{lon:.4f}_{lat:.4f}_{R:g}.png"
out.parent.mkdir(parents=True, exist_ok=True)
fig.savefig(out, dpi=150)
print(f"\\nwrote {out}")
