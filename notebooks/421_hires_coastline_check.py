# Where does a `resolve` coastline depart from OSM, and by how much?
#
# The first hires run reported a chord departure of 821.6 m while its NODES
# were within 47.97 m of the curve.  A node on the curve and a chord 821 m off
# it means the walk stepped across something -- and whether that something is
# a real inlet the base never had, or a defect in the resampler, decides
# whether the number is a result or a bug.
#
#   python notebooks/421_hires_coastline_check.py <refine output dir>
#
# Environment: FMESH_LAND (the OSM shoreline the run used).
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

MESH_EPSG = 32654
OUT = Path(sys.argv[1]).resolve()
LAND = Path(os.environ["FMESH_LAND"]).resolve()
mesh = read_fort14(next(OUT.glob("*.14")))
node_map = np.load(OUT / "node_map.npy")
base = read_fort14(Path(os.environ["FMESH_BASE14"])) if os.environ.get(
    "FMESH_BASE14") else None

is_new = np.ones(mesh.n_nodes, dtype=bool)
is_new[node_map[node_map >= 0]] = False
print(f"{mesh.n_nodes:,} nodes, {int(is_new.sum()):,} new")

e = np.sort(np.vstack([mesh.elements[:, [0, 1]], mesh.elements[:, [1, 2]],
                       mesh.elements[:, [2, 0]]]), axis=1)
u, c = np.unique(e, axis=0, return_counts=True)
bnd = u[c == 1]
new_bnd = bnd[is_new[bnd].any(axis=1)]
print(f"{len(bnd):,} boundary edges, {len(new_bnd):,} touching a new node")

land = gpd.read_file(LAND).to_crs(MESH_EPSG)
pts = shapely.MultiPoint(mesh.nodes[np.unique(new_bnd)])
reach = pts.buffer(3000.0)
rings = []
for g in land.geometry:
    for r in [g.exterior, *g.interiors]:
        if shapely.intersects(reach, r):
            rings.append(shapely.LineString(np.asarray(r.coords)))
curve = shapely.MultiLineString([np.asarray(ln.coords) for ln in rings])
print(f"{len(rings)} OSM ring(s) within 3 km")

# The curves the stretches were actually cut from, when the run saved them.
# `shapely.distance(MultiPoint, line)` is the CLOSEST APPROACH of the whole
# collection, not a distance per point -- it reported 0.00 m here, which is
# the nearest node of all of them and says nothing about the others.
cur = OUT / "coastline_curves.npz"
if cur.exists():
    z = np.load(cur)
    matched = shapely.MultiLineString([z[k] for k in z.files])
    print(f"{len(z.files)} matched source curve(s) saved by the run")
else:
    matched = None
    print("no coastline_curves.npz: measuring against every nearby ring only")

a, b = mesh.nodes[new_bnd[:, 0]], mesh.nodes[new_bnd[:, 1]]
f = np.linspace(0.0, 1.0, 21)[:, None, None]
samp = (a[None] + f * (b - a)[None])                       # (21, E, 2)
d = shapely.distance(shapely.points(samp.reshape(-1, 2)), curve).reshape(21, -1)
per_edge = d.max(axis=0)
order = np.argsort(per_edge)[::-1]
dm = (shapely.distance(shapely.points(samp.reshape(-1, 2)), matched).reshape(21, -1)
      if matched is not None else None)
per_edge_m = dm.max(axis=0) if dm is not None else None
print("\nworst chords (departure, length, midpoint lon/lat):")
to_ll = Transformer.from_crs(f"EPSG:{MESH_EPSG}", "EPSG:4326", always_xy=True)
for k in order[:10]:
    mid = 0.5 * (a[k] + b[k])
    lo, la = to_ll.transform(mid[0], mid[1])
    print(f"  {per_edge[k]:8.1f} m   edge {np.linalg.norm(b[k] - a[k]):7.1f} m   "
          f"{lo:.5f}, {la:.5f}   new={is_new[new_bnd[k]].tolist()}")
nd = shapely.distance(shapely.points(mesh.nodes[np.unique(new_bnd)]), curve)
print(f"\nnode departure from any ring: max {nd.max():.2f} m, "
      f"median {np.median(nd):.2f} m")
if matched is not None:
    ndm = shapely.distance(shapely.points(mesh.nodes[np.unique(new_bnd)]), matched)
    print(f"node departure from the MATCHED curve: max {ndm.max():.2f} m, "
          f"median {np.median(ndm):.2f} m")
    print(f"chord departure from the MATCHED curve: max {per_edge_m.max():.1f} m, "
          f"median {np.median(per_edge_m):.1f} m")
    lens = np.linalg.norm(b - a, axis=1)
    print("\ndeparture against the local chord length -- fidelity is bounded "
          "by element size, so this ratio is the number that travels:")
    for k in np.argsort(per_edge_m)[::-1][:6]:
        print(f"  chord {lens[k]:7.1f} m  departs {per_edge_m[k]:7.1f} m  "
              f"= {100 * per_edge_m[k] / lens[k]:5.1f} % of its own length")
print(f"chord departure: max {per_edge.max():.1f} m, p99 "
      f"{np.percentile(per_edge, 99):.1f} m, median {np.median(per_edge):.1f} m")
print(f"chords over 100 m from OSM: {int((per_edge > 100).sum())} of {len(per_edge)}")

# Where the violations the patch introduced actually are.  "Two QA failures"
# is a count; whether they sit ON the resolved coastline is the question that
# decides whether `resolve` caused them.
qa_path = next(OUT.glob("*_qa.json"), None)
if qa_path is not None:
    import json

    qa = json.loads(qa_path.read_text())
    bnd_nodes = set(np.unique(new_bnd).tolist())
    coast_line = shapely.MultiLineString(
        [shapely.LineString(mesh.nodes[[i, j]]) for i, j in new_bnd.tolist()])
    for chk in qa.get("checks", []):
        if chk.get("status") != "fail":
            continue
        ids = chk.get("offender_ids") or chk.get("offenders") or []
        here = []
        for off in ids[:400]:
            els = off.get("elements") or ([off["id"]] if off.get("kind") == "element"
                                          else [])
            if els:
                nn = mesh.elements[np.asarray(els, dtype=int)].ravel()
            elif off.get("kind") == "node":
                nn = np.array([off["id"]], dtype=int)
            else:
                continue
            c = mesh.nodes[nn].mean(axis=0)
            here.append((float(shapely.distance(shapely.Point(c), coast_line)),
                         bool(set(nn.tolist()) & bnd_nodes)))
        if here:
            d = np.array([x[0] for x in here])
            on = sum(x[1] for x in here)
            print(f"\n{chk['check_id']}: {len(here)} offender(s) placed; "
                  f"{on} touch the resolved coastline; distance to it "
                  f"min {d.min():.0f} m, median {np.median(d):.0f} m")

k = int(order[0])
mid = 0.5 * (a[k] + b[k])
half = max(2500.0, 3.0 * float(np.linalg.norm(b[k] - a[k])))
fig, ax = plt.subplots(figsize=(11, 9))
# The mesh, always: a coastline figure without the elements beside it cannot
# show whether the departure is a coastline problem or a spacing problem.
ax.triplot(mesh.nodes[:, 0], mesh.nodes[:, 1], mesh.elements,
           lw=0.25, color="0.75")
for ln in rings:
    xy = np.asarray(ln.coords)
    ax.plot(xy[:, 0], xy[:, 1], color="tab:blue", lw=1.2, zorder=3)
ax.plot([], [], color="tab:blue", lw=1.2, label="OSM shoreline")
for i, j in new_bnd:
    ax.plot(mesh.nodes[[i, j], 0], mesh.nodes[[i, j], 1],
            color="tab:red", lw=1.4, zorder=4)
ax.plot([], [], color="tab:red", lw=1.4, label="resolved coastline (mesh)")
ax.plot([a[k, 0], b[k, 0]], [a[k, 1], b[k, 1]], color="k", lw=3.0, zorder=5,
        label=f"worst chord: {per_edge[k]:.0f} m from OSM")
ax.set_xlim(mid[0] - half, mid[0] + half)
ax.set_ylim(mid[1] - half, mid[1] + half)
ax.set_aspect("equal")
ax.legend(loc="upper right", fontsize=9)
ax.set_title(f"{OUT.name}: resolved coastline against OSM")
fig.tight_layout()
fig.savefig(OUT / "coastline_vs_osm.png", dpi=150)
print(f"\nwrote {OUT / 'coastline_vs_osm.png'}")
