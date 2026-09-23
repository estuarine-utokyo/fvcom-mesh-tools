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

# Two panels, because the two questions are at two scales: the whole patch
# says whether the coastline was followed at all, and the declared region says
# whether it was followed at the size that was asked for.  The first figure of
# this run showed only the first, and its worst chord was 6 km from the core
# in a part of the transition the recipe never claimed to refine.
k = int(order[0])
mid = 0.5 * (a[k] + b[k])
import json as _json  # noqa: E402

rep = _json.loads((OUT / "report.json").read_text())
reg = rep["regions"][0]
rxy = np.asarray(reg["xy"], dtype=float)
rc_ = rxy.mean(axis=0)
rr = float(np.linalg.norm(rxy - rc_, axis=1).max())

bad = []
if qa_path is not None:
    for chk in qa.get("checks", []):
        if chk.get("status") != "fail" or chk["check_id"] != "c1_min_angle":
            continue
        # The same fallback the count above uses: a report may carry the
        # complete `offender_ids` or only the capped display list, and reading
        # one of the two left the figure with nothing to mark.
        for off in (chk.get("offender_ids") or chk.get("offenders") or [])[:50]:
            els = off.get("elements") or ([off["id"]] if off.get("kind") == "element"
                                          else [])
            if els:
                c = mesh.nodes[mesh.elements[np.asarray(els, int)].ravel()].mean(axis=0)
                bad.append(c)
                lo, la = to_ll.transform(c[0], c[1])
                print(f"  c1 offender at {lo:.5f}, {la:.5f} "
                      f"(element {els[0]})")

fig, axes = plt.subplots(1, 2, figsize=(17, 8.5))
for ax, (cx, cy, half, ttl) in zip(axes, [
        (mid[0], mid[1], max(3000.0, 3.0 * float(np.linalg.norm(b[k] - a[k]))),
         f"whole patch: worst chord {per_edge[k]:.0f} m from OSM"),
        (rc_[0], rc_[1], 2.2 * rr,
         f"the declared region ({reg['name']}, target {reg['target_h_m']:g} m)")]):
    # The mesh, always: a coastline figure without the elements beside it
    # cannot show whether a departure is the coastline or the spacing.
    ax.triplot(mesh.nodes[:, 0], mesh.nodes[:, 1], mesh.elements,
               lw=0.2, color="0.8")
    for ln in rings:
        q = np.asarray(ln.coords)
        ax.plot(q[:, 0], q[:, 1], color="tab:blue", lw=1.1, zorder=3)
    for i, j in new_bnd:
        ax.plot(mesh.nodes[[i, j], 0], mesh.nodes[[i, j], 1],
                color="tab:red", lw=1.3, zorder=4)
    ax.plot(np.append(rxy[:, 0], rxy[0, 0]), np.append(rxy[:, 1], rxy[0, 1]),
            color="tab:green", lw=1.6, ls="--", zorder=5)
    if bad:
        B = np.asarray(bad)
        ax.plot(B[:, 0], B[:, 1], "x", color="tab:orange", ms=11, mew=2.2,
                zorder=6)
    ax.set_xlim(cx - half, cx + half)
    ax.set_ylim(cy - half, cy + half)
    ax.set_aspect("equal")
    ax.set_title(ttl, fontsize=11)
axes[0].plot([a[k, 0], b[k, 0]], [a[k, 1], b[k, 1]], color="k", lw=3.0, zorder=7)
for ax, lab in zip(axes, [True, False]):
    ax.plot([], [], color="tab:blue", lw=1.1, label="OSM shoreline")
    ax.plot([], [], color="tab:red", lw=1.3, label="resolved coastline (mesh)")
    ax.plot([], [], color="tab:green", lw=1.6, ls="--", label="declared region")
    if bad:
        ax.plot([], [], "x", color="tab:orange", ms=9, mew=2,
                label=f"min-angle offender ({len(bad)})")
    ax.legend(loc="upper right", fontsize=8)
fig.suptitle(f"{OUT.name}: resolved coastline against OSM", fontsize=13)
fig.tight_layout()
fig.savefig(OUT / "coastline_vs_osm.png", dpi=150)
print(f"\nwrote {OUT / 'coastline_vs_osm.png'}")
