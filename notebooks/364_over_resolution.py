"""Over-resolution + thin-wall integrity check (owner 2026-07-13:
the channel-refinement experiment resolved narrow water the policy
never kept and punched through thin levees, and NO existing gate
caught either class -- the sample was better and the owner had to
find it by eye).

Two element classes are measured against the ORIGINAL OSM land and
the goto2023 sample:

STRAY over-resolution
    Our element lies outside the sample triangulation, its centroid
    is in ORIGINAL water narrower than 0.5 h (distance to original
    land), and it is NOT inside a kept-corridor tube or an applied
    manual edit. The policy never asked for this water -- target 0.

WALL crossings
    An element edge crosses original land over > 40 m of its
    length while both endpoints lie in original water: the mesh
    tunnels through a levee/wall thinner than the local cells.
    Intended corridor widenings are excluded the same way.

Both lists print with atlas refs and are drawn to
outputs/figures/over_resolution.png. Exit code 1 when STRAY or
WALL sites exist outside kept tubes (gate; run after 342).
"""

import json
import os
import sys
from pathlib import Path

import geopandas as gpd
import matplotlib
matplotlib.use("Agg")
import matplotlib.pyplot as plt
import numpy as np
import shapely
from pyproj import Transformer
from shapely.ops import unary_union
from shapely.strtree import STRtree

sys.path.insert(0, str(Path(__file__).resolve().parents[1] / "src"))
from fvcom_mesh_tools.gridref import TOKYO_BAY_GRID
from fvcom_mesh_tools.plotting import add_atlas_grid, use_readable_style

use_readable_style()
ROOT = Path(__file__).resolve().parents[1]
OUT = ROOT / "outputs/sample_repro"
FIG = ROOT / "outputs/figures"
H = 348.0


def load14(p):
    ls = Path(p).read_text().splitlines()
    ne, nn = map(int, ls[1].split()[:2])
    nod = np.array([ls[2 + i].split()[1:3] for i in range(nn)], float)
    tri = np.array([ls[2 + nn + i].split()[2:5]
                    for i in range(ne)], int) - 1
    return nod, tri


tr = Transformer.from_crs(32654, 4326, always_xy=True)
nod_s, tri_s = load14(OUT / "sample_original.14")
nod_o, tri_o = load14(OUT / "sample_repro_final.14")
lon_s, lat_s = tr.transform(nod_s[:, 0], nod_s[:, 1])
ps = np.column_stack([lon_s, lat_s])
lon_o, lat_o = tr.transform(nod_o[:, 0], nod_o[:, 1])
po = np.column_stack([lon_o, lat_o])
cent = po[tri_o].mean(axis=1)

stree = STRtree([shapely.Polygon(ps[t]) for t in tri_s])
cov = np.zeros(len(cent), bool)
cov[stree.query(shapely.points(cent[:, 0], cent[:, 1]),
                predicate="within")[0]] = True

land = unary_union(list(gpd.read_file(
    ROOT / "outputs/tb_varres_3r/land_osm_wide.shp").geometry))
lparts = list(getattr(land, "geoms", [land]))
ltree = STRtree(lparts)

# intended-resolution areas: kept corridor tubes + manual edits
tubes = []
wjson = OUT / "waterways.json"
if wjson.exists():
    for r in json.loads(wjson.read_text()):
        if r["action"] != "keep":
            continue
        pairs = r.get("arcs_done") or (
            [[r["arc"], r["widths_m"]]] if r.get("arc") else [])
        for a2, w2 in pairs:
            a2 = np.asarray(a2, float)
            w2 = np.asarray(w2, float)
            tubes.append(shapely.LineString(a2).buffer(
                0.65 * float(np.median(w2)) / 111e3))
_SR_EXCL = {s.strip() for s in os.environ.get(
    "SR_EDITS_EXCLUDE", "").split(",") if s.strip()}
for ed in sorted((ROOT / "recipes/edits/sample_repro").glob("*.json")):
    if ed.stem in _SR_EXCL:
        print(f"[ovr] edit {ed.stem}: EXCLUDED "
              f"(SR_EDITS_EXCLUDE)", flush=True)
        continue
    d = json.loads(ed.read_text())
    if d.get("type") in ("water_patch",) and d.get("polygon"):
        tubes.append(shapely.Polygon(d["polygon"]))
    elif d.get("arc"):
        a2 = np.asarray(d["arc"], float)
        w2 = float(np.max(d.get("width_m", 600.0))) \
            if np.ndim(d.get("width_m", 600.0)) else \
            float(d.get("width_m", 600.0))
        tubes.append(shapely.LineString(a2).buffer(0.8 * w2 / 111e3))
# choke widen-ops (finish stage): pushed banks are intended
wops_f = OUT / "widen_ops.json"
if wops_f.exists():
    from pyproj import Transformer as _Tr
    _tr2 = _Tr.from_crs(32654, 4326, always_xy=True)
    for op in json.loads(wops_f.read_text()):
        for o, n in zip(op["old"], op["new"]):
            ox, oy = _tr2.transform(*o)
            nx, ny = _tr2.transform(*n)
            tubes.append(shapely.LineString(
                [(ox, oy), (nx, ny)]).buffer(
                0.7 * op["h_loc"] / 111e3))
tube_u = unary_union(tubes) if tubes else shapely.Polygon()

# ---- STRAY over-resolution -------------------------------------
un = np.nonzero(~cov)[0]
upts = [shapely.Point(cent[i]) for i in un]
stray = []
for i, p in zip(un, upts):
    dl = lparts[ltree.nearest(p)].distance(p) * 111e3
    if dl >= 0.5 * H:
        continue                       # open-water class
    if tube_u.covers(p):
        continue                       # intended widening
    stray.append(i)

# ---- WALL crossings --------------------------------------------
# coastline-following edges touch the ORIGINAL land boundary over
# their whole length, so test against a 15 m ERODED land core:
# only edges genuinely tunnelling through a wall keep >25 m of
# intersection with the core.
ee = np.vstack([tri_o[:, [0, 1]], tri_o[:, [1, 2]], tri_o[:, [2, 0]]])
ee.sort(axis=1)
ee = np.unique(ee, axis=0)
segs = [shapely.LineString([po[a], po[b]]) for a, b in ee]
segtree = STRtree(segs)
er = [g.buffer(-15.0 / 111e3) for g in lparts]
er = [g for g in er if not g.is_empty]
hits = set()
for k in segtree.query(er, predicate="intersects").T:
    li, si = int(k[0]), int(k[1])
    a, b = ee[si]
    pa, pb = shapely.Point(po[a]), shapely.Point(po[b])
    if land.covers(pa) or land.covers(pb):
        continue                       # handled by 342 on-land
    ln = segs[si].intersection(er[li]).length * 111e3
    if ln < 25.0:
        continue
    mid = segs[si].interpolate(0.5, normalized=True)
    if tube_u.covers(mid):
        continue                       # intended widening
    hits.add(si)
wall = sorted(hits)

def _clusters(pts, ids):
    from collections import defaultdict
    cl = defaultdict(list)
    for i in ids:
        cl[(round(pts[i][0] / 0.008), round(pts[i][1] / 0.008))].append(i)
    return sorted(cl.values(), key=len, reverse=True)

print(f"[364] beyond-sample elements: {int((~cov).sum())}; "
      f"STRAY over-resolution: {len(stray)}; "
      f"WALL crossings: {len(wall)}", flush=True)
seg_mid = np.array([[segs[s].interpolate(0.5, normalized=True).x,
                     segs[s].interpolate(0.5, normalized=True).y]
                    for s in wall]) if wall else np.zeros((0, 2))
for name, ids, pts in (("STRAY", stray, cent),
                       ("WALL", list(range(len(wall))), seg_mid)):
    for c in _clusters(pts, ids)[:10]:
        p = pts[c].mean(axis=0) if len(c) > 1 else pts[c[0]]
        print(f"[364]   {name} x{len(c):3d} at "
              f"{TOKYO_BAY_GRID.point_to_subcell(p[0], p[1])} "
              f"({p[0]:.4f},{p[1]:.4f})", flush=True)

fig, ax = plt.subplots(figsize=(13, 15))
gpd.GeoSeries(lparts, crs="EPSG:4326").plot(
    ax=ax, color="0.88", edgecolor="0.65", linewidth=0.4)
ax.triplot(po[:, 0], po[:, 1], tri_o, color="steelblue",
           linewidth=0.25, alpha=0.6)
if stray:
    ax.scatter(cent[stray, 0], cent[stray, 1], s=46,
               facecolors="none", edgecolors="crimson",
               linewidths=1.6, label=f"stray over-resolution "
               f"({len(stray)} elems)")
if len(seg_mid):
    ax.scatter(seg_mid[:, 0], seg_mid[:, 1], s=70, marker="s",
               facecolors="none", edgecolors="darkorange",
               linewidths=1.8,
               label=f"wall crossings ({len(seg_mid)} edges)")
ax.set_xlim(139.60, 140.13)
ax.set_ylim(34.95, 35.75)
ax.set_aspect(1.0 / np.cos(np.deg2rad(35.35)))
add_atlas_grid(ax, crs="EPSG:4326")
ax.legend(loc="lower right")
ax.set_title("Over-resolution & thin-wall integrity vs ORIGINAL "
             "OSM land / goto2023 sample")
FIG.mkdir(parents=True, exist_ok=True)
fig.savefig(FIG / "over_resolution.png", dpi=150,
            bbox_inches="tight")
print(f"[364] figure -> {FIG / 'over_resolution.png'}", flush=True)

if stray or wall:
    print(f"[364] GATE FAIL: {len(stray)} stray + {len(wall)} "
          f"wall sites", flush=True)
    sys.exit(1)
print("[364] GATE PASS: no stray resolution, no wall crossings",
      flush=True)
