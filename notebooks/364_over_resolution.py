"""Over-resolution + thin-wall integrity check (owner 2026-07-13:
the channel-refinement experiment resolved narrow water the policy
never kept and punched through thin levees, and NO existing gate
caught either class -- the sample was better and the owner had to
find it by eye).

Measured against the ORIGINAL OSM land; the goto2023 sample is
reported for information only, never gated (owner 2026-09-20: the
sample is a reference, not ground truth, so water it omits is not a
defect).

NARROW water meshed
    The element's centroid sits in original water whose local width
    (2 x distance to the original coast) is under the element's own
    median edge length, and it is NOT inside a kept-corridor tube or
    an applied manual edit. w/h < 1 cannot carry one row and is
    advisory; w/h < 0.5 is severe and gates.

WALL crossings
    An element edge crosses original land over > 40 m of its
    length while both endpoints lie in original water: the mesh
    tunnels through a levee/wall thinner than the local cells.
    Intended corridor widenings are excluded the same way.

Both lists print with atlas refs and are drawn to
outputs/figures/over_resolution.png. Exit code 1 when severe narrow
water or WALL crossings exist outside kept tubes (gate; run after
342).
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

# ---- NARROW water meshed (was: STRAY over-resolution) ----------
# The old test called an element "stray" when goto2023 did not have
# it. goto2023 is a reference, not ground truth (owner 2026-09-20),
# so water it omits is not a defect by itself. What the policy does
# care about is water too narrow to carry a row of elements at the
# local size, meshed anyway outside a kept corridor: that is what
# produces one-cell-wide channels and bank-to-bank cells.
#
# The local WIDTH of the water must not be confused with the distance
# to the coast: every element lying against the shore of a wide bay is
# close to land while the water there is kilometres wide. Width is
# measured ACROSS the water: from the element centroid take d1 = the
# distance to the nearest land, then march along the ray pointing away
# from that land point until land is met again (d2, capped at 5 h).
# w = d1 + d2 is the channel width for a canal and "wider than the cap"
# for an element sitting on the shore of open water.
#
#   w/h < 1    cannot carry one row  -> advisory
#   w/h < 0.5  clearly too narrow    -> gate
from scipy.ndimage import distance_transform_edt  # noqa: E402

_ll = po[tri_o]
_dx = (_ll[:, [1, 2, 0], 0] - _ll[:, :, 0]) * 111e3 * np.cos(np.deg2rad(35.35))
_dy = (_ll[:, [1, 2, 0], 1] - _ll[:, :, 1]) * 111e3
h_elem = np.median(np.hypot(_dx, _dy), axis=1)

_PIX = 25.0                                   # raster step (m)
_cos = np.cos(np.deg2rad(35.35))
_x0, _y0 = po[:, 0].min() - 0.02, po[:, 1].min() - 0.02
_x1, _y1 = po[:, 0].max() + 0.02, po[:, 1].max() + 0.02
_nx = int((_x1 - _x0) * 111e3 * _cos / _PIX) + 1
_ny = int((_y1 - _y0) * 111e3 / _PIX) + 1
_gx = _x0 + (np.arange(_nx) + 0.5) * _PIX / (111e3 * _cos)
_gy = _y0 + (np.arange(_ny) + 0.5) * _PIX / 111e3
_XX, _YY = np.meshgrid(_gx, _gy)
_water = ~shapely.contains_xy(land, _XX, _YY)
_edt, _idx = distance_transform_edt(_water, sampling=_PIX, return_indices=True)
print(f"[364] width raster {_nx}x{_ny} at {_PIX:.0f} m", flush=True)


def _pix(lon, lat):
    return (int(round((lat - _gy[0]) * 111e3 / _PIX)),
            int(round((lon - _gx[0]) * 111e3 * _cos / _PIX)))


# SAMPLE THE ELEMENT, NOT ONE POINT (owner 2026-09-21, after the
# coastline fit): read from the centroid alone, the width flips
# categorically when a thin land finger pokes into the element -- three
# elements straddling a spit went from w/h = 5.1 to w/h = 0.30 because a
# 12 m node move put the centroid on the other side of the finger. The
# element is the same shape before and after; only the sample moved.
# Sample the centroid AND the three points halfway out to the vertices,
# all inside the triangle, and take the MEDIAN. A genuinely one-row
# channel reads narrow at every sample; a cell straddling a spit does not.
_samples = np.stack([cent] + [0.5 * (cent + po[tri_o[:, k]]) for k in range(3)], axis=1)


def _width_at(px, py, h):
    r, c = _pix(px, py)
    if not (0 <= r < _ny and 0 <= c < _nx) or not _water[r, c]:
        return np.nan                  # on land: 342's business
    d1 = float(_edt[r, c])
    lr, lc = int(_idx[0, r, c]), int(_idx[1, r, c])
    vr, vc = r - lr, c - lc            # away from the nearest land pixel
    norm = np.hypot(vr, vc)
    if norm == 0:
        return np.nan
    vr, vc = vr / norm, vc / norm
    cap = 5.0 * h
    d2 = cap
    for t in np.arange(_PIX, cap, _PIX):
        rr = int(round(r + vr * t / _PIX))
        cc = int(round(c + vc * t / _PIX))
        if not (0 <= rr < _ny and 0 <= cc < _nx):
            break
        if not _water[rr, cc]:
            d2 = float(t)
            break
    return d1 + d2


w_elem = np.full(len(cent), np.inf)
for i in range(len(cent)):
    ws = [_width_at(px, py, h_elem[i]) for px, py in _samples[i]]
    ws = [w for w in ws if np.isfinite(w)]
    if not ws:
        continue                       # every sample on land: 342's business
    w_elem[i] = float(np.median(ws))

narrow, narrow_wh = [], {}
for i in range(len(cent)):
    ratio = w_elem[i] / h_elem[i]
    if ratio >= 1.0:
        continue
    if tube_u.covers(shapely.Point(cent[i])):
        continue                       # intended widening / edit
    narrow.append(i)
    narrow_wh[i] = (round(float(w_elem[i]), 1), round(float(h_elem[i]), 1),
                    round(float(ratio), 3))
severe = [i for i in narrow if narrow_wh[i][2] < 0.5]
stray = narrow                         # keep the downstream name

# ---- WALL crossings --------------------------------------------
# An edge that FOLLOWS a slightly curved coastline overlaps the land
# over a long distance while barely entering it; the earlier test
# measured that overlap LENGTH and flagged 68 such edges, none of
# which tunnels through anything (owner 2026-09-20, measured: every
# one of them had a maximum penetration depth below 0.5 h). What a
# levee crossing actually looks like is
#   (a) the edge enters the land DEEPLY relative to the local cell
#       size, measured perpendicular to the coast, and
#   (b) it comes out on the far side: water on both sides of the land
#       strip it crosses, i.e. it connects two water areas through
#       land thinner than the local cells.
# Both are required here.
WALL_DEPTH_FRAC = 0.5      # penetration depth / local h
WALL_STRIP_FRAC = 1.0      # crossed land strip thickness / local h
ee = np.vstack([tri_o[:, [0, 1]], tri_o[:, [1, 2]], tri_o[:, [2, 0]]])
ee.sort(axis=1)
ee = np.unique(ee, axis=0)
segs = [shapely.LineString([po[a], po[b]]) for a, b in ee]
segtree = STRtree(segs)
land_bnd = land.boundary
hits = set()
wall_detail = {}
for k in segtree.query(lparts, predicate="intersects").T:
    li, si = int(k[0]), int(k[1])
    a, b = ee[si]
    pa, pb = shapely.Point(po[a]), shapely.Point(po[b])
    if land.covers(pa) or land.covers(pb):
        continue                       # handled by 342 on-land
    mid = segs[si].interpolate(0.5, normalized=True)
    if tube_u.covers(mid):
        continue                       # intended widening
    inter = segs[si].intersection(lparts[li])
    if inter.is_empty:
        continue
    pieces = list(inter.geoms) if hasattr(inter, "geoms") else [inter]
    h_loc = H                          # local target size (m)
    for piece in pieces:
        if piece.length <= 0:
            continue
        # (a) how deep inside the land does the edge run?
        depth = max(piece.interpolate(t, normalized=True).distance(land_bnd)
                    for t in np.linspace(0.0, 1.0, 21)) * 111e3
        if depth < WALL_DEPTH_FRAC * h_loc:
            continue
        # (b) does it come out into water again, through a thin strip?
        p0, p1 = piece.interpolate(0.0, normalized=True), \
            piece.interpolate(1.0, normalized=True)
        strip = p0.distance(p1) * 111e3           # land thickness on this line
        before = segs[si].project(p0) * 111e3
        after = (segs[si].length - segs[si].project(p1)) * 111e3
        if strip > WALL_STRIP_FRAC * h_loc:
            continue                   # not a thin wall: real land
        if min(before, after) < 1.0:
            continue                   # ends on the coast, no far side
        hits.add(si)
        wall_detail[si] = dict(depth_m=round(depth, 1), strip_m=round(strip, 1),
                               water_before_m=round(before, 1),
                               water_after_m=round(after, 1))
wall = sorted(hits)

def _clusters(pts, ids):
    from collections import defaultdict
    cl = defaultdict(list)
    for i in ids:
        cl[(round(pts[i][0] / 0.008), round(pts[i][1] / 0.008))].append(i)
    return sorted(cl.values(), key=len, reverse=True)

print(f"[364] beyond-sample elements (information only): "
      f"{int((~cov).sum())}; NARROW water meshed (w/h<1, outside kept "
      f"corridors): {len(narrow)}, of which severe (w/h<0.5): "
      f"{len(severe)}; WALL crossings: {len(wall)}", flush=True)
seg_mid = np.array([[segs[s].interpolate(0.5, normalized=True).x,
                     segs[s].interpolate(0.5, normalized=True).y]
                    for s in wall]) if wall else np.zeros((0, 2))
for name, ids, pts in (("NARROW", narrow, cent),
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
               linewidths=1.6, label=f"narrow water meshed, w/h<1 "
               f"({len(narrow)} elems)")
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

# Element / edge identities for the issue map (386), 0-indexed into the
# final mesh read above.
(OUT / "over_resolution.json").write_text(json.dumps({
    "narrow_elements": [int(i) for i in narrow],
    "narrow_detail": {str(int(i)): {"width_m": narrow_wh[i][0],
                                    "h_m": narrow_wh[i][1],
                                    "w_over_h": narrow_wh[i][2]} for i in narrow},
    "severe_elements": [int(i) for i in severe],
    "beyond_sample_elements": [int(i) for i in np.nonzero(~cov)[0]],
    "stray_elements": [int(i) for i in stray],
    "wall_edges": [[int(ee[s][0]), int(ee[s][1])] for s in wall],
    "wall_detail": {str(int(s)): wall_detail[s] for s in wall},
}, indent=1) + "\n")
print(f"[364] wrote {OUT / 'over_resolution.json'}", flush=True)

# The gate is now the two OBJECTIVE defects: a mesh edge tunnelling
# through a thin land wall, and water meshed at less than half the
# width one row needs. Plain w/h<1 water and "goto2023 does not have
# it" are reported, not gated (owner 2026-09-20).
if severe or wall:
    print(f"[364] GATE FAIL: {len(severe)} severe narrow (w/h<0.5) + "
          f"{len(wall)} wall sites", flush=True)
    sys.exit(1)
print(f"[364] GATE PASS: no wall crossings, no severe narrow water "
      f"({len(narrow)} advisory w/h<1 elements)", flush=True)
