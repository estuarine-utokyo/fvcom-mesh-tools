# Side-by-side comparison: goto2023 SAMPLE vs our reproduction
# (outputs/sample_repro). Full-bay + mouth views with atlas grid,
# OBC in red, and a banded edge-length profile comparison.
import os
import numpy as np
os.environ.setdefault("MPLBACKEND", "Agg")
import matplotlib.pyplot as plt
from pyproj import Transformer
from fvcom_mesh_tools.plotting import _add_coast, add_atlas_grid

DEG = 1.0 / 111e3
OUT = "outputs/figures"
os.makedirs(OUT, exist_ok=True)

# --- sample (EPSG:32654 -> lonlat) ---
gd = open(os.path.expanduser(
    '~/Github/TB-FVCOM/input/goto2023/grid/TokyoBay_grd.dat')).read().split('\n')
nn = int(gd[0].split('=')[1]); ne = int(gd[1].split('=')[1])
Ts = np.array([[int(w) for w in gd[2+i].split()[1:4]]
               for i in range(ne)]) - 1
Ps = np.array([[float(w) for w in gd[2+ne+i].split()[1:3]]
               for i in range(nn)])
obc_s = [int(l.split()[1]) - 1 for l in open(os.path.expanduser(
    '~/Github/TB-FVCOM/input/goto2023/grid/TokyoBay_obc.dat')
    ).read().strip().split('\n')[1:]]
tr = Transformer.from_crs("EPSG:32654", "EPSG:4326", always_xy=True)
lon_s, lat_s = tr.transform(Ps[:, 0], Ps[:, 1])

# --- ours (lonlat) ---
P = np.load("outputs/sample_repro/p.npy")
T = np.load("outputs/sample_repro/t.npy")
# open-boundary node strings from the fort.14
obc_o = []
with open("outputs/sample_repro/sample_repro.14") as f:
    lines = f.read().split("\n")
i = 2 + int(lines[1].split()[0]) + int(lines[1].split()[1])
nope = int(lines[i].split()[0]); i += 2
for _ in range(nope):
    n = int(lines[i].split()[0]); i += 1
    obc_o.extend(int(lines[i + k].split()[0]) - 1 for k in range(n))
    i += n
print(f"sample NP={nn:,} NE={ne:,} obc={len(obc_s)} | "
      f"ours NP={len(P):,} NE={len(T):,} obc={len(obc_o)}", flush=True)

COAST = (139.0, 34.5, 141.3, 36.2)
views = {
    "full": (139.57, 34.93, 140.15, 35.78, 0.35),
    "mouth": (139.60, 34.95, 140.05, 35.30, 0.55),
    "tama": (139.73, 35.48, 139.87, 35.58, 0.7),
}
for name, (x0, y0, x1, y1, lw) in views.items():
    fig, axes = plt.subplots(
        1, 2, figsize=(19, 9.5 * (y1 - y0) * 1.22 / (x1 - x0)))
    for ax, (lo, la, tt, ob, np_, ne_) in zip(axes, [
            (lon_s, lat_s, Ts, obc_s, nn, ne),
            (P[:, 0], P[:, 1], T, obc_o, len(P), len(T))]):
        _add_coast(ax, COAST, "EPSG:4326")
        ax.triplot(lo, la, tt, lw=lw, color="steelblue")
        if ob:
            ax.plot(np.asarray(lo)[ob], np.asarray(la)[ob], color="red",
                    lw=2.5, zorder=6, marker="o", ms=3,
                    label=f"OBC ({len(ob)} nodes)")
        ax.set_xlim(x0, x1); ax.set_ylim(y0, y1)
        add_atlas_grid(ax, crs="EPSG:4326")
        ax.set_aspect(1 / np.cos(np.deg2rad(0.5 * (y0 + y1))))
        ax.legend(loc="lower right")
        ax.set_title(f"NP={np_:,} NE={ne_:,}")
    axes[0].set_title("goto2023 SAMPLE  " + axes[0].get_title())
    axes[1].set_title("repro (faithful stack)  " + axes[1].get_title())
    fig.savefig(f"{OUT}/sample_cmp_{name}.png", dpi=200,
                bbox_inches="tight")
    print("saved", name, flush=True)

# --- banded edge-length profiles (length in m, band by latitude) ---
def edge_list(tri):
    e = np.vstack([tri[:, [0, 1]], tri[:, [1, 2]], tri[:, [2, 0]]])
    return np.unique(np.sort(e, axis=1), axis=0)

es = edge_list(Ts)
lat1 = 0.5 * (np.asarray(lat_s)[es[:, 0]] + np.asarray(lat_s)[es[:, 1]])
L1 = np.hypot(Ps[es[:, 0], 0] - Ps[es[:, 1], 0],
              Ps[es[:, 0], 1] - Ps[es[:, 1], 1])
eo = edge_list(T)
lat2 = 0.5 * (P[eo[:, 0], 1] + P[eo[:, 1], 1])
L2 = np.hypot((P[eo[:, 0], 0] - P[eo[:, 1], 0])
              * np.cos(np.deg2rad(lat2)),
              P[eo[:, 0], 1] - P[eo[:, 1], 1]) / DEG
bands = [(34.95, 35.10), (35.10, 35.20), (35.20, 35.30),
         (35.30, 35.40), (35.40, 35.55), (35.55, 35.70)]
print("band      sample p10/50/90      repro p10/50/90", flush=True)
for lo, hi in bands:
    m1 = (lat1 >= lo) & (lat1 < hi)
    m2 = (lat2 >= lo) & (lat2 < hi)
    a = np.percentile(L1[m1], [10, 50, 90]).round(0) if m1.sum() else []
    b = np.percentile(L2[m2], [10, 50, 90]).round(0) if m2.sum() else []
    print(f"{lo:.2f}-{hi:.2f}  {a}  {b}", flush=True)

# --- nodal h vs distance-to-own-coast: sample vs repro ---
from collections import defaultdict
import shapely
from pyproj import Transformer as _Tr

_tr_m = _Tr.from_crs("EPSG:4326", "EPSG:32654", always_xy=True)


def anatomy(pts_ll, tri, obc_nodes):
    xm, ym = _tr_m.transform(pts_ll[:, 0], pts_ll[:, 1])
    Pm = np.column_stack([xm, ym])
    ee = np.vstack([tri[:, [0, 1]], tri[:, [1, 2]], tri[:, [2, 0]]])
    ee = np.unique(np.sort(ee, axis=1), axis=0)
    Le = np.linalg.norm(Pm[ee[:, 0]] - Pm[ee[:, 1]], axis=1)
    hs = np.zeros(len(Pm)); hc = np.zeros(len(Pm))
    np.add.at(hs, ee[:, 0], Le); np.add.at(hc, ee[:, 0], 1)
    np.add.at(hs, ee[:, 1], Le); np.add.at(hc, ee[:, 1], 1)
    hn = hs / np.maximum(hc, 1)
    cnt2 = defaultdict(int)
    for a, b2, c2 in tri:
        for ed in ((a, b2), (b2, c2), (c2, a)):
            cnt2[tuple(sorted(ed))] += 1
    ob = set(int(v) for v in obc_nodes)
    ced = [ed for ed, k in cnt2.items() if k == 1
           and ed[0] not in ob and ed[1] not in ob]
    lines_c = shapely.MultiLineString(
        [[(Pm[a2, 0], Pm[a2, 1]), (Pm[b2, 0], Pm[b2, 1])]
         for a2, b2 in ced])
    dd = shapely.distance(shapely.points(Pm[:, 0], Pm[:, 1]), lines_c)
    cn = np.zeros(len(Pm), bool)
    cn[[v for ed in ced for v in ed]] = True
    return hn, dd, cn


hn1, dd1, cn1 = anatomy(np.column_stack([lon_s, lat_s]), Ts, obc_s)
hn2, dd2, cn2 = anatomy(P, T, obc_o)
print("\nnodal h by distance-to-coast   sample          repro", flush=True)
print(f"  coast nodes            "
      f"n={cn1.sum():5d} p50={np.percentile(hn1[cn1], 50):5.0f}   "
      f"n={cn2.sum():5d} p50={np.percentile(hn2[cn2], 50):5.0f}", flush=True)
for lo, hi in [(0, 300), (300, 700), (700, 1200), (1200, 2000),
               (2000, 3500), (3500, 6000), (6000, 99999)]:
    m1 = ~cn1 & (dd1 >= lo) & (dd1 < hi)
    m2 = ~cn2 & (dd2 >= lo) & (dd2 < hi)
    a = (f"n={m1.sum():5d} p50={np.percentile(hn1[m1], 50):5.0f}"
         if m1.sum() > 5 else "n=    -          ")
    b2 = (f"n={m2.sum():5d} p50={np.percentile(hn2[m2], 50):5.0f}"
          if m2.sum() > 5 else "n=    -          ")
    print(f"  d {lo:5d}-{hi:5d}          {a}   {b2}", flush=True)

# --- per-node Courant comparison at dt=18 (msh.CalcCFL port) ---
import sys
sys.path.insert(0, os.path.expanduser("~/Github/oceanmesh"))
from oceanmesh import calc_cfl

dep_s = np.loadtxt(os.path.expanduser(
    '~/Github/TB-FVCOM/input/goto2023/grid/TokyoBay_dep.dat'),
    skiprows=1)[:, 2]
dep_o = np.array([float(l.split()[3])
                  for l in lines[2:2 + int(lines[1].split()[1])]])
DT = 18.0
for tag, pts, tri, dep in [
        ("sample", np.column_stack([lon_s, lat_s]), Ts, dep_s),
        ("repro ", P, T, dep_o)]:
    cr = calc_cfl(pts, tri, dep, dt=DT)
    print(f"Cr(dt={DT:g}) {tag}: p50/p90/p99 = "
          f"{np.percentile(cr, [50, 90, 99]).round(3)} "
          f"max = {cr.max():.3f}  (n>0.5: {(cr > 0.5).sum()}, "
          f"n>0.6: {(cr > 0.6).sum()})", flush=True)
print("[cmp] done", flush=True)
