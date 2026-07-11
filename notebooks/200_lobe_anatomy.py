# Per-category crossing parity for NW-lobe points + windowed ring
# plot + matplotlib-Path even-odd referee.
import os, sys, logging
import numpy as np
logging.basicConfig(level=logging.WARNING, stream=sys.stdout)
os.environ.setdefault("MPLBACKEND", "Agg")
sys.path.insert(0, os.path.expanduser("~/Github/oceanmesh"))
import matplotlib.pyplot as plt
from pathlib import Path as P
from oceanmesh import DEM, Region, Shoreline
from oceanmesh import edges as om_edges
from oceanmesh.geometry import inpoly2

DSET = P(os.path.expanduser(
    "~/Github/OceanMesh2D/datasets/PostSandyNCEI"))
OUT = P("outputs/om2d_examples/ex2_ny")
DEG = 1.0 / 111e3
dem_probe = DEM(str(DSET / "PostSandyNCEI.nc"))
reg = Region(dem_probe.bbox, 4326)
shore = Shoreline(str(DSET / "PostSandyNCEI.shp"), reg.bbox,
                  30.0 * DEG)
pts = np.array([[-74.105, 40.780], [-74.100, 40.777],
                [-74.095, 40.770], [-74.110, 40.765]])

def cross_count(arr, q):
    a = np.asarray(arr, dtype=float)
    if a.size == 0:
        return 0
    e = om_edges.get_poly_edges(a)
    an = np.nan_to_num(a)
    p1 = an[e[:, 0]]; p2 = an[e[:, 1]]
    ca = p1[:, 1] > q[1]; cb = p2[:, 1] > q[1]
    cr = ca != cb
    with np.errstate(divide="ignore", invalid="ignore"):
        xi = p1[:, 0] + (q[1] - p1[:, 1]) * (p2[:, 0] - p1[:, 0]) \
            / (p2[:, 1] - p1[:, 1])
    return int((cr & (xi > q[0])).sum())

from matplotlib.path import Path as MPath
def mpl_evenodd(arr, q):
    a = np.asarray(arr, dtype=float)
    verts = []; codes = []
    idx = np.where(np.isnan(a[:, 0]))[0]
    start = 0
    for stop in list(idx) + [len(a)]:
        seg = a[start:stop]; start = stop + 1
        if len(seg) < 3:
            continue
        verts.extend(seg.tolist()); codes.extend(
            [MPath.MOVETO] + [MPath.LINETO] * (len(seg) - 1))
    return MPath(verts, codes).contains_point(q)

bb = np.asarray(shore.boubox)
ml = np.asarray(shore.mainland)
inn = np.asarray(shore.inner)
full = np.vstack((bb, ml, inn))
for q in pts:
    cb_ = cross_count(bb, q)
    cm_ = cross_count(ml, q)
    ci_ = cross_count(inn, q)
    tot = cb_ + cm_ + ci_
    ref = mpl_evenodd(full, q)
    print(f"[an] ({q[0]:.3f},{q[1]:.3f}): boubox={cb_} "
          f"mainland={cm_} inner={ci_} total={tot} "
          f"({'ODD=water' if tot % 2 else 'EVEN=land'}) "
          f"mpl_referee={'water' if ref else 'land'}", flush=True)

x0, x1, y0, y1 = -74.13, -74.07, 40.755, 40.795
fig, ax = plt.subplots(figsize=(13, 9))
for arr, c, lw in ((ml, "g", 1.0), (inn, "m", 0.7)):
    if arr.size:
        ax.plot(arr[:, 0], arr[:, 1], c, lw=lw)
ax.plot(pts[:, 0], pts[:, 1], "k*", ms=12)
ax.set_xlim(x0, x1); ax.set_ylim(y0, y1)
ax.set_aspect(1 / np.cos(np.deg2rad(40.775)))
ax.set_title("NW lobe rings: mainland g, inner m, probes *")
fig.savefig(OUT / "nw_lobe_rings.png", dpi=170, bbox_inches="tight")
print("[an] saved", flush=True)
