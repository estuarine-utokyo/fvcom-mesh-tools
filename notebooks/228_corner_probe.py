# SW-corner microscope: both raw meshes + our fd sign + ring fd.
import os, sys, logging
import numpy as np
logging.basicConfig(level=logging.WARNING, stream=sys.stdout)
os.environ.setdefault("MPLBACKEND", "Agg")
sys.path.insert(0, os.path.expanduser("~/Github/oceanmesh"))
import matplotlib.pyplot as plt
from pathlib import Path
from scipy.io import loadmat
import oceanmesh as om
from oceanmesh import DEM, Region, Shoreline
from oceanmesh.weirs import build_weir_geometry

OUT = Path("outputs/om2d_examples/jbay")
DSET = Path(os.path.expanduser(
    "~/Github/OceanMesh2D/datasets/PostSandyNCEI"))
WST = Path(os.path.expanduser(
    "~/Github/OceanMesh2D/datasets/weirs_struct.mat"))
DEG = 1.0 / 111e3
w = loadmat(str(WST), squeeze_me=True)["weirs"]
rings = []
for wi in np.atleast_1d(w):
    cl = np.column_stack([np.atleast_1d(wi["X"]),
                          np.atleast_1d(wi["Y"])])
    pf, ef, ib = build_weir_geometry(cl, float(wi["width"]),
                                     float(wi["min_ele"]))
    rings.append(np.vstack([pf, pf[0], [[np.nan, np.nan]]]))
reg = Region((-73.97, -73.75, 40.5, 40.68), 4326)
shore = Shoreline(str(DSET / "PostSandyNCEI.shp"), reg.bbox,
                  15.0 * DEG)
shore.inner = np.vstack(
    [np.asarray(shore.inner).reshape(-1, 2),
     np.array([[np.nan, np.nan]])] + rings)
sdf = om.signed_distance_function(shore)

# ring fd distribution (lonlat, no projection roundtrip)
ring = np.asarray(getattr(sdf, "boubox_ring"))
d = sdf.eval(ring)
print(f"[cn] ring pts={len(ring)}; fd<0: {(d<0).sum()} "
      f"fd==0: {(d==0).sum()} fd>0: {(d>0).sum()}; "
      f"fd p50={np.median(d):.2e} max={d.max():.2e}", flush=True)
# same through the tmerc roundtrip used in generate_mesh
from pyproj import Transformer
lon0, lat0 = -73.86, 40.59
tr = Transformer.from_crs(
    "EPSG:4326", f"+proj=tmerc +lon_0={lon0} +lat_0={lat0} "
    "+ellps=WGS84 +units=m", always_xy=True)
xp, yp = tr.transform(ring[:, 0], ring[:, 1])
xb, yb = tr.transform(xp, yp, direction="INVERSE")
rt = np.column_stack([xb, yb])
d2 = sdf.eval(rt)
err = np.hypot(rt[:, 0]-ring[:, 0], rt[:, 1]-ring[:, 1])
print(f"[cn] after tmerc roundtrip: fd<0: {(d2<0).sum()} "
      f"fd>0: {(d2>0).sum()}; roundtrip |dxy| p50={np.median(err):.2e} "
      f"max={err.max():.2e}", flush=True)

# corner figure
p2 = np.load(OUT / "p_raw_fresh_s2.npy")
t2 = np.load(OUT / "t_raw_fresh_s2.npy")
m = loadmat(os.path.expanduser(
    "~/Github/OceanMesh2D/Examples/Precleaned_grid.mat"))
pm = np.asarray(m["p"], float); tm = np.asarray(m["t"], int) - 1
x0, x1, y0, y1 = -73.972, -73.90, 40.498, 40.545
fig, axes = plt.subplots(1, 2, figsize=(17, 8))
for ax, (pp, tt, ttl) in zip(
    axes, [(p2, t2, "our raw (seed2)"), (pm, tm, "OM2D raw")]):
    ax.triplot(pp[:, 0], pp[:, 1], tt, lw=0.4, color="steelblue")
    ax.plot(pp[:, 0], pp[:, 1], "k.", ms=1)
    ax.set_xlim(x0, x1); ax.set_ylim(y0, y1)
    ax.set_aspect(1/np.cos(np.deg2rad(40.52)))
    ax.set_title(ttl)
fig.suptitle("SW offshore frame corner")
fig.savefig(OUT / "sw_corner.png", dpi=170, bbox_inches="tight")
print("[cn] saved", flush=True)
