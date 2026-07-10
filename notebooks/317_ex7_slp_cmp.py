# Isolate the Ex7 slope (slp) sizing: build OUR slope-only field
# (slope + hmin/max_el + gradation, no wl/fs/CFL — mirroring the
# MATLAB dump matlab_ex7_slponly.mat) and compare per depth band
# and on ridge windows.
import os, sys, logging, time
import numpy as np
import h5py
logging.basicConfig(level=logging.INFO, stream=sys.stdout,
                    format="%(levelname)s %(name)s: %(message)s")
sys.path.insert(0, os.path.expanduser("~/Github/oceanmesh"))
from pathlib import Path
import oceanmesh as om
from oceanmesh import DEM, Region

OM2D = Path(os.path.expanduser("~/Github/OceanMesh2D"))
DS = OM2D / "datasets"
SCR = Path(os.path.expanduser(
    "~/Github/fvcom-mesh-tools/outputs/om2d_examples/ex7glob"))
DEG = 1.0 / 111e3
t0 = time.time()

bbox = (-180.0, 180.0, -89.0, 90.0)
reg = Region(bbox, 4326)
dem = DEM(str(DS/"SRTM15+.nc"), bbox=reg)
# Example_7: no 'fl' -> edgefx.m defval=0 -> slope filter OFF
s = om.bathymetric_gradient_sizing_function(
    dem, slope_parameter=10,
    min_edge_length=4e3*DEG, max_edge_length=20e3*DEG,
    type_of_filter="none", grid_dx=4e3*DEG)
g, _ = om.finalize_sizing([s], dem=dem, hmin=4e3,
                          max_edge_length=20e3, gradation=0.25)
np.save(SCR/"g_slponly.npy", np.asarray(g.values, dtype=float))
print(f"[slp7] our slope-only field built +{time.time()-t0:.0f}s",
      flush=True)

with h5py.File(SCR/"matlab_ex7_slponly.mat", "r") as f:
    xg = np.asarray(f["xg"]).ravel()
    yg = np.asarray(f["yg"]).ravel()
    hh = np.asarray(f["hh"])
if hh.shape == (len(yg), len(xg)):
    hh = hh.T
step = max(1, len(xg)//1400)
xs = xg[::step]; ys = yg[::step]
X, Y = np.meshgrid(xs, ys, indexing="ij")
q = np.column_stack([X.ravel(), Y.ravel()])
ours = np.asarray(g.eval(q)).reshape(X.shape) / DEG
mlv = hh[::step, ::step]
z = np.asarray(dem.eval((np.clip(X, -179.999, 179.999),
                         np.clip(Y, -88.999, 89.999))))
r = ours / mlv
fin = np.isfinite(r) & (mlv > 0) & (z < 0)
print(f"[slp7] OCEAN ratio p10/50/90 = "
      f"{np.percentile(r[fin], [10, 50, 90]).round(3)}", flush=True)
for lo, hi, tag in [(0, 5e3, "<5k"), (5e3, 10e3, "5-10k"),
                    (10e3, 19e3, "10-19k"), (19e3, 1e9, ">=19k")]:
    m2 = fin & (mlv >= lo) & (mlv < hi)
    if m2.sum() > 100:
        print(f"[slp7]  ml {tag:7s}: p50={np.percentile(r[m2],50):.3f} "
              f"n={int(m2.sum()):,}", flush=True)
for lo, hi, tag in [(-11000, -4000, "deep"), (-4000, -2000, "mid"),
                    (-2000, -200, "shelf-slope"), (-200, 0, "shelf")]:
    m2 = fin & (z >= lo) & (z < hi)
    if m2.sum() > 100:
        print(f"[slp7]  depth {tag:11s}: p50={np.percentile(r[m2],50):.3f} "
              f"ml_p50={np.percentile(mlv[m2],50):,.0f} "
              f"ours_p50={np.percentile(ours[m2],50):,.0f}", flush=True)
# ridge windows
import matplotlib
matplotlib.use("Agg")
import matplotlib.pyplot as plt
wins = [("MAR 20-35S", (-20, -5, -35, -20)),
        ("EPR 5-20S", (-115, -100, -20, -5)),
        ("SWIR", (55, 70, -32, -22))]
fig, axs = plt.subplots(len(wins), 2, figsize=(14, 5.5*len(wins)))
for k, (nm, w) in enumerate(wins):
    m = (X >= w[0]) & (X <= w[1]) & (Y >= w[2]) & (Y <= w[3])
    for j, (fld, ti) in enumerate([(ours, "ours"), (mlv, "OM2D")]):
        ax = axs[k, j]
        im = ax.pcolormesh(xs[(xs >= w[0]) & (xs <= w[1])],
                           ys[(ys >= w[2]) & (ys <= w[3])],
                           fld[np.ix_((xs >= w[0]) & (xs <= w[1]),
                                      (ys >= w[2]) & (ys <= w[3]))].T,
                           vmin=4e3, vmax=20e3, cmap="turbo_r")
        ax.set_title(f"{nm} — {ti}", fontsize=9)
        fig.colorbar(im, ax=ax, shrink=0.8)
fig.savefig(SCR/"slp7_ridge_cmp.png", dpi=140, bbox_inches="tight")
print("[slp7] saved slp7_ridge_cmp.png", flush=True)
print(f"[slp7] done +{time.time()-t0:.0f}s", flush=True)
