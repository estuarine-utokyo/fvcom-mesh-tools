# Outer sizing field after nest smoothing: ours vs OM2D.
import os, sys, logging
import numpy as np
logging.basicConfig(level=logging.WARNING, stream=sys.stdout)
os.environ.setdefault("MPLBACKEND", "Agg")
sys.path.insert(0, os.path.expanduser("~/Github/oceanmesh"))
import matplotlib.pyplot as plt
from pathlib import Path
from scipy.io import loadmat
from oceanmesh.grid import Grid
from oceanmesh.edgefx import multiscale_sizing_function

DEG = 1.0/111e3
OUT = Path("outputs/om2d_examples/ex10")
def uniform_grid(bbox, h):
    g = Grid(bbox=bbox, dx=h, extrapolate=True, values=float(h),
             crs=4326)
    g.hmin = h
    g.build_interpolant()
    return g
g2 = uniform_grid((-1.0, 2.0, -1.0, 2.0), 1e4*DEG)
g1 = uniform_grid((0.0, 1.0, 0.0, 1.0), 1e3*DEG)
_, grids = multiscale_sizing_function([g2, g1], gradation=0.2)
ours = np.asarray(grids[0].values, float) / DEG
xv, yv = grids[0].create_vectors()

m = loadmat(str(OUT / "ml_outer_ef.mat"))
mx = m["xg"].ravel(); my = m["yg"].ravel()
mv = np.asarray(m["vals"], float)
if mv.shape == (len(my), len(mx)):
    mv = mv.T
mlm = mv * 111e3 if np.nanmedian(mv) < 1 else mv
print(f"[ef] ours lattice {ours.shape} ml {mlm.shape}", flush=True)

# cross sections through the box centre: W-E at y=0.5, S-N at x=0.5
iy = np.argmin(np.abs(yv - 0.5)); jx = np.argmin(np.abs(xv - 0.5))
miy = np.argmin(np.abs(my - 0.5)); mjx = np.argmin(np.abs(mx - 0.5))
print("[ef] W-E section y=0.5 (x, ours_m, ml_m):", flush=True)
for k in range(len(xv)):
    if -0.35 < xv[k] < 0.35 or 0.75 < xv[k] < 1.35:
        kk = np.argmin(np.abs(mx - xv[k]))
        print(f"  x={xv[k]:+.3f} ours={ours[k, iy]:7.0f} "
              f"ml={mlm[kk, miy]:7.0f}", flush=True)
fig, axes = plt.subplots(1, 2, figsize=(16, 7.5))
for ax, (Z, xs, ys, ttl) in zip(
    axes, [(ours, xv, yv, "ours (post smooth_outer port)"),
           (mlm, mx, my, "OM2D ef{1} (post smooth_outer)")]):
    X, Y = np.meshgrid(xs, ys, indexing="ij")
    pc = ax.pcolormesh(X, Y, Z, shading="auto", vmin=1e3, vmax=1e4)
    fig.colorbar(pc, ax=ax, shrink=0.8, label="h [m]")
    ax.plot([0,1,1,0,0],[0,0,1,1,0],"m-",lw=1.5)
    ax.set_aspect(1); ax.set_title(ttl)
fig.savefig(OUT / "outer_ef_cmp.png", dpi=150, bbox_inches="tight")
print("[ef] saved", flush=True)
