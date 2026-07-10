# Ex7: compare our global sizing (saved by 309: g_before_stereo.npy
# on the feature h0 lattice) against the MATLAB edgefx dump
# (matlab_ex7_fh.mat), and compute the ML field's expected-N with
# the same wet mask machinery as 309.
import os, sys, logging, time
import numpy as np
import h5py
logging.basicConfig(level=logging.WARNING, stream=sys.stdout)
sys.path.insert(0, os.path.expanduser("~/Github/oceanmesh"))
from pathlib import Path
import oceanmesh as om
from oceanmesh import Region, Shoreline, Grid

SCR = Path(os.path.expanduser(
    "~/Github/fvcom-mesh-tools/outputs/om2d_examples/ex7glob"))
DEG = 1.0 / 111e3
t0 = time.time()

with h5py.File(SCR/"matlab_ex7_fh.mat", "r") as f:
    xg = np.asarray(f["xg"]).ravel()
    yg = np.asarray(f["yg"]).ravel()
    hh = np.asarray(f["hh"])
if hh.shape == (len(yg), len(xg)):
    hh = hh.T
print(f"[fh7] ML lattice {hh.shape} x[{xg.min():.2f},{xg.max():.2f}] "
      f"y[{yg.min():.2f},{yg.max():.2f}] h p10/50/90 = "
      f"{np.percentile(hh,[10,50,90]).round(0)}", flush=True)

# our field: rebuild the Grid container for g_before/g_after
bbox = (-180.0, 180.0, -89.0, 90.0)
reg = Region(bbox, 4326)
sh = Shoreline(str(SCR/"gshhs_l1l6.shp"), reg.bbox, 4e3*DEG)
sdf = om.signed_distance_function(sh)
print(f"[fh7] shoreline +{time.time()-t0:.0f}s", flush=True)

# ML-lattice probes (subsample to ~1200 columns for the field maps)
step = max(1, len(xg)//1200)
xs = xg[::step]; ys = yg[::step]
X, Y = np.meshgrid(xs, ys, indexing="ij")
q = np.column_stack([X.ravel(), Y.ravel()])
mlv = hh[::step, ::step]

for tag in ("g_before_stereo", "g_after_stereo"):
    vals = np.load(SCR/f"{tag}.npy")
    # reconstruct the lattice this field lives on: feature h0 lattice
    nx, ny = vals.shape
    gx = np.linspace(bbox[0], bbox[1], nx)
    g = Grid(bbox=bbox, dx=float(gx[1]-gx[0]),
             dy=float((bbox[3]-bbox[2])/(ny-1)),
             extrapolate=True, values=vals, crs=4326)
    g.build_interpolant()
    ours = np.asarray(g.eval(q)).reshape(X.shape) / DEG  # deg -> m
    r = ours / mlv
    fin = np.isfinite(r) & (mlv > 0)
    print(f"[fh7] {tag}: ratio p10/50/90 = "
          f"{np.percentile(r[fin],[10,50,90]).round(3)}", flush=True)
    for lo, hi in ((0, 5e3), (5e3, 10e3), (10e3, 19e3), (19e3, 1e9)):
        m2 = fin & (mlv >= lo) & (mlv < hi)
        if m2.sum() > 100:
            print(f"[fh7]   ml in [{lo/1e3:.0f}k,{hi/1e3:.0f}k): "
                  f"p50={np.percentile(r[m2],50):.3f} n={int(m2.sum()):,}",
                  flush=True)

# expected-N of the ML field over the wet mask (ML lattice, full res)
Xf, Yf = np.meshgrid(xg, yg, indexing="ij")
qf = np.column_stack([Xf.ravel(), Yf.ravel()])
print(f"[fh7] evaluating sdf on {len(qf):,} pts...", flush=True)
d = np.asarray(sdf.eval(qf))
wet = (d < 0).reshape(Xf.shape)
dxm = float(np.mean(np.diff(xg))) * 111e3
dym = float(np.mean(np.diff(yg))) * 111e3
dA = dxm * dym * np.cos(np.deg2rad(np.clip(Yf, -89.9, 89.9)))
m = wet & np.isfinite(hh) & (hh > 0)
N = (2/np.sqrt(3)) * np.sum(dA[m] / hh[m]**2)
print(f"[fh7] ML field expected N = {N:,.0f} over {int(m.sum()):,} wet "
      f"cells (ours-before was 2,222,116; golden final NP 1,372,623)",
      flush=True)
print(f"[fh7] done +{time.time()-t0:.0f}s", flush=True)
