# Compare our Example_6b stage-1 BANDED sizing field against the
# MATLAB edgefx dump (matlab_ex6b_fh.mat), split by elevation band.
import os, sys, logging
import numpy as np
import h5py
logging.basicConfig(level=logging.WARNING, stream=sys.stdout)
sys.path.insert(0, os.path.expanduser("~/Github/oceanmesh"))
from pathlib import Path
from scipy.io import loadmat
import oceanmesh as om
from oceanmesh import DEM, Region, Shoreline

OM2D = Path(os.path.expanduser("~/Github/OceanMesh2D"))
DS = OM2D / "datasets"
OUT = Path("outputs/om2d_examples/ex6bfp")
DEG = 1.0 / 111e3
mm = loadmat(str(DS/"ECGC_Thalwegs.mat"), squeeze_me=False)
channels = [np.asarray(a, dtype=float)
            for a in mm["pts2"].ravel() if np.size(a) >= 4]
bbox = (-95.40, -94.4, 29.14, 30.09)
reg = Region(bbox, 4326)
sh = Shoreline(
    str(DS/"US_Medium_Shoreline/us_medium_shoreline_polygon.shp"),
    reg.bbox, 60.0*DEG)
sh.detect_inpoly_flip(str(DS/"GSHHS_shp/l/GSHHS_l_L1.shp"))
sdf = om.signed_distance_function(sh)
dem = DEM(str(DS/"galveston_13_mhw_2007.nc"), bbox=reg)
f = om.feature_sizing_function(sh, sdf, r=3,
                               max_edge_length=1e3*DEG,
                               lattice_anchor=(dem.bbox[0],
                                               dem.bbox[2]))
c = om.channel_sizing_function(dem, channels, ch=0.1,
    min_edge_length_channel=100.0, angle_of_reslope=60.0,
    min_edge_length=60.0, max_edge_length=1e3, dx=60.0*DEG)
max_el = np.array([[1e3, -np.inf, 0.0], [500.0, 0.0, np.inf]])
grade = np.array([[0.25, -np.inf, 0.0], [0.05, 0.0, np.inf]])
g, _ = om.finalize_sizing([f, c], dem=dem, shoreline=sh,
                          hmin=60.0, max_edge_length=max_el,
                          gradation=grade)
print("[fh6b] our sizing built", flush=True)

with h5py.File(str(OUT/"matlab_ex6b_fh.mat"), "r") as fml:
    xg = np.asarray(fml["xg"]).ravel()
    yg = np.asarray(fml["yg"]).ravel()
    vals = np.asarray(fml["hh"])
if vals.shape == (len(yg), len(xg)):
    vals = vals.T
step = max(1, len(xg)//1200)
xs = xg[::step]; ys = yg[::step]
X, Y = np.meshgrid(xs, ys, indexing="ij")
q = np.column_stack([X.ravel(), Y.ravel()])
ours = np.asarray(g.eval(q)).reshape(X.shape) / DEG
mlv = vals[::step, ::step]
# elevation on the same probes for band attribution
z = np.asarray(dem.eval((np.clip(X, dem.bbox[0]+1e-9, dem.bbox[1]-1e-9),
                         np.clip(Y, dem.bbox[2]+1e-9, dem.bbox[3]-1e-9))))
r = ours / mlv
fin = np.isfinite(r)
print(f"[fh6b] ratio p10/50/90 = {np.percentile(r[fin],10):.3f}/"
      f"{np.percentile(r[fin],50):.3f}/{np.percentile(r[fin],90):.3f}",
      flush=True)
for name, m in (("water z<0", fin & (z < 0)), ("land z>=0", fin & (z >= 0))):
    if m.sum() > 50:
        print(f"[fh6b] {name}: ratio p10/50/90 = "
              f"{np.percentile(r[m],10):.3f}/{np.percentile(r[m],50):.3f}/"
              f"{np.percentile(r[m],90):.3f}  (ml p50 "
              f"{np.percentile(mlv[m],50):.0f} m, ours "
              f"{np.percentile(ours[m],50):.0f} m, n={int(m.sum())})",
              flush=True)
for lo, hi in ((0,100),(100,300),(300,600),(600,1001)):
    m2 = fin & (mlv >= lo) & (mlv < hi)
    if m2.sum() > 50:
        print(f"[fh6b]  ml in [{lo},{hi}): p50={np.percentile(r[m2],50):.3f} "
              f"n={int(m2.sum())}", flush=True)
import matplotlib
matplotlib.use("Agg")
import matplotlib.pyplot as plt
for tag, cond in (("fine", r < 0.8), ("coarse", r > 1.25)):
    tail = fin & cond
    print(f"[fh6b] {tag}-tail cells: {int(tail.sum())} "
          f"({100*tail.sum()/fin.sum():.1f}%)", flush=True)
    fig, ax = plt.subplots(figsize=(9, 8))
    ax.scatter(X[tail], Y[tail], s=0.3,
               c="crimson" if tag == "fine" else "navy")
    ax.set_xlim(bbox[0], bbox[1]); ax.set_ylim(bbox[2], bbox[3])
    ax.set_aspect(1/np.cos(np.deg2rad(29.6)))
    ax.set_title(f"Ex6b banded fh: ours {'<0.8x' if tag=='fine' else '>1.25x'} ML")
    fig.savefig(OUT/f"fh_{tag}_tail_map.png", dpi=140, bbox_inches="tight")
print("[fh6b] done", flush=True)
