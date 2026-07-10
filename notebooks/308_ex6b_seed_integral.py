# Decisive split for the Ex6b stage-1 NP gap (-15%): expected seed
# count integral N ~ (2/sqrt(3)) * sum(dA / fh^2) over the WET
# domain, for OUR banded field vs the MATLAB dump on the SAME
# lattice/mask. If integrals agree, the gap is in the generator; if
# they differ ~15%, the fields differ where it matters (small fh).
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
print("[seed6b] our sizing built", flush=True)

with h5py.File(str(OUT/"matlab_ex6b_fh.mat"), "r") as fml:
    xg = np.asarray(fml["xg"]).ravel()
    yg = np.asarray(fml["yg"]).ravel()
    vals = np.asarray(fml["hh"])
if vals.shape == (len(yg), len(xg)):
    vals = vals.T
X, Y = np.meshgrid(xg, yg, indexing="ij")
q = np.column_stack([X.ravel(), Y.ravel()])
print(f"[seed6b] lattice {vals.shape}, evaluating sdf...", flush=True)
d = np.asarray(sdf.eval(q))
wet = (d < 0).reshape(X.shape)
ours = np.asarray(g.eval(q)).reshape(X.shape) / DEG
mlv = vals
cell_dx = float(np.mean(np.diff(xg))) * 111e3 * np.cos(np.deg2rad(29.6))
cell_dy = float(np.mean(np.diff(yg))) * 111e3
dA = cell_dx * cell_dy
for name, h in (("ours", ours), ("matlab", mlv)):
    m = wet & np.isfinite(h) & (h > 0)
    N = (2/np.sqrt(3)) * dA * np.sum(1.0/h[m]**2)
    print(f"[seed6b] {name}: expected N = {N:,.0f} over {int(m.sum()):,} wet cells",
          flush=True)
# where does the difference live? contribution-weighted ratio bands
m = wet & np.isfinite(ours) & np.isfinite(mlv) & (ours > 0) & (mlv > 0)
contrib_ml = 1.0/mlv[m]**2
r = ours[m]/mlv[m]
for lo, hi in ((0, 0.8), (0.8, 0.95), (0.95, 1.05), (1.05, 1.25), (1.25, 99)):
    b = (r >= lo) & (r < hi)
    print(f"[seed6b] ratio [{lo},{hi}): cells {int(b.sum()):,}, "
          f"ML-node-share {100*contrib_ml[b].sum()/contrib_ml.sum():.1f}%",
          flush=True)
print("[seed6b] done", flush=True)
