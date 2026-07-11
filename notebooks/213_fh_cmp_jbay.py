# JBAY sizing-field ratio: ours vs MATLAB edgefx dump.
import os, sys, logging
import numpy as np
import h5py
logging.basicConfig(level=logging.WARNING, stream=sys.stdout)
sys.path.insert(0, os.path.expanduser("~/Github/oceanmesh"))
from pathlib import Path
from scipy.io import loadmat
import oceanmesh as om
from oceanmesh import DEM, Region, Shoreline
from oceanmesh.weirs import build_weir_geometry

DSET = Path(os.path.expanduser(
    "~/Github/OceanMesh2D/datasets/PostSandyNCEI"))
WST = Path(os.path.expanduser(
    "~/Github/OceanMesh2D/datasets/weirs_struct.mat"))
OUT = Path("outputs/om2d_examples/jbay")
DEG = 1.0 / 111e3
bbox = (-73.97, -73.75, 40.5, 40.68)
w = loadmat(str(WST), squeeze_me=True)["weirs"]
rings = []
for wi in np.atleast_1d(w):
    cl = np.column_stack([np.atleast_1d(wi["X"]),
                          np.atleast_1d(wi["Y"])])
    pf, ef, ib = build_weir_geometry(cl, float(wi["width"]),
                                     float(wi["min_ele"]))
    rings.append(np.vstack([pf, pf[0], [[np.nan, np.nan]]]))
reg = Region(bbox, 4326)
shore = Shoreline(str(DSET / "PostSandyNCEI.shp"), reg.bbox,
                  15.0 * DEG)
shore.inner = np.vstack(
    [np.asarray(shore.inner).reshape(-1, 2),
     np.array([[np.nan, np.nan]])] + rings)
sdf = om.signed_distance_function(shore)
dem = DEM(str(DSET / "PostSandyNCEI.nc"), bbox=reg)
feat = om.feature_sizing_function(shore, sdf, r=3,
                                  max_edge_length=1e3 * DEG)
grid, _ = om.finalize_sizing(
    [feat], dem=dem, shoreline=shore, hmin=15.0,
    max_edge_length=1e3, gradation=0.15,
    courant={"timestep": 2.0},
)
try:
    with h5py.File(str(OUT.parent / "jbay" / "matlab_jbay_fh.mat"), "r") as f:
        xg = np.asarray(f["xg"]).ravel()
        yg = np.asarray(f["yg"]).ravel()
        vals = np.asarray(f["vals"])
except Exception:
    with h5py.File(str(Path("outputs/om2d_examples/ex2_ny/matlab_jbay_fh.mat")), "r") as f:
        xg = np.asarray(f["xg"]).ravel()
        yg = np.asarray(f["yg"]).ravel()
        vals = np.asarray(f["vals"])
if vals.shape == (len(yg), len(xg)):
    vals = vals.T
X, Y = np.meshgrid(xg, yg, indexing="ij")
ours = grid.eval(np.column_stack([X.ravel(), Y.ravel()])
                 ).reshape(vals.shape) / DEG
ml = vals / DEG if np.nanmedian(vals) < 1 else vals
ratio = ours / ml
m = np.isfinite(ratio)
print(f"[fhjb] ML lattice {vals.shape}", flush=True)
print(f"[fhjb] global ratio p10/50/90 = "
      f"{np.percentile(ratio[m],10):.3f}/"
      f"{np.percentile(ratio[m],50):.3f}/"
      f"{np.percentile(ratio[m],90):.3f}", flush=True)
print(f"[fhjb] ours p10/50/90 = "
      f"{np.percentile(ours[m],10):.0f}/"
      f"{np.percentile(ours[m],50):.0f}/"
      f"{np.percentile(ours[m],90):.0f} m", flush=True)
print(f"[fhjb] ml   p10/50/90 = "
      f"{np.percentile(ml[m],10):.0f}/"
      f"{np.percentile(ml[m],50):.0f}/"
      f"{np.percentile(ml[m],90):.0f} m", flush=True)
import matplotlib
matplotlib.use("Agg")
import matplotlib.pyplot as plt
fig, ax = plt.subplots(figsize=(11, 9))
pc = ax.pcolormesh(X, Y, np.clip(ratio, 0, 2), cmap="RdBu_r",
                   vmin=0, vmax=2, shading="auto")
fig.colorbar(pc, ax=ax, label="ours / MATLAB sizing ratio")
ax.set_aspect(1 / np.cos(np.deg2rad(40.59)))
fig.savefig(OUT / "fh_ratio.png", dpi=150, bbox_inches="tight")
print("[fhjb] saved", flush=True)
