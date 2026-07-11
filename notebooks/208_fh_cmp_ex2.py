# Sizing-field comparison at the Meadowlands: is their fh finer?
import os, sys, logging
import numpy as np
import h5py
logging.basicConfig(level=logging.WARNING, stream=sys.stdout)
sys.path.insert(0, os.path.expanduser("~/Github/oceanmesh"))
from pathlib import Path
import oceanmesh as om
from oceanmesh import DEM, Region, Shoreline

DSET = Path(os.path.expanduser(
    "~/Github/OceanMesh2D/datasets/PostSandyNCEI"))
OUT = Path("outputs/om2d_examples/ex2_ny")
DEG = 1.0 / 111e3
dem_probe = DEM(str(DSET / "PostSandyNCEI.nc"))
reg = Region(dem_probe.bbox, 4326)
shore = Shoreline(str(DSET / "PostSandyNCEI.shp"), reg.bbox,
                  30.0 * DEG)
sdf = om.signed_distance_function(shore)
dem = DEM(str(DSET / "PostSandyNCEI.nc"), bbox=reg)
feat = om.feature_sizing_function(shore, sdf, r=3,
                                  max_edge_length=1e3 * DEG)
grid, _ = om.finalize_sizing(
    [feat], dem=dem, shoreline=shore, hmin=30.0,
    max_edge_length=1e3, max_edge_length_nearshore=240.0,
    gradation=0.20, courant={"timestep": 2.0},
)
with h5py.File(str(OUT / "matlab_ex2_fh.mat"), "r") as f:
    xg = np.asarray(f["xg"]).ravel()
    yg = np.asarray(f["yg"]).ravel()
    vals = np.asarray(f["vals"])
if vals.shape == (len(yg), len(xg)):
    vals = vals.T
print(f"[fh2] ML lattice {vals.shape} "
      f"x[{xg.min():.4f},{xg.max():.4f}]", flush=True)
X, Y = np.meshgrid(xg, yg, indexing="ij")
q = np.column_stack([X.ravel(), Y.ravel()])
ours = grid.eval(q).reshape(vals.shape) / DEG  # meters
ml = vals * 1.0  # already meters? OM2D fh in meters for m_map...
# OM2D edgefx Values are in DEGREES (WGS84 grid) — convert
ml_m = vals / DEG if np.nanmedian(vals) < 1 else vals
ratio = ours / ml_m
WINS = {"global": (xg.min(), xg.max(), yg.min(), yg.max()),
        "lobe": (-74.118, -74.082, 40.7625, 40.790),
        "neck": (-74.104, -74.080, 40.752, 40.7725)}
for name, w in WINS.items():
    m = ((X > w[0]) & (X < w[1]) & (Y > w[2]) & (Y < w[3])
         & np.isfinite(ratio))
    r = ratio[m]
    o = ours[m]; g = ml_m[m]
    print(f"[fh2] {name}: ratio p10/50/90 = "
          f"{np.percentile(r,10):.3f}/{np.percentile(r,50):.3f}/"
          f"{np.percentile(r,90):.3f}  ours p50={np.percentile(o,50):.0f} m "
          f"ml p50={np.percentile(g,50):.0f} m", flush=True)
print("[fh2] done", flush=True)
