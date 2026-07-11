# Microscope: patches covered by OM2D raw but not ours — is our
# sizing coarser exactly there?
import os, sys, logging
import numpy as np
import h5py
logging.basicConfig(level=logging.WARNING, stream=sys.stdout)
os.environ.setdefault("MPLBACKEND", "Agg")
sys.path.insert(0, os.path.expanduser("~/Github/oceanmesh"))
import matplotlib.pyplot as plt
from pathlib import Path
from scipy.io import loadmat
from scipy.ndimage import label

OUT = Path("outputs/om2d_examples/jbay")
p2 = np.load(OUT / "p_raw_fresh_s2.npy")
t2 = np.load(OUT / "t_raw_fresh_s2.npy")
m = loadmat(os.path.expanduser(
    "~/Github/OceanMesh2D/Examples/Precleaned_grid.mat"))
pm = np.asarray(m["p"], float); tm = np.asarray(m["t"], int) - 1

x = np.linspace(-73.97, -73.75, 1100)
y = np.linspace(40.50, 40.68, 900)
X, Y = np.meshgrid(x, y)
from oceanmesh.cfl import _mesh_boundary_polygon
from oceanmesh import edges as om_edges
from oceanmesh.geometry import inpoly2
def cover(pp, tt):
    poly = _mesh_boundary_polygon(pp, tt)
    e = om_edges.get_poly_edges(poly)
    ins, _ = inpoly2(np.column_stack([X.ravel(), Y.ravel()]),
                     np.nan_to_num(poly), e)
    return ins.reshape(X.shape)
c_us = cover(p2, t2)
c_ml = cover(pm, tm)
miss = c_ml & ~c_us           # OM2D covers, we don't
both = c_ml & c_us
cell = (0.22/1099*111*np.cos(np.deg2rad(40.59))) * (0.18/899*111)
lab, n = label(miss)
sizes = np.bincount(lab.ravel())[1:]
big = np.argsort(sizes)[::-1][:10] + 1
print(f"[ms] missing total={miss.sum()*cell:.2f} km2 in {n} "
      f"patches; top10 areas km2: "
      f"{np.round(sizes[big-1]*cell, 3).tolist()}", flush=True)

# fh ratio inside missing patches vs common water
with h5py.File(str(OUT / "matlab_jbay_fh.mat"), "r") as f:
    xg = np.asarray(f["xg"]).ravel()
    yg = np.asarray(f["yg"]).ravel()
    vals = np.asarray(f["vals"])
if vals.shape == (len(yg), len(xg)):
    vals = vals.T
from scipy.interpolate import RegularGridInterpolator
Fml = RegularGridInterpolator((xg, yg), vals, bounds_error=False,
                              fill_value=np.nan)
# our fh on the lattice: rebuild grid quickly
from oceanmesh import DEM, Region, Shoreline
import oceanmesh as om
from oceanmesh.weirs import build_weir_geometry
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
dem = DEM(str(DSET / "PostSandyNCEI.nc"), bbox=reg)
feat = om.feature_sizing_function(shore, sdf, r=3,
                                  max_edge_length=1e3 * DEG)
grid, _ = om.finalize_sizing(
    [feat], dem=dem, shoreline=shore, hmin=15.0,
    max_edge_length=1e3, gradation=0.15,
    courant={"timestep": 2.0})
q = np.column_stack([X.ravel(), Y.ravel()])
ours = grid.eval(q).reshape(X.shape) / DEG
mlv = Fml(q).reshape(X.shape)
mlv = mlv / DEG if np.nanmedian(vals) < 1 else mlv
ratio = ours / mlv
for name, mask in (("missing-patches", miss), ("common-water", both)):
    r = ratio[mask & np.isfinite(ratio)]
    o = ours[mask & np.isfinite(ratio)]
    g = mlv[mask & np.isfinite(ratio)]
    print(f"[ms] {name}: fh ratio p10/50/90 = "
          f"{np.percentile(r,10):.3f}/{np.percentile(r,50):.3f}/"
          f"{np.percentile(r,90):.3f}  ours p50={np.percentile(o,50):.0f} m"
          f" ml p50={np.percentile(g,50):.0f} m", flush=True)
# vertex presence in top patches
from scipy.spatial import cKDTree
for k in big[:5]:
    msk = lab == k
    xs, ys = X[msk], Y[msk]
    x0, x1_, y0, y1_ = xs.min(), xs.max(), ys.min(), ys.max()
    nin_us = int(((p2[:, 0] > x0) & (p2[:, 0] < x1_)
                  & (p2[:, 1] > y0) & (p2[:, 1] < y1_)).sum())
    nin_ml = int(((pm[:, 0] > x0) & (pm[:, 0] < x1_)
                  & (pm[:, 1] > y0) & (pm[:, 1] < y1_)).sum())
    r = ratio[msk & np.isfinite(ratio)]
    print(f"[ms] patch#{k}: area={msk.sum()*cell:.3f} km2 "
          f"bbox=({x0:.4f},{y0:.4f}) our-verts={nin_us} "
          f"ml-verts={nin_ml} fh-ratio p50="
          f"{np.percentile(r,50):.3f}" if len(r) else
          f"[ms] patch#{k}: no ratio cells", flush=True)
print("[ms] done", flush=True)
