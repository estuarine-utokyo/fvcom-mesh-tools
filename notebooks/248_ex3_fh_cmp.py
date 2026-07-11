# Ex3 fh ratio per nest, banded by size class.
import os, sys, logging
import numpy as np
import h5py
logging.basicConfig(level=logging.WARNING, stream=sys.stdout)
sys.path.insert(0, os.path.expanduser("~/Github/oceanmesh"))
from pathlib import Path
import oceanmesh as om
from oceanmesh import DEM, Region, Shoreline

OM2D = Path(os.path.expanduser("~/Github/OceanMesh2D"))
DSET = OM2D / "datasets/PostSandyNCEI"
OUT = Path("outputs/om2d_examples/ex3full")
DEG = 1.0 / 111e3

def build(nest):
    if nest == 1:
        bbox = np.array([[-71.6, 42.7], [-64, 30], [-80, 24],
                         [-85, 38], [-71.6, 42.7]])
        reg = Region((-85.0, -64.0, 24.0, 42.7), 4326)
        sh = Shoreline(str(OM2D/"datasets/GSHHS_shp/f/GSHHS_f_L1.shp"),
                       bbox, 1e3*DEG)
        dem = DEM(str(OM2D/"datasets/SRTM15+.nc"), bbox=reg)
        sdf = om.signed_distance_function(sh)
        f = om.feature_sizing_function(sh, sdf, r=3,
                                       max_edge_length=50e3*DEG)
        w = om.wavelength_sizing_function(dem, wl=30,
                                          min_edgelength=1e3*DEG,
                                          max_edge_length=50e3*DEG)
        g, _ = om.finalize_sizing([f, w], dem=dem, shoreline=sh,
                                  hmin=1e3, max_edge_length=50e3,
                                  gradation=0.35,
                                  courant={"timestep": 0.0})
    else:
        bbox = np.array([[-74.25, 40.5], [-73.75, 40.55],
                         [-73.75, 41], [-74, 41], [-74.25, 40.5]])
        reg = Region((-74.25, -73.75, 40.5, 41.0), 4326)
        sh = Shoreline(str(DSET/"PostSandyNCEI.shp"), bbox, 30.0*DEG)
        dem = DEM(str(DSET/"PostSandyNCEI.nc"), bbox=reg)
        sdf = om.signed_distance_function(sh)
        f = om.feature_sizing_function(sh, sdf, r=3,
                                       max_edge_length=1e3*DEG)
        w = om.wavelength_sizing_function(dem, wl=30,
                                          min_edgelength=30*DEG,
                                          max_edge_length=1e3*DEG)
        g, _ = om.finalize_sizing([f, w], dem=dem, shoreline=sh,
                                  hmin=30.0, max_edge_length=1e3,
                                  max_edge_length_nearshore=240.0,
                                  gradation=0.35,
                                  courant={"timestep": 0.0})
    return g

for nest in (1, 2):
    g = build(nest)
    with h5py.File(str(OUT / f"ml_fh{nest}.mat"), "r") as f:
        xg = np.asarray(f["xg"]).ravel()
        yg = np.asarray(f["yg"]).ravel()
        vals = np.asarray(f["vals"])
    if vals.shape == (len(yg), len(xg)):
        vals = vals.T
    step = max(1, len(xg)//1500)
    xs = xg[::step]; ys = yg[::step]
    X, Y = np.meshgrid(xs, ys, indexing="ij")
    q = np.column_stack([X.ravel(), Y.ravel()])
    ours = np.asarray(g.eval(q)).reshape(X.shape) / DEG
    mlv = vals[::step, ::step]
    mlv = mlv / DEG if np.nanmedian(vals) < 1 else mlv
    r = ours / mlv
    fin = np.isfinite(r)
    print(f"[fh3] nest{nest}: ratio p10/50/90 = "
          f"{np.percentile(r[fin],10):.3f}/"
          f"{np.percentile(r[fin],50):.3f}/"
          f"{np.percentile(r[fin],90):.3f}", flush=True)
    for lo, hi in ((0, 100), (100, 500), (500, 2000),
                   (2000, 10000), (10000, 60000)):
        m = fin & (mlv >= lo) & (mlv < hi)
        if m.sum() > 50:
            print(f"[fh3]   ml∈[{lo},{hi}) m: ratio p50="
                  f"{np.percentile(r[m],50):.3f} n={int(m.sum())}",
                  flush=True)
print("[fh3] done", flush=True)
