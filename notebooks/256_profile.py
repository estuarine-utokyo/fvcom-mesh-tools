# Profile a representative pipeline: JBAY sizing + 15 mesh iters.
import os, sys, cProfile, pstats, io, time
import numpy as np
sys.path.insert(0, os.path.expanduser("~/Github/oceanmesh"))
import logging
logging.basicConfig(level=logging.WARNING)
from pathlib import Path
from scipy.io import loadmat
import oceanmesh as om
from oceanmesh import DEM, Region, Shoreline
from oceanmesh.weirs import build_weir_geometry

DSET = Path(os.path.expanduser(
    "~/Github/OceanMesh2D/datasets/PostSandyNCEI"))
WST = Path(os.path.expanduser(
    "~/Github/OceanMesh2D/datasets/weirs_struct.mat"))
DEG = 1.0 / 111e3
w = loadmat(str(WST), squeeze_me=True)["weirs"]
rings, pfix_all, egfix_all, wcfg = [], [], [], []
off = 0
for wi in np.atleast_1d(w):
    cl = np.column_stack([np.atleast_1d(wi["X"]),
                          np.atleast_1d(wi["Y"])])
    pf, ef, ib = build_weir_geometry(cl, float(wi["width"]),
                                     float(wi["min_ele"]))
    rings.append(np.vstack([pf, pf[0], [[np.nan, np.nan]]]))
    pfix_all.append(pf); egfix_all.append(ef + off); off += len(pf)
    wcfg.append({"crestline": cl, "min_ele_m": float(wi["min_ele"]),
                 "width_m": float(wi["width"])})
reg = Region((-73.97, -73.75, 40.5, 40.68), 4326)

pr = cProfile.Profile(); pr.enable()
shore = Shoreline(str(DSET/"PostSandyNCEI.shp"), reg.bbox, 15.0*DEG)
shore.inner = np.vstack([np.asarray(shore.inner).reshape(-1, 2),
                         np.array([[np.nan, np.nan]])] + rings)
sdf = om.signed_distance_function(shore)
dem = DEM(str(DSET/"PostSandyNCEI.nc"), bbox=reg)
feat = om.feature_sizing_function(shore, sdf, r=3,
                                  max_edge_length=1e3*DEG)
grid, _ = om.finalize_sizing([feat], dem=dem, shoreline=shore,
                             hmin=15.0, max_edge_length=1e3,
                             gradation=0.15,
                             courant={"timestep": 2.0},
                             weirs=wcfg)
pr.disable()
s = io.StringIO(); st = pstats.Stats(pr, stream=s)
st.sort_stats("cumulative").print_stats(18)
print("=== SIZING PROFILE ===")
print("\n".join(s.getvalue().splitlines()[4:28]), flush=True)

pfix = np.vstack(pfix_all); egfix = np.vstack(egfix_all)
pr2 = cProfile.Profile(); pr2.enable()
p, t = om.generate_mesh(sdf, grid, max_iter=15, seed=0,
                        pfix=pfix, egfix=egfix, cleanup="none")
pr2.disable()
s2 = io.StringIO(); st2 = pstats.Stats(pr2, stream=s2)
st2.sort_stats("cumulative").print_stats(20)
print("=== MESH-15ITER PROFILE ===")
print("\n".join(s2.getvalue().splitlines()[4:32]), flush=True)
