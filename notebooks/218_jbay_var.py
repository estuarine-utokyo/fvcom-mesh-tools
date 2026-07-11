# JBAY volume-variance attribution: seeds 1 & 2, raw + cleaned,
# with fort.14 export of raws for the MATLAB cross-clean.
import os, sys, logging, time
import numpy as np
logging.basicConfig(level=logging.INFO, stream=sys.stdout,
                    format="%(levelname)s %(name)s: %(message)s")
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
pfix_all, egfix_all, rings = [], [], []
off = 0
for wi in np.atleast_1d(w):
    cl = np.column_stack([np.atleast_1d(wi["X"]),
                          np.atleast_1d(wi["Y"])])
    pf, ef, ib = build_weir_geometry(cl, float(wi["width"]),
                                     float(wi["min_ele"]))
    pfix_all.append(pf); egfix_all.append(ef + off)
    off += len(pf)
    rings.append(np.vstack([pf, pf[0], [[np.nan, np.nan]]]))
pfix = np.vstack(pfix_all); egfix = np.vstack(egfix_all)

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
    courant={"timestep": 2.0})

def metrics(p, t, b=None):
    Re2 = 111.0 ** 2
    X, Y = p[:, 0], p[:, 1]
    x1, y1 = X[t[:, 0]], Y[t[:, 0]]
    x2, y2 = X[t[:, 1]], Y[t[:, 1]]
    x3, y3 = X[t[:, 2]], Y[t[:, 2]]
    pa = 0.5*np.abs((x2-x1)*(y3-y1)-(x3-x1)*(y2-y1))
    cosf = np.cos(np.deg2rad((y1+y2+y3)/3))
    A = float((pa*cosf).sum()*Re2)
    V = np.nan
    if b is not None:
        V = float((pa*cosf*(b[t].mean(axis=1)/1e3)).sum()*Re2)
    return A, V

for seed in (1, 2):
    pr, tr = om.generate_mesh(sdf, grid, max_iter=100, seed=seed,
                              pfix=pfix, egfix=egfix,
                              cleanup="none")
    Ar, _ = metrics(pr, tr)
    np.save(OUT / f"p_raw_s{seed}.npy", pr)
    np.save(OUT / f"t_raw_s{seed}.npy", tr)
    with open(OUT / f"raw_s{seed}.14", "w") as f:
        f.write(f"raw seed {seed}\n{len(tr)} {len(pr)}\n")
        for i, (x, y) in enumerate(pr, 1):
            f.write(f"{i} {x:.10f} {y:.10f} 0.0\n")
        for i, (a, b_, c) in enumerate(tr + 1, 1):
            f.write(f"{i} 3 {a} {b_} {c}\n")
        f.write("0\n0\n0\n0\n")
    pc, tc = om.generate_mesh(sdf, grid, max_iter=100, seed=seed,
                              pfix=pfix, egfix=egfix)
    bc = om.interp_bathymetry(pc, tc, dem, method="cell-averaging",
                              min_depth=1.0, nan_fill=True)
    Ac, Vc = metrics(pc, tc, bc)
    np.save(OUT / f"p_cln_s{seed}.npy", pc)
    np.save(OUT / f"t_cln_s{seed}.npy", tc)
    print(f"[var] seed{seed}: raw NP={len(pr):,} A={Ar:.2f} | "
          f"our-clean NP={len(pc):,} A={Ac:.2f} V={Vc:.4f}",
          flush=True)
print("[var] done", flush=True)
