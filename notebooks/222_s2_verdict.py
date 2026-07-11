# Verdict on the bad realization: same raw, both cleans, A & V.
import os, sys, logging
import numpy as np
logging.basicConfig(level=logging.WARNING, stream=sys.stdout)
sys.path.insert(0, os.path.expanduser("~/Github/oceanmesh"))
from pathlib import Path
from scipy.io import loadmat
import oceanmesh as om
from oceanmesh import DEM, Region

OUT = Path("outputs/om2d_examples/jbay")
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

reg = Region((-73.97, -73.75, 40.5, 40.68), 4326)
dem = DEM(os.path.expanduser(
    "~/Github/OceanMesh2D/datasets/PostSandyNCEI/PostSandyNCEI.nc"),
    bbox=reg)

pr = np.load(OUT / "p_raw_fresh_s2.npy")
tr = np.load(OUT / "t_raw_fresh_s2.npy")
Ar, _ = metrics(pr, tr)
print(f"[s2] raw: NP={len(pr):,} A={Ar:.2f}", flush=True)

p = np.load(OUT / "p.npy"); t = np.load(OUT / "t.npy")
b = np.load(OUT / "b.npy")
Ao, Vo = metrics(p, t, b)
print(f"[s2] our-clean: NP={len(p):,} A={Ao:.2f} V={Vo:.4f}",
      flush=True)

m = loadmat(str(OUT / "ml_clean_s2.mat"))
mp, mt = m["pc"], m["tc"].astype(int) - 1
mb = om.interp_bathymetry(mp, mt, dem, method="cell-averaging",
                          min_depth=1.0, nan_fill=True)
Am, Vm = metrics(mp, mt, mb)
print(f"[s2] ML-clean(same raw): NP={len(mp):,} A={Am:.2f} "
      f"V={Vm:.4f}", flush=True)
print(f"[s2] clean delta (ml-ours): A {Am-Ao:+.2f} km2, "
      f"V {Vm-Vo:+.4f} km3", flush=True)
