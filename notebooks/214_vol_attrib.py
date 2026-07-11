# Volume attribution: our CA depths on THEIR mesh vs their depths.
import os, sys, logging
import numpy as np
logging.basicConfig(level=logging.WARNING, stream=sys.stdout)
sys.path.insert(0, os.path.expanduser("~/Github/oceanmesh"))
from pathlib import Path
import oceanmesh as om
from oceanmesh import DEM, Region

f14 = Path(os.path.expanduser(
    "~/Github/OceanMesh2D/Examples/5_JBAY_w_weirs_mesh.14"))
with open(f14) as f:
    f.readline()
    ne, npn = map(int, f.readline().split()[:2])
    gp = np.empty((npn, 2)); gb = np.empty(npn)
    for i in range(npn):
        parts = f.readline().split()
        gp[i] = float(parts[1]), float(parts[2])
        gb[i] = float(parts[3])
    gt = np.empty((ne, 3), dtype=int)
    for i in range(ne):
        parts = f.readline().split()
        gt[i] = int(parts[2]), int(parts[3]), int(parts[4])
gt -= 1

def vol(p, t, b):
    Re2 = 111.0 ** 2
    X, Y = p[:, 0], p[:, 1]
    x1, y1 = X[t[:, 0]], Y[t[:, 0]]
    x2, y2 = X[t[:, 1]], Y[t[:, 1]]
    x3, y3 = X[t[:, 2]], Y[t[:, 2]]
    pa = 0.5 * np.abs((x2-x1)*(y3-y1) - (x3-x1)*(y2-y1))
    cosf = np.cos(np.deg2rad((y1+y2+y3)/3))
    bc = b[t].mean(axis=1) / 1e3
    return float((pa*cosf).sum()*Re2), float((pa*cosf*bc).sum()*Re2)

a_ml, v_ml = vol(gp, gt, gb)
print(f"[vol] THEIR mesh, THEIR depths: area={a_ml:.2f} "
      f"vol={v_ml:.4f}", flush=True)

reg = Region((-73.97, -73.75, 40.5, 40.68), 4326)
dem = DEM(os.path.expanduser(
    "~/Github/OceanMesh2D/datasets/PostSandyNCEI/PostSandyNCEI.nc"),
    bbox=reg)
b2 = om.interp_bathymetry(gp, gt, dem, method="cell-averaging",
                          min_depth=1.0, nan_fill=True)
a2, v2 = vol(gp, gt, b2)
print(f"[vol] THEIR mesh, OUR CA depths: vol={v2:.4f} "
      f"(delta {100*(v2-v_ml)/v_ml:+.2f}%)", flush=True)
db = b2 - gb
print(f"[vol] depth diff p10/50/90 = "
      f"{np.percentile(db,10):+.3f}/{np.percentile(db,50):+.3f}/"
      f"{np.percentile(db,90):+.3f} m", flush=True)

OUT = Path("outputs/om2d_examples/jbay")
p = np.load(OUT / "p.npy"); t = np.load(OUT / "t.npy")
b = np.load(OUT / "b.npy")
a3, v3 = vol(p, t, b)
print(f"[vol] OUR mesh, OUR depths: area={a3:.2f} vol={v3:.4f}",
      flush=True)
