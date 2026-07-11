# Coverage difference map: where does THEIR mesh cover water that
# OURS does not (and vice versa), weighted by depth?
import os, sys, logging
import numpy as np
logging.basicConfig(level=logging.WARNING, stream=sys.stdout)
os.environ.setdefault("MPLBACKEND", "Agg")
sys.path.insert(0, os.path.expanduser("~/Github/oceanmesh"))
import matplotlib.pyplot as plt
from pathlib import Path
from matplotlib.path import Path as MPath

OUT = Path("outputs/om2d_examples/jbay")
f14 = Path(os.path.expanduser(
    "~/Github/OceanMesh2D/Examples/5_JBAY_w_weirs_mesh.14"))
with open(f14) as f:
    f.readline()
    ne, npn = map(int, f.readline().split()[:2])
    gp = np.empty((npn, 2)); gb = np.empty(npn)
    for i in range(npn):
        q = f.readline().split()
        gp[i] = float(q[1]), float(q[2]); gb[i] = float(q[3])
    gt = np.empty((ne, 3), dtype=int)
    for i in range(ne):
        q = f.readline().split()
        gt[i] = int(q[2]), int(q[3]), int(q[4])
gt -= 1
p = np.load(OUT / "p.npy"); t = np.load(OUT / "t.npy")
b = np.load(OUT / "b.npy")

# lattice coverage via matplotlib triangulation finder
import matplotlib.tri as mtri
x = np.linspace(-73.97, -73.75, 1100)
y = np.linspace(40.50, 40.68, 900)
X, Y = np.meshgrid(x, y)
def cover(pp, tt):
    tr = mtri.Triangulation(pp[:, 0], pp[:, 1], tt)
    f = tr.get_trifinder()
    return f(X.ravel(), Y.ravel()).reshape(X.shape) >= 0
c_ml = cover(gp, gt)
c_us = cover(p, t)
only_ml = c_ml & ~c_us
only_us = c_us & ~c_ml
cell_km2 = (0.22/1099*111*np.cos(np.deg2rad(40.59))) * (0.18/899*111)
print(f"[cov] only-THEIRS: {only_ml.sum()*cell_km2:.2f} km2  "
      f"only-OURS: {only_us.sum()*cell_km2:.2f} km2", flush=True)
# depth-weight the difference using their nodal depths
from scipy.interpolate import NearestNDInterpolator
FZ = NearestNDInterpolator(gp, gb)
Z = FZ(np.column_stack([X.ravel(), Y.ravel()])).reshape(X.shape)
v_ml_only = (Z[only_ml].clip(min=0).sum() * cell_km2) / 1e3
v_us_only = (Z[only_us].clip(min=0).sum() * cell_km2) / 1e3
print(f"[cov] volume in only-THEIRS: {v_ml_only:.4f} km3  "
      f"only-OURS: {v_us_only:.4f} km3", flush=True)
fig, ax = plt.subplots(figsize=(13, 10))
img = np.zeros(X.shape)
img[only_ml] = 1; img[only_us] = -1
ax.pcolormesh(X, Y, img, cmap="bwr_r", vmin=-1, vmax=1,
              shading="auto")
ax.set_aspect(1/np.cos(np.deg2rad(40.59)))
ax.set_title("coverage diff: red=only THEIRS, blue=only OURS")
fig.savefig(OUT / "coverage_diff.png", dpi=160, bbox_inches="tight")
print("[cov] saved", flush=True)
