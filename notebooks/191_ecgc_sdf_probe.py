# Probe the ECGC polygon-boubox shoreline: plot classified rings
# and the SDF sign on a lattice over the striped SW region.
import os, sys, logging
import numpy as np
logging.basicConfig(level=logging.INFO, stream=sys.stdout)
os.environ.setdefault("MPLBACKEND", "Agg")
sys.path.insert(0, os.path.expanduser("~/Github/oceanmesh"))
import matplotlib.pyplot as plt
from pathlib import Path
import oceanmesh as om
from oceanmesh import Shoreline

OM2D = Path(os.path.expanduser("~/Github/OceanMesh2D"))
DEG = 1.0 / 111e3
bbox_poly = np.array([
    [-71.6, 42.7], [-64, 30], [-80, 24], [-85, 38], [-71.6, 42.7]])
shore = Shoreline(str(OM2D/"datasets/GSHHS_shp/f/GSHHS_f_L1.shp"),
                  bbox_poly, 1e3*DEG)
sdf = om.signed_distance_function(shore)
print("outer nan-parts:",
      int(np.isnan(np.asarray(shore.boubox)[:, 0]).sum()))
print("mainland pts:", len(shore.mainland),
      "inner pts:", len(shore.inner))

# closure check of every outer/mainland part
for label, arr in (("outer", np.asarray(shore.boubox)),
                   ("mainland", np.asarray(shore.mainland))):
    if len(arr) == 0:
        continue
    idx = np.where(np.isnan(arr[:, 0]))[0]
    start = 0
    nopen = 0
    for stop in list(idx) + [len(arr)]:
        seg = arr[start:stop]
        start = stop + 1
        if len(seg) < 3:
            continue
        if not np.allclose(seg[0], seg[-1]):
            nopen += 1
    print(f"[probe] {label}: {nopen} OPEN parts", flush=True)

# SDF sign lattice over the striped zone
x = np.linspace(-80.5, -76.5, 900)
y = np.linspace(24.2, 27.5, 700)
X, Y = np.meshgrid(x, y)
d = sdf.eval(np.column_stack([X.ravel(), Y.ravel()]))
D = d.reshape(X.shape)
fig, ax = plt.subplots(figsize=(12, 9))
ax.pcolormesh(X, Y, (D < 0), cmap="coolwarm_r", shading="auto")
for arr, c in ((shore.boubox, "k"), (shore.mainland, "g"),
               (shore.inner, "m")):
    a = np.asarray(arr)
    if len(a):
        ax.plot(a[:, 0], a[:, 1], c, lw=0.7)
ax.set_title("SDF sign (blue=inside) + rings: outer k, main g, inner m")
out = Path("outputs/om2d_examples/ecgc/sdf_probe.png")
fig.savefig(out, dpi=150, bbox_inches="tight")
print("[probe] saved", flush=True)
