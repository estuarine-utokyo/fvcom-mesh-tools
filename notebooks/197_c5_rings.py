# C5 (Hackensack Meadowlands) ring-level comparison: our
# classified shoreline vs MATLAB gdat rings + our SDF sign.
import os, sys, logging
import numpy as np
logging.basicConfig(level=logging.WARNING, stream=sys.stdout)
os.environ.setdefault("MPLBACKEND", "Agg")
sys.path.insert(0, os.path.expanduser("~/Github/oceanmesh"))
import matplotlib.pyplot as plt
from pathlib import Path
from scipy.io import loadmat
import oceanmesh as om
from oceanmesh import DEM, Region, Shoreline

DSET = Path(os.path.expanduser(
    "~/Github/OceanMesh2D/datasets/PostSandyNCEI"))
OUT = Path("outputs/om2d_examples/ex2_ny")
DEG = 1.0 / 111e3
dem_probe = DEM(str(DSET / "PostSandyNCEI.nc"))
x0, x1, y0, y1 = dem_probe.bbox
reg = Region((x0, x1, y0, y1), 4326)
shore = Shoreline(str(DSET / "PostSandyNCEI.shp"), reg.bbox,
                  30.0 * DEG)
sdf = om.signed_distance_function(shore)

m = loadmat(str(OUT / "matlab_gdat.mat"))
CELL = (-74.15, -74.10, 40.75, 40.80)  # C5

def clipped(a):
    a = np.asarray(a, dtype=float)
    if a.size == 0:
        return a.reshape(0, 2)
    return a

fig, axes = plt.subplots(1, 3, figsize=(21, 7.5))
# panel 1: our rings; panel 2: MATLAB rings; panel 3: our SDF sign
for ax, (rings, title) in zip(
    axes[:2],
    [((shore.boubox, shore.mainland, shore.inner),
      "improved oceanmesh rings"),
     ((m["outer"], m["mainland"], m["inner"]),
      "OceanMesh2D gdat rings")],
):
    for arr, c, lb in zip(rings, ("k", "g", "m"),
                          ("outer", "mainland", "inner")):
        a = clipped(arr)
        if len(a):
            ax.plot(a[:, 0], a[:, 1], c, lw=0.8, label=lb)
    ax.set_xlim(CELL[0], CELL[1]); ax.set_ylim(CELL[2], CELL[3])
    ax.set_aspect(1 / np.cos(np.deg2rad(40.77)))
    ax.legend(loc="upper right", fontsize=8)
    ax.set_title(title)
xg = np.linspace(CELL[0], CELL[1], 700)
yg = np.linspace(CELL[2], CELL[3], 700)
X, Y = np.meshgrid(xg, yg)
d = sdf.eval(np.column_stack([X.ravel(), Y.ravel()]))
axes[2].pcolormesh(X, Y, (d.reshape(X.shape) < 0),
                   cmap="coolwarm_r", shading="auto")
a = clipped(m["mainland"])
if len(a):
    axes[2].plot(a[:, 0], a[:, 1], "g", lw=0.5)
axes[2].set_xlim(CELL[0], CELL[1]); axes[2].set_ylim(CELL[2], CELL[3])
axes[2].set_aspect(1 / np.cos(np.deg2rad(40.77)))
axes[2].set_title("our SDF sign (blue=water) + ML mainland (g)")
fig.suptitle("C5 Hackensack Meadowlands: ring-level comparison")
fig.tight_layout()
fig.savefig(OUT / "c5_rings.png", dpi=170, bbox_inches="tight")
# ring counts within C5
for name, arr in (("ours.inner", shore.inner),
                  ("ml.inner", m["inner"]),
                  ("ours.mainland", shore.mainland),
                  ("ml.mainland", m["mainland"])):
    a = clipped(arr)
    if not len(a):
        print(f"[c5] {name}: EMPTY", flush=True)
        continue
    inbox = ((a[:, 0] > CELL[0]) & (a[:, 0] < CELL[1])
             & (a[:, 1] > CELL[2]) & (a[:, 1] < CELL[3]))
    nseg = int(np.isnan(a[:, 0]).sum())
    print(f"[c5] {name}: total-pts={len(a):,} pts-in-C5="
          f"{int(inbox.sum()):,} parts={nseg}", flush=True)
print("[c5] saved", flush=True)
