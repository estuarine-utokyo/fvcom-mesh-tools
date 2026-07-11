# Post-fix C5 verification: mesh zoom + ring recount.
import os, sys, logging
import numpy as np
logging.basicConfig(level=logging.WARNING, stream=sys.stdout)
os.environ.setdefault("MPLBACKEND", "Agg")
sys.path.insert(0, os.path.expanduser("~/Github/oceanmesh"))
import matplotlib.pyplot as plt
from pathlib import Path
from fvcom_mesh_tools.gridref import GridRef
from fvcom_mesh_tools.plotting import add_atlas_grid
from oceanmesh import DEM, Region, Shoreline

NY_GRID = GridRef(-74.25, 40.50, -73.75, 41.00, 0.05, 0.05)
DSET = Path(os.path.expanduser(
    "~/Github/OceanMesh2D/datasets/PostSandyNCEI"))
OUT = Path("outputs/om2d_examples/ex2_ny")
DEG = 1.0 / 111e3

# ring recount with the NEW classification
dem_probe = DEM(str(DSET / "PostSandyNCEI.nc"))
reg = Region(dem_probe.bbox, 4326)
shore = Shoreline(str(DSET / "PostSandyNCEI.shp"), reg.bbox,
                  30.0 * DEG)
CELL = (-74.15, -74.10, 40.75, 40.80)
for name, arr in (("inner", shore.inner),
                  ("mainland", shore.mainland)):
    a = np.asarray(arr, dtype=float)
    if a.size == 0:
        print(f"[c5b] ours.{name}: EMPTY", flush=True)
        continue
    inbox = ((a[:, 0] > CELL[0]) & (a[:, 0] < CELL[1])
             & (a[:, 1] > CELL[2]) & (a[:, 1] < CELL[3]))
    print(f"[c5b] ours.{name}: pts-in-C5={int(inbox.sum()):,}",
          flush=True)

p = np.load(OUT / "p.npy"); t = np.load(OUT / "t.npy")
f14 = Path(os.path.expanduser(
    "~/Github/OceanMesh2D/Examples/2_NY_mesh.14"))
with open(f14) as f:
    f.readline()
    ne, npn = map(int, f.readline().split()[:2])
    gp = np.empty((npn, 2))
    for i in range(npn):
        parts = f.readline().split()
        gp[i] = float(parts[1]), float(parts[2])
    gt = np.empty((ne, 3), dtype=int)
    for i in range(ne):
        parts = f.readline().split()
        gt[i] = int(parts[2]), int(parts[3]), int(parts[4])
gt -= 1
x0, x1, y0, y1 = -74.16, -74.06, 40.73, 40.81
asp = 1 / np.cos(np.deg2rad(40.77))
fig, axes = plt.subplots(1, 2, figsize=(17, 9))
for ax, (pp, tt, title, n_) in zip(
    axes,
    [(p, t, "improved oceanmesh (Python)", len(p)),
     (gp, gt, "OceanMesh2D (MATLAB, golden)", len(gp))],
):
    ax.triplot(pp[:, 0], pp[:, 1], tt, lw=0.3, color="steelblue")
    ax.set_xlim(x0, x1); ax.set_ylim(y0, y1)
    ax.set_aspect(asp)
    add_atlas_grid(ax, crs="EPSG:4326", grid=NY_GRID)
    ax.set_title(f"{title} NP={n_:,}")
fig.suptitle("Example_2_NY C4/C5 Hackensack Meadowlands (post-fix)")
fig.tight_layout()
fig.savefig(OUT / "cmp_c5_meadowlands.png", dpi=170,
            bbox_inches="tight")
print("[c5b] saved", flush=True)
