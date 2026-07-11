# Example_2_NY side-by-side: our translation vs MATLAB golden.
import os
import numpy as np
os.environ.setdefault("MPLBACKEND", "Agg")
import matplotlib.pyplot as plt
from pathlib import Path
from fvcom_mesh_tools.gridref import GridRef
from fvcom_mesh_tools.plotting import add_atlas_grid

# fixed atlas for the whole NY domain: cols A-J west->east,
# rows 1-10 north->south, 0.05 deg cells (stable across zooms)
NY_GRID = GridRef(-74.25, 40.50, -73.75, 41.00, 0.05, 0.05)

OUT = Path("outputs/om2d_examples/ex2_ny")
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
print(f"golden NP={npn:,} NE={ne:,}; ours NP={len(p):,} NT={len(t):,}")

views = {
    "full": (-74.26, -73.74, 40.49, 41.01, 0.15),
    "jamaica_bay": (-73.97, -73.73, 40.52, 40.70, 0.35),
    "upper_bay": (-74.12, -73.94, 40.60, 40.78, 0.35),
    "hudson_mouth": (-74.26, -74.02, 40.49, 40.70, 0.3),
}
for name, (x0, x1, y0, y1, lw) in views.items():
    asp = 1 / np.cos(np.deg2rad(0.5 * (y0 + y1)))
    h = 8.5 * (y1 - y0) * asp / (x1 - x0)
    fig, axes = plt.subplots(1, 2, figsize=(17, max(h, 5)))
    for ax, (pp, tt, title, n_) in zip(
        axes,
        [(p, t, "improved oceanmesh (Python)", len(p)),
         (gp, gt, "OceanMesh2D (MATLAB, golden)", len(gp))],
    ):
        ax.triplot(pp[:, 0], pp[:, 1], tt, lw=lw, color="steelblue")
        ax.set_xlim(x0, x1); ax.set_ylim(y0, y1)
        ax.set_aspect(asp)
        add_atlas_grid(ax, crs="EPSG:4326", grid=NY_GRID)
        ax.set_title(f"{title} NP={n_:,}")
    fig.suptitle(f"Example_2_NY (PostSandyNCEI, h0=30 m) — {name}")
    fig.tight_layout()
    fig.savefig(OUT / f"cmp_{name}.png", dpi=170, bbox_inches="tight")
    print("saved", name, flush=True)
