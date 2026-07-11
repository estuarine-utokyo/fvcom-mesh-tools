# Verdict: does MATLAB's clean keep the NE lobe when applied to
# OUR raw mesh?
import os
import numpy as np
os.environ.setdefault("MPLBACKEND", "Agg")
import matplotlib.pyplot as plt
from pathlib import Path
from scipy.io import loadmat
from fvcom_mesh_tools.gridref import GridRef
from fvcom_mesh_tools.plotting import add_atlas_grid

NY_GRID = GridRef(-74.25, 40.50, -73.75, 41.00, 0.05, 0.05)
OUT = Path("outputs/om2d_examples/ex2_ny")
m = loadmat(str(OUT / "matlab_cleaned_ours_ex2.mat"))
mp, mt = m["pc"], m["tc"].astype(int) - 1
p = np.load(OUT / "p.npy"); t = np.load(OUT / "t.npy")

LOBE = (-74.118, -74.082, 40.763, 40.790)
def lobe_nt(pp, tt):
    c = pp[tt].mean(axis=1)
    return int(((c[:, 0] > LOBE[0]) & (c[:, 0] < LOBE[1])
                & (c[:, 1] > LOBE[2]) & (c[:, 1] < LOBE[3])).sum())
print(f"[verdict] NE-lobe elements: ml-clean(our-raw)="
      f"{lobe_nt(mp, mt)} our-clean={lobe_nt(p, t)}", flush=True)

x0, x1, y0, y1 = -74.16, -74.06, 40.73, 40.81
asp = 1 / np.cos(np.deg2rad(40.77))
fig, axes = plt.subplots(1, 2, figsize=(17, 9))
for ax, (pp, tt, title) in zip(
    axes,
    [(mp, mt,
      f"MATLAB clean on OUR raw NP={len(mp):,}"),
     (p, t, f"our clean on our raw NP={len(p):,}")],
):
    ax.triplot(pp[:, 0], pp[:, 1], tt, lw=0.3, color="steelblue")
    ax.set_xlim(x0, x1); ax.set_ylim(y0, y1)
    ax.set_aspect(asp)
    add_atlas_grid(ax, crs="EPSG:4326", grid=NY_GRID)
    ax.set_title(title)
fig.suptitle("Cross-clean verdict on the C5/D5 lobe (same raw input)")
fig.tight_layout()
fig.savefig(OUT / "cmp_c5_cross_clean.png", dpi=170,
            bbox_inches="tight")
print("[verdict] saved", flush=True)
