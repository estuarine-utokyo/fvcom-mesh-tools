# Side-by-side ECGC comparison: our mesh vs the MATLAB TestECGC
# mesh (both post bound_courant_number).
import os, sys
import numpy as np
os.environ.setdefault("MPLBACKEND", "Agg")
import matplotlib.pyplot as plt
from pathlib import Path
from scipy.io import loadmat

OUT = Path("outputs/om2d_examples/ecgc")
p = np.load(OUT / "p.npy"); t = np.load(OUT / "t.npy")
m = loadmat(str(OUT / "matlab_ecgc.mat"))
gp, gt = m["pc"], m["tc"].astype(int) - 1

views = {
    "full": (-85.0, -64.0, 24.0, 42.7, 0.1),
    "newyork": (-74.6, -71.6, 40.0, 41.4, 0.4),
    "chesapeake": (-77.5, -74.8, 36.6, 39.6, 0.35),
    "hatteras": (-77.2, -74.8, 34.4, 36.4, 0.4),
}
for name, (x0, x1, y0, y1, lw) in views.items():
    asp = 1 / np.cos(np.deg2rad(0.5 * (y0 + y1)))
    h = 8 * (y1 - y0) * asp / (x1 - x0)
    fig, axes = plt.subplots(1, 2, figsize=(16, max(h, 4.5)))
    for ax, (pp, tt, title, np_) in zip(
        axes,
        [(p, t, "improved oceanmesh (Python)", len(p)),
         (gp, gt, "OceanMesh2D (MATLAB, golden)", len(gp))],
    ):
        ax.triplot(pp[:, 0], pp[:, 1], tt, lw=lw, color="steelblue")
        ax.set_xlim(x0, x1); ax.set_ylim(y0, y1)
        ax.set_aspect(asp)
        ax.set_title(f"{title} NP={np_:,}")
    fig.suptitle(f"TestECGC (Example_3 nest 1 + Courant bounds) — {name}")
    fig.tight_layout()
    fig.savefig(OUT / f"cmp_{name}.png", dpi=180, bbox_inches="tight")
    print("saved", name, flush=True)
