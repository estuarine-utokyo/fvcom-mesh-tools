import os, sys
import numpy as np
sys.path.insert(0, os.path.expanduser("~/Github/oceanmesh"))
from pathlib import Path
from scipy.io import loadmat
OUT = Path("outputs/om2d_examples/jbay")
corners = np.array([[-73.97, 40.5], [-73.75, 40.5],
                    [-73.75, 40.68], [-73.97, 40.68]])
DEG = 1.0/111e3
def census(name, p):
    for i, c in enumerate(corners):
        d = np.hypot(p[:, 0]-c[0], p[:, 1]-c[1]) / DEG
        n100 = int((d < 100).sum())
        print(f"[cc] {name} corner{i} ({c[0]},{c[1]}): "
              f"n<100m={n100} dmin={d.min():.1f} m", flush=True)
m = loadmat(os.path.expanduser(
    "~/Github/OceanMesh2D/Examples/Precleaned_grid.mat"))
census("ML-raw", np.asarray(m["p"], float))
census("ours-s2", np.load(OUT / "p_raw_fresh_s2.npy"))
census("ours-s1", np.load(OUT / "p_raw_s1.npy"))
print("[cc] done", flush=True)
