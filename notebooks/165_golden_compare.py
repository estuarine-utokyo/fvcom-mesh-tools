# Compare our Example_1_NZ mesh against the OM2D golden mesh
# produced on GENKAI (MATLAB R2024a, TestSanity passed).
import os, sys
import numpy as np
os.environ.setdefault("MPLBACKEND", "Agg")
sys.path.insert(0, os.path.expanduser("~/Github/oceanmesh"))
from pathlib import Path
from fvcom_mesh_tools.io import read_fort14
from scipy.spatial import cKDTree

DEG = 1.0 / 111e3
G = Path(os.path.expanduser("~/Github/OceanMesh2D/Examples"))
golden = None
for name in ("1_NZ_mesh.14", "1_NZ_mesh.grd", "fort.14"):
    if (G / name).exists():
        golden = read_fort14(str(G / name)); break
assert golden is not None, "golden fort.14 not found"
OUT = Path("outputs/om2d_examples/ex1_nz")
p = np.load(OUT / "p.npy"); t = np.load(OUT / "t.npy")
gp, gt = golden.nodes[:, :2], golden.elements

def nodal_res(p, t):
    e = np.vstack([t[:, [0, 1]], t[:, [1, 2]], t[:, [2, 0]]])
    e = np.unique(np.sort(e, axis=1), axis=0)
    L = np.linalg.norm(p[e[:, 0]] - p[e[:, 1]], axis=1) / DEG
    r = np.full(len(p), np.inf)
    np.minimum.at(r, e[:, 0], L); np.minimum.at(r, e[:, 1], L)
    return r, e, L

r_py, e_py, L_py = nodal_res(p, t)
r_gd, e_gd, L_gd = nodal_res(gp, gt)
print(f"[cmp] NP py={len(p):,} golden={len(gp):,}", flush=True)

# ratio of our resolution to golden's at golden nodes
tree = cKDTree(p)
_, nearest = tree.query(gp)
ratio = r_py[nearest] / r_gd
print(f"[cmp] res ratio (py/golden) p10/p50/p90 = "
      f"{np.percentile(ratio, [10, 50, 90]).round(2)}", flush=True)
# by distance to coast proxy: golden res value bands
for lo, hi, tag in [(0, 2e3, "fine<2k"), (2e3, 10e3, "2-10k"),
                    (10e3, 1e9, ">10k")]:
    s = (r_gd >= lo) & (r_gd < hi)
    if s.sum():
        print(f"[cmp] golden {tag:8s}: n={s.sum():6,} "
              f"ratio p50={np.percentile(ratio[s], 50):.2f}", flush=True)

# anisotropy: mean physical edge length by orientation (E-W vs N-S)
def aniso(p, e, L, tag):
    dx = np.abs(p[e[:, 0], 0] - p[e[:, 1], 0])
    dy = np.abs(p[e[:, 0], 1] - p[e[:, 1], 1])
    ew = dx > 2 * dy   # mostly E-W edges
    ns = dy > 2 * dx
    lat = 0.5 * (p[e[:, 0], 1] + p[e[:, 1], 1])
    Lm = L.copy()
    Lm[ew] = np.sqrt((dx[ew] * np.cos(np.deg2rad(lat[ew])))**2
                     + dy[ew]**2)[0:ew.sum()] / DEG if False else Lm[ew]
    # physical length: recompute with cos(lat) on the lon component
    dxm = dx * np.cos(np.deg2rad(lat))
    Lphys = np.sqrt(dxm**2 + dy**2) / DEG
    print(f"[cmp] {tag}: phys edge p50 EW={np.percentile(Lphys[ew],50):7.0f}"
          f"  NS={np.percentile(Lphys[ns],50):7.0f}  "
          f"EW/NS={np.percentile(Lphys[ew],50)/np.percentile(Lphys[ns],50):.3f}",
          flush=True)

aniso(p, e_py, L_py, "python")
aniso(gp, e_gd, L_gd, "golden")
