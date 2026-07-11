# 3-way raw coverage: our seed1 vs seed2 vs OM2D Precleaned.
import os, sys, logging
import numpy as np
logging.basicConfig(level=logging.WARNING, stream=sys.stdout)
os.environ.setdefault("MPLBACKEND", "Agg")
sys.path.insert(0, os.path.expanduser("~/Github/oceanmesh"))
import matplotlib.pyplot as plt
import matplotlib.tri as mtri
from pathlib import Path
from scipy.io import loadmat

OUT = Path("outputs/om2d_examples/jbay")
p1 = np.load(OUT / "p_raw_s1.npy"); t1 = np.load(OUT / "t_raw_s1.npy")
p2 = np.load(OUT / "p_raw_fresh_s2.npy")
t2 = np.load(OUT / "t_raw_fresh_s2.npy")
m = loadmat(os.path.expanduser(
    "~/Github/OceanMesh2D/Examples/Precleaned_grid.mat"))
pm = np.asarray(m["p"], float); tm = np.asarray(m["t"], int) - 1
print(f"[cov3] s1 NP={len(p1):,} s2 NP={len(p2):,} "
      f"ml NP={len(pm):,}", flush=True)

x = np.linspace(-73.97, -73.75, 1100)
y = np.linspace(40.50, 40.68, 900)
X, Y = np.meshgrid(x, y)
def cover(pp, tt):
    # even-odd over the mesh's own boundary polygon (referee-
    # verified inpoly2); robust to slivers that break trifinder
    from oceanmesh.cfl import _mesh_boundary_polygon
    from oceanmesh import edges as om_edges
    from oceanmesh.geometry import inpoly2
    poly = _mesh_boundary_polygon(pp, tt)
    e = om_edges.get_poly_edges(poly)
    ins, _ = inpoly2(np.column_stack([X.ravel(), Y.ravel()]),
                     np.nan_to_num(poly), e)
    return ins.reshape(X.shape)
c1, c2, cm = cover(p1, t1), cover(p2, t2), cover(pm, tm)
cell = (0.22/1099*111*np.cos(np.deg2rad(40.59))) * (0.18/899*111)
print(f"[cov3] area s1={c1.sum()*cell:.2f} s2={c2.sum()*cell:.2f} "
      f"ml={cm.sum()*cell:.2f}", flush=True)
print(f"[cov3] s1-not-s2: {(c1&~c2).sum()*cell:.2f} km2  "
      f"s2-not-s1: {(c2&~c1).sum()*cell:.2f}", flush=True)
print(f"[cov3] ml-not-s2: {(cm&~c2).sum()*cell:.2f} km2  "
      f"s2-not-ml: {(c2&~cm).sum()*cell:.2f}", flush=True)

fig, axes = plt.subplots(1, 2, figsize=(19, 9))
for ax, (a, b, ttl) in zip(
    axes,
    [(c1, c2, "red=only seed1 raw, blue=only seed2 raw"),
     (cm, c2, "red=only OM2D raw, blue=only seed2 raw")],
):
    img = np.zeros(X.shape)
    img[a & ~b] = 1; img[b & ~a] = -1
    ax.pcolormesh(X, Y, img, cmap="bwr_r", vmin=-1, vmax=1,
                  shading="auto")
    ax.set_aspect(1/np.cos(np.deg2rad(40.59)))
    ax.set_title(ttl)
fig.tight_layout()
fig.savefig(OUT / "raw_cov3.png", dpi=150, bbox_inches="tight")
print("[cov3] saved", flush=True)
