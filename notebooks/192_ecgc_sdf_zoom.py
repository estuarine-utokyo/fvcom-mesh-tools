# Zoom probe: SDF sign plaid block, numba vs cython inpoly, and
# outer-ring geometry audit (sliver parts, duplicate vertices).
import os, sys, logging
import numpy as np
logging.basicConfig(level=logging.WARNING, stream=sys.stdout)
os.environ.setdefault("MPLBACKEND", "Agg")
sys.path.insert(0, os.path.expanduser("~/Github/oceanmesh"))
import matplotlib.pyplot as plt
from pathlib import Path
from oceanmesh import Shoreline

OM2D = Path(os.path.expanduser("~/Github/OceanMesh2D"))
DEG = 1.0 / 111e3
bbox_poly = np.array([
    [-71.6, 42.7], [-64, 30], [-80, 24], [-85, 38], [-71.6, 42.7]])
shore = Shoreline(str(OM2D/"datasets/GSHHS_shp/f/GSHHS_f_L1.shp"),
                  bbox_poly, 1e3*DEG)

outer = np.asarray(shore.boubox)
inner = np.asarray(shore.inner)
# ring audit
idx = np.where(np.isnan(outer[:, 0]))[0]
start = 0
for k, stop in enumerate(list(idx) + [len(outer)]):
    seg = outer[start:stop]; start = stop + 1
    if len(seg) < 3:
        continue
    from oceanmesh.geodata import _poly_area
    a = abs(_poly_area(seg[:, 0], seg[:, 1]))
    d = np.sqrt((np.diff(seg[:, 0])**2 + np.diff(seg[:, 1])**2))
    ndup = int((d < 1e-12).sum())
    print(f"[zoom] outer part {k}: n={len(seg)} area={a:.3e} "
          f"dup-verts={ndup} bbox=({seg[:,0].min():.2f},"
          f"{seg[:,0].max():.2f},{seg[:,1].min():.2f},"
          f"{seg[:,1].max():.2f})", flush=True)

from oceanmesh import edges as om_edges
from oceanmesh.geometry import inpoly2
poly = np.vstack((outer, inner))
e = om_edges.get_poly_edges(poly)
x = np.linspace(-80.4, -78.0, 1200)
y = np.linspace(24.2, 26.2, 1000)
X, Y = np.meshgrid(x, y)
q = np.column_stack([X.ravel(), Y.ravel()])
res = {}
for meth in ("numba", "cython"):
    os.environ["OCEANMESH_INPOLY_METHOD"] = meth
    import importlib, oceanmesh.geometry.point_in_polygon as pip
    importlib.reload(pip)
    ins, _ = pip.inpoly2(q, np.nan_to_num(poly), e)
    res[meth] = ins
    print(f"[zoom] {meth}: inside={int(ins.sum())}", flush=True)
diff = res["numba"] ^ res["cython"]
print(f"[zoom] kernel mismatches: {int(diff.sum())}", flush=True)

fig, axes = plt.subplots(1, 2, figsize=(17, 7))
for ax, meth in zip(axes, ("numba", "cython")):
    ax.pcolormesh(X, Y, res[meth].reshape(X.shape),
                  cmap="coolwarm_r", shading="auto")
    ax.plot(outer[:, 0], outer[:, 1], "k.-", lw=0.6, ms=1.5)
    if len(inner):
        ax.plot(inner[:, 0], inner[:, 1], "m-", lw=0.5)
    ax.set_xlim(-80.4, -78.0); ax.set_ylim(24.2, 26.2)
    ax.set_title(f"inpoly[{meth}] inside(blue)")
out = Path("outputs/om2d_examples/ecgc/sdf_zoom.png")
fig.savefig(out, dpi=150, bbox_inches="tight")
print("[zoom] saved", flush=True)
