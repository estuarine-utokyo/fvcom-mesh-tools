# OM2D acceptance ladder, step 1: Example_1_NZ (1:1 translation).
# Ground truth from OceanMesh2D Tests/TestSanity.m:
#   NP = 5,968 +- 500; NT = 9,530 +- 1,500; min qual >= 0.25
#   reference qual (mean / lower-3rd / min) = 0.9325 / 0.7452 / 0.3548
# No DEM involved: isolates geodata chain + feature sizing +
# gradation + meshgen (+ auto-clean once P0-1 lands).
import logging
import os
import sys
import time
from pathlib import Path

import numpy as np

logging.basicConfig(level=logging.INFO, stream=sys.stdout,
                    format="%(levelname)s %(name)s: %(message)s")
os.environ.setdefault("MPLBACKEND", "Agg")

import oceanmesh as om  # noqa: E402
from oceanmesh import Shoreline  # noqa: E402

OM2D = Path(os.path.expanduser("~/Github/OceanMesh2D"))
GSHHS = OM2D / "datasets/GSHHS_shp/f/GSHHS_f_L1.shp"
OUT = Path("outputs/om2d_examples/ex1_nz")
OUT.mkdir(parents=True, exist_ok=True)
DEG = 1.0 / 111e3

bbox = (166.0, 176.0, -48.0, -40.0)
min_el = 1e3
max_el = 100e3
max_el_ns = 5e3
grade = 0.35
R = 3

t0 = time.time()
shore = Shoreline(str(GSHHS), bbox, min_el * DEG)
sdf = om.signed_distance_function(shore)
feat = om.feature_sizing_function(shore, sdf, r=R,
                                  max_edge_length=max_el * DEG)
grid, _ = om.finalize_sizing(
    [feat], shoreline=shore,
    hmin=min_el, max_edge_length=max_el,
    max_edge_length_nearshore=max_el_ns,
    gradation=grade,
)
print(f"[ex1] sizing done +{time.time()-t0:.0f}s", flush=True)
p, t = om.generate_mesh(sdf, grid, max_iter=100, seed=0)
print(f"[ex1] meshgen done +{time.time()-t0:.0f}s", flush=True)

from oceanmesh.fix_mesh import simp_qual  # noqa: E402

q = simp_qual(p, t)
print(f"[ex1] RAW    NP={len(p):,} NT={len(t):,} qual "
      f"mean/L3/min = {q.mean():.4f}/"
      f"{np.mean(np.sort(q)[:max(1, len(q)//3)]):.4f}/{q.min():.4f}",
      flush=True)

# targets
print("[ex1] TARGET NP=5,968+-500 NT=9,530+-1,500 "
      "qual 0.9325/0.7452/0.3548 (min>=0.25)", flush=True)
np.save(OUT / "p.npy", p)
np.save(OUT / "t.npy", t)

import matplotlib.pyplot as plt  # noqa: E402

fig, ax = plt.subplots(figsize=(9, 9))
ax.triplot(p[:, 0], p[:, 1], t, lw=0.25, color="steelblue")
ax.set_aspect(1 / np.cos(np.deg2rad(-44)))
ax.set_title(f"Example_1_NZ translation NP={len(p):,}")
fig.savefig(OUT / "ex1_nz.png", dpi=180, bbox_inches="tight")
fig, ax = plt.subplots(figsize=(9, 9))
ax.triplot(p[:, 0], p[:, 1], t, lw=0.4, color="steelblue")
ax.set_xlim(172, 176); ax.set_ylim(-42, -39)
ax.set_aspect(1 / np.cos(np.deg2rad(-40.5)))
ax.set_title("Example_1_NZ subdomain (as .m step 6)")
fig.savefig(OUT / "ex1_nz_sub.png", dpi=180, bbox_inches="tight")
print("[ex1] figures saved", flush=True)
