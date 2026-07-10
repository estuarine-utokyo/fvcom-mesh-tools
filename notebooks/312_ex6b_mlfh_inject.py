# Decisive A/B: mesh Ex6b stage-1 with the MATLAB banded fh field
# injected verbatim (matlab_ex6b_fh.mat). If initial NP ~60k (as
# MATLAB gets from this field), our generator is fine and the field
# comparison missed something; if ~50k, the generator's sampling of
# banded fields is at fault.
import os, sys, logging, time
import numpy as np
import h5py
logging.basicConfig(level=logging.INFO, stream=sys.stdout,
                    format="%(levelname)s %(name)s: %(message)s")
sys.path.insert(0, os.path.expanduser("~/Github/oceanmesh"))
from pathlib import Path
import oceanmesh as om
from oceanmesh import Region, Shoreline, Grid

OM2D = Path(os.path.expanduser("~/Github/OceanMesh2D"))
DS = OM2D / "datasets"
OUT = Path("outputs/om2d_examples/ex6bfp")
DEG = 1.0 / 111e3
t0 = time.time()

bbox = (-95.40, -94.4, 29.14, 30.09)
reg = Region(bbox, 4326)
sh = Shoreline(
    str(DS/"US_Medium_Shoreline/us_medium_shoreline_polygon.shp"),
    reg.bbox, 60.0*DEG)
sh.detect_inpoly_flip(str(DS/"GSHHS_shp/l/GSHHS_l_L1.shp"))
sdf = om.signed_distance_function(sh)

with h5py.File(OUT/"matlab_ex6b_fh.mat", "r") as f:
    xg = np.asarray(f["xg"]).ravel()
    yg = np.asarray(f["yg"]).ravel()
    hh = np.asarray(f["hh"])
if hh.shape == (len(yg), len(xg)):
    hh = hh.T
# ML values are metres; our Grid values must be degrees
g = Grid(bbox=(float(xg[0]), float(xg[-1]), float(yg[0]), float(yg[-1])),
         dx=float(np.mean(np.diff(xg))), dy=float(np.mean(np.diff(yg))),
         extrapolate=True, values=hh * DEG, crs=4326,
         hmin=60.0 * DEG)
g.build_interpolant()
print(f"[inject] ML fh grid {hh.shape} dx={g.dx:.2e} dy={g.dy:.2e} "
      f"+{time.time()-t0:.0f}s", flush=True)

p, t = om.generate_mesh(sdf, g, max_iter=100, seed=0)
print(f"[inject] final NP={len(p):,} NT={len(t):,} "
      f"(MATLAB got 49,564 from this field; ours-own-field 41,963) "
      f"+{time.time()-t0:.0f}s", flush=True)
np.save(OUT/"p_inject.npy", p); np.save(OUT/"t_inject.npy", t)
