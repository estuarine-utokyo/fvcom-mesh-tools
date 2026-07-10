# Smallest decisive test for the stereo improvement fix: Ex7
# pipeline at 10x coarser (h0=40 km). PASS = vertex count does NOT
# grow through iterations (pre-fix: 2.7M -> 6M at h0=4 km).
import os, sys, logging, time
import faulthandler
faulthandler.enable()
faulthandler.dump_traceback_later(900, repeat=True)
import numpy as np
logging.basicConfig(level=logging.INFO, stream=sys.stdout,
                    format="%(levelname)s %(name)s: %(message)s")
sys.path.insert(0, os.path.expanduser("~/Github/oceanmesh"))
from pathlib import Path
import oceanmesh as om
from oceanmesh import DEM, Region, Shoreline

OM2D = Path(os.path.expanduser("~/Github/OceanMesh2D"))
DS = OM2D / "datasets"
SCR = Path(os.path.expanduser(
    "~/Github/fvcom-mesh-tools/outputs/om2d_examples/ex7glob"))
DEG = 1.0 / 111e3
t0 = time.time()

bbox = (-180.0, 180.0, -89.0, 90.0)
reg = Region(bbox, 4326)
sh = Shoreline(str(SCR/"gshhs_l1l6.shp"), reg.bbox, 40e3*DEG)
sdf = om.signed_distance_function(sh)
dem = DEM(str(DS/"SRTM15+.nc"), bbox=reg)
print(f"[mini7] shoreline+dem +{time.time()-t0:.0f}s", flush=True)
f = om.feature_sizing_function(sh, sdf, r=3, max_edge_length=200e3*DEG)
w = om.wavelength_sizing_function(dem, wl=30, min_edgelength=40e3*DEG,
                                  max_edge_length=200e3*DEG,
                                  grid_dx=40e3*DEG)
g, dta = om.finalize_sizing([f, w], dem=dem, shoreline=sh,
                            hmin=40e3, max_edge_length=200e3,
                            gradation=0.25,
                            courant={"timestep": 0.0})
print(f"[mini7] sizing dt={dta:.2f} +{time.time()-t0:.0f}s", flush=True)

sh2 = Shoreline(str(SCR/"gshhs_l1l6_stereo.shp"), reg.bbox, 40e3*DEG,
                stereo=True)
dom = om.signed_distance_function(sh2)
p, t = om.generate_mesh(dom, g, stereo=True, max_iter=50, seed=0)
print(f"[mini7] mesh NP={len(p):,} NT={len(t):,} +{time.time()-t0:.0f}s",
      flush=True)
np.save(SCR/"p_mini.npy", p)
