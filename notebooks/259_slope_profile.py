import os, sys, cProfile, pstats, io, logging, time
import faulthandler
faulthandler.enable()
faulthandler.dump_traceback_later(120, repeat=True)
import numpy as np
logging.basicConfig(level=logging.INFO, stream=sys.stdout,
                    format="%(levelname)s %(name)s: %(message)s")
sys.path.insert(0, os.path.expanduser("~/Github/oceanmesh"))
from pathlib import Path
import oceanmesh as om
from oceanmesh import DEM, Region

OM2D = Path(os.path.expanduser("~/Github/OceanMesh2D"))
DEG = 1.0 / 111e3
reg = Region((-100.0, -53.0, 5.0, 52.5), 4326)
t0 = time.time()
dem = DEM(str(OM2D/"datasets/SRTM15+.nc"), bbox=reg)
print(f"[prof] DEM +{time.time()-t0:.0f}s", flush=True)
pr = cProfile.Profile(); pr.enable()
s_ = om.bathymetric_gradient_sizing_function(
    dem, slope_parameter=15, filter_quotient=50,
    min_edge_length=1000*DEG, max_edge_length=1000*DEG,
    type_of_filter="barotropic")
pr.disable()
print(f"[prof] slope done +{time.time()-t0:.0f}s", flush=True)
s = io.StringIO(); st = pstats.Stats(pr, stream=s)
st.sort_stats("cumulative").print_stats(22)
print("\n".join(s.getvalue().splitlines()[4:34]), flush=True)
