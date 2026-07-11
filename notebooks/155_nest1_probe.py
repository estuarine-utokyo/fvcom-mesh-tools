import faulthandler, logging, os, sys
faulthandler.enable()
faulthandler.dump_traceback_later(90, repeat=True)
logging.basicConfig(level=logging.INFO, stream=sys.stdout,
                    format="%(levelname)s %(name)s: %(message)s")
os.environ.setdefault("MPLBACKEND", "Agg")
from pathlib import Path
import numpy as np
import oceanmesh as om
from oceanmesh import DEM, Region, Shoreline

OM2D = Path(os.path.expanduser("~/Github/OceanMesh2D"))
exec(open("notebooks/150_tb_varres_3r_translation.py").read().split(
    "NESTS = [")[0].split("import oceanmesh as om")[1].split(
    "OUT = Path")[0] if False else "")
# minimal nest1 setup (same params as translation)
GSHHS = OM2D / "datasets/GSHHS_shp/f/GSHHS_f_L1.shp"
SRTM_PACIFIC = OM2D / "datasets/TokyoBay/dem/SRTM15_pacific_4min.nc"
DEG = 1.0 / 111e3
import numpy as np
x1 = [157.4634361,158.84262217,160.046858,161.07987092,161.94757837,162.65742188,163.21770287,163.63704393,163.92386584,164.08625707,164.1323057,164.06882521,163.9009549,163.63627547,163.31281843,162.89825466,162.38857843,161.80134206,161.13903484,160.4052709,159.60336662,158.7360147,157.80598566,156.81590455,155.7687459,154.66746366,153.5151905,152.3154099,151.0718658,149.78866673,148.4702309,147.12135316,145.74715292,144.35301788,142.94462003,141.52778842,140.10846886,138.69263968,137.2862313,135.89504337,134.52470652,133.18055169,131.86766327,130.59072034,129.35412832,128.16190388,127.01767503,125.92498131,124.88679107,123.90615375,122.98581257,122.12855282,121.60387397,115.64257669,158.59463098,167.27851682,157.4634361]
y1 = [51.36354136,50.35540912,49.27890879,48.14482943,46.96329793,45.74371019,44.4946652,43.2239559,41.93859546,40.64492778,39.34870289,38.0551143,36.76910533,35.56167205,34.36768078,33.12965026,31.91772707,30.7320848,29.57667543,28.45502315,27.37062929,26.32729387,25.32881846,24.37921688,23.48234609,22.64223232,21.86295783,21.14846084,20.50262818,19.92912188,19.4314301,19.01266254,18.67552551,18.42237308,18.25486888,18.17421142,18.18096748,18.27508101,18.45587251,18.7220402,19.07182143,19.5027674,20.01215433,20.59669911,21.25297733,21.97724182,22.76552456,23.61401919,24.51851667,25.47512814,26.47984282,27.52881366,28.22438211,39.405273,75.05163401,62.25276725,51.36354136]
bbox_01 = np.column_stack([x1, y1])
reg = Region((float(bbox_01[:,0].min()), float(bbox_01[:,0].max()),
              float(bbox_01[:,1].min()), float(bbox_01[:,1].max())), 4326)
shore = Shoreline(str(GSHHS), bbox_01, 10e3 * DEG)
sdf = om.signed_distance_function(shore)
dem = DEM(str(SRTM_PACIFIC), bbox=reg)
comps = [om.feature_sizing_function(shore, sdf, r=3, max_edge_length=50e3*DEG),
         om.wavelength_sizing_function(dem, wl=30, min_edgelength=10e3*DEG,
                                       max_edge_length=50e3*DEG)]
grid, dt_used = om.finalize_sizing(comps, dem=dem, hmin=10e3,
                                   max_edge_length=50e3, gradation=0.3,
                                   courant={"timestep": 0.0, "max": 0.5})
print("[probe] sizing ok, hmin(grid attr) =", grid.hmin, flush=True)
p, t = om.generate_mesh(sdf, grid, max_iter=5, seed=0)
print("[probe] nest1 alone: NP =", len(p), flush=True)
