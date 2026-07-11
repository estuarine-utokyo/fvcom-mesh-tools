import os, sys, logging
import numpy as np
logging.basicConfig(level=logging.WARNING)
os.environ.setdefault("MPLBACKEND", "Agg")
sys.path.insert(0, os.path.expanduser("~/Github/oceanmesh"))
from pathlib import Path
import oceanmesh as om
from oceanmesh import Shoreline
import scipy.spatial

OM2D = Path(os.path.expanduser("~/Github/OceanMesh2D"))
DEG = 1.0 / 111e3
bbox = (166.0, 176.0, -48.0, -40.0)
h0 = 1e3 * DEG
shore = Shoreline(str(OM2D / "datasets/GSHHS_shp/f/GSHHS_f_L1.shp"),
                  bbox, h0)
sdf = om.signed_distance_function(shore)
from oceanmesh.grid import Grid
grid_calc = Grid(bbox=shore.bbox, dx=h0, hmin=h0, values=0.0,
                 extrapolate=True, crs="EPSG:4326")
lon, lat = grid_calc.create_grid()
print(f"[probe] lattice shape lon={lon.shape} dx={grid_calc.dx:.6f} "
      f"dy={grid_calc.dy:.6f}", flush=True)
print(f"[probe] lon varies along axis "
      f"{0 if lon[1,0]!=lon[0,0] else 1}", flush=True)
qpts = np.column_stack((lon.flatten(), lat.flatten()))
d = sdf.eval(qpts).reshape(lon.shape)
water = d < 0
ddx, ddy = np.gradient(d, grid_calc.dx, grid_calc.dy)
gm = np.sqrt(ddx**2 + ddy**2)
print(f"[probe] |grad d| over water: p10/p50/p90 = "
      f"{np.percentile(gm[water], [10,50,90]).round(3)}", flush=True)
mask = (gm < 0.90) & (d < -0.5 * h0)
print(f"[probe] medial cells: {mask.sum():,} "
      f"({100*mask.sum()/water.sum():.1f}% of water)", flush=True)
mp = np.column_stack((lon[mask], lat[mask]))
# distance of medial points to coast (= |d| at those cells)
dd = np.abs(d[mask]) / DEG
print(f"[probe] medial |d| km: p10/p50/p90 = "
      f"{(np.percentile(dd, [10,50,90])/1e3).round(1)}", flush=True)
