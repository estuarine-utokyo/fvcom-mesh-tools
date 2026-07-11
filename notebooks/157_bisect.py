import os, sys, logging, time
import numpy as np
os.environ.setdefault("MPLBACKEND", "Agg")
logging.basicConfig(level=logging.WARNING)
import matplotlib.pyplot as plt
import oceanmesh as om
from oceanmesh import DEM, Region, Shoreline
from fvcom_mesh_tools.plotting import _add_coast, add_atlas_grid
from pathlib import Path

CFG = sys.argv[1]           # current | oldinpoly | noresample
if CFG == "oldinpoly":
    os.environ["OCEANMESH_INPOLY_METHOD"] = "cython"
if CFG == "noresample":
    os.environ["OCEANMESH_NO_RESAMPLE"] = "1"
DEG = 1.0/111e3
OM2D = Path(os.path.expanduser("~/Github/OceanMesh2D"))
x3 = [139.142156,139.142156,139.144724,139.152351,139.164815,139.18178,139.202839,139.227552,139.25548,139.286203,139.31934,139.354542,139.391505,139.429958,139.469663,139.510408,139.552007,139.594291,139.637105,139.680308,139.723765,139.767348,139.810932,139.854389,139.897592,139.940406,139.98269,140.024289,140.065034,140.104739,140.143192,140.180155,140.215357,140.248493,140.279217,140.307145,140.331858,140.352917,140.369882,140.382346,140.389973,140.392541,140.392541,139.142156]
y3 = [35.859592,35.225989,35.182486,35.139579,35.097818,35.057672,35.019512,34.983607,34.95014,34.919216,34.890891,34.865177,34.842062,34.82152,34.803513,34.788004,34.774954,34.764328,34.756096,34.750234,34.746724,34.745556,34.746724,34.750234,34.756096,34.764328,34.774954,34.788004,34.803513,34.82152,34.842062,34.865177,34.890891,34.919216,34.95014,34.983607,35.019512,35.057672,35.097818,35.139579,35.182486,35.225989,35.859592,35.859592]
bbox3 = np.column_stack([x3, y3])
reg = Region((139.142, 140.393, 34.745, 35.860), 4326)
t0 = time.time()
shore = Shoreline("outputs/tb_varres_3r/land_osm_wide.shp", bbox3, 100.0*DEG)
sdf = om.signed_distance_function(shore)
dem = DEM(str(OM2D/"datasets/TokyoBay/dem/SRTM15_kanto_15s.nc"), bbox=reg)
comps = [om.feature_sizing_function(shore, sdf, r=3, max_edge_length=500*DEG),
         om.wavelength_sizing_function(dem, wl=30, min_edgelength=100*DEG,
                                       max_edge_length=500*DEG)]
grid, _ = om.finalize_sizing(comps, dem=dem, hmin=100, max_edge_length=500,
                             gradation=0.1, courant=None)
p, t = om.generate_mesh(sdf, grid, max_iter=20, seed=0)
p, t = om.make_mesh_boundaries_traversable(p, t)
print(f"[{CFG}] NP={len(p):,} gen {time.time()-t0:.0f}s", flush=True)
fig, ax = plt.subplots(figsize=(10, 10))
_add_coast(ax, (139.0, 34.5, 141.3, 36.2), "EPSG:4326")
ax.triplot(p[:,0], p[:,1], t, lw=0.35, color="steelblue")
ax.set_xlim(139.60, 140.05); ax.set_ylim(34.95, 35.30)
ax.set_aspect(1/np.cos(np.deg2rad(35.12)))
add_atlas_grid(ax, crs="EPSG:4326")
ax.set_title(f"bisect {CFG}: nest3-only mouth")
fig.savefig(f"outputs/figures/bisect_{CFG}.png", dpi=200, bbox_inches="tight")
print(f"[{CFG}] saved", flush=True)
