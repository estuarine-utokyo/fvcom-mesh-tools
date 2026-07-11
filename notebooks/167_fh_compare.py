# Layered verification against the MATLAB edgefx dump:
#  (1) deterministic layer: our finalize_sizing grid vs OM2D's hh,
#      point-by-point on the SAME lattice
#  (2) generator layer: mesh from OM2D's OWN sizing values through
#      our generate_mesh -> compare with the golden mesh
import os, sys, logging
import numpy as np
logging.basicConfig(level=logging.WARNING)
os.environ.setdefault("MPLBACKEND", "Agg")
sys.path.insert(0, os.path.expanduser("~/Github/oceanmesh"))
from pathlib import Path
from scipy.io import loadmat
from scipy.interpolate import RegularGridInterpolator
import oceanmesh as om
from oceanmesh import Shoreline
from oceanmesh.grid import Grid

OM2D = Path(os.path.expanduser("~/Github/OceanMesh2D"))
DEG = 1.0 / 111e3
bbox = (166.0, 176.0, -48.0, -40.0)
m = loadmat(str(OM2D / "Examples/1_NZ_edgefx_dump.mat"))
xg, yg, hh = m["xg"].ravel(), m["yg"].ravel(), m["hh"]  # hh meters
print(f"[fh] MATLAB grid {hh.shape} "
      f"h range {hh.min():.0f}-{hh.max():.0f} m", flush=True)

# ---- (1) our sizing on their lattice
shore = Shoreline(str(OM2D / "datasets/GSHHS_shp/f/GSHHS_f_L1.shp"),
                  bbox, 1e3 * DEG)
sdf = om.signed_distance_function(shore)
feat = om.feature_sizing_function(shore, sdf, r=3,
                                  max_edge_length=100e3 * DEG)
grid, _ = om.finalize_sizing([feat], shoreline=shore, hmin=1e3,
                             max_edge_length=100e3,
                             max_edge_length_nearshore=5e3,
                             gradation=0.35)
X, Y = np.meshgrid(xg, yg, indexing="ij")
q = np.column_stack([X.ravel(), Y.ravel()])
ours_m = grid.eval(q).reshape(X.shape) / DEG  # -> meters
ratio = ours_m / hh
# water-side emphasis: where either is fine (< 20 km)
sel = (hh < 20e3) | (ours_m < 20e3)
print(f"[fh] ratio ours/OM2D  all: p10/p50/p90 = "
      f"{np.percentile(ratio, [10,50,90]).round(3)}", flush=True)
print(f"[fh] ratio (fine<20km): p10/p50/p90 = "
      f"{np.percentile(ratio[sel], [10,50,90]).round(3)}", flush=True)

# ---- (2) mesh from THEIR sizing through OUR generator
gm = Grid(bbox=bbox, dx=float(np.diff(xg).mean()),
          dy=float(np.diff(yg).mean()), hmin=1e3 * DEG, values=0.0,
          extrapolate=True, crs="EPSG:4326")
# fill our Grid object with THEIR values (meters -> degrees)
interp = RegularGridInterpolator((xg, yg), hh * DEG,
                                 bounds_error=False, fill_value=None)
lon, lat = gm.create_grid()
gm.values = interp(np.column_stack([lon.ravel(), lat.ravel()])
                   ).reshape(lon.shape)
gm.build_interpolant()
p, t = om.generate_mesh(sdf, gm, max_iter=100, seed=0)
from oceanmesh.fix_mesh import simp_qual
pp = p.copy(); pp[:, 0] *= np.cos(np.deg2rad(-44.0))
qq = simp_qual(pp, t)
print(f"[gen] from OM2D sizing: NP={len(p):,} NT={len(t):,} "
      f"qual mean/min = {qq.mean():.4f}/{qq.min():.4f} "
      f"(golden NP=6,017 mean 0.9344 min 0.4087)", flush=True)
