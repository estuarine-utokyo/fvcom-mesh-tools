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

OM2D = Path(os.path.expanduser("~/Github/OceanMesh2D"))
DEG = 1.0 / 111e3
bbox = (166.0, 176.0, -48.0, -40.0)
m = loadmat(str(OM2D / "Examples/1_NZ_edgefx_dump.mat"))
xg, yg, hh = m["xg"].ravel(), m["yg"].ravel(), m["hh"]
Fg = RegularGridInterpolator((xg, yg), hh, bounds_error=False,
                             fill_value=None)
shore = Shoreline(str(OM2D/"datasets/GSHHS_shp/f/GSHHS_f_L1.shp"),
                  bbox, 1e3*DEG)
sdf = om.signed_distance_function(shore)
feat = om.feature_sizing_function(shore, sdf, r=3,
                                  max_edge_length=100e3*DEG)
grid, _ = om.finalize_sizing([feat], shoreline=shore, hmin=1e3,
                             max_edge_length=100e3,
                             max_edge_length_nearshore=5e3,
                             gradation=0.35)
pts = np.array([[166.001,-39.999],[166.05,-40.05],[166.2,-40.2],
                [166.001,-47.999],[166.05,-47.95],[166.2,-47.8],
                [175.999,-47.999],[175.95,-47.95],[175.8,-47.8],
                [171.0,-44.0]])
ours = grid.eval(pts)/DEG/1e3
gold = Fg(pts)/1e3
for q, a, b in zip(pts, ours, gold):
    print(f"[corner-fh] ({q[0]:7.3f},{q[1]:7.3f}) ours={a:7.1f} km "
          f"OM2D={b:7.1f} km", flush=True)
# also the raw feature component at the corner
f0 = feat.eval(pts)/DEG/1e3
print("[corner-fh] feature-only ours km:", np.round(f0,1), flush=True)
