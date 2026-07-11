# Fiordland fan: sizing-field fault vs clean-erosion fault
import os, sys, logging
import numpy as np
logging.basicConfig(level=logging.WARNING)
os.environ.setdefault("MPLBACKEND", "Agg")
sys.path.insert(0, os.path.expanduser("~/Github/oceanmesh"))
from pathlib import Path
from scipy.io import loadmat
import oceanmesh as om
from oceanmesh import Shoreline

OM2D = Path(os.path.expanduser("~/Github/OceanMesh2D"))
DEG = 1.0 / 111e3
bbox = (166.0, 176.0, -48.0, -40.0)
BOX = (166.3, 166.85, -46.25, -45.3)   # the broken fan region

# (1) sizing transect through the fan at y=-45.7
m = loadmat(str(OM2D / "Examples/1_NZ_edgefx_dump.mat"))
xg, yg, hh = m["xg"].ravel(), m["yg"].ravel(), m["hh"]
from scipy.interpolate import RegularGridInterpolator
F_gold = RegularGridInterpolator((xg, yg), hh)
shore = Shoreline(str(OM2D / "datasets/GSHHS_shp/f/GSHHS_f_L1.shp"),
                  bbox, 1e3 * DEG)
sdf = om.signed_distance_function(shore)
feat = om.feature_sizing_function(shore, sdf, r=3,
                                  max_edge_length=100e3 * DEG)
grid, _ = om.finalize_sizing([feat], shoreline=shore, hmin=1e3,
                             max_edge_length=100e3,
                             max_edge_length_nearshore=5e3,
                             gradation=0.35)
xs = np.linspace(166.3, 167.0, 8)
pts = np.column_stack([xs, np.full(8, -45.7)])
ours = grid.eval(pts) / DEG
gold = F_gold(pts)
print("[probe] transect y=-45.7 (km):", flush=True)
print("[probe]   ours :", (ours/1e3).round(1), flush=True)
print("[probe]   OM2D :", (gold/1e3).round(1), flush=True)

# (2)(3) vertex counts in the fan box: final vs cleanup-none
OUT = Path("outputs/om2d_examples/ex1_nz")
p = np.load(OUT / "p.npy")
inbox = ((p[:,0]>=BOX[0])&(p[:,0]<=BOX[1])
         &(p[:,1]>=BOX[2])&(p[:,1]<=BOX[3]))
print(f"[probe] final mesh vertices in fan box: {inbox.sum()}",
      flush=True)
p2, t2 = om.generate_mesh(sdf, grid, max_iter=100, seed=0,
                          cleanup="none")
inbox2 = ((p2[:,0]>=BOX[0])&(p2[:,0]<=BOX[1])
          &(p2[:,1]>=BOX[2])&(p2[:,1]<=BOX[3]))
print(f"[probe] RAW (cleanup=none) vertices in fan box: "
      f"{inbox2.sum()}  (NP={len(p2):,})", flush=True)
