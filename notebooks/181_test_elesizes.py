# Ladder step 2a: port of OM2D Tests/TestEleSizes.m — NZ mesh at
# grades 0.15/0.25/0.35; >= 95% of nodal min-bar resolutions must
# lie within [min_el, max_el].
import os, sys, logging
import numpy as np
logging.basicConfig(level=logging.WARNING)
os.environ.setdefault("MPLBACKEND", "Agg")
sys.path.insert(0, os.path.expanduser("~/Github/oceanmesh"))
from pathlib import Path
import oceanmesh as om
from oceanmesh import Shoreline

OM2D = Path(os.path.expanduser("~/Github/OceanMesh2D"))
DEG = 1.0 / 111e3
bbox = (166.0, 176.0, -48.0, -40.0)
min_el, max_el, max_el_ns, R = 1e3, 100e3, 5e3, 3

shore = Shoreline(str(OM2D/"datasets/GSHHS_shp/f/GSHHS_f_L1.shp"),
                  bbox, min_el*DEG)
sdf = om.signed_distance_function(shore)
feat = om.feature_sizing_function(shore, sdf, r=R,
                                  max_edge_length=max_el*DEG)
ok_all = True
for grade in (0.15, 0.25, 0.35):
    import copy
    g = copy.deepcopy(feat)
    grid, _ = om.finalize_sizing([g], shoreline=shore, hmin=min_el,
                                 max_edge_length=max_el,
                                 max_edge_length_nearshore=max_el_ns,
                                 gradation=grade)
    p, t = om.generate_mesh(sdf, grid, max_iter=100, seed=0)
    e = np.vstack([t[:, [0,1]], t[:, [1,2]], t[:, [2,0]]])
    e = np.unique(np.sort(e, axis=1), axis=0)
    # GetBarLengths equivalent: physical bar lengths via tmerc-ish
    lat = 0.5*(p[e[:,0],1] + p[e[:,1],1])
    dx = (p[e[:,0],0]-p[e[:,1],0])*np.cos(np.deg2rad(lat))
    dyv = p[e[:,0],1]-p[e[:,1],1]
    L = np.sqrt(dx**2+dyv**2)/DEG
    reso = np.full(len(p), np.inf)
    np.minimum.at(reso, e[:,0], L)
    np.minimum.at(reso, e[:,1], L)
    reso = reso[np.isfinite(reso)]
    inb = 100*np.mean((reso > min_el) & (reso < max_el))
    ok = inb >= 95.0
    ok_all &= ok
    print(f"[elesz] grade {grade:.2f}: NP={len(p):,} in-bounds "
          f"{inb:.2f}% (>=95) {'PASS' if ok else 'FAIL'}",
          flush=True)
print(f"[elesz] TestEleSizes: {'PASSED' if ok_all else 'FAILED'}",
      flush=True)
