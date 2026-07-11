import os, sys, logging
import numpy as np
logging.basicConfig(level=logging.WARNING)
os.environ.setdefault("MPLBACKEND", "Agg")
sys.path.insert(0, os.path.expanduser("~/Github/oceanmesh"))
from pathlib import Path
import oceanmesh as om
from oceanmesh import Shoreline
from oceanmesh.clean import (_external_topology,
                             make_mesh_boundaries_traversable,
                             laplacian2)
from oceanmesh.fix_mesh import fix_mesh, simp_qual
from oceanmesh.mesh_improve import (collapse_thin_triangles,
                                    bound_connectivity)

OM2D = Path(os.path.expanduser("~/Github/OceanMesh2D"))
DEG = 1.0 / 111e3
bbox = (166.0, 176.0, -48.0, -40.0)
BOX = (166.3, 166.85, -46.25, -45.3)
CPP = np.cos(np.deg2rad(-44.0))

def nbox(p, tag, c=CPP):
    q = p.copy(); q[:, 0] = q[:, 0] / c   # back to lon for the box
    n = ((q[:,0]>=BOX[0])&(q[:,0]<=BOX[1])
         &(q[:,1]>=BOX[2])&(q[:,1]<=BOX[3])).sum()
    print(f"[step] {tag:22s}: NP={len(p):6,} box={n}", flush=True)

shore = Shoreline(str(OM2D/"datasets/GSHHS_shp/f/GSHHS_f_L1.shp"),
                  bbox, 1e3*DEG)
sdf = om.signed_distance_function(shore)
feat = om.feature_sizing_function(shore, sdf, r=3,
                                  max_edge_length=100e3*DEG)
grid, _ = om.finalize_sizing([feat], shoreline=shore, hmin=1e3,
                             max_edge_length=100e3,
                             max_edge_length_nearshore=5e3,
                             gradation=0.35)
p, t = om.generate_mesh(sdf, grid, max_iter=100, seed=0,
                        cleanup="none")
nbox(p, "raw (lonlat)", c=1.0)
# replicate the clean in the CPP frame as generate_mesh does
p[:, 0] *= CPP
for rnd in range(3):
    n_in = len(t)
    for _i in range(25):
        _, bv = _external_topology(p, t)
        touching = np.isin(t, bv).any(axis=1)
        q = simp_qual(p, t)
        bad = touching & (q < 0.25)
        if not bad.any():
            break
    # single-iteration counts only; apply once
        t = t[~bad]
        p, t, _ = fix_mesh(p, t, delete_unused=True)
    nbox(p, f"r{rnd} after db-loop")
    p, t = collapse_thin_triangles(p, t, min_qual=0.25)
    nbox(p, f"r{rnd} after collapse")
    p, t = make_mesh_boundaries_traversable(p, t,
                                            min_disconnected_area=0.25)
    nbox(p, f"r{rnd} after traversable")
    p, t = bound_connectivity(p, t, max_valence=9)
    nbox(p, f"r{rnd} after valence")
    p, t = laplacian2(p, t)
    nbox(p, f"r{rnd} after laplacian")
    if simp_qual(p, t).min() >= 0.25 or len(t) == n_in:
        break
