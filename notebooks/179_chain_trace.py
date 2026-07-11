# Step-by-step trace of our clean chain on ours_raw, mirroring the
# genuine MATLAB trace: db 101 -> collapse 0 -> MMBT 599+50 ->
# valence(9->8) -> smooth2d -> min 0.276
import os, sys, logging
import numpy as np
logging.basicConfig(level=logging.INFO, stream=sys.stdout,
                    format="%(levelname)s %(name)s: %(message)s")
os.environ.setdefault("MPLBACKEND", "Agg")
sys.path.insert(0, os.path.expanduser("~/Github/oceanmesh"))
from pathlib import Path
from oceanmesh.clean import (_external_topology,
                             make_mesh_boundaries_traversable)
from oceanmesh.fix_mesh import fix_mesh, simp_qual
from oceanmesh.mesh_improve import (collapse_thin_triangles,
                                    bound_connectivity)
from oceanmesh.smooth2d import smooth2d

OUT = Path("outputs/om2d_examples/ex1_nz")
p = np.load(OUT / "p_raw.npy"); t = np.load(OUT / "t_raw.npy")
CPP = np.cos(np.deg2rad(-44.0))
p[:, 0] *= CPP

def rep(tag, p, t):
    q = simp_qual(p, t)
    print(f"[trace] {tag:14s}: NP={len(p):,} NT={len(t):,} "
          f"min={q.min():.4f} n<0.25={int((q<0.25).sum())}",
          flush=True)

rep("raw", p, t)
# db-loop
total_del = 0
for _i in range(25):
    _, bv = _external_topology(p, t)
    touching = np.isin(t, bv).any(axis=1)
    q = simp_qual(p, t)
    bad = touching & (q < 0.25)
    if not bad.any():
        break
    total_del += int(bad.sum())
    t = t[~bad]
    p, t, _ = fix_mesh(p, t, delete_unused=True)
print(f"[trace] db-loop deleted {total_del} (their 101)", flush=True)
rep("after db", p, t)
p, t = collapse_thin_triangles(p, t, min_qual=0.25)
rep("after collapse", p, t)
p, t = make_mesh_boundaries_traversable(p, t,
                                        min_disconnected_area=0.25)
rep("after MMBT", p, t)
p, t = bound_connectivity(p, t, max_valence=9)
rep("after valence", p, t)
p, t = smooth2d(p, t, disp_every=8)
rep("after smooth2d", p, t)
