import os, sys, logging
import numpy as np
logging.basicConfig(level=logging.INFO, stream=sys.stdout,
                    format="%(levelname)s %(name)s: %(message)s")
os.environ.setdefault("MPLBACKEND", "Agg")
sys.path.insert(0, os.path.expanduser("~/Github/oceanmesh"))
from pathlib import Path
from oceanmesh.clean import om2d_default_clean
from oceanmesh.fix_mesh import simp_qual

OUT = Path("outputs/om2d_examples/ex1_nz")
p = np.load(OUT / "p_raw.npy"); t = np.load(OUT / "t_raw.npy")
CPP = np.cos(np.deg2rad(-44.0))
p[:, 0] *= CPP
p2, t2 = om2d_default_clean(p, t)
q2 = simp_qual(p2, t2)
print(f"[clean] OURS   NP={len(p2):,} NT={len(t2):,} "
      f"mean/min = {q2.mean():.4f}/{q2.min():.4f} "
      f"n<0.25={int((q2<0.25).sum())}", flush=True)
print("[clean] MATLAB target: NP=6,619 NT=10,473 mean 0.8877 "
      "min 0.2807", flush=True)

from oceanmesh.clean import _external_topology
bad = np.where(q2 < 0.25)[0]
_, bv = _external_topology(p2, t2)
bset = set(np.unique(np.asarray(bv).ravel()).astype(int).tolist())
val = np.bincount(t2.ravel(), minlength=len(p2))
for e in bad[:20]:
    tri = t2[e]
    onb = sum(1 for v in tri if int(v) in bset)
    L = np.linalg.norm(np.roll(p2[tri], -1, 0) - p2[tri], axis=1)
    Ls = np.sort(L / L.max())
    kind = "needle" if Ls[0] < 0.35 else "cap"
    print(f"[surv] q={q2[e]:.4f} bnd_verts={onb} {kind} "
          f"Lratio={np.round(Ls,2)} val={[int(val[v]) for v in tri]} "
          f"at {np.round(p2[tri].mean(0), 3)}", flush=True)
