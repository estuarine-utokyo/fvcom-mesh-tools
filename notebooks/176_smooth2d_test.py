# Acceptance for the smooth2d port: apply the om2d clean chain
# with ds=2 -> smooth2d to OUR RAW Example_1 mesh; target = the
# genuine MATLAB clean on the same input (6,619 / min 0.2807).
import os, sys, logging
import numpy as np
logging.basicConfig(level=logging.INFO, stream=sys.stdout,
                    format="%(levelname)s %(name)s: %(message)s")
os.environ.setdefault("MPLBACKEND", "Agg")
sys.path.insert(0, os.path.expanduser("~/Github/oceanmesh"))
from pathlib import Path
from oceanmesh.smooth2d import smooth2d
from oceanmesh.fix_mesh import simp_qual

OUT = Path("outputs/om2d_examples/ex1_nz")
p = np.load(OUT / "p_raw.npy"); t = np.load(OUT / "t_raw.npy")
CPP = np.cos(np.deg2rad(-44.0))
p[:, 0] *= CPP   # work in the physical frame like the loop does
q = simp_qual(p, t)
print(f"[s2d] RAW    NP={len(p):,} NT={len(t):,} "
      f"mean/min = {q.mean():.4f}/{q.min():.4f}", flush=True)
p2, t2 = smooth2d(p, t)
q2 = simp_qual(p2, t2)
print(f"[s2d] SMOOTH NP={len(p2):,} NT={len(t2):,} "
      f"mean/min = {q2.mean():.4f}/{q2.min():.4f} "
      f"n<0.25={int((q2<0.25).sum())}", flush=True)
print("[s2d] MATLAB target: NP=6,619 NT=10,473 mean 0.8877 "
      "min 0.2807", flush=True)
