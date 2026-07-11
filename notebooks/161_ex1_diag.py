# Where do the excess nodes live, and why does min-qual stay low?
import os, sys, logging
import numpy as np
logging.basicConfig(level=logging.WARNING)
os.environ.setdefault("MPLBACKEND", "Agg")
sys.path.insert(0, os.path.expanduser("~/Github/oceanmesh"))
from pathlib import Path
from oceanmesh.fix_mesh import simp_qual
from oceanmesh.clean import om2d_default_clean, _external_topology

OUT = Path("outputs/om2d_examples/ex1_nz")
p = np.load(OUT / "p.npy"); t = np.load(OUT / "t.npy")
DEG = 1.0 / 111e3

# nodal resolution = min connected bar length (meters)
e = np.vstack([t[:, [0, 1]], t[:, [1, 2]], t[:, [2, 0]]])
e = np.unique(np.sort(e, axis=1), axis=0)
L = np.linalg.norm(p[e[:, 0]] - p[e[:, 1]], axis=1) / DEG
res = np.full(len(p), np.inf)
np.minimum.at(res, e[:, 0], L); np.minimum.at(res, e[:, 1], L)

# distance to boundary (proxy for distance to coast)
_, bv = _external_topology(p, t)
from scipy.spatial import cKDTree
bvi = np.unique(np.asarray(bv).ravel()).astype(int)
dcoast = cKDTree(p[bvi]).query(p)[0] / DEG
for lo, hi in [(0, 2e3), (2e3, 5e3), (5e3, 15e3), (15e3, 50e3),
               (50e3, 1e9)]:
    s = (dcoast >= lo) & (dcoast < hi)
    if s.sum():
        print(f"[diag] d_coast {lo/1e3:5.0f}-{hi/1e3:5.0f} km: "
              f"n={s.sum():6,} res p10/p50/p90 = "
              f"{np.percentile(res[s], [10,50,90]).round(0)}", flush=True)

q = simp_qual(p, t)
worst = np.argsort(q)[:8]
bvset = set(bvi.tolist())
for w in worst:
    on_b = any(int(v) in bvset for v in t[w])
    print(f"[diag] worst q={q[w]:.4f} boundary={on_b} "
          f"at {p[t[w,0]].round(3)}", flush=True)

# re-clean probe: does a second default clean lift min qual?
p2, t2 = om2d_default_clean(p.copy(), t.copy())
q2 = simp_qual(p2, t2)
print(f"[diag] re-clean: NP {len(p)}->{len(p2)} min {q.min():.4f}->"
      f"{q2.min():.4f} mean {q.mean():.4f}->{q2.mean():.4f}", flush=True)
