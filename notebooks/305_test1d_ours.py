# Ladder: Tests/test_1d_original.m — our mesh1d port on the same
# 16-vertex polygon with fh = min(0.15*|x| + 0.01, 0.20), h0=0.01.
import os, sys
import numpy as np
sys.path.insert(0, os.path.expanduser("~/Github/oceanmesh"))
from pathlib import Path
from oceanmesh.mesh1d import mesh1d

OUT = Path("outputs/om2d_examples/test1d")
OUT.mkdir(parents=True, exist_ok=True)

poly = np.array([
    [0.2483, 0.4942], [0.1838, 0.6080], [0.2851, 0.7861],
    [0.4741, 0.8445], [0.5616, 0.8883], [0.6953, 0.7277],
    [0.6745, 0.3891], [0.8474, 0.3161], [0.9211, 0.6168],
    [0.9764, 0.6752], [0.9741, 0.1876], [0.7967, 0.0241],
    [0.6838, 0.0358], [0.3658, 0.1905], [0.3520, 0.345],
    [0.2483, 0.4942]])

minh, maxh, grade = 0.01, 0.20, 0.15
fh = lambda q: np.minimum(grade * np.abs(q[:, 0]) + minh, maxh)

p, t = mesh1d(poly, fh, minh)
print(f"[test1d] ours NP={len(p)} NE={len(t)}")
np.save(OUT/"p_ours.npy", p); np.save(OUT/"t_ours.npy", t)

# spacing-vs-target diagnostics
seg = p[t]
mid = seg.mean(axis=1)
ell = np.hypot(*(seg[:, 1] - seg[:, 0]).T)
ratio = ell / fh(mid)
print(f"[test1d] edge/fh ratio p5={np.percentile(ratio,5):.3f} "
      f"p50={np.percentile(ratio,50):.3f} p95={np.percentile(ratio,95):.3f}")
