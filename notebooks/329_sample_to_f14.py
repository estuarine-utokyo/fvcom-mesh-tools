# Convert the goto2023 sample (grd/dep/obc) to fort.14 (UTM54N) so
# fmesh-mesh-qa can gate the ORIGINAL sample as the benchmark.
import os
import numpy as np
from collections import defaultdict
import sys
sys.path.insert(0, os.path.expanduser("~/Github/oceanmesh"))
import oceanmesh as om

G = os.path.expanduser('~/Github/TB-FVCOM/input/goto2023/grid/')
gd = open(G + 'TokyoBay_grd.dat').read().split('\n')
nn = int(gd[0].split('=')[1]); ne = int(gd[1].split('=')[1])
T = np.array([[int(w) for w in gd[2 + i].split()[1:4]]
              for i in range(ne)]) - 1
P = np.array([[float(w) for w in gd[2 + ne + i].split()[1:3]]
              for i in range(nn)])
dep = np.loadtxt(G + 'TokyoBay_dep.dat', skiprows=1)[:, 2]
obc = [int(l.split()[1]) - 1 for l in
       open(G + 'TokyoBay_obc.dat').read().strip().split('\n')[1:]]

loops = om.boundary_loops(T)


def loop_area(lp):
    x, y = P[lp, 0], P[lp, 1]
    return 0.5 * float(np.dot(x, np.roll(y, -1))
                       - np.dot(y, np.roll(x, -1)))


outer = max(loops, key=lambda lp: abs(loop_area(lp)))
ring = np.asarray(outer)
i0 = int(np.where(ring == obc[0])[0][0])
ring = np.roll(ring, -i0)
if ring[1] != obc[1]:
    ring = np.roll(ring[::-1], 1)
assert list(ring[:len(obc)]) == list(obc), "obc not contiguous"
land = np.append(ring[len(obc) - 1:], ring[0])
bc = {"open": [np.asarray(obc)], "land": [land],
      "island": [lp for lp in loops if lp is not outer]}
om.write_fort14("outputs/sample_repro/sample_original.14", P, T,
                depth=dep, boundaries=bc)
print(f"wrote sample_original.14 NP={nn} NE={ne} "
      f"open={len(obc)} land={len(land)} "
      f"island={len(bc['island'])}", flush=True)
