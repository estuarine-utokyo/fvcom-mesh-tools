import os
import numpy as np
os.environ.setdefault("MPLBACKEND", "Agg")
from collections import defaultdict
from pyproj import Transformer
import shapely
from shapely.ops import unary_union

gd = open(os.path.expanduser(
    '~/Github/TB-FVCOM/input/goto2023/grid/TokyoBay_grd.dat')).read().split('\n')
nn = int(gd[0].split('=')[1]); ne = int(gd[1].split('=')[1])
T = np.array([[int(w) for w in gd[2+i].split()[1:4]] for i in range(ne)]) - 1
P = np.array([[float(w) for w in gd[2+ne+i].split()[1:3]] for i in range(nn)])
obc = set(int(l.split()[1])-1 for l in open(os.path.expanduser(
    '~/Github/TB-FVCOM/input/goto2023/grid/TokyoBay_obc.dat')).read().strip().split('\n')[1:])

cnt = defaultdict(int)
for a,b,c in T:
    for e in ((a,b),(b,c),(c,a)): cnt[tuple(sorted(e))] += 1
bedges = [e for e,k in cnt.items() if k==1]
bnodes = sorted({v for e in bedges for v in e})
coast_edges = [e for e in bedges if e[0] not in obc and e[1] not in obc]
L = np.array([np.linalg.norm(P[a]-P[b]) for a,b in coast_edges])
print(f"[anatomy] boundary nodes total={len(bnodes)} (of NP {nn}) "
      f"coast edges={len(coast_edges)}", flush=True)
print(f"[anatomy] coast edge length p10/p50/p90 = "
      f"{np.percentile(L,[10,50,90]).round(0)} m  min/max = "
      f"{L.min():.0f}/{L.max():.0f}", flush=True)

# deviation of sample coast nodes from OSM true-land boundary
from fvcom_mesh_tools.prep import fetch_true_land
g = fetch_true_land((139.55, 34.90, 140.20, 35.80),
                    min_water_area_deg2=5e-5)
gm = g.to_crs(32654)
coast_lines = unary_union(list(gm.geometry)).boundary
cn = [v for v in bnodes if v not in obc]
pts = shapely.points(P[cn,0], P[cn,1])
d = shapely.distance(pts, coast_lines)
print(f"[anatomy] sample coast-node deviation from OSM coast: "
      f"p50/p90/max = {np.percentile(d,50):.0f}/{np.percentile(d,90):.0f}"
      f"/{d.max():.0f} m", flush=True)
# node spacing vs deviation → is boundary a resampled version of a FINE line?
