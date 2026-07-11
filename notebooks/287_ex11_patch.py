# Ladder step 12b: Example_11 patch remesh + stitch.
# OM2D uses mesh2d for the patch; we use our generator with the
# patch boundary held fixed (pfix/egfix) — engine difference is
# by design and documented.
import os, sys, logging, time, pickle
import faulthandler
faulthandler.enable()
faulthandler.dump_traceback_later(900, repeat=True)
import numpy as np
logging.basicConfig(level=logging.INFO, stream=sys.stdout,
                    format="%(levelname)s %(name)s: %(message)s")
sys.path.insert(0, os.path.expanduser("~/Github/oceanmesh"))
from pathlib import Path
import oceanmesh as om
from oceanmesh.signed_distance_function import Domain
from oceanmesh.grid import Grid
from oceanmesh import edges as om_edges
from oceanmesh.geometry import inpoly2
import scipy.spatial

OUT = Path("outputs/om2d_examples/ex11rm")
p = np.load(OUT/"p_base.npy"); t = np.load(OUT/"t_base.npy")
hole = np.array([
    [172.7361, -44.0332], [172.4602, -44.3989],
    [171.2402, -45.2966], [171.4459, -46.6907],
    [174.5999, -45.1137], [173.5714, -43.6053],
    [173.0915, -44.1310], [172.7361, -44.0332]])
t0 = time.time()

# --- extract_subdomain: elements with centroid inside hole
ring = np.vstack([hole, [[np.nan, np.nan]]])
e_ = om_edges.get_poly_edges(ring)
cen = p[t].mean(axis=1)
ins, _ = inpoly2(cen, np.nan_to_num(ring), e_)
t_in = t[ins]; t_out = t[~ins]
print(f"[ex11p] patch elements {ins.sum():,} / {len(t):,}", flush=True)

# --- boundary chain of the patch (ordered loops)
from oceanmesh.clean import _external_topology
bedges, _ = _external_topology(p, t_in)
bedges = np.asarray(bedges, dtype=int)
# order edges into loops
adj = {}
for a, b in bedges:
    adj.setdefault(a, []).append(b)
    adj.setdefault(b, []).append(a)
unused = {tuple(sorted(ed)) for ed in bedges}
loops = []
while unused:
    e0 = unused.pop()
    loop = [e0[0], e0[1]]
    while True:
        cur, prev = loop[-1], loop[-2]
        nxts = [n for n in adj[cur]
                if tuple(sorted((cur, n))) in unused]
        if not nxts:
            break
        nxt = nxts[0]
        unused.discard(tuple(sorted((cur, nxt))))
        loop.append(nxt)
        if nxt == loop[0]:
            break
    loops.append(loop)
loops.sort(key=len, reverse=True)
print(f"[ex11p] boundary loops: {[len(l) for l in loops[:4]]}",
      flush=True)
main = np.asarray(loops[0])
if main[0] == main[-1]:
    main = main[:-1]
pfix = p[main]
n = len(pfix)
egfix = np.column_stack([np.arange(n), (np.arange(n)+1) % n])

# --- patch domain: parity vs the patch's own boundary polygon
bpoly = np.vstack([pfix, pfix[0], [[np.nan, np.nan]]])
be = om_edges.get_poly_edges(bpoly)
tree = scipy.spatial.cKDTree(pfix)

def _fd(q):
    q = np.asarray(q, dtype=float)
    d, _ = tree.query(q, k=1, workers=-1)
    ins_, _ = inpoly2(q, np.nan_to_num(bpoly), be)
    return np.where(ins_, -d, d)

bb = (float(pfix[:,0].min()), float(pfix[:,0].max()),
      float(pfix[:,1].min()), float(pfix[:,1].max()))
dom = Domain(bb, _fd)
dom.crs = "EPSG:4326"
dom.covering = _fd
dom.boubox_ring = pfix

with open(OUT/"fh.pkl", "rb") as fh_:
    gd = pickle.load(fh_)
g = Grid(bbox=gd["bbox"], dx=gd["dx"], dy=gd["dy"],
         values=gd["values"], extrapolate=True, crs="EPSG:4326")
g.hmin = 1e3/111e3
g.build_interpolant()

pp, tt = om.generate_mesh(dom, g, max_iter=50, seed=0,
                          pfix=pfix, egfix=egfix, cleanup="none")
print(f"[ex11p] patch NP={len(pp):,} NT={len(tt):,} "
      f"+{time.time()-t0:.0f}s", flush=True)

# --- stitch: m_w_hole + patch, merge coincident nodes
from oceanmesh.fix_mesh import fix_mesh
keep_idx = np.unique(t_out.reshape(-1))
remap = -np.ones(len(p), dtype=int)
remap[keep_idx] = np.arange(len(keep_idx))
p_out = p[keep_idx]; t_out2 = remap[t_out]
P = np.vstack([p_out, pp])
T = np.vstack([t_out2, tt + len(p_out)])
# weld nodes within 1e-8 deg
tr_ = scipy.spatial.cKDTree(P)
pairs = tr_.query_pairs(1e-8)
parent = np.arange(len(P))
def find(a):
    while parent[a] != a:
        parent[a] = parent[parent[a]]; a = parent[a]
    return a
for a, b in pairs:
    ra, rb = find(a), find(b)
    if ra != rb:
        parent[max(ra, rb)] = min(ra, rb)
root = np.array([find(i) for i in range(len(P))])
uniq, inv = np.unique(root, return_inverse=True)
P2 = P[uniq]; T2 = inv[root[T.reshape(-1)]].reshape(-1, 3)
P2, T2, _ = fix_mesh(P2, T2, dim=2, delete_unused=True)
print(f"[ex11p] stitched NP={len(P2):,} NT={len(T2):,}", flush=True)
np.save(OUT/"p_new.npy", P2); np.save(OUT/"t_new.npy", T2)
from oceanmesh.fix_mesh import simp_qual
from pyproj import Transformer
trx = Transformer.from_crs("EPSG:4326",
    "+proj=tmerc +lon_0=171 +lat_0=-44 +ellps=WGS84 +units=m",
    always_xy=True)
xx, yy = trx.transform(P2[:, 0], P2[:, 1])
q = simp_qual(np.column_stack([xx, yy]), T2)
print(f"[ex11p] qual min={q.min():.4f} mean={q.mean():.4f}", flush=True)
