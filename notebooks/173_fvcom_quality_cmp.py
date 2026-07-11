# FVCOM mesh-quality criteria (kickoff sec 7.1) applied to BOTH:
# improved oceanmesh vs OM2D golden, in the same projected frame.
#   C1 min interior angle >= 30 deg
#   C2 max interior angle <= 130 deg
#   C3 adjacent element area change (max-min)/max <= 0.5
#   C4 connecting elements per node <= 8 (interior nodes)
import os, sys
import numpy as np
os.environ.setdefault("MPLBACKEND", "Agg")
sys.path.insert(0, os.path.expanduser("~/Github/oceanmesh"))
from pathlib import Path
from pyproj import Transformer
from fvcom_mesh_tools.io import read_fort14

OM2D = Path(os.path.expanduser("~/Github/OceanMesh2D"))
OUT = Path("outputs/om2d_examples/ex1_nz")
tr = Transformer.from_crs(
    "EPSG:4326",
    "+proj=tmerc +lon_0=171 +lat_0=-44 +units=m +ellps=WGS84",
    always_xy=True,
)

def metrics(p, t, name):
    x, y = tr.transform(p[:, 0], p[:, 1])
    q = np.column_stack([x, y])
    a, b, c = q[t[:, 0]], q[t[:, 1]], q[t[:, 2]]
    def ang(u, v, w):
        d1 = v - u; d2 = w - u
        cosA = (d1 * d2).sum(1) / (
            np.linalg.norm(d1, axis=1) * np.linalg.norm(d2, axis=1)
        )
        return np.degrees(np.arccos(np.clip(cosA, -1, 1)))
    A0, A1, A2 = ang(a, b, c), ang(b, c, a), ang(c, a, b)
    mn = np.minimum(np.minimum(A0, A1), A2)
    mx = np.maximum(np.maximum(A0, A1), A2)
    area = 0.5 * np.abs(np.cross(b - a, c - a))
    # internal edge pairs
    e = np.vstack([t[:, [0, 1]], t[:, [1, 2]], t[:, [2, 0]]])
    eid = np.repeat(np.arange(len(t)), 3).reshape(3, -1).T.ravel()
    eid = np.tile(np.arange(len(t)), 3)
    e = np.vstack([t[:, [0, 1]], t[:, [1, 2]], t[:, [2, 0]]])
    key = np.sort(e, axis=1)
    order = np.lexsort((key[:, 1], key[:, 0]))
    key, eid = key[order], eid[order]
    same = (key[1:] == key[:-1]).all(1)
    pair = np.column_stack([eid[:-1][same], eid[1:][same]])
    ai, aj = area[pair[:, 0]], area[pair[:, 1]]
    big = np.maximum(ai, aj)
    ac = (big - np.minimum(ai, aj)) / np.maximum(big, 1e-30)
    # valence (interior nodes only)
    val = np.bincount(t.ravel(), minlength=len(p))
    bnd_e = key[:-1][~same] if len(same) else key
    # boundary edges appear once
    once = np.ones(len(key), bool)
    once[:-1][same] = False
    once[1:][same] = False
    bnd_nodes = np.unique(key[once])
    interior = np.ones(len(p), bool)
    interior[bnd_nodes] = False
    vi = val[interior]
    print(f"[fvcom] {name}")
    print(f"[fvcom]   NP={len(p):,} NT={len(t):,}")
    print(f"[fvcom]   C1 min angle : worst {mn.min():6.2f} deg  "
          f"violations(<30) {int((mn < 30).sum()):4d} "
          f"({100*(mn<30).mean():.2f}%)")
    print(f"[fvcom]   C2 max angle : worst {mx.max():6.2f} deg  "
          f"violations(>130) {int((mx > 130).sum()):4d} "
          f"({100*(mx>130).mean():.2f}%)")
    print(f"[fvcom]   C3 area chg  : worst {ac.max():6.3f}      "
          f"violations(>0.5) {int((ac > 0.5).sum()):4d} "
          f"({100*(ac>0.5).mean():.2f}%)")
    print(f"[fvcom]   C4 valence   : worst {int(vi.max()):2d}         "
          f"violations(>8) {int((vi > 8).sum()):4d}", flush=True)

p = np.load(OUT / "p.npy"); t = np.load(OUT / "t.npy")
metrics(p, t, "improved oceanmesh (Python)")
g = read_fort14(str(OM2D / "Examples/1_NZ_mesh.14"))
metrics(g.nodes[:, :2], np.asarray(g.elements, int), "OceanMesh2D golden (MATLAB)")
