# Careful A/B: numba kernel vs Cython kernel for inpoly2.
# Correctness on real + synthetic data (stat AND bnd outputs,
# incl. on-vertex and near-edge points, brute-force referee) and
# repeated timings at several scales.
import os, sys, time
import numpy as np
os.environ.setdefault("MPLBACKEND", "Agg")
sys.path.insert(0, os.path.expanduser("~/Github/oceanmesh"))
from pathlib import Path
import oceanmesh  # noqa: F401  (registers geometry module)
from oceanmesh.geometry import point_in_polygon as pip
from oceanmesh import edges as om_edges
from oceanmesh import Shoreline

DEG = 1.0 / 111e3
OM2D = Path(os.path.expanduser("~/Github/OceanMesh2D"))
rng = np.random.default_rng(11)

def brute(qpts, node, edge):
    x1 = node[edge[:, 0], 0]; y1 = node[edge[:, 0], 1]
    x2 = node[edge[:, 1], 0]; y2 = node[edge[:, 1], 1]
    out = np.zeros(len(qpts), bool)
    for i, (qx, qy) in enumerate(qpts):
        cond = (y1 > qy) != (y2 > qy)
        with np.errstate(divide="ignore", invalid="ignore"):
            xint = x1 + (qy - y1) * (x2 - x1) / (y2 - y1)
        out[i] = (np.count_nonzero(cond & (qx < xint)) % 2) == 1
    return out

def run_case(name, node, edge, bbox):
    n_q = 500_000
    q = np.column_stack([
        rng.uniform(bbox[0], bbox[1], n_q),
        rng.uniform(bbox[2], bbox[3], n_q),
    ])
    # adversarial points: exact vertices + near-edge offsets
    used = np.unique(edge)
    vsel = rng.choice(used, min(20_000, len(used)), replace=False)
    q_on = node[vsel]
    mids = 0.5 * (node[edge[:, 0]] + node[edge[:, 1]])
    msel = rng.choice(len(mids), min(20_000, len(mids)), replace=False)
    q_near = mids[msel] + rng.normal(0, 1e-9, (len(msel), 2))
    Q = np.vstack([q, q_on, q_near])

    res = {}
    for tag, fn in [("cython", pip.inpoly2_fast_wrapper
                     if hasattr(pip, "inpoly2_fast_wrapper") else None),
                    ]:
        pass
    # call kernels directly
    from oceanmesh.geometry.point_in_polygon import _inpoly_numba
    from oceanmesh.geometry.point_in_polygon_ import inpoly2_fast
    ts = {}
    outs = {}
    for tag, call in [
        ("numba", lambda: _inpoly_numba(Q, node, edge, 5e-14)),
        ("cython", lambda: inpoly2_fast(np.ascontiguousarray(Q), np.ascontiguousarray(node), np.ascontiguousarray(edge, dtype=np.int32), 5e-14)),
    ]:
        best = np.inf
        for _rep in range(3):
            t0 = time.time(); s, b = call(); dt = time.time() - t0
            best = min(best, dt)
        ts[tag] = best
        outs[tag] = (np.asarray(s, bool), np.asarray(b, bool))
    s_n, b_n = outs["numba"]; s_c, b_c = outs["cython"]
    mism_s = int((s_n != s_c).sum())
    mism_b = int((b_n != b_c).sum())
    # referee on stat disagreements + random sample (excl. bnd pts,
    # where crossing parity is ill-defined)
    ref_txt = "n/a"
    dis = np.where(s_n != s_c)[0]
    on_b = b_n | b_c
    dis = dis[~on_b[dis]][:150]
    sel = rng.choice(len(Q), 150, replace=False)
    sel = sel[~on_b[sel]]
    chk = np.unique(np.concatenate([dis, sel]))
    ref = brute(Q[chk], node, edge)
    agree_n = float((ref == s_n[chk]).mean())
    agree_c = float((ref == s_c[chk]).mean())
    print(f"[ab] {name}: nodes={len(node):,} edges={len(edge):,} "
          f"queries={len(Q):,}", flush=True)
    print(f"[ab]   stat mismatches: {mism_s}   "
          f"bnd mismatches: {mism_b} "
          f"(of which on-boundary pts: {int(on_b[(s_n!=s_c)].sum())})",
          flush=True)
    print(f"[ab]   referee agreement: numba {100*agree_n:.1f}% | "
          f"cython {100*agree_c:.1f}% (n={len(chk)})", flush=True)
    print(f"[ab]   time: numba {ts['numba']:.3f}s | "
          f"cython {ts['cython']:.3f}s | ratio "
          f"{ts['cython']/max(ts['numba'],1e-9):.0f}x", flush=True)

# case 1: NZ GSHHS shoreline (ladder step 1 data)
sh = Shoreline(str(OM2D / "datasets/GSHHS_shp/f/GSHHS_f_L1.shp"),
               (166.0, 176.0, -48.0, -40.0), 1e3 * DEG)
poly = np.vstack((sh.inner, sh.boubox))
run_case("NZ GSHHS h0=1km", np.nan_to_num(poly),
         om_edges.get_poly_edges(poly), (166, 176, -48, -40))

# case 2: Tokyo Bay OSM true-land (original incident data)
sh2 = Shoreline("outputs/tb_varres_3r/land_osm_wide.shp",
                np.column_stack([
                    [139.14, 139.14, 140.39, 140.39, 139.14],
                    [34.75, 35.86, 35.86, 34.75, 34.75]]),
                1e2 * DEG)
poly2 = np.vstack((sh2.inner, sh2.boubox))
run_case("TokyoBay OSM h0=100m", np.nan_to_num(poly2),
         om_edges.get_poly_edges(poly2), (139.14, 140.39, 34.75, 35.86))

# case 3: synthetic multi-ring + open chains + duplicates
segs = []
th = np.linspace(0, 2 * np.pi, 4001)
ring = np.column_stack([np.cos(th), np.sin(th)])
ring = np.insert(ring, 7, ring[7], axis=0)
segs.append(ring)
for cx, cy, r, n in ((0.3, 0.2, 0.2, 501), (-0.4, -0.3, 0.15, 301)):
    t2 = np.linspace(0, 2 * np.pi, n)
    segs.append(np.column_stack([cx + r * np.cos(t2),
                                 cy + r * np.sin(t2)]))
segs.append(np.column_stack([np.linspace(-0.9, 0.9, 800),
                             0.7 + 0.02 * rng.standard_normal(800)]))
parts = []
for s in segs:
    parts.append(s); parts.append(np.array([[np.nan, np.nan]]))
poly3 = np.vstack(parts)
run_case("synthetic multiring", np.nan_to_num(poly3),
         om_edges.get_poly_edges(poly3), (-1.2, 1.2, -1.2, 1.2))
