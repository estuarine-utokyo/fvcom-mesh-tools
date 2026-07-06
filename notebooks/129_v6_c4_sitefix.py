# PoC #129: single-pass per-site C4 fixer (user-approved option b).
# For each adjacent-element pair with area ratio > 2.0: free the
# INTERIOR nodes of the two-triangle patch (boundary nodes never
# move — whitelist), try a deterministic area-weighted centroid
# smooth first, then seed=42 stochastic perturbation with a 1-ring
# acceptance check (no flips, local C4 reduced, local min angle not
# degraded). One pass only; unfixed sites are reported.
import os
from collections import defaultdict
from pathlib import Path

import numpy as np

os.environ.setdefault("MPLBACKEND", "Agg")

from fvcom_mesh_tools.io import read_fort14, write_fort14  # noqa: E402

REPO = Path(__file__).resolve().parents[1]
F14 = REPO / "outputs" / "pipeline_v6r" / "tokyo_bay_v6_final.14"
THRESH = 2.0
rng = np.random.default_rng(42)

mesh = read_fort14(F14)
P = mesh.nodes.copy()
T = mesh.elements


def tri_areas(P, tris):
    a, b, c = P[tris[:, 0]], P[tris[:, 1]], P[tris[:, 2]]
    return 0.5 * np.cross(b - a, c - a)


def min_angles(P, tris):
    out = np.full(len(tris), 180.0)
    for k in range(3):
        u = P[tris[:, (k + 1) % 3]] - P[tris[:, k]]
        v = P[tris[:, (k + 2) % 3]] - P[tris[:, k]]
        cosv = (u * v).sum(1) / (
            np.linalg.norm(u, axis=1) * np.linalg.norm(v, axis=1)
        )
        ang = np.degrees(np.arccos(np.clip(cosv, -1, 1)))
        out = np.minimum(out, ang)
    return out


e2t = defaultdict(list)
for i, (a, b, c) in enumerate(T):
    for e in ((a, b), (b, c), (c, a)):
        e2t[tuple(sorted(e))].append(i)
bnd_nodes = set()
for e, ts in e2t.items():
    if len(ts) == 1:
        bnd_nodes.update(e)

areas = np.abs(tri_areas(P, T))
pairs = [
    (ts[0], ts[1])
    for e, ts in e2t.items()
    if len(ts) == 2
    and max(areas[ts[0]], areas[ts[1]])
    / min(areas[ts[0]], areas[ts[1]]) > THRESH
]
print(f"[fix] C4 sites: {len(pairs)}")

n2t = defaultdict(list)
for i, tri in enumerate(T):
    for v in tri:
        n2t[int(v)].append(i)


def local_metrics(P, ring):
    ring = np.asarray(sorted(ring))
    ar = tri_areas(P, T[ring])
    ok = (np.sign(ar) == sign0[ring]).all()
    aa = np.abs(ar)
    idx = {e: i for i, e in enumerate(ring)}
    worst = 1.0
    for e in ring:
        for ed in ((T[e, 0], T[e, 1]), (T[e, 1], T[e, 2]),
                   (T[e, 2], T[e, 0])):
            for o in e2t[tuple(sorted(ed))]:
                if o == e:
                    continue
                # cross-boundary pairs too: an outside neighbour's
                # area is unchanged but the RATIO changes
                a_o = (aa[idx[o]] if o in idx
                       else float(np.abs(tri_areas(P, T[[o]]))[0]))
                a_e = aa[idx[e]]
                r = max(a_e, a_o) / max(1e-30, min(a_e, a_o))
                worst = max(worst, r)
    return ok, worst, min_angles(P, T[ring]).min()


sign0 = np.sign(tri_areas(P, T))
fixed = 0
for ea, eb in pairs:
    patch_nodes = set(int(v) for v in T[ea]) | set(int(v) for v in T[eb])
    free = [v for v in patch_nodes if v not in bnd_nodes]
    ring = set()
    for v in patch_nodes:
        ring.update(n2t[v])
    slide_nodes = []
    if not free:
        # whitelist-compliant boundary op: slide nodes ALONG the
        # boundary polyline (never off it); OBC members excluded
        obc_set = set(
            int(v) for b in mesh.open_boundaries for v in b
        )
        bnd_nb = defaultdict(list)
        for e, ts in e2t.items():
            if len(ts) == 1:
                bnd_nb[e[0]].append(e[1])
                bnd_nb[e[1]].append(e[0])
        slide_nodes = [v for v in patch_nodes
                       if v not in obc_set and len(bnd_nb[v]) == 2]
        if not slide_nodes:
            print(f"[fix] site ({ea},{eb}): no movable node — skipped")
            continue
    _, r0, a0 = local_metrics(P, ring)
    best = None
    # deterministic: area-weighted centroid of each free node's ring
    P_try = P.copy()
    for v in (free or []):
        elems = n2t[v]
        cent = P[T[elems]].mean(axis=1)
        w = np.abs(tri_areas(P, T[elems]))
        P_try[v] = (cent * w[:, None]).sum(0) / w.sum()
    ok, r1, a1 = local_metrics(P_try, ring)
    if ok and r1 < min(r0, THRESH) and a1 >= min(a0, 30.0) - 1e-9:
        best = (P_try, r1, a1, "deterministic")
    if best is None:
        h = float(np.sqrt(np.abs(areas[ea]) + np.abs(areas[eb])))
        for _ in range(400):
            P_try = P.copy()
            if free:
                for v in free:
                    P_try[v] = P[v] + rng.normal(0.0, 0.10 * h, 2)
            else:
                for v in slide_nodes:
                    nb1, nb2 = bnd_nb[v]
                    tpar = rng.uniform(-0.45, 0.45)
                    tgt = nb1 if tpar < 0 else nb2
                    P_try[v] = P[v] + abs(tpar) * (P[tgt] - P[v])
            ok, r1, a1 = local_metrics(P_try, ring)
            if (ok and r1 < THRESH
                    and a1 >= min(a0, 30.0) - 1e-9
                    and (best is None or r1 < best[1])):
                best = (P_try, r1, a1,
                        "stochastic" if free else "boundary-slide")
    if best is None:
        print(f"[fix] site ({ea},{eb}): r0={r0:.2f} NOT fixed "
              "(reported as residual)")
        continue
    P, r1, a1, how = best
    areas = np.abs(tri_areas(P, T))
    fixed += 1
    print(f"[fix] site ({ea},{eb}): {r0:.2f} -> {r1:.2f} "
          f"(minang {a0:.1f} -> {a1:.1f}) [{how}]")

print(f"[fix] fixed {fixed}/{len(pairs)} sites (single pass)")
if fixed:
    import shutil

    shutil.copy(F14, F14.with_suffix(".14.bak129"))
    mesh.nodes = P
    write_fort14(mesh, F14)
    print(f"[fix] wrote {F14} (backup .bak129)")
