"""Stage-2 executor: apply one operator per planned patch with a
cross-halo acceptance check; revert the patch on failure (per-patch
atomicity). Single pass — no re-detection loop (design rule)."""

from __future__ import annotations

from collections import defaultdict

import numpy as np

from .detect import THRESH, _angles


def _local_ok(P, T, ring, e2t, sign0, area_of, th):
    ring = np.asarray(sorted(ring))
    a, b, c = P[T[ring, 0]], P[T[ring, 1]], P[T[ring, 2]]
    ar = 0.5 * np.cross(b - a, c - a)
    if (np.sign(ar) != sign0[ring]).any() or (ar == 0).any():
        return False, np.inf, 0.0
    aa = np.abs(ar)
    idx = {int(e): i for i, e in enumerate(ring)}
    worst = 1.0
    for e in ring:
        for k in range(3):
            ed = tuple(sorted((int(T[e, k]), int(T[e, (k + 1) % 3]))))
            for o in e2t[ed]:
                if o == e:
                    continue
                a_o = aa[idx[o]] if o in idx else area_of(o)
                r = max(aa[idx[e]], a_o) / max(
                    1e-30, min(aa[idx[e]], a_o))
                worst = max(worst, r)
    ang = _angles(P, T[ring])
    return True, worst, float(ang.min())


def execute_patches(points, cells, patches, obc_nodes=None,
                    log=print, n_trials=400):
    """Apply operators patch-by-patch. Returns (points, ledger):
    points possibly modified (topology never changes in micro/bound
    ops), ledger = per-patch outcome records."""
    P = np.asarray(points, float).copy()
    T = np.asarray(cells, int)
    th = dict(THRESH)
    obc = set(int(v) for v in (obc_nodes or []))

    e2t = defaultdict(list)
    for i, (a, b, c) in enumerate(T):
        for e in ((a, b), (b, c), (c, a)):
            e2t[tuple(sorted(e))].append(i)
    bnd_nb = defaultdict(list)
    bnd_nodes = set()
    for e, ts in e2t.items():
        if len(ts) == 1:
            bnd_nodes.update(e)
            bnd_nb[e[0]].append(e[1])
            bnd_nb[e[1]].append(e[0])
    n2t = defaultdict(list)
    for i, tri in enumerate(T):
        for v in tri:
            n2t[int(v)].append(i)

    a, b, c = P[T[:, 0]], P[T[:, 1]], P[T[:, 2]]
    sign0 = np.sign(0.5 * np.cross(b - a, c - a))

    def area_of(e):
        aa = P[T[e]]
        return abs(0.5 * np.cross(aa[1] - aa[0], aa[2] - aa[0]))

    ledger = []
    for q in patches:
        if "note" in q:
            ledger.append(q)
            continue
        rec = dict(q)
        rec.pop("elements", None)
        if q["class"] == "obc-locked":
            rec["outcome"] = "skipped (obc-locked)"
            ledger.append(rec)
            continue
        seeds = q["seed_elements"]
        patch_nodes = {int(v) for ie in seeds for v in T[ie]}
        ring = set(q.get("elements") or
                   {ie for v in patch_nodes for ie in n2t[v]})
        free = [v for v in patch_nodes
                if v not in bnd_nodes and v not in obc]
        slide = [v for v in patch_nodes
                 if v in bnd_nodes and v not in obc
                 and len(bnd_nb[v]) == 2]
        ok0, r0, a0 = _local_ok(P, T, ring, e2t, sign0, area_of, th)
        rng = np.random.default_rng(42 + q["id"])
        best = None
        h = float(np.sqrt(sum(area_of(e) for e in seeds)
                          / max(len(seeds), 1)))

        def accept(P_try):
            ok, r1, a1 = _local_ok(P_try, T, ring, e2t, sign0,
                                   area_of, th)
            return (ok and r1 < th["c4_ratio"]
                    and a1 >= min(a0, th["c1_min_deg"]) - 1e-9
                    and r1 <= r0 + 1e-9), r1, a1

        # 1) deterministic centroid smooth (interior nodes only)
        if free:
            P_try = P.copy()
            for v in free:
                elems = n2t[v]
                cent = P[T[elems]].mean(axis=1)
                w = np.array([area_of(e) for e in elems])
                P_try[v] = (cent * w[:, None]).sum(0) / w.sum()
            good, r1, a1 = accept(P_try)
            if good:
                best = (P_try, r1, a1, "deterministic")
        # 2) stochastic (free interior) / boundary slide
        if best is None and (free or slide):
            for _ in range(n_trials):
                P_try = P.copy()
                # free interior nodes and boundary-slide nodes are
                # perturbed TOGETHER: 'bound' patches often need
                # both (interior-only misses coast-pair C4 sites)
                for v in free:
                    P_try[v] = P[v] + rng.normal(0, 0.10 * h, 2)
                for v in slide:
                    n1, n2 = bnd_nb[v]
                    tpar = rng.uniform(-0.45, 0.45)
                    tgt = n1 if tpar < 0 else n2
                    P_try[v] = P[v] + abs(tpar) * (P[tgt] - P[v])
                good, r1, a1 = accept(P_try)
                if good and (best is None or r1 < best[1]):
                    best = (P_try, r1, a1,
                            "stochastic+slide" if (free and slide)
                            else ("stochastic" if free
                                  else "boundary-slide"))
        if best is None:
            rec["outcome"] = f"unfixed (r0={r0:.2f}, a0={a0:.1f})"
        else:
            P, r1, a1, how = best
            rec["outcome"] = (f"fixed [{how}] r {r0:.2f}->{r1:.2f} "
                              f"minang {a0:.1f}->{a1:.1f}")
        log(f"[finishing] patch {q['id']} ({q['class']}): "
            f"{rec['outcome']}")
        ledger.append(rec)
    return P, ledger
