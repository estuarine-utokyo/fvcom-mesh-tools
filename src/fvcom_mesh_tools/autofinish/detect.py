"""Stage-2 detectors: QA-class violations as (check, element/node
sites). Pure analysis; thresholds mirror the QA gates."""

from __future__ import annotations

from collections import defaultdict

import numpy as np

THRESH = {
    "c1_min_deg": 30.0,
    "c2_max_deg": 130.0,
    "c4_ratio": 2.0,
    "c5_valence": 8,
}


def _angles(P, T):
    out = np.empty((len(T), 3))
    for k in range(3):
        u = P[T[:, (k + 1) % 3]] - P[T[:, k]]
        v = P[T[:, (k + 2) % 3]] - P[T[:, k]]
        c = (u * v).sum(1) / np.maximum(
            np.linalg.norm(u, axis=1) * np.linalg.norm(v, axis=1),
            1e-30,
        )
        out[:, k] = np.degrees(np.arccos(np.clip(c, -1, 1)))
    return out


def detect_violations(points, cells, thresholds=None):
    """Returns dict check_id -> {"elements": [...]} plus shared
    topology info reused by the planner."""
    th = dict(THRESH)
    th.update(thresholds or {})
    P = np.asarray(points, float)
    T = np.asarray(cells, int)

    e2t = defaultdict(list)
    for i, (a, b, c) in enumerate(T):
        for e in ((a, b), (b, c), (c, a)):
            e2t[tuple(sorted(e))].append(i)
    bnd_nodes = set()
    bnd_valence = defaultdict(int)
    for e, ts in e2t.items():
        if len(ts) == 1:
            bnd_nodes.update(e)
            bnd_valence[e[0]] += 1
            bnd_valence[e[1]] += 1

    ang = _angles(P, T)
    a, b, c = P[T[:, 0]], P[T[:, 1]], P[T[:, 2]]
    area = 0.5 * np.abs(np.cross(b - a, c - a))

    out = {}
    out["c1"] = {"elements": np.where(
        ang.min(axis=1) < th["c1_min_deg"])[0].tolist()}
    out["c2"] = {"elements": np.where(
        ang.max(axis=1) > th["c2_max_deg"])[0].tolist()}
    c4 = set()
    for e, ts in e2t.items():
        if len(ts) == 2:
            r = max(area[ts[0]], area[ts[1]]) / max(
                1e-30, min(area[ts[0]], area[ts[1]]))
            if r > th["c4_ratio"]:
                c4.update(ts)
    out["c4"] = {"elements": sorted(c4)}
    val = np.bincount(T.ravel(), minlength=len(P))
    bad_v = [int(v) for v in np.where(val > th["c5_valence"])[0]
             if v not in bnd_nodes]
    out["c5"] = {"elements": sorted(
        {ie for v in bad_v for ie in np.where((T == v).any(1))[0]})}
    # pinch: boundary node incident to >2 boundary edges
    pinch = [v for v, k in bnd_valence.items() if k > 2]
    out["pinch"] = {"elements": sorted(
        {ie for v in pinch for ie in np.where((T == v).any(1))[0]})}
    out["_topo"] = {"e2t": e2t, "bnd_nodes": bnd_nodes,
                    "area": area, "angles": ang}
    return out
