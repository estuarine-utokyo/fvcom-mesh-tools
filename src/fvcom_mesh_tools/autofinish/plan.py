"""Stage-2 planner: cluster violation elements into k-ring patches
and classify them (micro / patch / bound). Analysis only."""

from __future__ import annotations

import json
from collections import defaultdict
from pathlib import Path

import numpy as np


def _kring(T, seeds, k, n2t):
    cur = set(seeds)
    for _ in range(k):
        nodes = {int(v) for ie in cur for v in T[ie]}
        cur |= {ie for v in nodes for ie in n2t[v]}
    return cur


def plan_patches(points, cells, detections, kring=2,
                 obc_nodes=None, max_patches=50):
    """Cluster all detected elements -> merged patches with class
    and per-patch metrics. Returns list of patch dicts sorted by
    severity (worst first)."""
    P = np.asarray(points, float)
    T = np.asarray(cells, int)
    topo = detections["_topo"]
    bnd = topo["bnd_nodes"]
    obc = set(int(v) for v in (obc_nodes or []))

    n2t = defaultdict(list)
    for i, tri in enumerate(T):
        for v in tri:
            n2t[int(v)].append(i)

    seed_checks = defaultdict(set)
    for chk, d in detections.items():
        if chk.startswith("_"):
            continue
        for ie in d["elements"]:
            seed_checks[ie].add(chk)
    seeds = sorted(seed_checks)
    if not seeds:
        return []

    # union-find over seeds whose k-rings overlap
    rings = {s: _kring(T, [s], kring, n2t) for s in seeds}
    parent = {s: s for s in seeds}

    def find(x):
        while parent[x] != x:
            parent[x] = parent[parent[x]]
            x = parent[x]
        return x

    elem_owner = {}
    for s in seeds:
        for ie in rings[s]:
            if ie in elem_owner:
                ra, rb = find(elem_owner[ie]), find(s)
                if ra != rb:
                    parent[rb] = ra
            else:
                elem_owner[ie] = s

    groups = defaultdict(list)
    for s in seeds:
        groups[find(s)].append(s)

    patches = []
    for gid, members in groups.items():
        elems = set()
        checks = set()
        for s in members:
            elems |= rings[s]
            checks |= seed_checks[s]
        nodes = {int(v) for ie in elems for v in T[ie]}
        seed_elems = sorted(set(members))
        touches_obc = bool(nodes & obc)
        touches_bnd = bool(nodes & bnd)
        if touches_obc:
            cls = "obc-locked"
        elif touches_bnd:
            cls = "bound"
        elif len(seed_elems) <= 2:
            cls = "micro"
        else:
            cls = "patch"
        ang = topo["angles"]
        sel = np.asarray(sorted(elems))
        cen = P[T[sel]].mean(axis=(0, 1))
        patches.append({
            "id": len(patches),
            "class": cls,
            "checks": sorted(checks),
            "seed_elements": seed_elems,
            "n_elements": len(elems),
            "elements": sorted(int(e) for e in elems),
            "centroid_xy": [float(cen[0]), float(cen[1])],
            "worst_min_angle": float(ang[sel].min()),
            "severity": float(
                max(0.0, 30.0 - ang[sel].min()) + 10 * len(checks)
            ),
        })
    patches.sort(key=lambda q: -q["severity"])
    for i, q in enumerate(patches):
        q["id"] = i
    if len(patches) > max_patches:
        dropped = len(patches) - max_patches
        patches = patches[:max_patches]
        patches.append({"note": f"{dropped} patches dropped by "
                        f"max_patches={max_patches} cap"})
    return patches


def write_ledger(patches, out_path):
    out_path = Path(out_path)
    slim = [{k: v for k, v in q.items() if k != "elements"}
            for q in patches]
    out_path.write_text(json.dumps(slim, indent=2))
    return out_path
