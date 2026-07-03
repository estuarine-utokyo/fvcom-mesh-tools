"""Quality-gated LOCAL open-boundary perpendicularity fixer.

Promoted from PoC #59e. The global first-ring aligner
(:func:`align_open_boundary_first_ring`) moves the ring of EVERY OBC
node, so alternating it with a quality finisher (``phase_h_finish``)
ping-pongs: each pass undoes the other's work on already-passing
nodes. This fixer instead moves only the first-ring interior partners
of the *violating* nodes, one node at a time, accepting a move iff

* the target OBC node's best-edge deviation from perpendicular drops
  to ``dev_target`` or below;
* every element of the moved node's 1-ring keeps ``min_angle`` /
  ``max_angle`` and a positive area, and every internal edge touching
  the ring keeps the adjacent-area-change bound;
* no OBC node adjacent to the moved node newly violates ``dev_max``
  (already-violating neighbours must not get worse).

Candidate positions walk damped steps toward the exact perpendicular
target (edge length preserved) plus seeded Gaussian jitter, so the
result is bit-reproducible. On the Tokyo-Bay #59d residual this
fixed 8/8 violating nodes in a single pass without breaking any
C1/C2/C4/C5 gate.

The module is deliberately self-contained (numpy + Fort14Mesh only)
to stay import-cycle-free from the QA and Phase-H layers.
"""

from __future__ import annotations

from typing import Any

import numpy as np

from fvcom_mesh_tools.algorithms.perpendicularity import boundary_tangents
from fvcom_mesh_tools.io import Fort14Mesh

DEFAULT_DEV_MAX: float = 20.0
DEFAULT_DEV_TARGET: float = 19.0
_MIN_ANG_EPS = 1e-9


def _edge_arrays(elements: np.ndarray, n_nodes: int):
    """Unique undirected edges, per-node neighbour lists, boundary-node
    mask, node->elements map, and element->edge-neighbour map."""
    raw = np.vstack([
        elements[:, [0, 1]], elements[:, [1, 2]], elements[:, [2, 0]],
    ])
    raw.sort(axis=1)
    codes = raw[:, 0].astype(np.int64) * n_nodes + raw[:, 1]
    elem_of = np.tile(np.arange(elements.shape[0], dtype=np.int64), 3)
    order = np.argsort(codes, kind="stable")
    codes_s, elem_s = codes[order], elem_of[order]
    uniq, first_idx, counts = np.unique(
        codes_s, return_index=True, return_counts=True,
    )
    uv = np.column_stack([uniq // n_nodes, uniq % n_nodes]).astype(np.int64)

    boundary_nodes = np.zeros(n_nodes, dtype=bool)
    buv = uv[counts == 1]
    if buv.size:
        boundary_nodes[buv.ravel()] = True

    adj: dict[int, list[int]] = {}
    for u, v in uv:
        adj.setdefault(int(u), []).append(int(v))
        adj.setdefault(int(v), []).append(int(u))

    rows = elements.ravel()
    cols = np.repeat(np.arange(elements.shape[0]), 3)
    o = np.argsort(rows, kind="stable")
    rows, cols = rows[o], cols[o]
    starts = np.searchsorted(rows, np.arange(n_nodes + 1))
    n2e = {v: cols[starts[v]:starts[v + 1]] for v in range(n_nodes)
           if starts[v] < starts[v + 1]}

    e2nbr: dict[int, list[int]] = {}
    internal_mask = counts == 2
    fi = first_idx[internal_mask]
    for a, b in zip(elem_s[fi], elem_s[fi + 1]):
        e2nbr.setdefault(int(a), []).append(int(b))
        e2nbr.setdefault(int(b), []).append(int(a))
    return adj, boundary_nodes, n2e, e2nbr


def _tri_quality(pts: np.ndarray):
    """Vectorised (min_angle_deg, max_angle_deg, twice_signed_area) for
    ``pts`` shaped (R, 3, 2)."""
    p0, p1, p2 = pts[:, 0], pts[:, 1], pts[:, 2]
    twice = (
        (p1[:, 0] - p0[:, 0]) * (p2[:, 1] - p0[:, 1])
        - (p1[:, 1] - p0[:, 1]) * (p2[:, 0] - p0[:, 0])
    )
    e0 = np.linalg.norm(p1 - p0, axis=1)
    e1 = np.linalg.norm(p2 - p1, axis=1)
    e2 = np.linalg.norm(p0 - p2, axis=1)

    def _ang(opp, ea, eb):
        denom = 2.0 * ea * eb
        cos = np.where(denom > 0, (ea**2 + eb**2 - opp**2)
                       / np.where(denom > 0, denom, 1.0), 1.0)
        return np.degrees(np.arccos(np.clip(cos, -1.0, 1.0)))

    a0, a1, a2 = _ang(e1, e2, e0), _ang(e2, e0, e1), _ang(e0, e1, e2)
    return (
        np.minimum(np.minimum(a0, a1), a2),
        np.maximum(np.maximum(a0, a1), a2),
        twice,
    )


def _dev_of(nodes, v, tangent, nbrs, in_seg):
    """Best-edge deviation (deg) of OBC node ``v``; inf without an
    incident edge leaving the segment."""
    best = np.inf
    for i in nbrs:
        if in_seg[i]:
            continue
        vec = nodes[i] - nodes[v]
        nrm = np.hypot(vec[0], vec[1])
        if nrm == 0:
            continue
        cosang = abs(float(vec @ tangent)) / nrm
        best = min(best, 90.0 - np.degrees(np.arccos(np.clip(cosang, 0.0, 1.0))))
    return best


def _local_ok(nodes, elements, ring, e2nbr, areas_all, i, p_new,
              *, min_angle, max_angle, max_area_change):
    tri = elements[ring]
    pts = nodes[tri].copy()
    pts[tri == i] = p_new
    min_ang, max_ang, twice = _tri_quality(pts)
    if (
        (twice <= 0).any()
        or min_ang.min() < min_angle - _MIN_ANG_EPS
        or max_ang.max() > max_angle + _MIN_ANG_EPS
    ):
        return False, None
    new_areas = 0.5 * np.abs(twice)
    ring_pos = {int(e): k for k, e in enumerate(ring)}
    for k, e in enumerate(ring):
        a_e = new_areas[k]
        for nb in e2nbr.get(int(e), ()):
            a_nb = new_areas[ring_pos[nb]] if nb in ring_pos else areas_all[nb]
            hi = max(a_e, a_nb)
            if hi > 0 and (hi - min(a_e, a_nb)) / hi > max_area_change + 1e-12:
                return False, None
    return True, new_areas


def align_open_boundary_local(
    mesh: Fort14Mesh,
    *,
    dev_max: float = DEFAULT_DEV_MAX,
    dev_target: float = DEFAULT_DEV_TARGET,
    min_angle: float = 30.0,
    max_angle: float = 130.0,
    max_area_change: float = 0.5,
    seed: int = 4242,
    max_outer: int = 6,
    w_steps: tuple[float, ...] = (1.0, 0.8, 0.6, 0.4, 0.2),
    n_jitter: int = 40,
    jitter_sigma_frac: float = 0.05,
) -> tuple[Fort14Mesh, dict[str, Any]]:
    """Fix per-node open-boundary perpendicularity locally, without
    regressing the quality gates. Returns ``(new_mesh, info)``;
    ``info["remaining"]`` lists OBC nodes still above ``dev_max``
    (empty on full success).
    """
    nodes = mesh.nodes.copy()
    elements = mesh.elements
    n_nodes = mesh.n_nodes
    rng = np.random.default_rng(seed)
    passes: list[dict[str, int]] = []
    unresolved: list[int] = []

    for _outer in range(max_outer):
        adj, boundary_nodes, n2e, e2nbr = _edge_arrays(elements, n_nodes)
        p0_, p1_, p2_ = (nodes[elements[:, k]] for k in range(3))
        areas_all = 0.5 * np.abs(
            (p1_[:, 0] - p0_[:, 0]) * (p2_[:, 1] - p0_[:, 1])
            - (p1_[:, 1] - p0_[:, 1]) * (p2_[:, 0] - p0_[:, 0])
        )
        n_viol_total = 0
        n_accept = 0
        unresolved = []
        for seg in mesh.open_boundaries:
            seg = np.asarray(seg, dtype=np.int64)
            if seg.size < 2:
                continue
            in_seg = np.zeros(n_nodes, dtype=bool)
            in_seg[seg] = True
            tang = boundary_tangents(nodes[seg])
            tangent_of = {int(v): tang[k] for k, v in enumerate(seg)}
            dev_before = {
                int(v): _dev_of(nodes, int(v), tangent_of[int(v)],
                                adj.get(int(v), ()), in_seg)
                for v in seg
            }
            viol = [v for v, d in dev_before.items() if d > dev_max]
            n_viol_total += len(viol)
            for v in viol:
                t_v = tangent_of[v]
                n_hat = np.array([-t_v[1], t_v[0]])
                cands = [i for i in adj.get(v, ())
                         if not in_seg[i] and not boundary_nodes[i]]

                def _dev_edge(i, v=v, t_v=t_v):
                    vec = nodes[i] - nodes[v]
                    nrm = max(np.hypot(vec[0], vec[1]), 1e-30)
                    cosang = abs(float(vec @ t_v)) / nrm
                    return 90.0 - np.degrees(
                        np.arccos(np.clip(cosang, 0.0, 1.0)))

                cands.sort(key=_dev_edge)
                accepted = False
                for i in cands:
                    edge = nodes[i] - nodes[v]
                    length = np.hypot(edge[0], edge[1])
                    if length == 0:
                        continue
                    side = np.sign(float(edge @ n_hat)) or 1.0
                    target = nodes[v] + length * side * n_hat
                    ring = n2e[i]
                    watch = [int(w) for w in adj.get(i, ())
                             if in_seg[w] and int(w) != v]
                    p_old = nodes[i].copy()
                    trials = [(w, np.zeros(2)) for w in w_steps]
                    trials += [
                        (w, rng.normal(0.0, jitter_sigma_frac * length, 2))
                        for w in w_steps
                        for _ in range(max(1, n_jitter // len(w_steps)))
                    ]
                    for w, jit in trials:
                        p_new = (1.0 - w) * p_old + w * target + jit
                        nodes[i] = p_new
                        dev_v = _dev_of(nodes, v, t_v, adj.get(v, ()), in_seg)
                        ok = dev_v <= dev_target
                        if ok:
                            for wn in watch:
                                d_new = _dev_of(nodes, wn, tangent_of[wn],
                                                adj.get(wn, ()), in_seg)
                                d_old = dev_before[wn]
                                if (d_old <= dev_max and d_new > dev_max) or (
                                    d_old > dev_max and d_new > d_old + 1e-9
                                ):
                                    ok = False
                                    break
                        nodes[i] = p_old
                        if not ok:
                            continue
                        ok_q, new_areas = _local_ok(
                            nodes, elements, ring, e2nbr, areas_all, i, p_new,
                            min_angle=min_angle, max_angle=max_angle,
                            max_area_change=max_area_change,
                        )
                        if not ok_q:
                            continue
                        nodes[i] = p_new
                        areas_all[ring] = new_areas
                        dev_before[v] = dev_v
                        for wn in watch:
                            dev_before[wn] = _dev_of(
                                nodes, wn, tangent_of[wn],
                                adj.get(wn, ()), in_seg,
                            )
                        n_accept += 1
                        accepted = True
                        break
                    if accepted:
                        break
                if not accepted:
                    unresolved.append(int(v))
        passes.append({"violations": n_viol_total, "accepted": n_accept})
        if n_viol_total == 0 or n_accept == 0:
            break

    out = Fort14Mesh(
        title=mesh.title,
        nodes=nodes,
        depths=mesh.depths.copy(),
        elements=mesh.elements.copy(),
        open_boundaries=[np.asarray(s).copy() for s in mesh.open_boundaries],
        land_boundaries=[(ib, np.asarray(s).copy())
                         for ib, s in mesh.land_boundaries],
    )
    info: dict[str, Any] = {
        "passes": passes,
        "accepted_total": int(sum(p["accepted"] for p in passes)),
        "remaining": unresolved,
        "dev_max": float(dev_max),
        "seed": int(seed),
    }
    return out, info


__all__ = ["align_open_boundary_local"]
