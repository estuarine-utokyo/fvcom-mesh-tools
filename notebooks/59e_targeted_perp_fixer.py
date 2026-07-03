"""PoC #59e — quality-gated LOCAL perpendicularity fixer.

#59d converged 20/21 gates; the residual is ``obc_perpendicularity``
on 8 Uraga-transect nodes (worst 29.0°). The remaining obstacle is
that ``align_open_boundary_first_ring`` moves the first ring of EVERY
OBC node, so each global perpfix undoes finish-stage quality work on
the ~42 already-passing nodes and vice versa (damped oscillation).

This PoC moves only the first-ring interior partners of the VIOLATING
nodes, one node at a time, with strict local acceptance (the
stochastic-local-fixer recipe applied to perpendicularity):

* new best-edge deviation of the target OBC node <= threshold;
* every element of the moved node's 1-ring keeps C1 >= 30°,
  C2 <= 130°, positive area; every internal edge touching the ring
  keeps C4 <= 0.5;
* no OBC node adjacent to the moved node newly violates the
  perpendicularity threshold.

Candidate positions: damped steps toward the exact perpendicular
target (edge length preserved) plus seeded Gaussian jitter.

Output: outputs/59e_gate_passed.14 (+ _qa.json).
"""

from __future__ import annotations

import json
import time
from pathlib import Path

import numpy as np

from fvcom_mesh_tools.algorithms.perpendicularity import boundary_tangents
from fvcom_mesh_tools.io import read_fort14, write_fort14
from fvcom_mesh_tools.mesh_clean_phase_h import _inline_quality
from fvcom_mesh_tools.qa import _edge_topology, format_report, run_qa

REPO = Path(__file__).resolve().parents[1]
SRC = REPO / "outputs" / "59d_gate_passed.14"
OUT = REPO / "outputs" / "59e_gate_passed.14"
QA_JSON = REPO / "outputs" / "59e_gate_passed_qa.json"

DEV_MAX = 20.0          # QA gate
DEV_TARGET = 19.0       # accept below this (1 deg headroom)
MIN_ANG = 30.0 - 1e-9
MAX_ANG = 130.0 + 1e-9
MAX_AC = 0.5 + 1e-12
MAX_OUTER = 6
W_STEPS = (1.0, 0.8, 0.6, 0.4, 0.2)
N_JITTER = 40
JITTER_SIGMA_FRAC = 0.05
SEED = 4242


def _build_ctx(mesh):
    """Per-pass topology context."""
    n = mesh.n_nodes
    topo = _edge_topology(mesh.elements, n)
    boundary_nodes = np.zeros(n, dtype=bool)
    buv = topo.uv[topo.counts == 1]
    boundary_nodes[buv.ravel()] = True

    # node -> incident elements
    rows = mesh.elements.ravel()
    cols = np.repeat(np.arange(mesh.n_elements), 3)
    order = np.argsort(rows, kind="stable")
    rows, cols = rows[order], cols[order]
    starts = np.searchsorted(rows, np.arange(n + 1))
    n2e = {v: cols[starts[v]:starts[v + 1]] for v in range(n)
           if starts[v] < starts[v + 1]}

    # element -> edge-neighbour elements
    e2nbr: dict[int, list[int]] = {}
    for a, b in topo.internal_pair:
        e2nbr.setdefault(int(a), []).append(int(b))
        e2nbr.setdefault(int(b), []).append(int(a))

    # node -> mesh-edge neighbours
    adj: dict[int, list[int]] = {}
    for u, v in topo.uv:
        adj.setdefault(int(u), []).append(int(v))
        adj.setdefault(int(v), []).append(int(u))
    return topo, boundary_nodes, n2e, e2nbr, adj


def _areas(nodes, elements):
    p0, p1, p2 = nodes[elements[:, 0]], nodes[elements[:, 1]], nodes[elements[:, 2]]
    return 0.5 * np.abs(
        (p1[:, 0] - p0[:, 0]) * (p2[:, 1] - p0[:, 1])
        - (p1[:, 1] - p0[:, 1]) * (p2[:, 0] - p0[:, 0])
    )


def _dev_of(nodes, v, tangent, nbrs, in_seg):
    """Best-edge deviation (deg) of OBC node v; inf if no interior edge."""
    best = np.inf
    for i in nbrs:
        if in_seg[i]:
            continue
        vec = nodes[i] - nodes[v]
        nrm = np.hypot(vec[0], vec[1])
        if nrm == 0:
            continue
        cosang = abs(float(vec @ tangent)) / nrm
        dev = 90.0 - np.degrees(np.arccos(np.clip(cosang, 0.0, 1.0)))
        best = min(best, dev)
    return best


def _local_ok(nodes, elements, ring, e2nbr, areas_all, i, p_new):
    """Quality gates on i's 1-ring with i moved to p_new."""
    tri = elements[ring]
    pts = nodes[tri].copy()          # (R, 3, 2)
    pts[tri == i] = p_new
    _alpha, min_ang, max_ang, twice = _inline_quality(
        pts[:, 0], pts[:, 1], pts[:, 2],
    )
    if (twice <= 0).any() or min_ang.min() < MIN_ANG or max_ang.max() > MAX_ANG:
        return False, None
    new_areas = 0.5 * np.abs(twice)
    ring_pos = {int(e): k for k, e in enumerate(ring)}
    for k, e in enumerate(ring):
        a_e = new_areas[k]
        for nb in e2nbr.get(int(e), ()):
            a_nb = new_areas[ring_pos[nb]] if nb in ring_pos else areas_all[nb]
            hi = max(a_e, a_nb)
            if hi > 0 and (hi - min(a_e, a_nb)) / hi > MAX_AC:
                return False, None
    return True, new_areas


def main() -> int:
    t0 = time.perf_counter()
    mesh = read_fort14(SRC)
    rng = np.random.default_rng(SEED)
    nodes = mesh.nodes  # mutated in place
    n_accept_total = 0

    for outer in range(MAX_OUTER):
        topo, boundary_nodes, n2e, e2nbr, adj = _build_ctx(mesh)
        areas_all = _areas(nodes, mesh.elements)
        seg = np.asarray(mesh.open_boundaries[0], dtype=np.int64)
        in_seg = np.zeros(mesh.n_nodes, dtype=bool)
        in_seg[seg] = True
        tang = boundary_tangents(nodes[seg])
        tangent_of = {int(v): tang[k] for k, v in enumerate(seg)}

        dev_before = {
            int(v): _dev_of(nodes, int(v), tangent_of[int(v)], adj[int(v)], in_seg)
            for v in seg
        }
        viol = [v for v, d in dev_before.items() if d > DEV_MAX]
        print(f"[59e] pass {outer}: {len(viol)} violating OBC nodes", flush=True)
        if not viol:
            break

        n_accept = 0
        for v in viol:
            t_v = tangent_of[v]
            n_hat = np.array([-t_v[1], t_v[0]])
            cands = [i for i in adj[v]
                     if not in_seg[i] and not boundary_nodes[i]]
            # Try the currently-best edge partner first.
            cands.sort(key=lambda i: abs(
                90.0 - np.degrees(np.arccos(np.clip(
                    abs(float((nodes[i] - nodes[v]) @ t_v))
                    / max(np.hypot(*(nodes[i] - nodes[v])), 1e-30),
                    0.0, 1.0,
                )))
            ))
            accepted = False
            for i in cands:
                edge = nodes[i] - nodes[v]
                length = np.hypot(edge[0], edge[1])
                side = np.sign(float(edge @ n_hat)) or 1.0
                target = nodes[v] + length * side * n_hat
                ring = n2e[i]
                # OBC nodes whose deviation must not newly violate:
                watch = [int(w) for w in adj[i]
                         if in_seg[w] and int(w) != v]
                p_old = nodes[i].copy()
                trials = [(w, np.zeros(2)) for w in W_STEPS]
                trials += [
                    (w, rng.normal(0.0, JITTER_SIGMA_FRAC * length, 2))
                    for w in W_STEPS for _ in range(N_JITTER // len(W_STEPS))
                ]
                for w, jit in trials:
                    p_new = (1.0 - w) * p_old + w * target + jit
                    nodes[i] = p_new
                    dev_v = _dev_of(nodes, v, t_v, adj[v], in_seg)
                    ok = dev_v <= DEV_TARGET
                    if ok:
                        for wn in watch:
                            d_new = _dev_of(
                                nodes, wn, tangent_of[wn], adj[wn], in_seg,
                            )
                            d_old = dev_before[wn]
                            # Passing neighbours must stay passing;
                            # already-violating ones must not worsen.
                            if (d_old <= DEV_MAX and d_new > DEV_MAX) or (
                                d_old > DEV_MAX and d_new > d_old + 1e-9
                            ):
                                ok = False
                                break
                    nodes[i] = p_old
                    if not ok:
                        continue
                    ok_q, new_areas = _local_ok(
                        nodes, mesh.elements, ring, e2nbr, areas_all, i, p_new,
                    )
                    if not ok_q:
                        continue
                    nodes[i] = p_new
                    areas_all[ring] = new_areas
                    dev_before[v] = dev_v
                    for wn in watch:
                        dev_before[wn] = _dev_of(
                            nodes, wn, tangent_of[wn], adj[wn], in_seg,
                        )
                    n_accept += 1
                    accepted = True
                    break
                if accepted:
                    break
            if not accepted:
                print(f"[59e]   node {v}: no acceptable move found", flush=True)
        n_accept_total += n_accept
        print(f"[59e] pass {outer}: {n_accept} accepted moves", flush=True)
        if n_accept == 0:
            break

    report = run_qa(mesh, name=OUT.name, path=OUT)
    mesh.title = "PoC 59e UTM54N Tokyo Bay, locally perp-fixed Uraga OBC"
    write_fort14(mesh, OUT)
    QA_JSON.write_text(
        json.dumps(report.to_dict(), indent=2, ensure_ascii=False), encoding="utf-8",
    )
    print(format_report(report, lang="ja"), flush=True)
    print(f"[59e] total accepted moves: {n_accept_total}", flush=True)
    print(f"[59e] wrote {OUT}", flush=True)
    print(f"[59e] wall: {time.perf_counter() - t0:.1f} s", flush=True)
    return 0 if report.passed else 1


if __name__ == "__main__":
    raise SystemExit(main())
