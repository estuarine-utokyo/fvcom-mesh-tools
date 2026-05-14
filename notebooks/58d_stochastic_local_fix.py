"""PoC #58d: stochastic local fixer on PoC #57 Stage 1 output.

User-suggested approach (2026-05-14): mimic SMS manual editing via
stochastic local perturbation with rejection sampling. Pseudo-code:

  rng = numpy.random.default_rng(seed=42)
  for each fail element / fail edge in priority order:
      for try_i in range(max_tries):
          op = rng.choice(["move_vertex", "insert_vertex"], p=[0.7, 0.3])
          apply op with random parameters
          if local patch (fail element + 1-ring buddies) now satisfies
                  C1, C2, C4, C5 AND no flipped triangle:
              accept; break
          else:
              revert
      else:
          mark fail as "stuck" and move on

The seed-fixed RNG gives reproducibility. The local-only check
mirrors what a human does manually — they don't recompute the whole
mesh's metrics after every drag of a node, they look at the elements
touching the moved vertex.

This PoC is move-only first (operator weight 1.0 for "move",
0.0 for "insert") to validate the framework on the C1=2 + C4=48
residual. If move-only doesn't cover the full residual, the next
iteration adds insert.

Operator definitions:

  move_vertex:
    1. pick a vertex uniformly from the fail element's 3 vertices
    2. interior nodes: Gaussian perturbation scaled to mean local
       edge length × N(0, 0.3)
    3. boundary nodes: 1-D Gaussian along the boundary tangent line,
       projected onto the coastline if a projector is supplied,
       clamped so the parameter stays in [0.05, 0.95] of the segment

  insert_vertex (not used in this PoC, reserved for follow-up):
    1. sample barycentric (u, v, w) uniformly on the open 2-simplex
       inside the fail element
    2. add as a new node, split the element into 3 sub-triangles
       (and update the topology bookkeeping)

Local patch = the 1-ring of the moved vertex v, i.e. every element
containing v. This is **exactly** the set of elements whose
geometry changes when v moves, so a local quality check on this
patch is a complete check of the move's consequences (no new fail
can sneak in elsewhere — triangles not containing v are unchanged).

Outputs:
    outputs/58d_stochastic_local_fix.14
    outputs/58d_summary.{txt,json}
"""
from __future__ import annotations

import json
import time
from collections import defaultdict
from pathlib import Path

import numpy as np

from fvcom_mesh_tools.algorithms.quality import min_interior_angle
from fvcom_mesh_tools.diagnostics import node_valence
from fvcom_mesh_tools.io import Fort14Mesh, read_fort14, write_fort14
from fvcom_mesh_tools.mesh_clean_phase_h import (
    BOUNDARY_TANGENT_T_MIN,
    BOUNDARY_TANGENT_T_MAX,
    _boundary_node_mask,
    _boundary_topology,
    _edge_use_counts,
    _per_edge_area_change,
    _per_element_quality,
    build_coastline_projector,
)

REPO = Path(__file__).resolve().parent.parent
INPUT = REPO / "outputs" / "57_phase_h_pass_f_only.14"
COASTLINE = (
    REPO / "data" / "coastline" / "tokyo_bay"
    / "MLIT_C23" / "C23-06_TOKYOBAY.shp"
)
OUT_DIR = REPO / "outputs"
OUTPUT = OUT_DIR / "58d_stochastic_local_fix.14"
SUMMARY_TXT = OUT_DIR / "58d_summary.txt"
SUMMARY_JSON = OUT_DIR / "58d_summary.json"

SEED = 42
MIN_ANGLE_TARGET = 30.0
MAX_ANGLE_TARGET = 130.0
AREA_RATIO_TARGET = 0.5
MAX_VALENCE = 8

MAX_TRIES_PER_FAIL = 500
PERTURBATION_SIGMA = 0.30   # fraction of mean local edge length
MAX_OUTER_PASSES = 5        # restart from worst residual N times
OPERATOR_MOVE_WEIGHT = 1.0   # PoC #58d is move-only


def _max_interior_angle(mesh: Fort14Mesh) -> np.ndarray:
    p0 = mesh.nodes[mesh.elements[:, 0]]
    p1 = mesh.nodes[mesh.elements[:, 1]]
    p2 = mesh.nodes[mesh.elements[:, 2]]
    e0 = np.linalg.norm(p1 - p0, axis=1)
    e1 = np.linalg.norm(p2 - p1, axis=1)
    e2 = np.linalg.norm(p0 - p2, axis=1)

    def _ang(opp, ea, eb):
        denom = 2.0 * ea * eb
        safe = np.where(denom > 0, denom, 1.0)
        cos = np.where(
            denom > 0, (ea * ea + eb * eb - opp * opp) / safe, 0.0,
        )
        return np.arccos(np.clip(cos, -1.0, 1.0))

    return np.degrees(
        np.maximum(
            np.maximum(_ang(e1, e2, e0), _ang(e2, e0, e1)),
            _ang(e0, e1, e2),
        ),
    )


def _global_metrics(mesh: Fort14Mesh) -> dict:
    m = min_interior_angle(mesh)
    M = _max_interior_angle(mesh)
    _u, _p, ac = _per_edge_area_change(mesh.nodes, mesh.elements)
    val = node_valence(mesh.elements, mesh.n_nodes)
    return {
        "NP": int(mesh.n_nodes),
        "NE": int(mesh.n_elements),
        "C1": int((m < MIN_ANGLE_TARGET).sum()),
        "C2": int((M > MAX_ANGLE_TARGET).sum()),
        "C4": int((ac > AREA_RATIO_TARGET).sum()),
        "C5": int((val > MAX_VALENCE).sum()),
    }


def _vertex_ring_patch(
    n2e: dict[int, np.ndarray], v: int,
) -> np.ndarray:
    """Return the element IDs in the 1-ring of vertex ``v`` — i.e.,
    every element containing ``v``. This is the complete set of
    elements whose geometry changes when ``v`` moves."""
    ring = n2e.get(int(v))
    if ring is None:
        return np.empty(0, dtype=np.int64)
    return np.asarray(ring, dtype=np.int64)


def _affected_internal_edge_pairs(
    mesh: Fort14Mesh, patch_eids: np.ndarray,
    edge_uses: dict[tuple[int, int], list[int]],
) -> np.ndarray:
    """Return ``(M, 2)`` array of (e_i, e_j) element pairs sharing an
    internal edge such that at least one of e_i, e_j is in
    ``patch_eids``."""
    patch_set = set(int(e) for e in patch_eids)
    pairs: list[tuple[int, int]] = []
    seen: set[tuple[int, int]] = set()
    for e in patch_eids:
        tri = mesh.elements[int(e)]
        for k in range(3):
            a = int(tri[k])
            b = int(tri[(k + 1) % 3])
            key = (min(a, b), max(a, b))
            if key in seen:
                continue
            buds = edge_uses.get(key, [])
            if len(buds) != 2:
                continue
            seen.add(key)
            if int(buds[0]) in patch_set or int(buds[1]) in patch_set:
                pairs.append((int(buds[0]), int(buds[1])))
    if not pairs:
        return np.empty((0, 2), dtype=np.int64)
    return np.asarray(pairs, dtype=np.int64)


def _local_quality_ok(
    mesh: Fort14Mesh, patch_eids: np.ndarray,
    affected_pairs: np.ndarray,
    valence: np.ndarray,
) -> tuple[bool, dict]:
    """Check C1/C2/C4/C5 + signed area on the local patch. Returns
    ``(ok, diagnostics)``."""
    if patch_eids.size == 0:
        return True, {"reason": "empty_patch"}
    block = mesh.elements[patch_eids]
    p0 = mesh.nodes[block[:, 0]]
    p1 = mesh.nodes[block[:, 1]]
    p2 = mesh.nodes[block[:, 2]]
    cross = (
        (p1[:, 0] - p0[:, 0]) * (p2[:, 1] - p0[:, 1])
        - (p1[:, 1] - p0[:, 1]) * (p2[:, 0] - p0[:, 0])
    )
    if (cross <= 0).any():
        return False, {"reason": "flipped_triangle"}
    _alpha, min_ang, max_ang = _per_element_quality(mesh.nodes, block)
    c1 = int((min_ang < MIN_ANGLE_TARGET).sum())
    c2 = int((max_ang > MAX_ANGLE_TARGET).sum())
    if affected_pairs.shape[0]:
        pair_elems = mesh.elements[affected_pairs.ravel()].reshape(
            affected_pairs.shape[0], 2, 3,
        )
        pi = mesh.nodes[pair_elems[:, 0]]
        pj = mesh.nodes[pair_elems[:, 1]]
        ai = np.abs(0.5 * (
            (pi[:, 1, 0] - pi[:, 0, 0]) * (pi[:, 2, 1] - pi[:, 0, 1])
            - (pi[:, 1, 1] - pi[:, 0, 1]) * (pi[:, 2, 0] - pi[:, 0, 0])
        ))
        aj = np.abs(0.5 * (
            (pj[:, 1, 0] - pj[:, 0, 0]) * (pj[:, 2, 1] - pj[:, 0, 1])
            - (pj[:, 1, 1] - pj[:, 0, 1]) * (pj[:, 2, 0] - pj[:, 0, 0])
        ))
        larger = np.maximum(ai, aj)
        smaller = np.minimum(ai, aj)
        ratio = (larger - smaller) / np.maximum(larger, 1e-30)
        c4 = int((ratio > AREA_RATIO_TARGET).sum())
    else:
        c4 = 0
    # C5 is topology-invariant under move; check anyway for completeness.
    patch_verts = np.unique(block.ravel())
    c5 = int((valence[patch_verts] > MAX_VALENCE).sum())
    ok = (c1 == 0 and c2 == 0 and c4 == 0 and c5 == 0)
    return ok, {
        "c1_local": c1, "c2_local": c2, "c4_local": c4, "c5_local": c5,
        "reason": "ok" if ok else "criterion_fail",
    }


def _propose_move(
    mesh: Fort14Mesh, eid: int, rng,
    *,
    boundary_node_mask: np.ndarray,
    boundary_prev: np.ndarray,
    boundary_next: np.ndarray,
    coastline_projector,
    n2e: dict[int, np.ndarray],
) -> tuple[int, np.ndarray, np.ndarray] | None:
    """Propose a random node move. Returns
    ``(vertex_id, new_pos, ring_eids)`` or ``None`` if no candidate
    move can be constructed. ``ring_eids`` is the 1-ring of the
    chosen vertex — the patch on which to validate the move."""
    tri = mesh.elements[int(eid)]
    v = int(rng.choice(tri))
    ring_eids = _vertex_ring_patch(n2e, v)
    if ring_eids.size == 0:
        return None

    # Local mean edge length (over the vertex's 1-ring) as scale.
    block = mesh.elements[ring_eids]
    pa = mesh.nodes[block[:, 0]]
    pb = mesh.nodes[block[:, 1]]
    pc = mesh.nodes[block[:, 2]]
    edge_lens = np.concatenate([
        np.linalg.norm(pb - pa, axis=1),
        np.linalg.norm(pc - pb, axis=1),
        np.linalg.norm(pa - pc, axis=1),
    ])
    local_h = float(edge_lens.mean())
    if local_h <= 0:
        return None

    sigma = PERTURBATION_SIGMA * local_h
    if boundary_node_mask[v]:
        pv = int(boundary_prev[v])
        nv = int(boundary_next[v])
        if pv < 0 or nv < 0:
            return None  # corner, skip
        a_pos = mesh.nodes[pv]
        b_pos = mesh.nodes[nv]
        ab = b_pos - a_pos
        ab_len = float(np.linalg.norm(ab))
        if ab_len < 1e-20:
            return None
        # Current parameter along the tangent line.
        t_now = float(((mesh.nodes[v] - a_pos) @ ab) / (ab @ ab))
        # Gaussian perturbation in the parameter; scale by σ / ab_len so
        # the spatial sigma is ~PERTURBATION_SIGMA × local_h.
        delta_t = float(rng.standard_normal()) * (sigma / ab_len)
        t_new = max(
            BOUNDARY_TANGENT_T_MIN,
            min(BOUNDARY_TANGENT_T_MAX, t_now + delta_t),
        )
        new_pos = a_pos + t_new * ab
        if coastline_projector is not None:
            snapped = coastline_projector(new_pos)
            if snapped is not None:
                new_pos = snapped
    else:
        delta = rng.standard_normal(2) * sigma
        new_pos = mesh.nodes[v] + delta
    return v, np.asarray(new_pos, dtype=mesh.nodes.dtype), ring_eids


def _process_fail_element(
    mesh: Fort14Mesh, eid: int, rng,
    *,
    boundary_node_mask: np.ndarray,
    boundary_prev: np.ndarray,
    boundary_next: np.ndarray,
    coastline_projector,
    valence: np.ndarray,
    n2e: dict[int, np.ndarray],
    edge_uses: dict[tuple[int, int], list[int]],
) -> tuple[bool, int]:
    """Try up to MAX_TRIES_PER_FAIL random moves on a vertex of
    element ``eid``. Each try recomputes the patch as the 1-ring of
    the proposed vertex (the set of elements whose geometry actually
    changes), and accepts iff that patch passes every criterion.

    Returns ``(fixed, n_tries_used)``. ``fixed`` is True if the fail
    element no longer fails after some try; False if no accepted move
    cleared it within ``MAX_TRIES_PER_FAIL``.
    """
    # Pre-check: if the fail element already passes (could happen if
    # an earlier try in this outer pass already cleared it via a
    # neighbour's vertex move), report success with 0 tries.
    def _eid_still_fails() -> bool:
        block = mesh.elements[[int(eid)]]
        _a, m, M = _per_element_quality(mesh.nodes, block)
        if float(m[0]) < MIN_ANGLE_TARGET:
            return True
        if float(M[0]) > MAX_ANGLE_TARGET:
            return True
        # Check C4 on the two internal edges of eid against their
        # buddies.
        tri = mesh.elements[int(eid)]
        for k in range(3):
            a = int(tri[k])
            b = int(tri[(k + 1) % 3])
            key = (min(a, b), max(a, b))
            buds = edge_uses.get(key, [])
            if len(buds) != 2:
                continue
            be = buds[0] if buds[1] == int(eid) else buds[1]
            block2 = mesh.elements[[int(eid), int(be)]]
            p0 = mesh.nodes[block2[:, 0]]
            p1 = mesh.nodes[block2[:, 1]]
            p2 = mesh.nodes[block2[:, 2]]
            areas = np.abs(0.5 * (
                (p1[:, 0] - p0[:, 0]) * (p2[:, 1] - p0[:, 1])
                - (p1[:, 1] - p0[:, 1]) * (p2[:, 0] - p0[:, 0])
            ))
            larger = areas.max()
            smaller = areas.min()
            if (larger - smaller) / max(larger, 1e-30) > AREA_RATIO_TARGET:
                return True
        return False

    if not _eid_still_fails():
        return True, 0
    for try_i in range(1, MAX_TRIES_PER_FAIL + 1):
        proposal = _propose_move(
            mesh, eid, rng,
            boundary_node_mask=boundary_node_mask,
            boundary_prev=boundary_prev,
            boundary_next=boundary_next,
            coastline_projector=coastline_projector,
            n2e=n2e,
        )
        if proposal is None:
            continue
        v, new_pos, ring_eids = proposal
        affected_pairs = _affected_internal_edge_pairs(
            mesh, ring_eids, edge_uses,
        )
        old_pos = mesh.nodes[v].copy()
        mesh.nodes[v] = new_pos
        ok, _diag = _local_quality_ok(
            mesh, ring_eids, affected_pairs, valence,
        )
        if ok and not _eid_still_fails():
            return True, try_i
        mesh.nodes[v] = old_pos
    return False, MAX_TRIES_PER_FAIL


def main() -> int:
    if not INPUT.exists():
        raise SystemExit(f"input missing: {INPUT}")
    OUT_DIR.mkdir(parents=True, exist_ok=True)

    mesh = read_fort14(INPUT)
    before = _global_metrics(mesh)
    print(
        f"[58d] input : NP={before['NP']:,} NE={before['NE']:,} "
        f"C1={before['C1']} C2={before['C2']} "
        f"C4={before['C4']} C5={before['C5']}",
        flush=True,
    )

    projector = build_coastline_projector(
        [COASTLINE],
        max_snap_distance_m=500.0,
        mean_latitude_deg=float(mesh.nodes[:, 1].mean()),
    )
    rng = np.random.default_rng(SEED)

    fix_records: list[dict] = []
    t_total = time.perf_counter()

    for outer in range(1, MAX_OUTER_PASSES + 1):
        m = min_interior_angle(mesh)
        M = _max_interior_angle(mesh)
        _u, _p, ac = _per_edge_area_change(mesh.nodes, mesh.elements)
        c1_fail = m < MIN_ANGLE_TARGET
        c2_fail = M > MAX_ANGLE_TARGET
        # Map C4 fail edges -> their incident elements.
        edge_uv, elem_pair, _ac_full = _per_edge_area_change(
            mesh.nodes, mesh.elements,
        )
        c4_edge_fail = _ac_full > AREA_RATIO_TARGET
        c4_fail_elems = set()
        for pair in elem_pair[c4_edge_fail]:
            c4_fail_elems.add(int(pair[0]))
            c4_fail_elems.add(int(pair[1]))

        fail_set = (
            set(int(e) for e in np.where(c1_fail)[0])
            | set(int(e) for e in np.where(c2_fail)[0])
            | c4_fail_elems
        )
        if not fail_set:
            print(f"[58d] outer pass {outer}: all fails cleared",
                  flush=True)
            break

        # Sort fail elements deterministically (by id) for reproducibility.
        fail_list = sorted(fail_set)
        print(
            f"[58d] outer pass {outer}: {len(fail_list)} fail elements"
            f" (C1={int(c1_fail.sum())}, C2={int(c2_fail.sum())}, "
            f"C4_elems={len(c4_fail_elems)})",
            flush=True,
        )

        bnd_node = _boundary_node_mask(mesh)
        bnd_prev, bnd_next, _ = _boundary_topology(mesh)
        valence_arr = node_valence(mesh.elements, mesh.n_nodes)
        edge_uses = _edge_use_counts(mesh.elements)
        n2e_map: dict[int, np.ndarray] = defaultdict(
            lambda: np.empty(0, dtype=np.int64),
        )
        # Build node-to-element map.
        from collections import defaultdict as _defaultdict
        tmp: dict[int, list[int]] = _defaultdict(list)
        for k, tri in enumerate(mesh.elements):
            tmp[int(tri[0])].append(k)
            tmp[int(tri[1])].append(k)
            tmp[int(tri[2])].append(k)
        n2e_map = {
            int(v): np.asarray(es, dtype=np.int64)
            for v, es in tmp.items()
        }

        outer_fixed = 0
        outer_stuck = 0
        outer_tries_sum = 0
        for eid in fail_list:
            fixed, ntries = _process_fail_element(
                mesh, eid, rng,
                boundary_node_mask=bnd_node,
                boundary_prev=bnd_prev,
                boundary_next=bnd_next,
                coastline_projector=projector,
                valence=valence_arr,
                n2e=n2e_map,
                edge_uses=edge_uses,
            )
            fix_records.append({
                "outer": outer, "eid": int(eid),
                "fixed": bool(fixed), "n_tries": int(ntries),
            })
            outer_tries_sum += ntries
            if fixed:
                outer_fixed += 1
            else:
                outer_stuck += 1

        avg_tries = (
            outer_tries_sum / max(1, outer_fixed + outer_stuck)
        )
        print(
            f"[58d] outer pass {outer}: fixed={outer_fixed} "
            f"stuck={outer_stuck} avg_tries={avg_tries:.1f}",
            flush=True,
        )
        if outer_fixed == 0:
            print(f"[58d] no progress in outer pass {outer}; "
                  "terminating", flush=True)
            break

    wall = time.perf_counter() - t_total
    after = _global_metrics(mesh)
    write_fort14(mesh, OUTPUT)
    print(
        f"[58d] output: NP={after['NP']:,} NE={after['NE']:,} "
        f"C1={after['C1']} C2={after['C2']} "
        f"C4={after['C4']} C5={after['C5']} "
        f"(wall {wall:.1f} s)",
        flush=True,
    )

    payload = {
        "input": str(INPUT.resolve()),
        "output": str(OUTPUT.resolve()),
        "seed": SEED,
        "max_tries_per_fail": MAX_TRIES_PER_FAIL,
        "perturbation_sigma": PERTURBATION_SIGMA,
        "max_outer_passes": MAX_OUTER_PASSES,
        "before": before,
        "after": after,
        "delta": {k: after[k] - before[k] for k in before},
        "wall_seconds": wall,
        "n_records": len(fix_records),
        "n_fixed_total": sum(1 for r in fix_records if r["fixed"]),
        "n_stuck_total": sum(1 for r in fix_records if not r["fixed"]),
    }
    SUMMARY_JSON.write_text(
        json.dumps(payload, indent=2, default=str), encoding="utf-8",
    )

    total_b = before["C1"] + before["C2"] + before["C4"] + before["C5"]
    total_a = after["C1"] + after["C2"] + after["C4"] + after["C5"]
    lines = [
        "PoC #58d — stochastic local fixer (move-only, seed=42)",
        f"input : {INPUT.name}",
        f"output: {OUTPUT.name}",
        "",
        f"  {'stage':<28} | {'NP':>6} | {'NE':>6} | "
        f"{'C1':>4} | {'C2':>4} | {'C4':>4} | {'C5':>3} | total",
        "  " + "-" * 82,
        f"  {'PoC #57 Stage 1 (input)':<28} | "
        f"{before['NP']:>6,} | {before['NE']:>6,} | "
        f"{before['C1']:>4} | {before['C2']:>4} | "
        f"{before['C4']:>4} | {before['C5']:>3} | {total_b}",
        f"  {'PoC #58d (after fixer)':<28} | "
        f"{after['NP']:>6,} | {after['NE']:>6,} | "
        f"{after['C1']:>4} | {after['C2']:>4} | "
        f"{after['C4']:>4} | {after['C5']:>3} | {total_a}",
        "",
        f"  wall          : {wall:.2f} s",
        f"  n_fix_records : {len(fix_records)}",
        f"  n_fixed       : {payload['n_fixed_total']}",
        f"  n_stuck       : {payload['n_stuck_total']}",
        f"  delta total   : {total_a - total_b:+d}",
    ]
    SUMMARY_TXT.write_text("\n".join(lines) + "\n", encoding="utf-8")
    print()
    print("\n".join(lines))
    print()
    print(f"wrote {OUTPUT}")
    print(f"wrote {SUMMARY_TXT}")
    print(f"wrote {SUMMARY_JSON}")
    return 0


if __name__ == "__main__":
    raise SystemExit(main())
