"""PoC #44: 2-step lookahead dry-run on Phase H v3 residual.

PoC #43 left ~11k elements at ``alpha < 0.95 ∨ min_angle < 20°``. PoC
#42 hypothesised that these residuals are "fundamental constraints of
greedy local edits": every 1-step move that would lift the element
drags a neighbour below the threshold, so the strict 1-ring penalty
gate rejects it.

The natural escape is **2-step lookahead** — apply op1 even if it
raises the local penalty (force the move), then immediately apply op2
that exploits op1's new local geometry. Accept the pair iff the
union-penalty over (op1 ∪ op2) affected nodes drops vs the initial
mesh. This PoC measures, before implementing the driver, how much of
the residual is "fixable by some (op1, op2)" pair.

Method (per fail element ``E`` on outputs/43_phase_h_v3_optimized.14):
  1. Sanity: 1-step strict gate on E.  We expect this to fail for
     ~all residuals because v3 already exhausted it; the counter
     guards against drift between PoC runs.
  2. For each ``op1`` in OP1_INVENTORY and each local variant of op1
     applied to E (3 edges / 3 vertices / etc.), call
     ``_apply_<op>(force=True)``.  Validity is still enforced
     (signed-area, boundary-mask, segment-topology).
  3. On the resulting m1, identify candidate ``e2`` elements: the
     elements of m1 that overlap op1's "affected" region (the small
     block op1 modified).  Restrict op2 candidates further to those
     that still fail the per-element gate on m1.
  4. For each (e2, op2 in OP2_INVENTORY × local variant), apply with
     ``force=True``.  Compute union penalty in initial mesh and m2
     over the set of nodes touched by op1 or op2.  Accept the pair
     iff union_after + 1e-12 < union_before strictly.
  5. Record the first accepting pair (op1, op2); break.

Outputs (read-only — no mesh is written):
    outputs/44_phase_h_2step_summary.txt
    outputs/44_phase_h_2step_summary.json
"""
from __future__ import annotations

import json
import time
from collections import Counter
from pathlib import Path

import numpy as np

from fvcom_mesh_tools.io import Fort14Mesh, read_fort14
from fvcom_mesh_tools.mesh_clean_phase_h import (
    _apply_edge_split_boundary,
    _apply_edge_split_interior,
    _apply_edge_swap,
    _apply_smooth_node,
    _apply_vertex_remove,
    _boundary_node_mask,
    _boundary_topology,
    _edge_use_counts,
    _is_fail,
    _node_to_elements,
    _penalty,
    _per_element_quality,
)

REPO = Path(__file__).resolve().parent.parent
INPUT = REPO / "outputs" / "43_phase_h_v3_optimized.14"
OUT_DIR = REPO / "outputs"
SUMMARY_TXT = OUT_DIR / "44_phase_h_2step_summary.txt"
SUMMARY_JSON = OUT_DIR / "44_phase_h_2step_summary.json"

ALPHA_TARGET = 0.95
MIN_ANGLE_TARGET = 20.0

# op1 sweeps the full Phase H v3 inventory. op2 is restricted to the
# two cheapest operators because a smooth or swap recoupment is the
# canonical "barrier-crossing" closer; splits/removes as op2 would
# stack two topology changes which is rarely productive and triples
# the search cost. If the dry-run says "0 % fixable" we will widen
# op2_inventory before declaring lookahead a dead end.
OP1_INVENTORY = (
    "smooth_node",
    "edge_swap",
    "edge_split_interior",
    "edge_split_boundary",
    "vertex_remove",
)
OP2_INVENTORY = ("smooth_node", "edge_swap")

def _max_fail_env() -> int | None:
    """``POC44_MAX_FAIL`` overrides the cap (smoke-test/local debug)."""
    import os
    raw = os.environ.get("POC44_MAX_FAIL")
    if raw is None or raw == "":
        return None
    return int(raw)


# Default: process all fail elements. Override via ``POC44_MAX_FAIL=N``
# environment variable for smoke-tests.
MAX_FAIL_ELEMENTS: int | None = _max_fail_env()


def _ctx_for_mesh(mesh: Fort14Mesh) -> dict:
    n2e = _node_to_elements(mesh.elements, mesh.n_nodes)
    eu = _edge_use_counts(mesh.elements)
    bnd_node = _boundary_node_mask(mesh)
    bnd_edges = {k for k, v in eu.items() if len(v) == 1}
    _bp, _bn, e2s = _boundary_topology(mesh)
    return {
        "n2e": n2e,
        "eu": eu,
        "bnd_node": bnd_node,
        "bnd_edges": bnd_edges,
        "e2s": e2s,
    }


def _affected_nodes(info: dict) -> set[int]:
    """Node IDs whose 1-ring may have been modified by an op.

    Returned IDs reference the *resulting* mesh (i.e. they are valid
    in m1 after op1 applied). New-node IDs from splits are included.
    """
    op = info["operator"]
    if op == "smooth_node":
        return {int(info["vertex"])}
    if op == "edge_swap":
        return {int(info["edge"][0]), int(info["edge"][1])}
    if op == "edge_split_interior":
        return {int(info["edge"][0]), int(info["edge"][1]), int(info["new_node"])}
    if op == "edge_split_boundary":
        return {int(info["edge"][0]), int(info["edge"][1]), int(info["new_node"])}
    if op == "vertex_remove":
        return {int(info["vertex"])}
    raise ValueError(f"unknown op: {op}")


def _affected_elements_in_m1(mesh: Fort14Mesh, info: dict) -> list[int]:
    """Element IDs in m1 that overlap op1's modified region.

    For pure-displacement ops (smooth_node) this is the 1-ring of the
    moved vertex.  For topology ops the relevant elements are either
    the modified pair (edge_swap) or the newly-appended block (splits,
    vertex_remove).
    """
    op = info["operator"]
    if op == "smooth_node":
        return [int(x) for x in info["affected_elements"]]
    if op == "edge_swap":
        return [int(x) for x in info["elements_modified"]]
    if op == "edge_split_interior":
        # 4 new triangles at the end of m1.elements after vstack.
        ne = int(mesh.n_elements)
        return list(range(ne - 4, ne))
    if op == "edge_split_boundary":
        ne = int(mesh.n_elements)
        return list(range(ne - 2, ne))
    if op == "vertex_remove":
        ne = int(mesh.n_elements)
        n_new = int(info["n_new_elements"])
        return list(range(ne - n_new, ne))
    raise ValueError(f"unknown op: {op}")


def _union_penalty(mesh: Fort14Mesh, nodes: set[int]) -> float:
    """Sum penalty over all elements in ``mesh`` that contain any node
    in ``nodes``. Out-of-range IDs (e.g. new-node IDs in the initial
    mesh) are silently skipped.
    """
    if not nodes:
        return 0.0
    n2e = _node_to_elements(mesh.elements, mesh.n_nodes)
    affected_eids: set[int] = set()
    for v in nodes:
        if 0 <= v < mesh.n_nodes:
            for e in n2e.get(int(v), ()):
                affected_eids.add(int(e))
    if not affected_eids:
        return 0.0
    eids = np.fromiter(affected_eids, dtype=np.int64)
    block = mesh.elements[eids]
    a, m = _per_element_quality(mesh.nodes, block)
    return float(_penalty(
        a, m, alpha_target=ALPHA_TARGET, min_angle_target=MIN_ANGLE_TARGET,
    ).sum())


def _try_op_variants(
    mesh: Fort14Mesh,
    op_name: str,
    elem_id: int,
    *,
    force: bool,
    ctx: dict,
):
    """Yield ``(new_mesh, info)`` for each local variant of op_name
    applied to elem_id with the given ``force`` flag.
    """
    if op_name == "smooth_node":
        for v in mesh.elements[elem_id]:
            ring = ctx["n2e"].get(int(v))
            if ring is None:
                continue
            out = _apply_smooth_node(
                mesh, int(v), ring,
                alpha_target=ALPHA_TARGET,
                min_angle_target=MIN_ANGLE_TARGET,
                boundary_node_mask=ctx["bnd_node"],
                force=force,
            )
            if out is not None:
                yield out
    elif op_name == "edge_swap":
        for k in range(3):
            out = _apply_edge_swap(
                mesh, int(elem_id), k,
                alpha_target=ALPHA_TARGET,
                min_angle_target=MIN_ANGLE_TARGET,
                edge_uses=ctx["eu"],
                boundary_edge_keys=ctx["bnd_edges"],
                force=force,
            )
            if out is not None:
                yield out
    elif op_name == "edge_split_interior":
        for k in range(3):
            out = _apply_edge_split_interior(
                mesh, int(elem_id), k,
                alpha_target=ALPHA_TARGET,
                min_angle_target=MIN_ANGLE_TARGET,
                edge_uses=ctx["eu"],
                boundary_edge_keys=ctx["bnd_edges"],
                force=force,
            )
            if out is not None:
                yield out
    elif op_name == "edge_split_boundary":
        for k in range(3):
            out = _apply_edge_split_boundary(
                mesh, int(elem_id), k,
                alpha_target=ALPHA_TARGET,
                min_angle_target=MIN_ANGLE_TARGET,
                edge_uses=ctx["eu"],
                edge_to_segment=ctx["e2s"],
                force=force,
            )
            if out is not None:
                yield out
    elif op_name == "vertex_remove":
        for v in mesh.elements[elem_id]:
            ring = ctx["n2e"].get(int(v))
            if ring is None:
                continue
            out = _apply_vertex_remove(
                mesh, int(v), ring,
                alpha_target=ALPHA_TARGET,
                min_angle_target=MIN_ANGLE_TARGET,
                boundary_node_mask=ctx["bnd_node"],
                force=force,
            )
            if out is not None:
                yield out
    else:
        raise ValueError(f"unknown op: {op_name}")


def _has_one_step_fix(mesh: Fort14Mesh, elem_id: int, ctx: dict) -> dict | None:
    """Return the first op-info that 1-step strict-gate accepts, else None."""
    for op_name in OP1_INVENTORY:
        for _m_new, info in _try_op_variants(
            mesh, op_name, elem_id, force=False, ctx=ctx,
        ):
            return info
    return None


def _try_two_step(
    mesh: Fort14Mesh, elem_id: int, ctx: dict,
) -> tuple[str, str] | None:
    """Search for a successful (op1, op2) pair on ``elem_id``.

    Returns ``(op1_name, op2_name)`` on the first accepting pair (in
    inventory order); ``None`` if every pair is rejected.
    """
    for op1_name in OP1_INVENTORY:
        for m1, info1 in _try_op_variants(
            mesh, op1_name, elem_id, force=True, ctx=ctx,
        ):
            affected_1_nodes = _affected_nodes(info1)
            affected_1_eids = _affected_elements_in_m1(m1, info1)
            if not affected_1_eids:
                continue
            # Quality of affected block in m1.
            block = m1.elements[affected_1_eids]
            a_blk, m_blk = _per_element_quality(m1.nodes, block)
            fail_blk = _is_fail(
                a_blk, m_blk,
                alpha_target=ALPHA_TARGET, min_angle_target=MIN_ANGLE_TARGET,
            )
            candidate_e2 = [
                eid for eid, f in zip(affected_1_eids, fail_blk) if f
            ]
            if not candidate_e2:
                continue
            ctx_m1 = _ctx_for_mesh(m1)
            for e2 in candidate_e2:
                for op2_name in OP2_INVENTORY:
                    for m2, info2 in _try_op_variants(
                        m1, op2_name, int(e2), force=True, ctx=ctx_m1,
                    ):
                        affected_2_nodes = _affected_nodes(info2)
                        union = affected_1_nodes | affected_2_nodes
                        pen_before = _union_penalty(mesh, union)
                        pen_after = _union_penalty(m2, union)
                        if pen_after + 1e-12 < pen_before:
                            return (op1_name, op2_name)
    return None


def main() -> None:
    print(f"[PoC #44] loading {INPUT}", flush=True)
    mesh = read_fort14(INPUT)
    print(f"  NP={mesh.n_nodes:,}  NE={mesh.n_elements:,}", flush=True)

    ctx = _ctx_for_mesh(mesh)
    alpha, min_ang = _per_element_quality(mesh.nodes, mesh.elements)
    fail = _is_fail(
        alpha, min_ang,
        alpha_target=ALPHA_TARGET, min_angle_target=MIN_ANGLE_TARGET,
    )
    fail_eids = np.where(fail)[0]
    print(f"  fail elements: {fail_eids.size:,} "
          f"(α<{ALPHA_TARGET} ∨ min_ang<{MIN_ANGLE_TARGET}°)", flush=True)

    if MAX_FAIL_ELEMENTS is not None and MAX_FAIL_ELEMENTS < fail_eids.size:
        rng = np.random.default_rng(42)
        sampled = np.sort(rng.choice(
            fail_eids, size=MAX_FAIL_ELEMENTS, replace=False,
        ))
        print(
            f"  RANDOM-SAMPLED {MAX_FAIL_ELEMENTS:,} of {fail_eids.size:,} "
            f"(seed=42); rate extrapolates to full residual",
            flush=True,
        )
        fail_eids = sampled

    one_step = 0
    two_step = 0
    unfixable = 0
    op1_hist: Counter = Counter()
    pair_hist: Counter = Counter()
    sample_residual: list[dict] = []

    OUT_DIR.mkdir(parents=True, exist_ok=True)

    def _write_summary(n_processed: int, elapsed_sec: float, partial: bool) -> None:
        total = max(one_step + two_step + unfixable, 1)
        pair_lines = "\n".join(
            f"  {p[0]:>22s} + {p[1]:<14s} : {c:>5d}"
            for p, c in pair_hist.most_common(10)
        ) or "  (none)"
        status = "PARTIAL (job in progress / interrupted)" if partial else "FINAL"
        SUMMARY_TXT.write_text(
            f"PoC #44 — 2-step lookahead dry-run on {INPUT.name}\n"
            f"status           = {status}\n"
            f"NP={mesh.n_nodes:,}  NE={mesh.n_elements:,}  "
            f"fail_total={fail_eids.size:,}  processed={n_processed}  "
            f"wall={elapsed_sec:.0f}s\n"
            f"\n"
            f"alpha_target     = {ALPHA_TARGET}\n"
            f"min_angle_target = {MIN_ANGLE_TARGET}°\n"
            f"op1 inventory    = {OP1_INVENTORY}\n"
            f"op2 inventory    = {OP2_INVENTORY}\n"
            f"\n"
            f"1-step fixable : {one_step:6d} ({one_step/total:.1%})\n"
            f"2-step fixable : {two_step:6d} ({two_step/total:.1%})\n"
            f"unfixable      : {unfixable:6d} ({unfixable/total:.1%})\n"
            f"\n"
            f"op-pair histogram (top 10):\n{pair_lines}\n"
        )
        SUMMARY_JSON.write_text(json.dumps({
            "status": status,
            "input": str(INPUT),
            "n_nodes": int(mesh.n_nodes),
            "n_elements": int(mesh.n_elements),
            "n_fail_total": int(fail_eids.size),
            "n_processed": int(n_processed),
            "alpha_target": ALPHA_TARGET,
            "min_angle_target": MIN_ANGLE_TARGET,
            "op1_inventory": list(OP1_INVENTORY),
            "op2_inventory": list(OP2_INVENTORY),
            "wall_seconds": float(elapsed_sec),
            "one_step_fixable": int(one_step),
            "two_step_fixable": int(two_step),
            "unfixable": int(unfixable),
            "op1_histogram": dict(op1_hist),
            "op_pair_histogram": {
                f"{k[0]}+{k[1]}": int(v) for k, v in pair_hist.items()
            },
            "sample_residual": sample_residual,
        }, indent=2))

    CHECKPOINT_EVERY = 100
    t0 = time.time()
    for i, eid in enumerate(fail_eids):
        eid = int(eid)
        info_1 = _has_one_step_fix(mesh, eid, ctx)
        if info_1 is not None:
            one_step += 1
            op1_hist[info_1["operator"]] += 1
        else:
            pair = _try_two_step(mesh, eid, ctx)
            if pair is not None:
                two_step += 1
                pair_hist[pair] += 1
            else:
                unfixable += 1
                if len(sample_residual) < 20:
                    sample_residual.append({
                        "elem_id": eid,
                        "alpha": float(alpha[eid]),
                        "min_angle_deg": float(min_ang[eid]),
                    })
        if (i + 1) % CHECKPOINT_EVERY == 0:
            elapsed = time.time() - t0
            rate = (i + 1) / elapsed
            eta_min = (fail_eids.size - i - 1) / rate / 60.0
            print(
                f"  [{i+1:>6d}/{fail_eids.size}] "
                f"1-step={one_step}  2-step={two_step}  unfix={unfixable}  "
                f"rate={rate:.2f}/s  eta={eta_min:.1f}min",
                flush=True,
            )
            _write_summary(i + 1, elapsed, partial=True)

    elapsed = time.time() - t0
    total = max(one_step + two_step + unfixable, 1)
    print(f"[PoC #44] done in {elapsed:.0f}s", flush=True)
    print(f"  1-step fixable  : {one_step:6d}  ({one_step/total:.1%})", flush=True)
    print(f"  2-step fixable  : {two_step:6d}  ({two_step/total:.1%})", flush=True)
    print(f"  unfixable       : {unfixable:6d}  ({unfixable/total:.1%})", flush=True)
    print(f"  op1 hist (1-step): {dict(op1_hist.most_common())}", flush=True)
    print(f"  op-pair hist top 10: "
          f"{dict(pair_hist.most_common(10))}", flush=True)

    _write_summary(fail_eids.size, elapsed, partial=False)
    print(f"  wrote {SUMMARY_TXT}", flush=True)
    print(f"  wrote {SUMMARY_JSON}", flush=True)


if __name__ == "__main__":
    main()
