"""PoC #44b: 2-step lookahead dry-run with the v4.1
``target_exits_fail`` gate.

PoC #44 measured "fixable" pairs under the v4 ``union_penalty`` gate
(any strict drop in the union penalty over op1 ∪ op2 affected nodes).
PoC #45 showed that gate is too permissive — accepts of O(10⁻³)
penalty drops do not lift the target out of fail status, and the
iterative driver thrashes. PoC #46 (initial bug-version) then
showed the alternative "fixed-by-elimination" interpretation of
``target_exits_fail`` is catastrophic. The corrected v4.1 gate is:

    accept iff target E exists in m_after AND
            alpha(E) >= alpha_target ∧ min_angle(E) >= min_angle_target

This PoC re-runs the PoC #44 measurement under that strict gate to
quote a defensible theoretical ceiling for v4.1 *before* PoC #46's
rerun result arrives. Sample is the same 1,000-element random
subset of the v3 residual (seed=42, n=1000) for direct comparison
to PoC #44.

Op1/op2 inventory matches the v4.1 driver default
(``DEFAULT_LOOKAHEAD_OP1_INVENTORY = ("smooth_node",)`` and
``DEFAULT_LOOKAHEAD_OP2_INVENTORY = ("smooth_node",)``). Destructive
ops (``vertex_remove``, ``edge_split_*``, ``edge_swap``) erase E's
vertex set so the strict gate rejects them by construction; running
them would only burn compute.

Outputs (read-only — no mesh is written):
    outputs/44b_phase_h_target_exits_summary.txt
    outputs/44b_phase_h_target_exits_summary.json
"""
from __future__ import annotations

import json
import os
import time
from collections import Counter
from pathlib import Path

import numpy as np

from fvcom_mesh_tools.io import read_fort14
from fvcom_mesh_tools.mesh_clean_phase_h import (
    DEFAULT_LOOKAHEAD_OP1_INVENTORY,
    DEFAULT_LOOKAHEAD_OP2_INVENTORY,
    _ctx_for_lookahead,
    _is_fail,
    _iter_op_candidates,
    _per_element_quality,
    _try_lookahead_pair,
)

REPO = Path(__file__).resolve().parent.parent
INPUT = REPO / "outputs" / "43_phase_h_v3_optimized.14"
OUT_DIR = REPO / "outputs"
SUMMARY_TXT = OUT_DIR / "44b_phase_h_target_exits_summary.txt"
SUMMARY_JSON = OUT_DIR / "44b_phase_h_target_exits_summary.json"

ALPHA_TARGET = 0.95
MIN_ANGLE_TARGET = 20.0

# Inventory matching v4.1 driver default — keeps the dry-run rate
# comparable to what the driver can in principle realise.
OP1_INVENTORY = DEFAULT_LOOKAHEAD_OP1_INVENTORY
OP2_INVENTORY = DEFAULT_LOOKAHEAD_OP2_INVENTORY
GATE = "target_exits_fail"

# Op inventory for the 1-step sanity check — full v3 inventory.
ONE_STEP_INVENTORY: tuple[str, ...] = (
    "smooth_node",
    "edge_swap",
    "edge_split_interior",
    "edge_split_boundary",
    "vertex_remove",
)


def _max_fail_env() -> int | None:
    raw = os.environ.get("POC44B_MAX_FAIL")
    if raw is None or raw == "":
        return None
    return int(raw)


MAX_FAIL_ELEMENTS: int | None = _max_fail_env()


def _has_one_step_fix(mesh, eid: int, ctx: dict) -> str | None:
    """Return the operator name that 1-step strict-gate accepts on
    ``eid``, or ``None``. Sanity check — for the v3 residual we
    expect this to fail almost always.
    """
    for op_name in ONE_STEP_INVENTORY:
        for _m_new, info in _iter_op_candidates(
            mesh, eid, op_name, force=False, ctx=ctx,
            alpha_target=ALPHA_TARGET, min_angle_target=MIN_ANGLE_TARGET,
        ):
            return info["operator"]
    return None


def main() -> None:
    print(f"[PoC #44b] loading {INPUT}", flush=True)
    mesh = read_fort14(INPUT)
    print(f"  NP={mesh.n_nodes:,}  NE={mesh.n_elements:,}", flush=True)

    ctx = _ctx_for_lookahead(mesh)
    alpha, min_ang = _per_element_quality(mesh.nodes, mesh.elements)
    fail = _is_fail(
        alpha, min_ang,
        alpha_target=ALPHA_TARGET, min_angle_target=MIN_ANGLE_TARGET,
    )
    fail_eids = np.where(fail)[0]
    print(
        f"  fail elements: {fail_eids.size:,} "
        f"(α<{ALPHA_TARGET} ∨ min_ang<{MIN_ANGLE_TARGET}°)",
        flush=True,
    )

    if MAX_FAIL_ELEMENTS is not None and MAX_FAIL_ELEMENTS < fail_eids.size:
        rng = np.random.default_rng(42)
        sampled = np.sort(rng.choice(
            fail_eids, size=MAX_FAIL_ELEMENTS, replace=False,
        ))
        print(
            f"  RANDOM-SAMPLED {MAX_FAIL_ELEMENTS:,} of {fail_eids.size:,} "
            f"(seed=42); shares the PoC #44 sample for A/B comparison",
            flush=True,
        )
        fail_eids = sampled

    one_step = 0
    op1_only = 0
    two_step = 0
    unfixable = 0
    one_step_hist: Counter = Counter()
    pair_hist: Counter = Counter()
    sample_residual: list[dict] = []

    OUT_DIR.mkdir(parents=True, exist_ok=True)

    def _write_summary(n_processed: int, elapsed_sec: float, partial: bool) -> None:
        total = max(one_step + op1_only + two_step + unfixable, 1)
        pair_lines = "\n".join(
            f"  {label:<28s} : {c:>5d}"
            for label, c in pair_hist.most_common(10)
        ) or "  (none)"
        status = "PARTIAL (job in progress / interrupted)" if partial else "FINAL"
        SUMMARY_TXT.write_text(
            f"PoC #44b — target_exits_fail dry-run on {INPUT.name}\n"
            f"status           = {status}\n"
            f"NP={mesh.n_nodes:,}  NE={mesh.n_elements:,}  "
            f"fail_total={fail_eids.size:,}  processed={n_processed}  "
            f"wall={elapsed_sec:.0f}s\n"
            f"\n"
            f"alpha_target     = {ALPHA_TARGET}\n"
            f"min_angle_target = {MIN_ANGLE_TARGET}°\n"
            f"op1 inventory    = {OP1_INVENTORY}\n"
            f"op2 inventory    = {OP2_INVENTORY}\n"
            f"gate             = {GATE}\n"
            f"\n"
            f"1-step fixable     : {one_step:6d} ({one_step/total:.1%})\n"
            f"op1-only fixable   : {op1_only:6d} ({op1_only/total:.1%})\n"
            f"2-step (op1+op2)   : {two_step:6d} ({two_step/total:.1%})\n"
            f"unfixable          : {unfixable:6d} ({unfixable/total:.1%})\n"
            f"  total fixable    : "
            f"{one_step + op1_only + two_step:6d} "
            f"({(one_step + op1_only + two_step)/total:.1%})\n"
            f"\n"
            f"accepted-pair histogram:\n{pair_lines}\n"
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
            "gate": GATE,
            "wall_seconds": float(elapsed_sec),
            "one_step_fixable": int(one_step),
            "op1_only_fixable": int(op1_only),
            "two_step_fixable": int(two_step),
            "unfixable": int(unfixable),
            "total_fixable_share": (
                (one_step + op1_only + two_step) / total
            ),
            "one_step_op_histogram": dict(one_step_hist),
            "accepted_pair_histogram": dict(pair_hist),
            "sample_residual": sample_residual,
        }, indent=2))

    CHECKPOINT_EVERY = 100
    t0 = time.time()
    for i, eid in enumerate(fail_eids):
        eid = int(eid)
        op1 = _has_one_step_fix(mesh, eid, ctx)
        if op1 is not None:
            one_step += 1
            one_step_hist[op1] += 1
        else:
            applied = _try_lookahead_pair(
                mesh, eid, ctx,
                alpha_target=ALPHA_TARGET,
                min_angle_target=MIN_ANGLE_TARGET,
                op1_inventory=OP1_INVENTORY,
                op2_inventory=OP2_INVENTORY,
                coastline_projector=None,
                gate=GATE,
            )
            if applied is not None:
                _new_mesh, pair_label = applied
                pair_hist[pair_label] += 1
                if pair_label.endswith("+none"):
                    op1_only += 1
                else:
                    two_step += 1
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
                f"1-step={one_step}  op1-only={op1_only}  "
                f"2-step={two_step}  unfix={unfixable}  "
                f"rate={rate:.2f}/s  eta={eta_min:.1f}min",
                flush=True,
            )
            _write_summary(i + 1, elapsed, partial=True)

    elapsed = time.time() - t0
    total = max(one_step + op1_only + two_step + unfixable, 1)
    print(f"[PoC #44b] done in {elapsed:.0f}s", flush=True)
    print(
        f"  1-step fixable     : {one_step:6d} ({one_step/total:.1%})",
        flush=True,
    )
    print(
        f"  op1-only fixable   : {op1_only:6d} ({op1_only/total:.1%})",
        flush=True,
    )
    print(
        f"  2-step fixable     : {two_step:6d} ({two_step/total:.1%})",
        flush=True,
    )
    print(
        f"  unfixable          : {unfixable:6d} ({unfixable/total:.1%})",
        flush=True,
    )
    print(
        f"  total fixable      : "
        f"{one_step + op1_only + two_step:6d} "
        f"({(one_step + op1_only + two_step)/total:.1%})",
        flush=True,
    )
    print(
        f"  pair histogram     : {dict(pair_hist.most_common(10))}",
        flush=True,
    )

    _write_summary(fail_eids.size, elapsed, partial=False)
    print(f"  wrote {SUMMARY_TXT}", flush=True)
    print(f"  wrote {SUMMARY_JSON}", flush=True)


if __name__ == "__main__":
    main()
