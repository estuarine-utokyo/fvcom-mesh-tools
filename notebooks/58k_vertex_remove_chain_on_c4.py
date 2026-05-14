"""PoC #58k: chain vertex_remove + stochastic cleanup on residual C4.

PoC #58j cleared the seed-invariant C1 fail (elem 12914 via
vertex_remove of its valence-7 interior node + stochastic cleanup),
leaving 2 marginal C4 fails on the PoC #58e residual. PoC #58k
applies the same recipe to those 2 surviving C4 edges: for each
C4 fail edge ``(u, v)`` between elements ``e_i`` and ``e_j``, try
``_apply_vertex_remove`` on every interior (non-boundary) vertex
of either element, run the stochastic local fixer on the result,
and keep the variant that strictly decreases the global FVCOM
violation count.

If both C4 edges have all-boundary incident elements (no interior
vertex), vertex_remove cannot help and we report the residual as
structurally beyond local topology.  Otherwise we expect to drive
the total residual to 0 or 1, completing the SMS phase-out goal
on this mesh.

Outputs:
    outputs/58k_chained.14
    outputs/58k_summary.{txt,json}
"""
from __future__ import annotations

import importlib.util
import json
import time
from collections import defaultdict
from pathlib import Path

import numpy as np

from fvcom_mesh_tools.algorithms.quality import min_interior_angle
from fvcom_mesh_tools.diagnostics import node_valence
from fvcom_mesh_tools.io import Fort14Mesh, read_fort14, write_fort14
from fvcom_mesh_tools.mesh_clean_phase_h import (
    _apply_vertex_remove,
    _boundary_node_mask,
    _boundary_topology,
    _edge_use_counts,
    _node_to_elements,
    _per_edge_area_change,
    build_coastline_projector,
)

REPO = Path(__file__).resolve().parent.parent
INPUT = REPO / "outputs" / "58j_chained.14"
COASTLINE = (
    REPO / "data" / "coastline" / "tokyo_bay"
    / "MLIT_C23" / "C23-06_TOKYOBAY.shp"
)
FIXER_SCRIPT = REPO / "notebooks" / "58d_stochastic_local_fix.py"
OUT_DIR = REPO / "outputs"
OUTPUT = OUT_DIR / "58k_chained.14"
SUMMARY_TXT = OUT_DIR / "58k_summary.txt"
SUMMARY_JSON = OUT_DIR / "58k_summary.json"

MIN_ANGLE_TARGET = 30.0
MAX_ANGLE_TARGET = 130.0
AREA_RATIO_TARGET = 0.5
MAX_VALENCE = 8
ALPHA_TARGET = 0.95
SEED = 42

# Import PoC #58d fixer (numeric-prefix filename precludes a plain import).
_spec = importlib.util.spec_from_file_location("fixer58d", str(FIXER_SCRIPT))
fixer58d = importlib.util.module_from_spec(_spec)
_spec.loader.exec_module(fixer58d)


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


def _metrics(mesh: Fort14Mesh) -> dict:
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


def _total(m: dict) -> int:
    return m["C1"] + m["C2"] + m["C4"] + m["C5"]


def _clone(mesh: Fort14Mesh) -> Fort14Mesh:
    return Fort14Mesh(
        title=mesh.title,
        nodes=mesh.nodes.copy(),
        depths=mesh.depths.copy(),
        elements=mesh.elements.copy(),
        open_boundaries=[np.asarray(s).copy() for s in mesh.open_boundaries],
        land_boundaries=[
            (int(ib), np.asarray(s).copy())
            for ib, s in mesh.land_boundaries
        ],
    )


def _run_stochastic_fixer(
    mesh: Fort14Mesh, projector,
) -> tuple[Fort14Mesh, dict]:
    """Run PoC #58d fixer on mesh and return updated mesh + stats."""
    rng = np.random.default_rng(SEED)
    fix_records: list[dict] = []
    t0 = time.perf_counter()
    for outer in range(1, fixer58d.MAX_OUTER_PASSES + 1):
        m = min_interior_angle(mesh)
        M = _max_interior_angle(mesh)
        _u, _p, ac = _per_edge_area_change(mesh.nodes, mesh.elements)
        c1_fail = m < MIN_ANGLE_TARGET
        c2_fail = M > MAX_ANGLE_TARGET
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
            break
        fail_list = sorted(fail_set)
        bnd_node = _boundary_node_mask(mesh)
        bnd_prev, bnd_next, _ = _boundary_topology(mesh)
        valence_arr = node_valence(mesh.elements, mesh.n_nodes)
        edge_uses = _edge_use_counts(mesh.elements)
        tmp: dict[int, list[int]] = defaultdict(list)
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
        for eid in fail_list:
            fixed, ntries = fixer58d._process_fail_element(
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
            if fixed:
                outer_fixed += 1
            else:
                outer_stuck += 1
        if outer_fixed == 0:
            break
    wall = time.perf_counter() - t0
    return mesh, {
        "wall_seconds": float(wall),
        "n_records": len(fix_records),
        "n_fixed": sum(1 for r in fix_records if r["fixed"]),
        "n_stuck": sum(1 for r in fix_records if not r["fixed"]),
    }


def _attempt_vertex_remove_and_clean(
    base_mesh: Fort14Mesh, vertex_id: int, projector,
    *, label: str,
) -> tuple[Fort14Mesh, dict, dict] | None:
    """Clone ``base_mesh``, force-apply ``_apply_vertex_remove`` on
    ``vertex_id``, run the stochastic cleanup, and return
    ``(new_mesh, metrics, info)``. Returns ``None`` if vertex_remove
    rejects (boundary / degenerate / Delaunay failure)."""
    trial = _clone(base_mesh)
    bnd_node = _boundary_node_mask(trial)
    n2e_map = _node_to_elements(trial.elements, trial.n_nodes)
    ring = n2e_map.get(int(vertex_id))
    if ring is None or ring.size < 3:
        return None
    out = _apply_vertex_remove(
        trial, int(vertex_id), ring,
        alpha_target=ALPHA_TARGET,
        min_angle_target=MIN_ANGLE_TARGET,
        max_angle_target=MAX_ANGLE_TARGET,
        boundary_node_mask=bnd_node,
        force=True,
    )
    if out is None:
        return None
    trial, vr_info = out
    m_after_vr = _metrics(trial)
    trial, fix_stats = _run_stochastic_fixer(trial, projector)
    m_final = _metrics(trial)
    return trial, m_final, {
        "label": label,
        "vertex_removed": int(vertex_id),
        "valence": int(vr_info.get("rim_size", -1)),
        "metrics_after_vr": m_after_vr,
        "metrics_after_clean": m_final,
        "fix_stats": fix_stats,
    }


def main() -> int:
    if not INPUT.exists():
        raise SystemExit(f"input missing: {INPUT}")
    if not FIXER_SCRIPT.exists():
        raise SystemExit(f"fixer script missing: {FIXER_SCRIPT}")
    OUT_DIR.mkdir(parents=True, exist_ok=True)

    mesh = read_fort14(INPUT)
    before = _metrics(mesh)
    print(
        f"[58k] input : NP={before['NP']:,} NE={before['NE']:,} "
        f"C1={before['C1']} C2={before['C2']} "
        f"C4={before['C4']} C5={before['C5']} (total {_total(before)})",
        flush=True,
    )

    projector = build_coastline_projector(
        [COASTLINE],
        max_snap_distance_m=500.0,
        mean_latitude_deg=float(mesh.nodes[:, 1].mean()),
    )

    # Iterate: each round attempts vertex_remove on every interior
    # candidate vertex of every current C4 fail pair, keeps the
    # trial that strictly drops the global total, and re-evaluates
    # the residual.  Stops when no candidate improves.
    trial_records: list[dict] = []
    iteration = 0
    while True:
        iteration += 1
        cur_metrics = _metrics(mesh)
        cur_total = _total(cur_metrics)
        if cur_metrics["C4"] == 0 and cur_metrics["C1"] == 0:
            print(f"[58k] iteration {iteration}: zero residual reached",
                  flush=True)
            break

        edge_uv, elem_pair, ac = _per_edge_area_change(
            mesh.nodes, mesh.elements,
        )
        c4_mask = ac > AREA_RATIO_TARGET
        bnd_node = _boundary_node_mask(mesh)

        candidates: set[int] = set()
        for idx in np.where(c4_mask)[0]:
            for eid in elem_pair[idx]:
                tri = mesh.elements[int(eid)]
                for v in tri:
                    if not bnd_node[int(v)]:
                        candidates.add(int(v))
        if not candidates:
            print(
                f"[58k] iteration {iteration}: no interior vertex on any "
                f"C4 fail pair (all_bnd cluster) — terminating",
                flush=True,
            )
            break

        # Order by valence descending (most central vertex first).
        valence_arr = node_valence(mesh.elements, mesh.n_nodes)
        ordered = sorted(
            candidates,
            key=lambda v: -int(valence_arr[v]),
        )
        print(
            f"[58k] iteration {iteration}: cur total={cur_total} "
            f"(C1={cur_metrics['C1']} C4={cur_metrics['C4']}) "
            f"trying {len(ordered)} interior vertices",
            flush=True,
        )

        best: tuple[int, dict, Fort14Mesh, dict] | None = None
        for v in ordered:
            label = f"vremove(v={v}, valence={int(valence_arr[v])})"
            result = _attempt_vertex_remove_and_clean(
                mesh, v, projector, label=label,
            )
            if result is None:
                trial_records.append({
                    "iteration": iteration, "vertex": int(v),
                    "result": "rejected_or_degenerate",
                })
                continue
            cand_mesh, cand_metrics, info = result
            cand_total = _total(cand_metrics)
            trial_records.append({
                "iteration": iteration, "vertex": int(v),
                "result": "applied",
                "metrics": cand_metrics,
                "total": cand_total,
                "info": info,
            })
            print(
                f"[58k]   {label}: after-clean total={cand_total} "
                f"({cand_metrics}) "
                f"fixer_wall={info['fix_stats']['wall_seconds']:.1f}s",
                flush=True,
            )
            if cand_total < cur_total and (
                best is None or cand_total < best[1]["total"]
            ):
                best = (
                    v,
                    {"total": cand_total, "metrics": cand_metrics},
                    cand_mesh,
                    info,
                )
                # Greedy: accept the first improving candidate, since
                # subsequent vertex_removes on the now-changed mesh
                # may behave differently.
                break

        if best is None:
            print(
                f"[58k] iteration {iteration}: no improving candidate "
                "— terminating",
                flush=True,
            )
            break
        v_best, best_score, mesh, _info = best
        print(
            f"[58k] iteration {iteration}: ACCEPTED vremove(v={v_best}) "
            f"→ total {cur_total} → {best_score['total']}",
            flush=True,
        )

    after = _metrics(mesh)
    write_fort14(mesh, OUTPUT)
    print(
        f"[58k] output: NP={after['NP']:,} NE={after['NE']:,} "
        f"C1={after['C1']} C2={after['C2']} "
        f"C4={after['C4']} C5={after['C5']} (total {_total(after)})",
        flush=True,
    )

    payload = {
        "input": str(INPUT.resolve()),
        "output": str(OUTPUT.resolve()),
        "seed": SEED,
        "before": before,
        "after": after,
        "n_iterations": int(iteration),
        "trial_records": trial_records,
    }
    SUMMARY_JSON.write_text(
        json.dumps(payload, indent=2, default=str), encoding="utf-8",
    )

    lines = [
        "PoC #58k — vertex_remove + stochastic chain on C4 residual",
        f"input : {INPUT.name}",
        f"output: {OUTPUT.name}",
        "",
        f"  {'stage':<28} | {'NP':>6} | {'NE':>6} | "
        f"{'C1':>4} | {'C2':>4} | {'C4':>4} | {'C5':>3} | total",
        "  " + "-" * 82,
        f"  {'PoC #58j (input)':<28} | "
        f"{before['NP']:>6,} | {before['NE']:>6,} | "
        f"{before['C1']:>4} | {before['C2']:>4} | "
        f"{before['C4']:>4} | {before['C5']:>3} | {_total(before)}",
        f"  {'PoC #58k (after chain)':<28} | "
        f"{after['NP']:>6,} | {after['NE']:>6,} | "
        f"{after['C1']:>4} | {after['C2']:>4} | "
        f"{after['C4']:>4} | {after['C5']:>3} | {_total(after)}",
        "",
        f"  iterations  : {iteration}",
        f"  trial_count : {len(trial_records)}",
        f"  delta total : {_total(after) - _total(before):+d}",
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
