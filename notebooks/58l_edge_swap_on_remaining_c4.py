"""PoC #58l: edge_swap + stochastic cleanup on remaining C4 fail.

PoC #58k took the residual to 1 violation — a single C4 fail edge
whose adjacent pair of elements has all 4 vertices on the
boundary, so vertex_remove is not applicable.  The shared edge
between the two elements is internal (by definition: 2 elements
share it), and ``_apply_edge_swap(force=True)`` can swap its
diagonal to rebalance the two areas without touching any vertex.

Pipeline:

  1. Read ``outputs/58k_chained.14`` (1 C4 fail).
  2. Locate the C4 fail edge ``(u, v)`` and its element pair
     ``(e_i, e_j)``.
  3. Call ``_apply_edge_swap(force=True)`` on the edge.
  4. Validate globally: total residual must not regress.
  5. If applied + no regression, run the PoC #58d stochastic
     cleanup to absorb any incidental fails created by the swap.
  6. Write output.

If the swap is rejected (invalid block or regression), we report
1 C4 as the structural floor for the auto-pipeline and recommend
either coastline editing or accepting the residual.

Outputs:
    outputs/58l_chained.14
    outputs/58l_summary.{txt,json}
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
    _apply_edge_swap,
    _boundary_node_mask,
    _boundary_topology,
    _edge_use_counts,
    _per_edge_area_change,
    build_coastline_projector,
)

REPO = Path(__file__).resolve().parent.parent
INPUT = REPO / "outputs" / "58k_chained.14"
COASTLINE = (
    REPO / "data" / "coastline" / "tokyo_bay"
    / "MLIT_C23" / "C23-06_TOKYOBAY.shp"
)
FIXER_SCRIPT = REPO / "notebooks" / "58d_stochastic_local_fix.py"
OUT_DIR = REPO / "outputs"
OUTPUT = OUT_DIR / "58l_chained.14"
SUMMARY_TXT = OUT_DIR / "58l_summary.txt"
SUMMARY_JSON = OUT_DIR / "58l_summary.json"

MIN_ANGLE_TARGET = 30.0
MAX_ANGLE_TARGET = 130.0
AREA_RATIO_TARGET = 0.5
MAX_VALENCE = 8
ALPHA_TARGET = 0.95
SEED = 42

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


def _run_stochastic_fixer(
    mesh: Fort14Mesh, projector,
) -> tuple[Fort14Mesh, dict]:
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
        if outer_fixed == 0:
            break
    return mesh, {
        "wall_seconds": float(time.perf_counter() - t0),
        "n_records": len(fix_records),
        "n_fixed": sum(1 for r in fix_records if r["fixed"]),
    }


def main() -> int:
    if not INPUT.exists():
        raise SystemExit(f"input missing: {INPUT}")
    OUT_DIR.mkdir(parents=True, exist_ok=True)

    mesh = read_fort14(INPUT)
    before = _metrics(mesh)
    print(
        f"[58l] input : NP={before['NP']:,} NE={before['NE']:,} "
        f"C1={before['C1']} C2={before['C2']} "
        f"C4={before['C4']} C5={before['C5']} (total {_total(before)})",
        flush=True,
    )

    projector = build_coastline_projector(
        [COASTLINE],
        max_snap_distance_m=500.0,
        mean_latitude_deg=float(mesh.nodes[:, 1].mean()),
    )

    # Locate the C4 fail edge.
    edge_uv, elem_pair, ac = _per_edge_area_change(
        mesh.nodes, mesh.elements,
    )
    c4_mask = ac > AREA_RATIO_TARGET
    fail_indices = np.where(c4_mask)[0]
    if not fail_indices.size:
        print("[58l] no C4 fails — already at zero", flush=True)
        after = before
    else:
        # Iterate over each C4 fail edge and try the swap.
        edge_uses_map = _edge_use_counts(mesh.elements)
        boundary_edge_keys = {
            k for k, v in edge_uses_map.items() if len(v) == 1
        }
        applied = False
        for idx in fail_indices:
            u, v = int(edge_uv[idx, 0]), int(edge_uv[idx, 1])
            e_i, e_j = int(elem_pair[idx, 0]), int(elem_pair[idx, 1])
            print(
                f"[58l] trying edge_swap on edge ({u}, {v}) between "
                f"elements {e_i} and {e_j} (ratio {ac[idx]:.3f})",
                flush=True,
            )
            # edge_local index of (u, v) in e_i.
            tri = mesh.elements[e_i]
            edge_local = -1
            for k in range(3):
                a = int(tri[k])
                b = int(tri[(k + 1) % 3])
                if {a, b} == {u, v}:
                    edge_local = k
                    break
            if edge_local < 0:
                print(
                    f"[58l]   could not locate edge_local for e_i={e_i}; "
                    "skip",
                    flush=True,
                )
                continue
            out = _apply_edge_swap(
                mesh, e_i, edge_local,
                alpha_target=ALPHA_TARGET,
                min_angle_target=MIN_ANGLE_TARGET,
                max_angle_target=MAX_ANGLE_TARGET,
                edge_uses=edge_uses_map,
                boundary_edge_keys=boundary_edge_keys,
                force=True,
            )
            if out is None:
                print(
                    f"[58l]   edge_swap rejected (boundary edge / "
                    "validity)",
                    flush=True,
                )
                continue
            candidate, info = out
            print(
                f"[58l]   edge_swap applied: penalty {info['penalty_before']:.3e} → "
                f"{info['penalty_after']:.3e}",
                flush=True,
            )
            cand_m = _metrics(candidate)
            print(
                f"[58l]   post-swap metrics: C1={cand_m['C1']} C2={cand_m['C2']} "
                f"C4={cand_m['C4']} C5={cand_m['C5']} "
                f"(total {_total(cand_m)})",
                flush=True,
            )
            if _total(cand_m) > _total(before):
                print(
                    f"[58l]   total regressed {_total(before)} → "
                    f"{_total(cand_m)} — running stochastic cleanup "
                    "to absorb",
                    flush=True,
                )
            # Run cleanup regardless — never hurts.
            candidate, fix_stats = _run_stochastic_fixer(
                candidate, projector,
            )
            final_m = _metrics(candidate)
            print(
                f"[58l]   after cleanup: C1={final_m['C1']} "
                f"C2={final_m['C2']} C4={final_m['C4']} "
                f"C5={final_m['C5']} (total {_total(final_m)}, "
                f"wall {fix_stats['wall_seconds']:.1f}s)",
                flush=True,
            )
            if _total(final_m) < _total(before):
                mesh = candidate
                applied = True
                print(
                    f"[58l]   ACCEPTED: total {_total(before)} → "
                    f"{_total(final_m)}",
                    flush=True,
                )
                break
        if not applied:
            print(
                "[58l] no swap candidate improved the total — "
                "1 C4 fail is structurally beyond local topology",
                flush=True,
            )

    after = _metrics(mesh)
    write_fort14(mesh, OUTPUT)
    print(
        f"[58l] output: NP={after['NP']:,} NE={after['NE']:,} "
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
    }
    SUMMARY_JSON.write_text(
        json.dumps(payload, indent=2, default=str), encoding="utf-8",
    )

    lines = [
        "PoC #58l — edge_swap + stochastic cleanup on remaining C4",
        f"input : {INPUT.name}",
        f"output: {OUTPUT.name}",
        "",
        f"  {'stage':<28} | {'NP':>6} | {'NE':>6} | "
        f"{'C1':>4} | {'C2':>4} | {'C4':>4} | {'C5':>3} | total",
        "  " + "-" * 82,
        f"  {'PoC #58k (input)':<28} | "
        f"{before['NP']:>6,} | {before['NE']:>6,} | "
        f"{before['C1']:>4} | {before['C2']:>4} | "
        f"{before['C4']:>4} | {before['C5']:>3} | {_total(before)}",
        f"  {'PoC #58l (after chain)':<28} | "
        f"{after['NP']:>6,} | {after['NE']:>6,} | "
        f"{after['C1']:>4} | {after['C2']:>4} | "
        f"{after['C4']:>4} | {after['C5']:>3} | {_total(after)}",
        "",
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
