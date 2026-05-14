"""PoC #58j: chain vertex_remove + stochastic cleanup.

PoC #58i confirmed that ``_apply_vertex_remove`` on elem 12914's
interior node (#46234, valence 7) DOES clear elem 12914 — but its
Delaunay retriangulation of the 7-vertex rim creates 2 new C1
fails and 1 new C4 fail elsewhere, taking the global count from
4 → 6 (rejected by the count-comparison gate).

The hypothesis here: those 2 new C1 + 1 new C4 are marginal
geometry artefacts of the local Delaunay, not deep structural
fails, and the stochastic local fixer should absorb most of them.
If the chain ends up below 4 we have a true zero-residual or
near-zero finishing step.

Pipeline:

  1. Read ``outputs/58e_stochastic_with_insert.14`` (4 violations).
  2. ``_apply_vertex_remove(eid=12914 interior=46234, force=True)``,
     ACCEPTING the +2 regression unconditionally.
  3. Run the move-only stochastic fixer (seed=42, identical
     parameters to PoC #58d) on the resulting mesh.
  4. Report the final residual.

The fixer functions are imported from PoC #58d's notebook script
via ``importlib`` so we don't duplicate ~500 lines.

Outputs:
    outputs/58j_chained.14
    outputs/58j_summary.{txt,json}
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
INPUT = REPO / "outputs" / "58e_stochastic_with_insert.14"
COASTLINE = (
    REPO / "data" / "coastline" / "tokyo_bay"
    / "MLIT_C23" / "C23-06_TOKYOBAY.shp"
)
FIXER_SCRIPT = REPO / "notebooks" / "58d_stochastic_local_fix.py"
OUT_DIR = REPO / "outputs"
OUTPUT = OUT_DIR / "58j_chained.14"
SUMMARY_TXT = OUT_DIR / "58j_summary.txt"
SUMMARY_JSON = OUT_DIR / "58j_summary.json"

MIN_ANGLE_TARGET = 30.0
MAX_ANGLE_TARGET = 130.0
AREA_RATIO_TARGET = 0.5
MAX_VALENCE = 8
ALPHA_TARGET = 0.95
SEED = 42


# ---------------------------------------------------------------------------
# Import PoC #58d fixer functions via importlib (numeric-prefix filename
# precludes a plain ``import``).
# ---------------------------------------------------------------------------
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


def _print_metrics(label: str, mesh: Fort14Mesh) -> dict:
    m = _metrics(mesh)
    t = m["C1"] + m["C2"] + m["C4"] + m["C5"]
    print(
        f"[58j] {label:<24}: NP={m['NP']:,} NE={m['NE']:,} "
        f"C1={m['C1']} C2={m['C2']} C4={m['C4']} C5={m['C5']} "
        f"(total {t})",
        flush=True,
    )
    return m


def _run_stochastic_fixer(mesh: Fort14Mesh, projector) -> tuple[Fort14Mesh, dict]:
    """Run the PoC #58d move-only stochastic fixer on ``mesh``.
    Returns ``(updated_mesh, stats)``."""
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
            print(f"[58j] outer {outer}: all cleared", flush=True)
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
        print(
            f"[58j] outer {outer}: {len(fail_list)} fails, "
            f"fixed={outer_fixed} stuck={outer_stuck}",
            flush=True,
        )
        if outer_fixed == 0:
            print(f"[58j] no progress in outer {outer}; terminating",
                  flush=True)
            break

    wall = time.perf_counter() - t0
    stats = {
        "wall_seconds": float(wall),
        "n_records": len(fix_records),
        "n_fixed": sum(1 for r in fix_records if r["fixed"]),
        "n_stuck": sum(1 for r in fix_records if not r["fixed"]),
    }
    return mesh, stats


def main() -> int:
    if not INPUT.exists():
        raise SystemExit(f"input missing: {INPUT}")
    if not FIXER_SCRIPT.exists():
        raise SystemExit(f"fixer script missing: {FIXER_SCRIPT}")
    OUT_DIR.mkdir(parents=True, exist_ok=True)

    mesh = read_fort14(INPUT)
    before = _print_metrics("input", mesh)

    projector = build_coastline_projector(
        [COASTLINE],
        max_snap_distance_m=500.0,
        mean_latitude_deg=float(mesh.nodes[:, 1].mean()),
    )

    # ====================================================================
    # Stage 1: force-accept vertex_remove on elem 12914 interior vertex.
    # ====================================================================
    print("[58j] stage 1: vertex_remove on elem 12914", flush=True)
    m_arr = min_interior_angle(mesh)
    c1_eids = sorted(int(e) for e in np.where(m_arr < MIN_ANGLE_TARGET)[0])
    if not c1_eids:
        print("[58j] no C1 fails to target; skipping stage 1",
              flush=True)
        stage1_metrics = before
        stage1_info: dict = {}
    else:
        bnd_node = _boundary_node_mask(mesh)
        n2e_map = _node_to_elements(mesh.elements, mesh.n_nodes)
        valence_arr = node_valence(mesh.elements, mesh.n_nodes)
        eid = c1_eids[0]
        tri = mesh.elements[int(eid)]
        interior_verts = [
            int(v) for v in tri if not bnd_node[int(v)]
        ]
        if not interior_verts:
            print(f"[58j] elem {eid}: no interior vertex; abort",
                  flush=True)
            stage1_metrics = before
            stage1_info = {}
        else:
            interior_verts.sort(key=lambda v: -int(valence_arr[v]))
            chosen_v = interior_verts[0]
            ring = n2e_map[chosen_v]
            print(
                f"[58j]   targeting elem {eid} via interior vertex "
                f"{chosen_v} (valence {int(valence_arr[chosen_v])}, "
                f"ring size {int(ring.size)})",
                flush=True,
            )
            out = _apply_vertex_remove(
                mesh, chosen_v, ring,
                alpha_target=ALPHA_TARGET,
                min_angle_target=MIN_ANGLE_TARGET,
                max_angle_target=MAX_ANGLE_TARGET,
                boundary_node_mask=bnd_node,
                force=True,
            )
            if out is None:
                print(
                    "[58j]   vertex_remove returned None — skip",
                    flush=True,
                )
                stage1_metrics = before
                stage1_info = {"applied": False}
            else:
                mesh, info = out
                stage1_info = {
                    "applied": True,
                    "vertex": int(chosen_v),
                    "rim_size": int(info["rim_size"]),
                    "n_new_elements": int(info["n_new_elements"]),
                }
                stage1_metrics = _print_metrics("stage 1 done", mesh)

    # ====================================================================
    # Stage 2: stochastic cleanup (PoC #58d fixer).
    # ====================================================================
    print("[58j] stage 2: stochastic cleanup (move-only, seed=42)",
          flush=True)
    mesh, fixer_stats = _run_stochastic_fixer(mesh, projector)
    after = _print_metrics("stage 2 done", mesh)
    write_fort14(mesh, OUTPUT)

    total_b = before["C1"] + before["C2"] + before["C4"] + before["C5"]
    total_s1 = (
        stage1_metrics["C1"] + stage1_metrics["C2"]
        + stage1_metrics["C4"] + stage1_metrics["C5"]
    )
    total_a = after["C1"] + after["C2"] + after["C4"] + after["C5"]

    payload = {
        "input": str(INPUT.resolve()),
        "output": str(OUTPUT.resolve()),
        "seed": SEED,
        "before": before,
        "stage1": stage1_metrics,
        "stage1_info": stage1_info,
        "after": after,
        "fixer_stats": fixer_stats,
    }
    SUMMARY_JSON.write_text(
        json.dumps(payload, indent=2, default=str), encoding="utf-8",
    )

    lines = [
        "PoC #58j — chained vertex_remove + stochastic cleanup",
        f"input : {INPUT.name}",
        f"output: {OUTPUT.name}",
        "",
        f"  {'stage':<28} | {'NP':>6} | {'NE':>6} | "
        f"{'C1':>4} | {'C2':>4} | {'C4':>4} | {'C5':>3} | total",
        "  " + "-" * 82,
        f"  {'PoC #58e (input)':<28} | "
        f"{before['NP']:>6,} | {before['NE']:>6,} | "
        f"{before['C1']:>4} | {before['C2']:>4} | "
        f"{before['C4']:>4} | {before['C5']:>3} | {total_b}",
        f"  {'stage 1: vertex_remove':<28} | "
        f"{stage1_metrics['NP']:>6,} | {stage1_metrics['NE']:>6,} | "
        f"{stage1_metrics['C1']:>4} | {stage1_metrics['C2']:>4} | "
        f"{stage1_metrics['C4']:>4} | {stage1_metrics['C5']:>3} | "
        f"{total_s1}",
        f"  {'stage 2: stochastic clean':<28} | "
        f"{after['NP']:>6,} | {after['NE']:>6,} | "
        f"{after['C1']:>4} | {after['C2']:>4} | "
        f"{after['C4']:>4} | {after['C5']:>3} | {total_a}",
        "",
        f"  delta input → final  : {total_a - total_b:+d}",
        f"  fixer wall (stage 2) : {fixer_stats['wall_seconds']:.2f} s",
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
