"""PoC #58i: targeted vertex_remove on stuck C1 fail elem 12914.

PoC #58g identified elem 12914 (verts [21227, 1847, 46234],
min_ang = 28.0°, valences = [2, 2, 7]) as the structurally stuck
C1 residual that the stochastic local fixer cannot clear across
all 5 tested seeds — moving any vertex of the triangle leaves the
thin angle in place because two of the three vertices are boundary-
locked.

SMS manual practice would: delete the **interior obtuse vertex**
(here node 46234, valence 7) and let Delaunay retriangulate the
7-element 1-ring polygon into 5 sub-triangles. ``mesh_clean_phase_h``
already exposes :func:`_apply_vertex_remove` as a single-element
variant of the Stage 2 medial-axis re-mesh, which does exactly
this.  We just call it directly on the right node — bypassing the
existing Pass B / C dispatching that would never propose this
vertex because the per-element fail enumeration starts from the
element, not the obtuse vertex.

Pipeline:

  1. Read ``outputs/58e_stochastic_with_insert.14`` (1 C1 + 3 C4).
  2. Locate every current C1 fail element.  For each, identify the
     interior vertex (the one not on the boundary mask) — that is
     the candidate for ``vertex_remove``.
  3. Call ``_apply_vertex_remove(mesh, v, ring, force=True)`` —
     ``force=True`` bypasses the local penalty gate (penalty is
     based on alpha + min/max angle, which won't strictly drop
     when the residual is marginal).
  4. Validate the new mesh **globally**: every FVCOM criterion
     (C1 / C2 / C4 / C5) must be ≤ its current value.  If yes,
     keep; if no, revert.
  5. Run the stochastic local fixer (PoC #58e style) over any new
     residual to absorb perturbations.

Outputs:
    outputs/58i_vertex_remove_result.14
    outputs/58i_summary.{txt,json}
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
    _apply_vertex_remove,
    _boundary_node_mask,
    _node_to_elements,
    _per_edge_area_change,
    _per_element_quality,
)

REPO = Path(__file__).resolve().parent.parent
INPUT = REPO / "outputs" / "58e_stochastic_with_insert.14"
OUT_DIR = REPO / "outputs"
OUTPUT = OUT_DIR / "58i_vertex_remove_result.14"
SUMMARY_TXT = OUT_DIR / "58i_summary.txt"
SUMMARY_JSON = OUT_DIR / "58i_summary.json"

MIN_ANGLE_TARGET = 30.0
MAX_ANGLE_TARGET = 130.0
AREA_RATIO_TARGET = 0.5
MAX_VALENCE = 8
ALPHA_TARGET = 0.95


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


def main() -> int:
    if not INPUT.exists():
        raise SystemExit(f"input missing: {INPUT}")
    OUT_DIR.mkdir(parents=True, exist_ok=True)

    mesh = read_fort14(INPUT)
    before = _metrics(mesh)
    total_b = before["C1"] + before["C2"] + before["C4"] + before["C5"]
    print(
        f"[58i] input : NP={before['NP']:,} NE={before['NE']:,} "
        f"C1={before['C1']} C2={before['C2']} "
        f"C4={before['C4']} C5={before['C5']} (total {total_b})",
        flush=True,
    )

    # Identify current C1 fails.
    m = min_interior_angle(mesh)
    c1_mask = m < MIN_ANGLE_TARGET
    c1_eids = sorted(int(e) for e in np.where(c1_mask)[0])
    print(f"[58i] {len(c1_eids)} C1 fail elements: {c1_eids}",
          flush=True)

    bnd_node = _boundary_node_mask(mesh)
    records: list[dict] = []
    n_accepted = 0

    for eid in c1_eids:
        tri = mesh.elements[int(eid)]
        verts = [int(v) for v in tri]
        # The interior vertex is the one NOT on the boundary mask.
        interior_verts = [v for v in verts if not bnd_node[v]]
        if not interior_verts:
            print(
                f"[58i]   elem {eid}: all three vertices boundary — "
                "vertex_remove cannot help here",
                flush=True,
            )
            records.append({
                "eid": int(eid), "interior_vertex": None,
                "result": "all_boundary_skip",
            })
            continue
        # If multiple interior vertices, pick the one with highest
        # valence (most likely to be the "extra" node SMS would remove).
        n2e_map = _node_to_elements(mesh.elements, mesh.n_nodes)
        valence_arr = node_valence(mesh.elements, mesh.n_nodes)
        interior_verts.sort(key=lambda v: -int(valence_arr[v]))
        chosen_v = int(interior_verts[0])
        ring = n2e_map.get(chosen_v)
        print(
            f"[58i]   elem {eid}: interior vertex = {chosen_v} "
            f"(valence {int(valence_arr[chosen_v])}, "
            f"ring size = {0 if ring is None else int(ring.size)})",
            flush=True,
        )

        if ring is None or ring.size < 3:
            records.append({
                "eid": int(eid), "interior_vertex": chosen_v,
                "result": "ring_too_small",
            })
            continue

        t0 = time.perf_counter()
        out = _apply_vertex_remove(
            mesh, chosen_v, ring,
            alpha_target=ALPHA_TARGET,
            min_angle_target=MIN_ANGLE_TARGET,
            max_angle_target=MAX_ANGLE_TARGET,
            boundary_node_mask=bnd_node,
            force=True,
        )
        dt = time.perf_counter() - t0
        if out is None:
            print(
                f"[58i]   elem {eid}: _apply_vertex_remove returned "
                "None (boundary / degenerate / Delaunay failed)",
                flush=True,
            )
            records.append({
                "eid": int(eid), "interior_vertex": chosen_v,
                "result": "vertex_remove_returned_none",
            })
            continue
        candidate, info = out
        cand_metrics = _metrics(candidate)
        cand_total = (
            cand_metrics["C1"] + cand_metrics["C2"]
            + cand_metrics["C4"] + cand_metrics["C5"]
        )
        cur_total = (
            before["C1"] + before["C2"] + before["C4"] + before["C5"]
            if n_accepted == 0
            else _metrics(mesh)["C1"] + _metrics(mesh)["C2"]
                 + _metrics(mesh)["C4"] + _metrics(mesh)["C5"]
        )
        cur_metrics = _metrics(mesh)

        # Acceptance: every criterion must not regress, AND total
        # strictly decreases.
        regressed = False
        for k in ("C1", "C2", "C4", "C5"):
            if cand_metrics[k] > cur_metrics[k]:
                regressed = True
                break
        if regressed or cand_total >= cur_total:
            print(
                f"[58i]   elem {eid}: vertex_remove rejected — "
                f"after {cand_metrics} vs before {cur_metrics} "
                f"(total {cand_total} vs {cur_total})",
                flush=True,
            )
            records.append({
                "eid": int(eid), "interior_vertex": chosen_v,
                "result": "regression_reject",
                "after_metrics": cand_metrics,
                "wall_seconds": dt,
            })
            continue

        # Accept.
        mesh = candidate
        n_accepted += 1
        print(
            f"[58i]   elem {eid}: vertex_remove ACCEPTED — "
            f"after {cand_metrics} (wall {dt:.3f} s, "
            f"rim_size={info['rim_size']}, "
            f"new_elements={info['n_new_elements']})",
            flush=True,
        )
        records.append({
            "eid": int(eid), "interior_vertex": chosen_v,
            "result": "accepted",
            "after_metrics": cand_metrics,
            "wall_seconds": dt,
            "rim_size": int(info["rim_size"]),
            "removed_elements": list(info["removed_elements"]),
            "n_new_elements": int(info["n_new_elements"]),
        })

    after = _metrics(mesh)
    total_a = after["C1"] + after["C2"] + after["C4"] + after["C5"]
    write_fort14(mesh, OUTPUT)
    print(
        f"[58i] output: NP={after['NP']:,} NE={after['NE']:,} "
        f"C1={after['C1']} C2={after['C2']} "
        f"C4={after['C4']} C5={after['C5']} (total {total_a})",
        flush=True,
    )

    payload = {
        "input": str(INPUT.resolve()),
        "output": str(OUTPUT.resolve()),
        "before": before,
        "after": after,
        "delta": {k: after[k] - before[k] for k in before},
        "n_c1_targets": len(c1_eids),
        "n_accepted": n_accepted,
        "records": records,
    }
    SUMMARY_JSON.write_text(
        json.dumps(payload, indent=2, default=str), encoding="utf-8",
    )

    lines = [
        "PoC #58i — targeted vertex_remove on stuck C1 fails",
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
        f"  {'PoC #58i (after fixer)':<28} | "
        f"{after['NP']:>6,} | {after['NE']:>6,} | "
        f"{after['C1']:>4} | {after['C2']:>4} | "
        f"{after['C4']:>4} | {after['C5']:>3} | {total_a}",
        "",
        f"  C1 targets attempted : {len(c1_eids)}",
        f"  n_accepted           : {n_accepted}",
        f"  delta total          : {total_a - total_b:+d}",
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
