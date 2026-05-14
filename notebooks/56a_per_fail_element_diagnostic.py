"""PoC #56a: per-fail-element diagnostic on PoC #54c residuals.

PoC #54d identified the residual 77 violations after PoC #54c as
100 % boundary-driven (C1 8 singletons + C4 69 mostly boundary-
endpoint). Before designing Pass F we need to know, per fail
element, exactly:

  * which 3 vertices it owns (positions, lat/lon)
  * which edges are on the coastline (boundary), which are interior
  * its 1-ring neighbours (adjacent elements)
  * its angle distribution (min / max / which vertex)
  * what each candidate operator would change locally:
      - edge_split_boundary on each boundary edge
      - edge_split_interior on the longest interior edge
      - edge_swap on each interior edge
  * for each candidate, simulate the FVCOM fail count in the
    affected block BEFORE vs AFTER (count-comparison gate)

This is read-only and runs on the login node in < 30 s.

The output table tells us, per fail element, whether at least ONE
operator under a count-comparison gate could reduce its local
fail count — i.e., whether Pass F has a chance.

Output:
    outputs/56a_per_fail_diagnostic.txt
    outputs/56a_per_fail_diagnostic.json
"""
from __future__ import annotations

import json
from pathlib import Path

import numpy as np

from fvcom_mesh_tools.algorithms.quality import min_interior_angle
from fvcom_mesh_tools.diagnostics import node_valence
from fvcom_mesh_tools.io import read_fort14
from fvcom_mesh_tools.mesh_clean_phase_h import (
    DEFAULT_MAX_VALENCE,
    _apply_edge_split_boundary,
    _apply_edge_split_interior,
    _apply_edge_swap,
    _boundary_topology,
    _edge_use_counts,
    _per_edge_area_change,
    _per_element_quality,
    build_coastline_projector,
)

REPO = Path(__file__).resolve().parent.parent
INPUT = REPO / "outputs" / "54c_phase_h_optimized.14"
COASTLINE = (
    REPO / "data" / "coastline" / "tokyo_bay"
    / "MLIT_C23" / "C23-06_TOKYOBAY.shp"
)
OUT_DIR = REPO / "outputs"
SUMMARY_TXT = OUT_DIR / "56a_per_fail_diagnostic.txt"
SUMMARY_JSON = OUT_DIR / "56a_per_fail_diagnostic.json"

MIN_ANGLE_TARGET = 30.0
MAX_ANGLE_TARGET = 130.0
AREA_RATIO_TARGET = 0.5
MAX_VALENCE = DEFAULT_MAX_VALENCE


def _max_interior_angle(mesh) -> np.ndarray:
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


def _global_counts(mesh):
    m = min_interior_angle(mesh)
    M = _max_interior_angle(mesh)
    _u, _p, ac = _per_edge_area_change(mesh.nodes, mesh.elements)
    val = node_valence(mesh.elements, mesh.n_nodes)
    return {
        "C1": int((m < MIN_ANGLE_TARGET).sum()),
        "C2": int((M > MAX_ANGLE_TARGET).sum()),
        "C4": int((ac > AREA_RATIO_TARGET).sum()),
        "C5": int((val > MAX_VALENCE).sum()),
    }


def _try_operator(label, mesh, before_counts, op):
    """Run ``op(force=True)`` and report the global FVCOM count
    delta. ``op`` is a thunk returning ``(new_mesh, info) | None``.
    """
    out = op()
    if out is None:
        return {"op": label, "result": "rejected_validity"}
    new_mesh, info = out
    after_counts = _global_counts(new_mesh)
    return {
        "op": label,
        "result": "applied",
        "before": before_counts,
        "after": after_counts,
        "delta": {
            k: after_counts[k] - before_counts[k]
            for k in ("C1", "C2", "C4", "C5")
        },
        "info_keys": list(info.keys()),
    }


def main() -> int:
    if not INPUT.exists():
        raise SystemExit(f"input missing: {INPUT}")

    mesh = read_fort14(INPUT)
    print(f"[56a] input: NP={mesh.n_nodes:,} NE={mesh.n_elements:,}",
          flush=True)

    # Global compliance and per-element angle metrics.
    m = min_interior_angle(mesh)
    M = _max_interior_angle(mesh)
    a_q, min_ang, max_ang = _per_element_quality(mesh.nodes, mesh.elements)
    edge_uv, elem_pair, ac = _per_edge_area_change(
        mesh.nodes, mesh.elements,
    )
    val = node_valence(mesh.elements, mesh.n_nodes)

    c1_mask = m < MIN_ANGLE_TARGET
    c2_mask = M > MAX_ANGLE_TARGET
    c4_edge_mask = ac > AREA_RATIO_TARGET

    print(
        f"[56a] global: C1={int(c1_mask.sum())} "
        f"C2={int(c2_mask.sum())} "
        f"C4={int(c4_edge_mask.sum())} "
        f"C5={int((val > 8).sum())}",
        flush=True,
    )

    edge_uses = _edge_use_counts(mesh.elements)
    boundary_edge_keys = {k for k, v in edge_uses.items() if len(v) == 1}
    _bp, _bn, edge_to_segment = _boundary_topology(mesh)

    projector = build_coastline_projector(
        [COASTLINE],
        max_snap_distance_m=500.0,
        mean_latitude_deg=float(mesh.nodes[:, 1].mean()),
    )

    global_before = _global_counts(mesh)

    # Compute set of C4 fail elements (elements touching any C4 fail edge).
    c4_fail_edges = edge_uv[c4_edge_mask]
    c4_fail_pairs = elem_pair[c4_edge_mask]
    c4_fail_elem_set = set(int(e) for e in c4_fail_pairs.ravel())

    c1_eids = sorted(int(e) for e in np.where(c1_mask)[0])
    print(
        f"[56a] inspecting {len(c1_eids)} C1 fail elements", flush=True,
    )

    per_fail_records = []
    for eid in c1_eids:
        verts = [int(v) for v in mesh.elements[eid]]
        coords = [
            [float(mesh.nodes[v, 0]), float(mesh.nodes[v, 1])]
            for v in verts
        ]
        # Edge classification on this element
        edge_records = []
        for el in range(3):
            a = verts[el]
            b = verts[(el + 1) % 3]
            k = (min(a, b), max(a, b))
            buds = edge_uses.get(k, [])
            edge_records.append({
                "local": el,
                "uv": [a, b],
                "boundary": k in boundary_edge_keys,
                "n_incident": len(buds),
                "length_deg": float(
                    np.linalg.norm(mesh.nodes[a] - mesh.nodes[b]),
                ),
            })
        rec = {
            "elem_id": eid,
            "vertices": verts,
            "coords_lonlat": coords,
            "min_angle": float(min_ang[eid]),
            "max_angle": float(max_ang[eid]),
            "alpha": float(a_q[eid]),
            "valence_per_vertex": [int(val[v]) for v in verts],
            "edges": edge_records,
            "operator_trials": [],
        }

        # Try each candidate operator with force=True and measure
        # the global FVCOM count delta.
        # 1. edge_split_boundary on each boundary edge
        for el_rec in edge_records:
            if not el_rec["boundary"]:
                continue
            el = el_rec["local"]
            rec["operator_trials"].append(_try_operator(
                f"edge_split_boundary(el={el})",
                mesh, global_before,
                lambda mesh=mesh, eid=eid, el=el: _apply_edge_split_boundary(
                    mesh, eid, el,
                    alpha_target=0.95,
                    min_angle_target=MIN_ANGLE_TARGET,
                    max_angle_target=MAX_ANGLE_TARGET,
                    edge_uses=edge_uses,
                    edge_to_segment=edge_to_segment,
                    coastline_projector=projector,
                    force=True,
                ),
            ))
        # 2. edge_split_interior on each interior edge
        for el_rec in edge_records:
            if el_rec["boundary"]:
                continue
            el = el_rec["local"]
            rec["operator_trials"].append(_try_operator(
                f"edge_split_interior(el={el})",
                mesh, global_before,
                lambda mesh=mesh, eid=eid, el=el: _apply_edge_split_interior(
                    mesh, eid, el,
                    alpha_target=0.95,
                    min_angle_target=MIN_ANGLE_TARGET,
                    max_angle_target=MAX_ANGLE_TARGET,
                    edge_uses=edge_uses,
                    boundary_edge_keys=boundary_edge_keys,
                    force=True,
                ),
            ))
        # 3. edge_swap on each interior edge
        for el_rec in edge_records:
            if el_rec["boundary"]:
                continue
            el = el_rec["local"]
            rec["operator_trials"].append(_try_operator(
                f"edge_swap(el={el})",
                mesh, global_before,
                lambda mesh=mesh, eid=eid, el=el: _apply_edge_swap(
                    mesh, eid, el,
                    alpha_target=0.95,
                    min_angle_target=MIN_ANGLE_TARGET,
                    max_angle_target=MAX_ANGLE_TARGET,
                    edge_uses=edge_uses,
                    boundary_edge_keys=boundary_edge_keys,
                    force=True,
                ),
            ))

        # Identify the BEST candidate: largest C1+C2+C4 reduction
        # under the count-comparison gate (no C5 / flipped regressions).
        best = None
        for trial in rec["operator_trials"]:
            if trial["result"] != "applied":
                continue
            d = trial["delta"]
            if d["C5"] > 0:
                continue  # C5 regression
            score = d["C1"] + d["C2"] + d["C4"]
            if best is None or score < best[0]:
                best = (score, trial)
        rec["best_op"] = best[1] if best else None
        rec["best_delta_combined_c1_c2_c4"] = best[0] if best else None

        per_fail_records.append(rec)

    # Summary tally
    n_fixable = sum(
        1 for r in per_fail_records
        if r["best_delta_combined_c1_c2_c4"] is not None
        and r["best_delta_combined_c1_c2_c4"] < 0
    )
    n_neutral = sum(
        1 for r in per_fail_records
        if r["best_delta_combined_c1_c2_c4"] is not None
        and r["best_delta_combined_c1_c2_c4"] == 0
    )
    n_no_op_valid = sum(
        1 for r in per_fail_records
        if r["best_op"] is None
    )

    payload = {
        "input": str(INPUT.resolve()),
        "global_counts": global_before,
        "n_c1_fails_examined": len(c1_eids),
        "n_fixable_under_count_gate": n_fixable,
        "n_neutral_under_count_gate": n_neutral,
        "n_no_valid_operator": n_no_op_valid,
        "per_fail": per_fail_records,
    }
    SUMMARY_JSON.write_text(
        json.dumps(payload, indent=2, default=str), encoding="utf-8",
    )

    # Text report
    lines = [
        "PoC #56a — per-fail-element diagnostic on PoC #54c residuals",
        f"input: {INPUT.name}",
        f"global counts: {global_before}",
        "",
        f"  C1 fail elements: {len(c1_eids)}",
        f"  fixable by some operator (count gate, no C5 regression): "
        f"{n_fixable}",
        f"  neutral (count unchanged): {n_neutral}",
        f"  no valid operator at all: {n_no_op_valid}",
        "",
        "Per-element details:",
    ]
    for rec in per_fail_records:
        bnd_count = sum(1 for er in rec["edges"] if er["boundary"])
        lines.append(
            f"  --- elem {rec['elem_id']} | "
            f"verts={rec['vertices']} | "
            f"min_ang={rec['min_angle']:.1f}° | "
            f"max_ang={rec['max_angle']:.1f}° | "
            f"alpha={rec['alpha']:.3f} | "
            f"bnd_edges={bnd_count}/3 | "
            f"valences={rec['valence_per_vertex']}"
        )
        for er in rec["edges"]:
            tag = "BND" if er["boundary"] else "INT"
            lines.append(
                f"    edge {er['local']}: {tag} {er['uv']} "
                f"len={er['length_deg']:.4e}°"
            )
        for trial in rec["operator_trials"]:
            if trial["result"] == "rejected_validity":
                lines.append(f"    [{trial['op']}] rejected (validity)")
            else:
                d = trial["delta"]
                lines.append(
                    f"    [{trial['op']}] applied  "
                    f"ΔC1={d['C1']:+}  ΔC2={d['C2']:+}  "
                    f"ΔC4={d['C4']:+}  ΔC5={d['C5']:+}"
                )
        if rec["best_op"] is not None:
            lines.append(
                f"    BEST: {rec['best_op']['op']}  combined Δ="
                f"{rec['best_delta_combined_c1_c2_c4']:+}"
            )
        else:
            lines.append("    BEST: (no valid operator)")
    SUMMARY_TXT.write_text("\n".join(lines) + "\n", encoding="utf-8")
    print()
    print(f"wrote {SUMMARY_TXT}")
    print(f"wrote {SUMMARY_JSON}")
    return 0


if __name__ == "__main__":
    raise SystemExit(main())
