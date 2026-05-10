"""PoC #40: Stage 1 of Phase H — residual quality + operator dry-run.

The previous CLIs (`fmesh-mesh-pipeline`) check whether a mesh
satisfies *aggregate* thresholds (alpha_mean, frac<20°, max_valence,
...). The "professional FVCOM mesh" SMS workflow demands a stricter
*per-element* gate: every triangle should reach
``alpha >= alpha_target`` and ``min_angle >= min_angle_target``.

Phase H — the planned per-element greedy optimizer — would visit
each fail element and try a sequence of local-edit operators
(smooth_node, edge_swap, edge_split, edge_collapse, vertex_remove)
until one improves quality without making any neighbour worse.

This PoC measures, *before implementing Phase H*, two things:

  (a) **Residual size**: how many elements on the pipeline-rung-1
      output of PoC #19 fail the per-element gate?
  (b) **Operator coverage**: for each fail element, dry-run the two
      simplest operators (smooth_node on each of 3 vertices,
      edge_swap on each of 3 internal edges) and record whether
      *any* dry-run improves the element's quality. The fraction
      "fixable by smooth + swap alone" sets the lower bound on what
      Phase H gains from those two operators; the residual is the
      candidate territory for split / collapse / cluster-level
      re-mesh operators.

The dry-runs are read-only (each tentative coordinate / topology
change is reverted before continuing). No mesh is written; this PoC
produces only a JSON / TXT summary plus a sample of the residual
elements with full per-element context.

Outputs:
    outputs/40_phase_h_dry_run_summary.txt
    outputs/40_phase_h_dry_run_summary.json
"""
from __future__ import annotations

import json
from collections import Counter
from pathlib import Path

import numpy as np
from scipy.sparse.csgraph import connected_components

from fvcom_mesh_tools.algorithms.quality import (
    alpha_quality,
    min_interior_angle,
)
from fvcom_mesh_tools.diagnostics import face_face_adjacency
from fvcom_mesh_tools.io import Fort14Mesh, read_fort14

REPO = Path(__file__).resolve().parent.parent
INPUT = REPO / "outputs" / "33_pipeline_passing.14"
OUT_DIR = REPO / "outputs"
SUMMARY_TXT = OUT_DIR / "40_phase_h_dry_run_summary.txt"
SUMMARY_JSON = OUT_DIR / "40_phase_h_dry_run_summary.json"

ALPHA_TARGET = 0.95
MIN_ANGLE_TARGET = 20.0
N_SAMPLE_RESIDUAL = 20  # how many "unfixable" elements to dump into JSON


# ---------------------------------------------------------------------------
# Per-element penalty / quality helpers (read-only, vectorised where useful)
# ---------------------------------------------------------------------------


def _per_element_alpha_minangle(
    nodes: np.ndarray, elements: np.ndarray,
) -> tuple[np.ndarray, np.ndarray]:
    """``alpha`` and ``min_angle_deg`` arrays of shape ``(NE,)``."""
    mesh = Fort14Mesh(
        title="dryrun", nodes=nodes,
        depths=np.zeros(nodes.shape[0]),
        elements=elements,
        open_boundaries=[], land_boundaries=[],
    )
    return alpha_quality(mesh), min_interior_angle(mesh)


def _is_fail(alpha: np.ndarray, min_ang: np.ndarray) -> np.ndarray:
    return (alpha < ALPHA_TARGET) | (min_ang < MIN_ANGLE_TARGET)


def _penalty(alpha: np.ndarray, min_ang: np.ndarray) -> np.ndarray:
    """Element-level penalty used to compare 'before/after' a dry-run."""
    a_pen = np.maximum(0.0, ALPHA_TARGET - alpha) ** 2
    g_pen = np.maximum(0.0, MIN_ANGLE_TARGET - min_ang) ** 2 / 100.0
    return a_pen + g_pen


# ---------------------------------------------------------------------------
# Dry-run smooth_node
# ---------------------------------------------------------------------------


def _try_smooth_node(
    mesh: Fort14Mesh, vertex_id: int, ring_elem_ids: np.ndarray,
    boundary_node_mask: np.ndarray,
) -> dict:
    """Move ``vertex_id`` to the centroid of its 1-ring neighbour nodes
    (excluding itself), recompute alpha + min_angle on the affected
    elements, and report the per-element penalty before / after.
    """
    if boundary_node_mask[vertex_id]:
        return {"applicable": False, "reason": "boundary node"}

    neighbours = np.unique(mesh.elements[ring_elem_ids].ravel())
    neighbours = neighbours[neighbours != vertex_id]
    if neighbours.size == 0:
        return {"applicable": False, "reason": "no neighbours"}

    new_xy = mesh.nodes[neighbours].mean(axis=0)
    nodes_proposed = mesh.nodes.copy()
    nodes_proposed[vertex_id] = new_xy

    alpha_before, ang_before = _per_element_alpha_minangle(
        mesh.nodes, mesh.elements[ring_elem_ids],
    )
    alpha_after, ang_after = _per_element_alpha_minangle(
        nodes_proposed, mesh.elements[ring_elem_ids],
    )
    pen_before = _penalty(alpha_before, ang_before).sum()
    pen_after = _penalty(alpha_after, ang_after).sum()

    p0 = nodes_proposed[mesh.elements[ring_elem_ids, 0]]
    p1 = nodes_proposed[mesh.elements[ring_elem_ids, 1]]
    p2 = nodes_proposed[mesh.elements[ring_elem_ids, 2]]
    cross = (
        (p1[:, 0] - p0[:, 0]) * (p2[:, 1] - p0[:, 1])
        - (p1[:, 1] - p0[:, 1]) * (p2[:, 0] - p0[:, 0])
    )
    if (cross <= 0).any():
        return {"applicable": False, "reason": "would flip"}

    return {
        "applicable": True,
        "penalty_before": float(pen_before),
        "penalty_after": float(pen_after),
        "improves": pen_after + 1e-12 < pen_before,
        "alpha_after_min": float(alpha_after.min()),
        "ang_after_min": float(ang_after.min()),
    }


# ---------------------------------------------------------------------------
# Dry-run edge_swap
# ---------------------------------------------------------------------------


def _edge_buddy(
    elements: np.ndarray, full_adj, elem_id: int, edge_local: int,
) -> tuple[int, int] | None:
    """Return ``(buddy_elem_id, buddy_local_third_vertex_id)`` or
    ``None`` if the edge is on the boundary.
    """
    nbrs_csr = full_adj[elem_id]
    nbr_ids = nbrs_csr.indices
    a = elements[elem_id, edge_local]
    b = elements[elem_id, (edge_local + 1) % 3]
    target = {int(a), int(b)}
    for nb in nbr_ids:
        nb_set = set(int(x) for x in elements[nb])
        if target.issubset(nb_set):
            third = (nb_set - target).pop()
            return int(nb), int(third)
    return None


def _try_edge_swap(
    mesh: Fort14Mesh, full_adj, elem_id: int, edge_local: int,
    boundary_edge_set: set,
) -> dict:
    """Swap the shared edge between ``elem_id`` and its buddy across
    ``edge_local``. Returns dry-run penalty before / after.
    """
    a = int(mesh.elements[elem_id, edge_local])
    b = int(mesh.elements[elem_id, (edge_local + 1) % 3])
    if (min(a, b), max(a, b)) in boundary_edge_set:
        return {"applicable": False, "reason": "boundary edge"}

    buddy = _edge_buddy(mesh.elements, full_adj, elem_id, edge_local)
    if buddy is None:
        return {"applicable": False, "reason": "no buddy element"}
    buddy_id, fourth = buddy
    third = int(mesh.elements[elem_id, (edge_local + 2) % 3])

    # Two-element block before the swap.
    block_before = mesh.elements[[elem_id, buddy_id]]
    # After swapping (a, b) -> (third, fourth), reform two triangles:
    #   tri1 = (a, third, fourth)
    #   tri2 = (b, fourth, third)
    block_after = np.array(
        [[a, third, fourth], [b, fourth, third]], dtype=mesh.elements.dtype,
    )

    alpha_before, ang_before = _per_element_alpha_minangle(
        mesh.nodes, block_before,
    )
    alpha_after, ang_after = _per_element_alpha_minangle(
        mesh.nodes, block_after,
    )
    pen_before = _penalty(alpha_before, ang_before).sum()
    pen_after = _penalty(alpha_after, ang_after).sum()

    p0 = mesh.nodes[block_after[:, 0]]
    p1 = mesh.nodes[block_after[:, 1]]
    p2 = mesh.nodes[block_after[:, 2]]
    cross = (
        (p1[:, 0] - p0[:, 0]) * (p2[:, 1] - p0[:, 1])
        - (p1[:, 1] - p0[:, 1]) * (p2[:, 0] - p0[:, 0])
    )
    if (cross <= 0).any():
        return {"applicable": False, "reason": "would flip"}

    return {
        "applicable": True,
        "penalty_before": float(pen_before),
        "penalty_after": float(pen_after),
        "improves": pen_after + 1e-12 < pen_before,
        "alpha_after_min": float(alpha_after.min()),
        "ang_after_min": float(ang_after.min()),
    }


# ---------------------------------------------------------------------------
# Boundary masks
# ---------------------------------------------------------------------------


def _boundary_node_mask(mesh: Fort14Mesh) -> np.ndarray:
    mask = np.zeros(mesh.n_nodes, dtype=bool)
    for seg in mesh.open_boundaries:
        mask[np.asarray(seg, dtype=np.int64)] = True
    for _ib, seg in mesh.land_boundaries:
        mask[np.asarray(seg, dtype=np.int64)] = True
    return mask


def _boundary_edge_set(mesh: Fort14Mesh) -> set:
    """Edges used by exactly one element are on the geometric boundary."""
    e = np.vstack([
        mesh.elements[:, [0, 1]],
        mesh.elements[:, [1, 2]],
        mesh.elements[:, [2, 0]],
    ])
    e_sorted = np.sort(e, axis=1)
    keys = (
        e_sorted[:, 0].astype(np.int64) << 32
    ) | e_sorted[:, 1].astype(np.int64)
    uniq, counts = np.unique(keys, return_counts=True)
    one_use = uniq[counts == 1]
    out = set()
    for k in one_use:
        a = int(k >> 32)
        b = int(k & 0xFFFFFFFF)
        out.add((min(a, b), max(a, b)))
    return out


# ---------------------------------------------------------------------------
# Driver
# ---------------------------------------------------------------------------


def _classify_element(
    elem_id: int, mesh: Fort14Mesh, full_adj, ring_by_node,
    boundary_node_mask, boundary_edge_set,
) -> dict:
    """Run smooth + swap dry-runs on ``elem_id``. Return per-operator
    result + best operator label."""
    record = {
        "elem_id": int(elem_id),
        "vertices": [int(v) for v in mesh.elements[elem_id]],
        "boundary_touch": bool(
            boundary_node_mask[mesh.elements[elem_id]].any()
        ),
        "smooth_attempts": [],
        "swap_attempts": [],
    }

    for k_local, v in enumerate(mesh.elements[elem_id]):
        ring = ring_by_node.get(int(v), np.array([], dtype=np.int64))
        if ring.size == 0:
            record["smooth_attempts"].append(
                {"vertex_local": k_local, "applicable": False,
                 "reason": "empty 1-ring"}
            )
            continue
        out = _try_smooth_node(
            mesh, int(v), ring, boundary_node_mask,
        )
        out["vertex_local"] = k_local
        record["smooth_attempts"].append(out)

    for k_edge in range(3):
        out = _try_edge_swap(
            mesh, full_adj, int(elem_id), k_edge, boundary_edge_set,
        )
        out["edge_local"] = k_edge
        record["swap_attempts"].append(out)

    smooth_improves = any(a.get("improves", False)
                          for a in record["smooth_attempts"])
    swap_improves = any(a.get("improves", False)
                        for a in record["swap_attempts"])
    if smooth_improves and swap_improves:
        record["best_operator"] = "either"
    elif smooth_improves:
        record["best_operator"] = "smooth"
    elif swap_improves:
        record["best_operator"] = "swap"
    else:
        record["best_operator"] = "none"
    return record


def _ring_by_node(elements: np.ndarray, n_nodes: int
                  ) -> dict[int, np.ndarray]:
    """Map node-id → array of incident element ids."""
    rows = elements.ravel()
    cols = np.tile(np.arange(elements.shape[0]), 3)
    order = np.argsort(rows)
    rows = rows[order]
    cols = cols[order]
    boundaries = np.searchsorted(rows, np.arange(n_nodes + 1))
    out: dict[int, np.ndarray] = {}
    for n in range(n_nodes):
        s, e = boundaries[n], boundaries[n + 1]
        if s < e:
            out[n] = cols[s:e]
    return out


def main() -> int:
    OUT_DIR.mkdir(parents=True, exist_ok=True)
    if not INPUT.exists():
        raise SystemExit(f"input missing: {INPUT}")

    mesh = read_fort14(INPUT)
    print(f"[40] input: NP={mesh.n_nodes:,}  NE={mesh.n_elements:,}")
    alpha = alpha_quality(mesh)
    ang = min_interior_angle(mesh)
    fail_mask = _is_fail(alpha, ang)
    fail_idx = np.where(fail_mask)[0]
    print(
        f"[40] residual: {fail_idx.size:,} fail elements "
        f"({fail_idx.size / mesh.n_elements:.4%}) at thresholds "
        f"alpha>={ALPHA_TARGET}  min_angle>={MIN_ANGLE_TARGET}°"
    )

    if fail_idx.size == 0:
        SUMMARY_TXT.write_text(
            f"PoC #40: 0 fail elements at the gate "
            f"(alpha>={ALPHA_TARGET}, min_angle>={MIN_ANGLE_TARGET}°).\n",
            encoding="utf-8",
        )
        SUMMARY_JSON.write_text(json.dumps({"n_fail": 0}, indent=2),
                                encoding="utf-8")
        return 0

    full_adj = face_face_adjacency(mesh.elements)
    ring_by_node = _ring_by_node(mesh.elements, mesh.n_nodes)
    bnd_node = _boundary_node_mask(mesh)
    bnd_edge = _boundary_edge_set(mesh)
    print(
        f"[40] boundary nodes: {int(bnd_node.sum()):,}  "
        f"boundary edges: {len(bnd_edge):,}"
    )

    # Optional clustering of fail elements via face-face adjacency restricted
    # to the fail subgraph — diagnoses "isolated sliver" vs "boundary band".
    sub = full_adj[fail_idx][:, fail_idx]
    n_fail_comp, _ = connected_components(sub, directed=False, return_labels=True)

    print("[40] dry-run smooth + swap on every fail element ...")
    classifications: list[dict] = []
    bins: Counter = Counter()
    for n, eid in enumerate(fail_idx, start=1):
        rec = _classify_element(
            int(eid), mesh, full_adj, ring_by_node, bnd_node, bnd_edge,
        )
        rec["alpha"] = float(alpha[eid])
        rec["min_angle_deg"] = float(ang[eid])
        classifications.append(rec)
        bins[rec["best_operator"]] += 1
        if n % 50 == 0:
            print(f"[40]   processed {n}/{fail_idx.size}")

    n_total = len(classifications)
    n_smooth = bins["smooth"]
    n_swap = bins["swap"]
    n_either = bins["either"]
    n_none = bins["none"]
    n_fixable = n_smooth + n_swap + n_either

    n_boundary_touching_fail = sum(
        1 for r in classifications if r["boundary_touch"]
    )
    n_fail_unfixable_on_boundary = sum(
        1 for r in classifications
        if r["boundary_touch"] and r["best_operator"] == "none"
    )

    out_lines = [
        "PoC #40 Phase H Stage 1 — residual + dry-run analysis",
        f"input: {INPUT.name}",
        f"thresholds: alpha>={ALPHA_TARGET}  min_angle>={MIN_ANGLE_TARGET}°",
        "",
        f"  NE                              : {mesh.n_elements:,}",
        f"  fail elements                   : {n_total:,}  "
        f"({n_total / mesh.n_elements:.4%})",
        f"  fail components (face-face adj) : {int(n_fail_comp):,}",
        f"  fail elements touching boundary : {n_boundary_touching_fail:,}",
        "",
        "  per-element operator coverage",
        f"    fixable by smooth only        : {n_smooth:,}",
        f"    fixable by swap only          : {n_swap:,}",
        f"    fixable by both               : {n_either:,}",
        f"    unfixable by smooth or swap   : {n_none:,}",
        f"    fixable total                 : {n_fixable:,}  "
        f"({n_fixable / max(n_total, 1):.2%})",
        "",
        f"  boundary-touching unfixable     : "
        f"{n_fail_unfixable_on_boundary:,}",
    ]
    print("\n".join(out_lines))

    # Sample of unfixable elements for manual inspection.
    unfixable = [r for r in classifications if r["best_operator"] == "none"]
    sample = unfixable[:N_SAMPLE_RESIDUAL]
    payload = {
        "thresholds": {"alpha": ALPHA_TARGET, "min_angle_deg": MIN_ANGLE_TARGET},
        "input": str(INPUT.resolve()),
        "input_NP": int(mesh.n_nodes),
        "input_NE": int(mesh.n_elements),
        "n_fail": int(n_total),
        "n_fail_components": int(n_fail_comp),
        "n_boundary_touching_fail": int(n_boundary_touching_fail),
        "operator_coverage": {
            "smooth_only": int(n_smooth),
            "swap_only": int(n_swap),
            "either": int(n_either),
            "none": int(n_none),
            "fixable_total": int(n_fixable),
            "fixable_fraction": float(n_fixable / max(n_total, 1)),
        },
        "unfixable_sample": sample,
        "all_classifications": classifications,
    }
    SUMMARY_JSON.write_text(json.dumps(payload, indent=2, default=str),
                            encoding="utf-8")
    SUMMARY_TXT.write_text("\n".join(out_lines) + "\n", encoding="utf-8")
    print(f"\nwrote {SUMMARY_TXT}")
    print(f"wrote {SUMMARY_JSON}")
    return 0


if __name__ == "__main__":
    raise SystemExit(main())
