"""PoC #58a: per-fail diagnostic on PoC #57 Stage 1 output (Pass F only).

PoC #57 Stage 1 dropped C4 from 68 to 48 (Pass F clearing 20 in
three sweeps) while leaving C1 = 2 untouched. Before designing the
next operator (C4-aware proposal, Pass G C1-aware smoothing, or
relaxed gate) we need to know, on the new input
``57_phase_h_pass_f_only.14``:

  * whether the 2 remaining C1 fails are still the same elements as
    PoC #56c (Pass F is C4-only so they should be, but verify)
  * how the 48 C4 fails now distribute spatially — has Pass F's
    smoothing shifted the cluster topology, or is the residual a
    cleaner subset of the original 68?
  * for each cluster, whether the residual is now boundary-corner
    bound (where Pass F's tangent path skips), making the next move
    necessarily a topology-changing one

This is read-only and runs on the login node in < 30 s.

Outputs:
    outputs/58a_residual_diagnostic.txt
    outputs/58a_residual_diagnostic.json
"""
from __future__ import annotations

import json
from collections import defaultdict, deque
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
INPUT = REPO / "outputs" / "57_phase_h_pass_f_only.14"
COASTLINE = (
    REPO / "data" / "coastline" / "tokyo_bay"
    / "MLIT_C23" / "C23-06_TOKYOBAY.shp"
)
OUT_DIR = REPO / "outputs"
SUMMARY_TXT = OUT_DIR / "58a_residual_diagnostic.txt"
SUMMARY_JSON = OUT_DIR / "58a_residual_diagnostic.json"

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
        "C5": int((val > 8).sum()),
    }


def _try_operator(label, mesh, before_counts, op):
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


def _face_face_adjacency(elements: np.ndarray) -> dict[int, list[int]]:
    """Return ``adj[e] = [e', e'', ...]`` listing every element that
    shares an edge with element ``e``.
    """
    edge_to_elems: dict[tuple[int, int], list[int]] = defaultdict(list)
    for ei, tri in enumerate(elements):
        for el in range(3):
            a = int(tri[el])
            b = int(tri[(el + 1) % 3])
            edge_to_elems[(min(a, b), max(a, b))].append(ei)
    adj: dict[int, list[int]] = defaultdict(list)
    for elems in edge_to_elems.values():
        if len(elems) == 2:
            e1, e2 = elems
            adj[e1].append(e2)
            adj[e2].append(e1)
    return adj


def _cluster_fail_elements(
    fail_elem_set: set[int],
    adj: dict[int, list[int]],
) -> list[list[int]]:
    """Group fail elements into connected components via face-face
    adjacency restricted to fail elements only.
    """
    seen: set[int] = set()
    clusters: list[list[int]] = []
    for start in fail_elem_set:
        if start in seen:
            continue
        comp: list[int] = []
        q = deque([start])
        seen.add(start)
        while q:
            e = q.popleft()
            comp.append(e)
            for nb in adj.get(e, []):
                if nb in fail_elem_set and nb not in seen:
                    seen.add(nb)
                    q.append(nb)
        clusters.append(comp)
    clusters.sort(key=len, reverse=True)
    return clusters


def main() -> int:
    if not INPUT.exists():
        raise SystemExit(f"input missing: {INPUT}")

    mesh = read_fort14(INPUT)
    print(
        f"[58a] input: NP={mesh.n_nodes:,} NE={mesh.n_elements:,}",
        flush=True,
    )

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
        f"[58a] global: C1={int(c1_mask.sum())} "
        f"C2={int(c2_mask.sum())} "
        f"C4={int(c4_edge_mask.sum())} "
        f"C5={int((val > 8).sum())}",
        flush=True,
    )

    edge_uses = _edge_use_counts(mesh.elements)
    boundary_edge_keys = {k for k, v in edge_uses.items() if len(v) == 1}
    boundary_node_set: set[int] = set()
    for (a, b) in boundary_edge_keys:
        boundary_node_set.add(a)
        boundary_node_set.add(b)
    _bp, _bn, edge_to_segment = _boundary_topology(mesh)

    projector = build_coastline_projector(
        [COASTLINE],
        max_snap_distance_m=500.0,
        mean_latitude_deg=float(mesh.nodes[:, 1].mean()),
    )

    global_before = _global_counts(mesh)

    # ====================================================================
    # 1. Per-element diagnostic on the 2 remaining C1 fails
    # ====================================================================
    c1_eids = sorted(int(e) for e in np.where(c1_mask)[0])
    print(
        f"[58a] inspecting {len(c1_eids)} C1 fail elements", flush=True,
    )
    c1_records = []
    for eid in c1_eids:
        verts = [int(v) for v in mesh.elements[eid]]
        coords = [
            [float(mesh.nodes[v, 0]), float(mesh.nodes[v, 1])]
            for v in verts
        ]
        edge_records = []
        for el in range(3):
            a = verts[el]
            b = verts[(el + 1) % 3]
            k = (min(a, b), max(a, b))
            edge_records.append({
                "local": el,
                "uv": [a, b],
                "boundary": k in boundary_edge_keys,
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
        best = None
        for trial in rec["operator_trials"]:
            if trial["result"] != "applied":
                continue
            d = trial["delta"]
            if d["C5"] > 0:
                continue
            score = d["C1"] + d["C2"] + d["C4"]
            if best is None or score < best[0]:
                best = (score, trial)
        rec["best_op"] = best[1] if best else None
        rec["best_delta_combined_c1_c2_c4"] = best[0] if best else None
        c1_records.append(rec)

    # ====================================================================
    # 2. C4 cluster analysis
    # ====================================================================
    print(
        f"[58a] analyzing {int(c4_edge_mask.sum())} C4 fail edges",
        flush=True,
    )
    c4_fail_edges = edge_uv[c4_edge_mask]
    c4_fail_pairs = elem_pair[c4_edge_mask]
    c4_fail_ratios = ac[c4_edge_mask]

    # Build full face-face adjacency
    adj = _face_face_adjacency(mesh.elements)

    # Set of fail elements: any element on either side of a C4 fail edge
    c4_fail_elem_set: set[int] = set()
    for pair in c4_fail_pairs:
        c4_fail_elem_set.add(int(pair[0]))
        c4_fail_elem_set.add(int(pair[1]))

    clusters = _cluster_fail_elements(c4_fail_elem_set, adj)
    cluster_sizes = [len(c) for c in clusters]

    # Boundary-touching characterization per cluster
    cluster_records = []
    for ci, comp in enumerate(clusters[:20]):  # top 20 clusters
        comp_verts = set()
        for e in comp:
            for v in mesh.elements[e]:
                comp_verts.add(int(v))
        boundary_verts_in_comp = comp_verts & boundary_node_set
        # Find boundary edges within this cluster
        bnd_edges_in_comp = 0
        for e in comp:
            tri = mesh.elements[e]
            for el in range(3):
                a = int(tri[el])
                b = int(tri[(el + 1) % 3])
                k = (min(a, b), max(a, b))
                if k in boundary_edge_keys:
                    bnd_edges_in_comp += 1
        # C4 fail edges incident on this cluster
        c4_in_comp = 0
        worst_ratio = 0.0
        for ei, pair in enumerate(c4_fail_pairs):
            if int(pair[0]) in set(comp) or int(pair[1]) in set(comp):
                c4_in_comp += 1
                worst_ratio = max(worst_ratio, float(c4_fail_ratios[ei]))
        comp_xy = np.array([
            [float(mesh.nodes[mesh.elements[e][0], 0]),
             float(mesh.nodes[mesh.elements[e][0], 1])]
            for e in comp
        ])
        cluster_records.append({
            "cluster_id": ci,
            "n_elements": len(comp),
            "n_vertices": len(comp_verts),
            "n_boundary_vertices": len(boundary_verts_in_comp),
            "n_boundary_edges": bnd_edges_in_comp,
            "n_c4_fail_edges": c4_in_comp,
            "worst_area_ratio": worst_ratio,
            "centroid_lonlat": [
                float(comp_xy[:, 0].mean()),
                float(comp_xy[:, 1].mean()),
            ],
            "all_touch_boundary": len(boundary_verts_in_comp) == len(comp_verts),
            "any_touch_boundary": len(boundary_verts_in_comp) > 0,
        })

    # Distribution stats
    from collections import Counter
    size_hist = Counter(cluster_sizes)
    n_clusters = len(clusters)
    n_boundary_touching_clusters = sum(
        1 for cr in cluster_records if cr["any_touch_boundary"]
    )
    largest_cluster = max(cluster_sizes) if cluster_sizes else 0
    pct_in_large_clusters = (
        sum(s for s in cluster_sizes if s >= 4) / max(sum(cluster_sizes), 1)
    )

    # ====================================================================
    # 3. Write outputs
    # ====================================================================
    payload = {
        "input": str(INPUT.resolve()),
        "global_counts": global_before,
        "c1_per_fail": c1_records,
        "c4_n_fail_edges": int(c4_edge_mask.sum()),
        "c4_n_fail_elements": len(c4_fail_elem_set),
        "c4_n_clusters": n_clusters,
        "c4_cluster_size_histogram": dict(sorted(size_hist.items())),
        "c4_largest_cluster": largest_cluster,
        "c4_n_boundary_touching_clusters_top20": n_boundary_touching_clusters,
        "c4_top20_clusters": cluster_records,
        "c4_pct_elements_in_clusters_ge4": pct_in_large_clusters,
    }
    SUMMARY_JSON.write_text(
        json.dumps(payload, indent=2, default=str), encoding="utf-8",
    )

    # Text report
    lines = [
        "PoC #58a — residual diagnostic on PoC #57 Stage 1 output",
        f"input: {INPUT.name}",
        f"global counts: {global_before}",
        "",
        "Section 1 — C1 fail per-element details",
        "-" * 50,
    ]
    for rec in c1_records:
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

    lines.extend([
        "",
        "Section 2 — C4 fail cluster analysis",
        "-" * 50,
        f"  C4 fail edges       : {int(c4_edge_mask.sum())}",
        f"  C4 fail elements    : {len(c4_fail_elem_set)} "
        f"(elements on either side of any C4 fail edge)",
        f"  connected components: {n_clusters}",
        f"  largest cluster     : {largest_cluster} elements",
        f"  cluster size hist   : "
        f"{dict(sorted(size_hist.items()))}",
        f"  % elements in clusters of size >= 4: "
        f"{pct_in_large_clusters * 100:.1f}%",
        "",
        "  Top-20 clusters (by element count):",
        f"    {'cid':>3} | {'n_el':>4} | {'n_v':>4} | {'n_bv':>4} | "
        f"{'n_be':>4} | {'n_c4':>4} | {'worst_r':>7} | "
        f"{'all_bnd':>7} | centroid_lonlat",
        "    " + "-" * 90,
    ])
    for cr in cluster_records:
        lines.append(
            f"    {cr['cluster_id']:>3} | "
            f"{cr['n_elements']:>4} | "
            f"{cr['n_vertices']:>4} | "
            f"{cr['n_boundary_vertices']:>4} | "
            f"{cr['n_boundary_edges']:>4} | "
            f"{cr['n_c4_fail_edges']:>4} | "
            f"{cr['worst_area_ratio']:>7.3f} | "
            f"{str(cr['all_touch_boundary']):>7} | "
            f"({cr['centroid_lonlat'][0]:.4f}, "
            f"{cr['centroid_lonlat'][1]:.4f})"
        )

    lines.extend([
        "",
        "Notes:",
        "  n_el   = elements in cluster",
        "  n_v    = unique vertices in cluster",
        "  n_bv   = vertices on boundary",
        "  n_be   = boundary edges within cluster",
        "  n_c4   = C4 fail edges incident on cluster",
        "  worst_r = largest area_change ratio in cluster",
        "  all_bnd = every vertex is on boundary?",
    ])

    SUMMARY_TXT.write_text("\n".join(lines) + "\n", encoding="utf-8")
    print()
    print("\n".join(lines))
    print()
    print(f"wrote {SUMMARY_TXT}")
    print(f"wrote {SUMMARY_JSON}")
    return 0


if __name__ == "__main__":
    raise SystemExit(main())
