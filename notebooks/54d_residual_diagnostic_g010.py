"""PoC #54d: diagnose residual FVCOM violations on PoC #54c output.

PoC #54c (g=0.10 + thin_chain=none + Phase H A+B+E) left:
    C1 (min_ang<30°)        : 8 fails
    C2 (max_ang>130°)       : 0 ✓
    C4 (area_change>0.5)    : 69 fails
    C5 (valence>8)          : 0 ✓
    n_flipped               : 0 ✓

The PoC #52b baseline (g=0.15, thin_chain=widen) was 90 violations
(9 C1 / 1 C2 / 80 C4 / 0 C5). The g=0.10 + thin_chain=none switch
removed 13 violations net AND drove C2 to zero. To decide whether
the next attack should be coastline-aware C1 operators or even
tighter grading (g=0.05) or a different mesh-size strategy
altogether, we need to know the *structure* of the new residual:

  - C1 fails: same boundary-pinned pattern as PoC #53, or new
    interior clusters from the tighter mesh?
  - C4 fails: still barely-failing (ac in 0.50-0.55) and
    boundary-adjacent, or have they shifted toward extreme ratios
    in interior?

This is read-only and runs on the login node in under 30 s.

Outputs:
    outputs/54d_diagnostic_summary.txt
    outputs/54d_diagnostic_summary.json
"""
from __future__ import annotations

import json
from collections import Counter
from pathlib import Path

import numpy as np

from fvcom_mesh_tools.algorithms.quality import min_interior_angle
from fvcom_mesh_tools.diagnostics import (
    face_face_adjacency,
    node_valence,
)
from fvcom_mesh_tools.io import read_fort14
from fvcom_mesh_tools.mesh_clean_phase_h import _per_edge_area_change

REPO = Path(__file__).resolve().parent.parent
INPUT = REPO / "outputs" / "54c_phase_h_optimized.14"
OUT_DIR = REPO / "outputs"
SUMMARY_TXT = OUT_DIR / "54d_diagnostic_summary.txt"
SUMMARY_JSON = OUT_DIR / "54d_diagnostic_summary.json"

MIN_ANGLE_TARGET = 30.0
MAX_ANGLE_TARGET = 130.0
AREA_RATIO_TARGET = 0.5


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


def _boundary_node_set(mesh) -> set[int]:
    s: set[int] = set()
    for seg in mesh.open_boundaries:
        s.update(int(x) for x in np.asarray(seg))
    for _ib, seg in mesh.land_boundaries:
        s.update(int(x) for x in np.asarray(seg))
    return s


def _elem_centroid(nodes, elements, eid):
    verts = elements[eid]
    return nodes[verts].mean(axis=0)


def _elem_edge_lengths(nodes, elements, eid):
    verts = elements[eid]
    p = nodes[verts]
    e0 = np.linalg.norm(p[1] - p[0])
    e1 = np.linalg.norm(p[2] - p[1])
    e2 = np.linalg.norm(p[0] - p[2])
    return np.array([e0, e1, e2])


def _cluster_face_adj(elements, fail_mask):
    """Connected components of fail elements under face-face adjacency."""
    from scipy.sparse.csgraph import (  # noqa: PLC0415
        connected_components,
    )
    fail_eids = np.where(fail_mask)[0]
    if fail_eids.size == 0:
        return []
    adj = face_face_adjacency(elements)
    fail_sub = adj[fail_mask][:, fail_mask]
    _n, labels = connected_components(
        fail_sub, directed=False, return_labels=True,
    )
    clusters = []
    for cid in range(int(labels.max()) + 1 if labels.size else 0):
        members = fail_eids[labels == cid]
        clusters.append(members.astype(np.int64))
    clusters.sort(key=len, reverse=True)
    return clusters


def _cluster_summary(nodes, elements, clusters, bnd_nodes):
    out = []
    for i, c in enumerate(clusters):
        bnd_count = 0
        for eid in c:
            verts = elements[int(eid)]
            if any(int(v) in bnd_nodes for v in verts):
                bnd_count += 1
        centroids = np.array(
            [_elem_centroid(nodes, elements, int(e)) for e in c],
        )
        edge_lens = np.concatenate(
            [_elem_edge_lengths(nodes, elements, int(e)) for e in c],
        )
        out.append({
            "cluster_index": i,
            "size": int(c.size),
            "boundary_touching_elems": bnd_count,
            "all_boundary": bool(bnd_count == int(c.size)),
            "lon_range": [
                float(centroids[:, 0].min()),
                float(centroids[:, 0].max()),
            ],
            "lat_range": [
                float(centroids[:, 1].min()),
                float(centroids[:, 1].max()),
            ],
            "center_lonlat": [
                float(centroids[:, 0].mean()),
                float(centroids[:, 1].mean()),
            ],
            "edge_len_median_deg": float(np.median(edge_lens)),
            "elements": [int(e) for e in c],
        })
    return out


def main() -> int:
    if not INPUT.exists():
        raise SystemExit(f"input missing: {INPUT}")
    mesh = read_fort14(INPUT)
    nodes = mesh.nodes
    elements = mesh.elements
    NE = elements.shape[0]
    print(
        f"[54d] input: NP={mesh.n_nodes:,}  NE={NE:,}",
        flush=True,
    )

    # Compute per-element angle stats
    m = min_interior_angle(mesh)
    M = _max_interior_angle(mesh)
    c1_mask = m < MIN_ANGLE_TARGET
    c2_mask = M > MAX_ANGLE_TARGET
    n_c1 = int(c1_mask.sum())
    n_c2 = int(c2_mask.sum())

    # C4 per internal edge
    edge_uv, elem_pair, ac = _per_edge_area_change(nodes, elements)
    c4_edge_mask = ac > AREA_RATIO_TARGET
    n_c4 = int(c4_edge_mask.sum())
    c4_edges = edge_uv[c4_edge_mask]
    c4_pairs = elem_pair[c4_edge_mask]
    c4_acs = ac[c4_edge_mask]

    # C5
    val = node_valence(elements, mesh.n_nodes)
    n_c5 = int((val > 8).sum())
    max_val = int(val.max())

    bnd_nodes = _boundary_node_set(mesh)

    # Cluster C1 (face-face adjacency over c1 fail elements)
    c1_clusters = _cluster_face_adj(elements, c1_mask)
    c1_clu = _cluster_summary(nodes, elements, c1_clusters, bnd_nodes)

    # Cluster C2
    c2_clusters = _cluster_face_adj(elements, c2_mask)
    c2_clu = _cluster_summary(nodes, elements, c2_clusters, bnd_nodes)

    # Cluster C4: elements that touch any C4 fail edge
    c4_touching = np.unique(c4_pairs.ravel())
    c4_elem_mask = np.zeros(NE, dtype=bool)
    c4_elem_mask[c4_touching] = True
    c4_clusters = _cluster_face_adj(elements, c4_elem_mask)
    c4_clu = _cluster_summary(nodes, elements, c4_clusters, bnd_nodes)

    # C4 edge-level details
    c4_edge_records = []
    for i in range(c4_edges.shape[0]):
        u, v = int(c4_edges[i, 0]), int(c4_edges[i, 1])
        e_i, e_j = int(c4_pairs[i, 0]), int(c4_pairs[i, 1])
        # Larger / smaller triangle area
        from fvcom_mesh_tools.mesh_clean_phase_h import (  # noqa: PLC0415
            _signed_areas,
        )
        areas = np.abs(_signed_areas(nodes, elements[[e_i, e_j]]))
        big, small = float(areas.max()), float(areas.min())
        ratio = big / max(small, 1e-30)
        on_bnd = (u in bnd_nodes) or (v in bnd_nodes)
        midpoint = 0.5 * (nodes[u] + nodes[v])
        c4_edge_records.append({
            "edge": [u, v],
            "elem_pair": [e_i, e_j],
            "area_change": float(c4_acs[i]),
            "area_ratio_LS": float(ratio),
            "endpoint_on_boundary": bool(on_bnd),
            "midpoint_lonlat": [float(midpoint[0]), float(midpoint[1])],
        })

    # Histograms
    c1_cluster_size_hist = Counter(c.size for c in c1_clusters)
    c2_cluster_size_hist = Counter(c.size for c in c2_clusters)
    c4_cluster_size_hist = Counter(c.size for c in c4_clusters)
    c4_ratio_buckets = Counter()
    for r in c4_edge_records:
        v = r["area_change"]
        if v < 0.55:
            c4_ratio_buckets["0.50-0.55"] += 1
        elif v < 0.60:
            c4_ratio_buckets["0.55-0.60"] += 1
        elif v < 0.70:
            c4_ratio_buckets["0.60-0.70"] += 1
        elif v < 0.80:
            c4_ratio_buckets["0.70-0.80"] += 1
        else:
            c4_ratio_buckets["0.80-1.00"] += 1

    # Boundary share
    n_c1_bnd_clusters = sum(
        1 for cs in c1_clu if cs["boundary_touching_elems"] > 0
    )
    n_c2_bnd_clusters = sum(
        1 for cs in c2_clu if cs["boundary_touching_elems"] > 0
    )
    n_c4_bnd_clusters = sum(
        1 for cs in c4_clu if cs["boundary_touching_elems"] > 0
    )
    n_c4_bnd_edges = sum(
        1 for r in c4_edge_records if r["endpoint_on_boundary"]
    )

    # Assemble output
    payload = {
        "input": str(INPUT.resolve()),
        "n_nodes": int(mesh.n_nodes),
        "n_elements": int(NE),
        "counts": {
            "C1_min_ang_lt_30": n_c1,
            "C2_max_ang_gt_130": n_c2,
            "C4_area_change_gt_0_5": n_c4,
            "C5_valence_gt_8": n_c5,
            "max_valence": max_val,
        },
        "C1": {
            "n_fails": n_c1,
            "n_clusters": len(c1_clu),
            "cluster_size_histogram": dict(c1_cluster_size_hist),
            "n_clusters_touching_boundary": n_c1_bnd_clusters,
            "clusters": c1_clu,
        },
        "C2": {
            "n_fails": n_c2,
            "n_clusters": len(c2_clu),
            "cluster_size_histogram": dict(c2_cluster_size_hist),
            "n_clusters_touching_boundary": n_c2_bnd_clusters,
            "clusters": c2_clu,
        },
        "C4": {
            "n_fails": n_c4,
            "n_clusters": len(c4_clu),
            "cluster_size_histogram": dict(c4_cluster_size_hist),
            "n_clusters_touching_boundary": n_c4_bnd_clusters,
            "n_edges_touching_boundary": n_c4_bnd_edges,
            "area_change_histogram": dict(c4_ratio_buckets),
            "clusters": c4_clu,
            "edges": c4_edge_records,
        },
    }
    SUMMARY_JSON.write_text(
        json.dumps(payload, indent=2, default=str), encoding="utf-8",
    )

    # ASCII summary
    lines = [
        "PoC #54d — residual FVCOM violation diagnostic (g=0.10 + Phase H)",
        f"input: {INPUT.name}",
        f"NP={mesh.n_nodes:,}  NE={NE:,}",
        "",
        "Residual counts:",
        f"  C1 (min_ang < 30°)     : {n_c1}",
        f"  C2 (max_ang > 130°)    : {n_c2}",
        f"  C4 (area_change > 0.5) : {n_c4}",
        f"  C5 (valence > 8)       : {n_c5}  (max_valence = {max_val})",
        "",
    ]

    def _emit_clusters(label, clusters, hist, bnd_share):
        lines.append(f"--- {label} ---")
        lines.append(
            f"  total clusters: {len(clusters)}  (boundary-touching: {bnd_share})"
        )
        lines.append("  cluster size histogram:")
        for size in sorted(hist, reverse=True):
            lines.append(f"    size {size}: {hist[size]} clusters")
        lines.append("  top clusters (largest 5):")
        for cs in clusters[:5]:
            lines.append(
                f"    cluster #{cs['cluster_index']:>3}  size={cs['size']:>3}  "
                f"bnd_elems={cs['boundary_touching_elems']:>3}  "
                f"center=({cs['center_lonlat'][0]:.4f}, "
                f"{cs['center_lonlat'][1]:.4f})  "
                f"h_med={cs['edge_len_median_deg']:.4e}°"
            )

    _emit_clusters("C1", c1_clu, c1_cluster_size_hist, n_c1_bnd_clusters)
    lines.append("")
    _emit_clusters("C2", c2_clu, c2_cluster_size_hist, n_c2_bnd_clusters)
    lines.append("")
    _emit_clusters(
        "C4 (triangles touching C4 fail edge)", c4_clu,
        c4_cluster_size_hist, n_c4_bnd_clusters,
    )
    lines.append("")
    lines.append("--- C4 edge details ---")
    lines.append("  area_change histogram:")
    for bucket in sorted(c4_ratio_buckets):
        lines.append(f"    {bucket}: {c4_ratio_buckets[bucket]}")
    lines.append(
        f"  edges touching boundary: {n_c4_bnd_edges} / {n_c4}"
    )
    # Top-10 worst C4 fails
    worst = sorted(
        c4_edge_records, key=lambda r: -r["area_change"],
    )[:10]
    lines.append("  worst 10 C4 fails:")
    for r in worst:
        u, v = r["edge"]
        lines.append(
            f"    edge ({u:>5}, {v:>5})  ac={r['area_change']:.3f}  "
            f"ratio L:S={r['area_ratio_LS']:.2f}  "
            f"bnd={'Y' if r['endpoint_on_boundary'] else 'N'}  "
            f"mid=({r['midpoint_lonlat'][0]:.4f}, "
            f"{r['midpoint_lonlat'][1]:.4f})"
        )

    SUMMARY_TXT.write_text("\n".join(lines) + "\n", encoding="utf-8")
    print()
    print("\n".join(lines))
    print()
    print(f"wrote {SUMMARY_TXT}")
    print(f"wrote {SUMMARY_JSON}")
    return 0


if __name__ == "__main__":
    raise SystemExit(main())
