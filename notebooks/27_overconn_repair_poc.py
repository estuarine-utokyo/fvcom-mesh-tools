"""PoC #27: greedy edge-swap repair of over-connected nodes.

Targets the last detector still flagged after `fmesh-mesh-clean` Phase
A+B+C: nodes with valence > MAX_NBR_ELEM (FVCOM's per-node element
neighbour cap). The PoC explores whether a Lawson-style edge flip,
scored by reduction of the per-edge "valence excess" (sum over the
four affected nodes of `max(0, valence - max_nbr_elem)`), can drive
the valence distribution below the threshold without destroying
triangle quality.

Per interior edge `(i, j)` shared by triangles T1=(i,j,k) and
T2=(i,j,m), a flip turns the diagonal from `i-j` into `k-m`:

    valence(i): -1   valence(j): -1
    valence(k): +1   valence(m): +1

So a flip strictly reduces `valence(i) + valence(j)` and increases
`valence(k) + valence(m)`. We accept the flip if it does not invert a
triangle, does not introduce a new over-connection, optionally keeps
the worst min-angle of the two triangles above a floor, and reduces
the total excess. Non-conflicting flips (no two share a triangle in a
pass) are applied greedily by descending score; multiple passes run
until convergence or `max_iters`.

Cases:

    (a) synthetic 12-wedge fan            — interior node valence 12,
        target 8. Pure topology test, no boundary noise.
    (b) outputs/19_tokyo_bay_oceanmesh_cleaned.14 — 3 over-connected
        nodes, max valence 9. Easy real-world case.
    (c) outputs/16_tokyo_bay_with_rivers.14       — 440 over-connected
        nodes, max valence 26. The hard OCSMesh+gmsh case from PoC #25.

Outputs:
    outputs/27_overconn_repair_summary.txt
"""

from __future__ import annotations

from dataclasses import replace
from pathlib import Path

import numpy as np

from fvcom_mesh_tools.algorithms import alpha_quality, min_interior_angle
from fvcom_mesh_tools.algorithms.edge_swap import (
    _interior_edge_pairs,
    _min_angle,
    _signed_area,
)
from fvcom_mesh_tools.io import Fort14Mesh, read_fort14

REPO_ROOT = Path(__file__).resolve().parent.parent
OUT_DIR = REPO_ROOT / "outputs"
SUMMARY_TXT = OUT_DIR / "27_overconn_repair_summary.txt"

CASE_PATHS: list[tuple[str, Path | None]] = [
    ("synthetic_fan_v12", None),
    (
        "19_tokyo_bay_cleaned",
        OUT_DIR / "19_tokyo_bay_oceanmesh_cleaned.14",
    ),
    (
        "16_tokyo_bay_with_rivers",
        OUT_DIR / "16_tokyo_bay_with_rivers.14",
    ),
]

MAX_NBR_ELEM = 8
MAX_ITERS = 100
MIN_ANGLE_FLOOR_DEG = 20.0


# ---------------------------------------------------------------------------
# Algorithm
# ---------------------------------------------------------------------------


def _node_valence(elements: np.ndarray, n_nodes: int) -> np.ndarray:
    counts = np.zeros(n_nodes, dtype=np.int64)
    np.add.at(counts, elements.ravel(), 1)
    return counts


def _excess(val: np.ndarray, cap: int) -> np.ndarray:
    return np.maximum(0, val - cap)


def swap_edges_for_valence(
    mesh: Fort14Mesh,
    *,
    max_nbr_elem: int = MAX_NBR_ELEM,
    max_iters: int = MAX_ITERS,
    min_angle_floor_deg: float = MIN_ANGLE_FLOOR_DEG,
) -> tuple[Fort14Mesh, dict]:
    """Greedy Lawson-style flips that reduce node valence excess.

    A flip is accepted when:

        * the diagonal is interior (boundary edges are excluded from
          the candidate list automatically by ``_interior_edge_pairs``);
        * neither resulting triangle is inverted (signed area > 0);
        * worst-of-pair min-angle does not drop below
          ``min(min_angle_floor_deg, current_worst)`` — avoids creating
          triangles worse than the FVCOM 20° threshold unless the input
          itself was below it;
        * total per-edge excess strictly decreases (excess defined as
          ``sum(max(0, valence - cap))`` over the four nodes the flip
          touches).

    Non-conflicting flips (no shared triangle within one pass) are
    applied greedily by descending score; the loop stops when no
    over-connected node remains, or no candidate improves the score.
    """
    elements = mesh.elements.copy()
    nodes = mesh.nodes
    n_nodes = mesh.n_nodes
    history: list[int] = []
    initial_excess = int(_excess(_node_valence(elements, n_nodes), max_nbr_elem).sum())

    for _ in range(max_iters):
        val = _node_valence(elements, n_nodes)
        excess = _excess(val, max_nbr_elem)
        if excess.sum() == 0:
            break

        edges, tri_pairs = _interior_edge_pairs(elements)
        if edges.size == 0:
            break

        i, j = edges[:, 0], edges[:, 1]
        T1 = elements[tri_pairs[:, 0]]
        T2 = elements[tri_pairs[:, 1]]
        k = T1.sum(axis=1) - i - j
        m = T2.sum(axis=1) - i - j

        cand1 = np.stack([i, m, k], axis=1)
        cand2 = np.stack([j, k, m], axis=1)
        for cand in (cand1, cand2):
            sa = _signed_area(nodes, cand)
            flip = sa < 0
            if flip.any():
                cand[flip] = cand[flip][:, [0, 2, 1]]
        sa1 = _signed_area(nodes, cand1)
        sa2 = _signed_area(nodes, cand2)
        no_invert = (sa1 > 0) & (sa2 > 0)

        v_i, v_j, v_k, v_m = val[i], val[j], val[k], val[m]
        cap = max_nbr_elem
        excess_before = (
            np.maximum(0, v_i - cap)
            + np.maximum(0, v_j - cap)
            + np.maximum(0, v_k - cap)
            + np.maximum(0, v_m - cap)
        )
        excess_after = (
            np.maximum(0, v_i - 1 - cap)
            + np.maximum(0, v_j - 1 - cap)
            + np.maximum(0, v_k + 1 - cap)
            + np.maximum(0, v_m + 1 - cap)
        )
        score = excess_before - excess_after

        q_before = np.minimum(_min_angle(nodes, T1), _min_angle(nodes, T2))
        q_after = np.minimum(_min_angle(nodes, cand1), _min_angle(nodes, cand2))
        quality_ok = q_after >= np.minimum(q_before, min_angle_floor_deg) - 1e-9

        improves = no_invert & quality_ok & (score > 0)
        if not improves.any():
            break

        order = np.argsort(-score, kind="stable")
        used = np.zeros(elements.shape[0], dtype=bool)
        applied = 0
        for k_idx in order:
            if not improves[k_idx]:
                break
            t1 = int(tri_pairs[k_idx, 0])
            t2 = int(tri_pairs[k_idx, 1])
            if used[t1] or used[t2]:
                continue
            elements[t1] = cand1[k_idx]
            elements[t2] = cand2[k_idx]
            used[t1] = True
            used[t2] = True
            applied += 1
        history.append(applied)
        if applied == 0:
            break

    out_mesh = replace(mesh, elements=elements)
    final_val = _node_valence(elements, n_nodes)
    info = {
        "max_nbr_elem": int(max_nbr_elem),
        "swaps_per_iter": history,
        "total_swaps": int(sum(history)),
        "iterations_run": len(history),
        "initial_total_excess": int(initial_excess),
        "final_total_excess": int(_excess(final_val, max_nbr_elem).sum()),
        "max_valence_before": int(_node_valence(mesh.elements, n_nodes).max())
        if n_nodes else 0,
        "max_valence_after": int(final_val.max()) if n_nodes else 0,
        "n_overconn_before": int((_node_valence(mesh.elements, n_nodes) > max_nbr_elem).sum()),
        "n_overconn_after": int((final_val > max_nbr_elem).sum()),
    }
    return out_mesh, info


# ---------------------------------------------------------------------------
# Test cases
# ---------------------------------------------------------------------------


def synthetic_fan(n_wedges: int) -> Fort14Mesh:
    """Regular fan of ``n_wedges`` triangles around a central node.

    Centre node valence = n_wedges. Outer ring is land boundary.
    """
    centre = np.array([[0.0, 0.0]])
    angles = np.linspace(0, 2 * np.pi, n_wedges + 1)[:-1]
    rim = np.column_stack([np.cos(angles), np.sin(angles)])
    nodes = np.vstack([centre, rim])
    elements = np.array(
        [[0, 1 + i, 1 + (i + 1) % n_wedges] for i in range(n_wedges)],
        dtype=np.int64,
    )
    rim_chain = np.concatenate([np.arange(1, n_wedges + 1), [1]])
    return Fort14Mesh(
        title="fan",
        nodes=nodes,
        depths=np.zeros(len(nodes)),
        elements=elements,
        open_boundaries=[],
        land_boundaries=[(0, rim_chain.astype(np.int64))],
    )


# ---------------------------------------------------------------------------
# Driver
# ---------------------------------------------------------------------------


def _summary_lines(name: str, before: Fort14Mesh, after: Fort14Mesh, info: dict) -> list[str]:
    val_b = _node_valence(before.elements, before.n_nodes)
    val_a = _node_valence(after.elements, after.n_nodes)
    q_b = alpha_quality(before)
    q_a = alpha_quality(after)
    a_b = min_interior_angle(before)
    a_a = min_interior_angle(after)
    return [
        f"=== {name} ===",
        f"NP={before.n_nodes:,}  NE={before.n_elements:,}",
        f"max valence : {int(val_b.max()):2d} -> {int(val_a.max()):2d}",
        f"n over_connected (>8) : {int((val_b > 8).sum()):,} -> {int((val_a > 8).sum()):,}",
        f"total excess: {info['initial_total_excess']} -> {info['final_total_excess']}",
        f"swaps per iter: {info['swaps_per_iter'][:20]}"
        + (" ..." if len(info["swaps_per_iter"]) > 20 else ""),
        f"total swaps : {info['total_swaps']:,}  "
        f"(iterations {info['iterations_run']})",
        f"alpha mean  : {q_b.mean():.4f} -> {q_a.mean():.4f}",
        f"alpha p10   : {np.percentile(q_b, 10):.4f} -> "
        f"{np.percentile(q_a, 10):.4f}",
        f"min-angle p50: {np.percentile(a_b, 50):.2f} -> "
        f"{np.percentile(a_a, 50):.2f}  deg",
        f"frac<20 deg : {(a_b < 20).mean() * 100:.2f} -> "
        f"{(a_a < 20).mean() * 100:.2f} %",
    ]


FLOOR_VARIANTS: list[tuple[str, float]] = [
    ("conservative_floor20", 20.0),
    ("relaxed_floor10",      10.0),
    ("aggressive_floor0",     0.0),
]


def main() -> None:
    OUT_DIR.mkdir(parents=True, exist_ok=True)
    sections: list[str] = []
    for name, path in CASE_PATHS:
        if name == "synthetic_fan_v12":
            mesh = synthetic_fan(12)
        else:
            if path is None or not path.exists():
                sections.append(f"=== {name} ===\nSKIP (missing {path})")
                continue
            mesh = read_fort14(path)
        print(
            f"[27] === case={name}  "
            f"NP={mesh.n_nodes:,}  NE={mesh.n_elements:,} ==="
        )
        for variant, floor in FLOOR_VARIANTS:
            repaired, info = swap_edges_for_valence(
                mesh,
                max_nbr_elem=MAX_NBR_ELEM,
                max_iters=MAX_ITERS,
                min_angle_floor_deg=floor,
            )
            heading = f"{name} / {variant} (min_angle_floor={floor:.0f}°)"
            section = "\n".join(_summary_lines(heading, mesh, repaired, info))
            print(section)
            print()
            sections.append(section)

    summary = "\n\n".join(sections) + "\n"
    SUMMARY_TXT.write_text(summary, encoding="utf-8")
    print(f"[27] wrote {SUMMARY_TXT}")


if __name__ == "__main__":
    main()
