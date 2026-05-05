"""PoC #3: does MeshKernel orthogonalization improve FVCOM open-boundary perpendicularity?

Runs MeshKernel's mesh2d_compute_orthogonalization (defaults) on
tb_futtsu20220311.14, retrieves the modified geometry, and compares the
two metrics from PoC #2 between before and after.

Outputs:
    outputs/03_orthogonalize_summary.txt
    outputs/03_perp_before_after_hist.png
    outputs/03_open_bdy_xy.png
    outputs/03_tb_futtsu20220311_orthogonalized.14
"""

from __future__ import annotations

import time
from pathlib import Path

import matplotlib

matplotlib.use("Agg")

import matplotlib.pyplot as plt  # noqa: E402
import numpy as np  # noqa: E402
from meshkernel import (  # noqa: E402
    GeometryList,
    Mesh2d,
    MeshKernel,
    OrthogonalizationParameters,
    ProjectionType,
    ProjectToLandBoundaryOption,
)

from fvcom_mesh_tools.io import Fort14Mesh, read_fort14, write_fort14  # noqa: E402

REPO_ROOT = Path(__file__).resolve().parent.parent
MESH_PATH = REPO_ROOT / "data" / "mesh" / "reference" / "tokyo_bay" / "tb_futtsu20220311.14"
OUT_DIR = REPO_ROOT / "outputs"
SUMMARY_TXT = OUT_DIR / "03_orthogonalize_summary.txt"
PERP_HIST_BA = OUT_DIR / "03_perp_before_after_hist.png"
OPEN_BDY_XY = OUT_DIR / "03_open_bdy_xy.png"
ORTHO_OUT_F14 = OUT_DIR / "03_tb_futtsu20220311_orthogonalized.14"

MK_MISSING = -999.0


def unique_edges(elements: np.ndarray) -> np.ndarray:
    e = np.vstack([
        elements[:, [0, 1]], elements[:, [1, 2]], elements[:, [2, 0]],
    ])
    e.sort(axis=1)
    return np.unique(e, axis=0)


def to_mesh2d(mesh: Fort14Mesh) -> Mesh2d:
    edges = unique_edges(mesh.elements)
    return Mesh2d(
        node_x=np.ascontiguousarray(mesh.nodes[:, 0], dtype=np.float64),
        node_y=np.ascontiguousarray(mesh.nodes[:, 1], dtype=np.float64),
        edge_nodes=edges.flatten().astype(np.int32),
        face_nodes=mesh.elements.flatten().astype(np.int32),
        nodes_per_face=np.full(mesh.n_elements, 3, dtype=np.int32),
    )


def mk_orthogonality(mesh: Fort14Mesh) -> tuple[np.ndarray, np.ndarray]:
    mk = MeshKernel(projection=ProjectionType.SPHERICAL)
    mk.mesh2d_set(to_mesh2d(mesh))
    geom = mk.mesh2d_get_orthogonality()
    values = np.asarray(geom.values, dtype=np.float64)
    return values, values != MK_MISSING


def open_bdy_perp(mesh: Fort14Mesh) -> np.ndarray:
    if not mesh.open_boundaries:
        return np.array([])
    bdy = np.asarray(mesh.open_boundaries[0], dtype=np.int64)
    bdy_xy = mesh.nodes[bdy]

    tangents = np.empty_like(bdy_xy)
    tangents[1:-1] = bdy_xy[2:] - bdy_xy[:-2]
    tangents[0] = bdy_xy[1] - bdy_xy[0]
    tangents[-1] = bdy_xy[-1] - bdy_xy[-2]
    norms = np.linalg.norm(tangents, axis=1, keepdims=True)
    tangents = tangents / np.where(norms == 0, 1.0, norms)

    inv_map = np.full(mesh.n_nodes, -1, dtype=np.int64)
    inv_map[bdy] = np.arange(len(bdy))

    edges = unique_edges(mesh.elements)
    a_in = inv_map[edges[:, 0]] >= 0
    b_in = inv_map[edges[:, 1]] >= 0
    incident = a_in ^ b_in
    inc = edges[incident]
    inc_a_in = a_in[incident]
    bdy_node = np.where(inc_a_in, inc[:, 0], inc[:, 1])
    int_node = np.where(inc_a_in, inc[:, 1], inc[:, 0])

    edge_vec = mesh.nodes[int_node] - mesh.nodes[bdy_node]
    edge_norms = np.linalg.norm(edge_vec, axis=1, keepdims=True)
    edge_vec = edge_vec / np.where(edge_norms == 0, 1.0, edge_norms)
    edge_tangent = tangents[inv_map[bdy_node]]
    cos_angle = np.clip((edge_vec * edge_tangent).sum(axis=1), -1.0, 1.0)
    angle = np.degrees(np.arccos(np.abs(cos_angle)))
    return 90.0 - angle


def orthogonalize(mesh: Fort14Mesh) -> tuple[Fort14Mesh, dict]:
    """Run MK orthogonalization with default parameters and return a Fort14Mesh
    with updated node coordinates (topology and boundaries are preserved).
    """
    mk = MeshKernel(projection=ProjectionType.SPHERICAL)
    mk.mesh2d_set(to_mesh2d(mesh))

    params = OrthogonalizationParameters()
    t0 = time.perf_counter()
    mk.mesh2d_compute_orthogonalization(
        ProjectToLandBoundaryOption.DO_NOT_PROJECT_TO_LANDBOUNDARY,
        params,
        GeometryList(),
    )
    elapsed = time.perf_counter() - t0

    new_m2d = mk.mesh2d_get()
    if len(new_m2d.node_x) != mesh.n_nodes:
        raise RuntimeError(
            f"node count changed during orthogonalization: "
            f"{mesh.n_nodes} -> {len(new_m2d.node_x)}"
        )

    new_nodes = np.column_stack([new_m2d.node_x, new_m2d.node_y]).astype(np.float64)
    after = Fort14Mesh(
        title=mesh.title + " (MK orthogonalized)",
        nodes=new_nodes,
        depths=mesh.depths.copy(),
        elements=mesh.elements.copy(),
        open_boundaries=[a.copy() for a in mesh.open_boundaries],
        land_boundaries=[(ib, a.copy()) for (ib, a) in mesh.land_boundaries],
    )
    info = {
        "elapsed_sec": elapsed,
        "params": vars(params),
    }
    return after, info


def _stats(name: str, x: np.ndarray) -> str:
    if x.size == 0:
        return f"{name}: empty"
    return (
        f"{name}: n={x.size:,}  min={np.nanmin(x):.4f}  "
        f"p50={np.nanpercentile(x, 50):.4f}  p95={np.nanpercentile(x, 95):.4f}  "
        f"max={np.nanmax(x):.4f}  mean={np.nanmean(x):.4f}"
    )


def plot_hist(perp_b: np.ndarray, perp_a: np.ndarray, png: Path) -> None:
    fig, ax = plt.subplots(figsize=(7, 4), dpi=120)
    bins = np.linspace(0, 90, 60)
    ax.hist(perp_b, bins=bins, alpha=0.6, label="before", color="tab:red")
    ax.hist(perp_a, bins=bins, alpha=0.6, label="after", color="tab:blue")
    ax.set_xlabel("perp deviation from 90 deg (deg)")
    ax.set_ylabel("count")
    ax.set_title("FVCOM open-boundary perpendicularity: before vs after MK orthogonalize")
    ax.legend()
    ax.grid(True, lw=0.3, color="0.85")
    fig.tight_layout()
    fig.savefig(png, bbox_inches="tight")
    plt.close(fig)


def plot_open_bdy_xy(before: Fort14Mesh, after: Fort14Mesh, png: Path) -> None:
    if not before.open_boundaries:
        return
    bdy_idx = before.open_boundaries[0]
    fig, ax = plt.subplots(figsize=(8, 7), dpi=120)
    ax.scatter(
        before.nodes[:, 0], before.nodes[:, 1],
        s=0.2, color="0.9", linewidths=0,
    )
    ax.plot(
        before.nodes[bdy_idx, 0], before.nodes[bdy_idx, 1], "-",
        color="tab:red", lw=1.2, label="before",
    )
    ax.plot(
        after.nodes[bdy_idx, 0], after.nodes[bdy_idx, 1], "-",
        color="tab:blue", lw=1.0, label="after",
    )
    ax.set_aspect("equal")
    ax.set_xlabel("lon (deg)")
    ax.set_ylabel("lat (deg)")
    ax.set_title("Open boundary: before vs after MK orthogonalization")
    ax.legend()
    ax.grid(True, lw=0.3, color="0.9")
    fig.tight_layout()
    fig.savefig(png, bbox_inches="tight")
    plt.close(fig)


def main() -> None:
    if not MESH_PATH.exists():
        raise SystemExit(f"mesh not found: {MESH_PATH}")
    OUT_DIR.mkdir(parents=True, exist_ok=True)

    print(f"[03] reading {MESH_PATH}")
    before = read_fort14(MESH_PATH)
    print(f"[03] NP={before.n_nodes:,}  NE={before.n_elements:,}")

    print("[03] measuring metrics on the original mesh ...")
    mk_b, mk_b_v = mk_orthogonality(before)
    perp_b = open_bdy_perp(before)

    print("[03] running MK mesh2d_compute_orthogonalization (defaults) ...")
    after, info = orthogonalize(before)
    print(f"[03] orthogonalization elapsed: {info['elapsed_sec']:.1f}s")

    print("[03] measuring metrics on the orthogonalized mesh ...")
    mk_a, mk_a_v = mk_orthogonality(after)
    perp_a = open_bdy_perp(after)

    disp = np.linalg.norm(after.nodes - before.nodes, axis=1)

    lines = [
        f"file:  {MESH_PATH}",
        f"NP={before.n_nodes:,}  NE={before.n_elements:,}",
        f"orthogonalization elapsed: {info['elapsed_sec']:.1f}s",
        f"params: {info['params']}",
        "",
        "[Per-node displacement (degrees on the lon/lat sphere)]",
        _stats("  ||after - before||", disp),
        "",
        "[MK orthogonality (valid edges only, 0 = orthogonal)]",
        _stats("  before", mk_b[mk_b_v]),
        _stats("  after ", mk_a[mk_a_v]),
        "",
        "[FVCOM open-boundary perpendicularity (deg from 90 deg)]",
        _stats("  before per-edge", perp_b),
        _stats("  after  per-edge", perp_a),
    ]
    summary = "\n".join(lines)
    print(summary)
    SUMMARY_TXT.write_text(summary + "\n", encoding="utf-8")
    print(f"[03] wrote {SUMMARY_TXT}")

    write_fort14(after, ORTHO_OUT_F14)
    print(f"[03] wrote {ORTHO_OUT_F14}")

    plot_hist(perp_b, perp_a, PERP_HIST_BA)
    print(f"[03] wrote {PERP_HIST_BA}")

    plot_open_bdy_xy(before, after, OPEN_BDY_XY)
    print(f"[03] wrote {OPEN_BDY_XY}")


if __name__ == "__main__":
    main()
