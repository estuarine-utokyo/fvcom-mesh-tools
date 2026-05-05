"""PoC #6: parity comparison between OCSMesh-generated mesh and reference.

Compares the minimal OCSMesh output (PoC #5) against the legacy reference
mesh ``tb_futtsu20220311.14`` produced by OceanMesh2D + manual editing.
The aim is to quantify *what is missing* from the Python pipeline so we
can scope Phase 4 / next-iteration work, not to rank the meshes.

Per-mesh metrics:
    - NP, NE, bounding box (lon/lat)
    - Edge-length distribution (great-circle metres)
    - Triangle alpha-quality (1 = equilateral, 0 = degenerate)
    - Triangle minimum-interior-angle (deg)
    - Boundary segment counts (open / land) and total boundary nodes

Outputs:
    outputs/06_parity_summary.txt
    outputs/06_parity_meshes_side_by_side.png
    outputs/06_parity_edge_length_hist.png
    outputs/06_parity_quality_hist.png
"""

from __future__ import annotations

from pathlib import Path

import matplotlib

matplotlib.use("Agg")

import matplotlib.pyplot as plt  # noqa: E402
import numpy as np  # noqa: E402

from fvcom_mesh_tools.algorithms import unique_edges  # noqa: E402
from fvcom_mesh_tools.io import Fort14Mesh, read_fort14  # noqa: E402

REPO_ROOT = Path(__file__).resolve().parent.parent
REF_MESH = REPO_ROOT / "data" / "mesh" / "reference" / "tokyo_bay" / "tb_futtsu20220311.14"
OCS_MESH = REPO_ROOT / "outputs" / "05_tokyo_bay_minimal.14"

OUT_DIR = REPO_ROOT / "outputs"
SUMMARY_TXT = OUT_DIR / "06_parity_summary.txt"
SIDE_BY_SIDE_PNG = OUT_DIR / "06_parity_meshes_side_by_side.png"
EDGE_LEN_PNG = OUT_DIR / "06_parity_edge_length_hist.png"
QUALITY_PNG = OUT_DIR / "06_parity_quality_hist.png"

EARTH_R_M = 6_371_000.0


def haversine_m(p0: np.ndarray, p1: np.ndarray) -> np.ndarray:
    """Great-circle distance in metres between matched lon/lat point arrays."""
    lon0 = np.deg2rad(p0[:, 0])
    lat0 = np.deg2rad(p0[:, 1])
    lon1 = np.deg2rad(p1[:, 0])
    lat1 = np.deg2rad(p1[:, 1])
    dlat = lat1 - lat0
    dlon = lon1 - lon0
    a = np.sin(dlat / 2) ** 2 + np.cos(lat0) * np.cos(lat1) * np.sin(dlon / 2) ** 2
    return 2 * EARTH_R_M * np.arcsin(np.sqrt(a))


def edge_lengths_m(mesh: Fort14Mesh) -> np.ndarray:
    e = unique_edges(mesh.elements)
    return haversine_m(mesh.nodes[e[:, 0]], mesh.nodes[e[:, 1]])


def alpha_quality(mesh: Fort14Mesh) -> np.ndarray:
    """Per-triangle quality alpha = 4*sqrt(3)*A / (l1^2 + l2^2 + l3^2)."""
    p0 = mesh.nodes[mesh.elements[:, 0]]
    p1 = mesh.nodes[mesh.elements[:, 1]]
    p2 = mesh.nodes[mesh.elements[:, 2]]
    l01 = haversine_m(p0, p1)
    l12 = haversine_m(p1, p2)
    l20 = haversine_m(p2, p0)
    # Use the same haversine-derived lengths for area-via-Heron, so the
    # quality stays scale-consistent with the edge-length metric.
    s = 0.5 * (l01 + l12 + l20)
    inner = np.maximum(s * (s - l01) * (s - l12) * (s - l20), 0.0)
    area = np.sqrt(inner)
    denom = l01 ** 2 + l12 ** 2 + l20 ** 2
    return 4.0 * np.sqrt(3.0) * area / np.where(denom == 0, 1.0, denom)


def min_interior_angle_deg(mesh: Fort14Mesh) -> np.ndarray:
    p0 = mesh.nodes[mesh.elements[:, 0]]
    p1 = mesh.nodes[mesh.elements[:, 1]]
    p2 = mesh.nodes[mesh.elements[:, 2]]
    l01 = haversine_m(p0, p1)
    l12 = haversine_m(p1, p2)
    l20 = haversine_m(p2, p0)
    a, b, c = l12, l20, l01
    # Law of cosines per vertex.
    def _ang(opp, e1, e2):
        cos = (e1 ** 2 + e2 ** 2 - opp ** 2) / np.where(
            e1 * e2 == 0, 1.0, 2.0 * e1 * e2
        )
        return np.degrees(np.arccos(np.clip(cos, -1.0, 1.0)))
    A = _ang(a, b, c)
    B = _ang(b, c, a)
    C = _ang(c, a, b)
    return np.minimum(np.minimum(A, B), C)


def _stats(name: str, x: np.ndarray, fmt: str = ".4f") -> str:
    return (
        f"{name}: n={x.size:,}  "
        f"min={np.nanmin(x):{fmt}}  "
        f"p05={np.nanpercentile(x, 5):{fmt}}  "
        f"p50={np.nanpercentile(x, 50):{fmt}}  "
        f"p95={np.nanpercentile(x, 95):{fmt}}  "
        f"max={np.nanmax(x):{fmt}}  "
        f"mean={np.nanmean(x):{fmt}}"
    )


def plot_meshes_side_by_side(
    ref: Fort14Mesh, ocs: Fort14Mesh, png: Path
) -> None:
    fig, axes = plt.subplots(1, 2, figsize=(14, 7), dpi=120, sharex=True, sharey=True)
    for ax, mesh, title in (
        (axes[0], ref, f"Reference (NP={ref.n_nodes:,}, NE={ref.n_elements:,})"),
        (axes[1], ocs, f"OCSMesh minimal (NP={ocs.n_nodes:,}, NE={ocs.n_elements:,})"),
    ):
        ax.triplot(
            mesh.nodes[:, 0], mesh.nodes[:, 1], mesh.elements,
            color="0.4", lw=0.15,
        )
        ax.set_aspect("equal")
        ax.set_title(title)
        ax.set_xlabel("lon (deg)")
        ax.grid(True, lw=0.3, color="0.9")
    axes[0].set_ylabel("lat (deg)")
    fig.suptitle("Tokyo Bay: reference vs OCSMesh minimal pipeline", y=1.02)
    fig.tight_layout()
    fig.savefig(png, bbox_inches="tight")
    plt.close(fig)


def plot_edge_len_hist(ref_e: np.ndarray, ocs_e: np.ndarray, png: Path) -> None:
    fig, ax = plt.subplots(figsize=(8, 4.5), dpi=120)
    bins = np.logspace(
        np.log10(max(min(ref_e.min(), ocs_e.min()), 1.0)),
        np.log10(max(ref_e.max(), ocs_e.max()) * 1.05),
        60,
    )
    ax.hist(ref_e, bins=bins, alpha=0.6, label="reference", color="tab:gray")
    ax.hist(ocs_e, bins=bins, alpha=0.6, label="OCSMesh minimal", color="tab:blue")
    ax.set_xscale("log")
    ax.set_xlabel("edge length (m, great circle)")
    ax.set_ylabel("count")
    ax.set_title("Edge-length distribution")
    ax.legend()
    ax.grid(True, which="both", lw=0.3, color="0.9")
    fig.tight_layout()
    fig.savefig(png, bbox_inches="tight")
    plt.close(fig)


def plot_quality_hist(ref_q: np.ndarray, ocs_q: np.ndarray, png: Path) -> None:
    fig, ax = plt.subplots(figsize=(8, 4.5), dpi=120)
    bins = np.linspace(0.0, 1.0, 41)
    ax.hist(ref_q, bins=bins, alpha=0.6, label="reference", color="tab:gray")
    ax.hist(ocs_q, bins=bins, alpha=0.6, label="OCSMesh minimal", color="tab:blue")
    ax.set_xlabel("triangle alpha-quality (1 = equilateral)")
    ax.set_ylabel("count")
    ax.set_title("Triangle quality distribution")
    ax.legend()
    ax.grid(True, lw=0.3, color="0.9")
    fig.tight_layout()
    fig.savefig(png, bbox_inches="tight")
    plt.close(fig)


def _boundary_summary(mesh: Fort14Mesh) -> str:
    n_open_segs = len(mesh.open_boundaries)
    n_open_nodes = sum(len(b) for b in mesh.open_boundaries)
    n_land_segs = len(mesh.land_boundaries)
    n_land_nodes = sum(len(ids) for _, ids in mesh.land_boundaries)
    ibtypes = sorted({ib for ib, _ in mesh.land_boundaries})
    return (
        f"  open : {n_open_segs} segments, {n_open_nodes:,} nodes\n"
        f"  land : {n_land_segs} segments, {n_land_nodes:,} nodes  "
        f"(ibtypes={ibtypes})"
    )


def main() -> None:
    if not REF_MESH.exists():
        raise SystemExit(f"reference mesh not found: {REF_MESH}")
    if not OCS_MESH.exists():
        raise SystemExit(
            f"OCSMesh output not found: {OCS_MESH}\n"
            f"Run notebooks/05_ocsmesh_minimal.py first."
        )
    OUT_DIR.mkdir(parents=True, exist_ok=True)

    print(f"[06] reading {REF_MESH}")
    ref = read_fort14(REF_MESH)
    print(f"[06] reading {OCS_MESH}")
    ocs = read_fort14(OCS_MESH)

    ref_e = edge_lengths_m(ref)
    ocs_e = edge_lengths_m(ocs)
    ref_q = alpha_quality(ref)
    ocs_q = alpha_quality(ocs)
    ref_a = min_interior_angle_deg(ref)
    ocs_a = min_interior_angle_deg(ocs)

    lines = [
        f"reference : {REF_MESH}",
        f"OCSMesh   : {OCS_MESH}",
        "",
        "[Mesh size]",
        f"  reference : NP={ref.n_nodes:,}  NE={ref.n_elements:,}",
        f"  OCSMesh   : NP={ocs.n_nodes:,}  NE={ocs.n_elements:,}",
        f"  ratio (OCS / ref) : NP x{ocs.n_nodes / ref.n_nodes:.3f}  "
        f"NE x{ocs.n_elements / ref.n_elements:.3f}",
        "",
        "[Bounding box (lon/lat deg)]",
        f"  reference : "
        f"x[{ref.bbox[0]:.4f}, {ref.bbox[2]:.4f}]  "
        f"y[{ref.bbox[1]:.4f}, {ref.bbox[3]:.4f}]",
        f"  OCSMesh   : "
        f"x[{ocs.bbox[0]:.4f}, {ocs.bbox[2]:.4f}]  "
        f"y[{ocs.bbox[1]:.4f}, {ocs.bbox[3]:.4f}]",
        "",
        "[Edge length (m, haversine)]",
        _stats("  reference", ref_e, ".1f"),
        _stats("  OCSMesh  ", ocs_e, ".1f"),
        "",
        "[Triangle alpha-quality (1 = equilateral)]",
        _stats("  reference", ref_q),
        _stats("  OCSMesh  ", ocs_q),
        f"  reference frac alpha < 0.3 : {(ref_q < 0.3).mean() * 100:.2f} %",
        f"  OCSMesh   frac alpha < 0.3 : {(ocs_q < 0.3).mean() * 100:.2f} %",
        "",
        "[Triangle minimum interior angle (deg)]",
        _stats("  reference", ref_a, ".2f"),
        _stats("  OCSMesh  ", ocs_a, ".2f"),
        f"  reference frac min-angle < 20 deg : {(ref_a < 20).mean() * 100:.2f} %",
        f"  OCSMesh   frac min-angle < 20 deg : {(ocs_a < 20).mean() * 100:.2f} %",
        "",
        "[Boundary structure]",
        "  reference :",
        _boundary_summary(ref),
        "  OCSMesh   :",
        _boundary_summary(ocs),
    ]
    summary = "\n".join(lines)
    print(summary)
    SUMMARY_TXT.write_text(summary + "\n", encoding="utf-8")
    print(f"[06] wrote {SUMMARY_TXT}")

    plot_meshes_side_by_side(ref, ocs, SIDE_BY_SIDE_PNG)
    print(f"[06] wrote {SIDE_BY_SIDE_PNG}")

    plot_edge_len_hist(ref_e, ocs_e, EDGE_LEN_PNG)
    print(f"[06] wrote {EDGE_LEN_PNG}")

    plot_quality_hist(ref_q, ocs_q, QUALITY_PNG)
    print(f"[06] wrote {QUALITY_PNG}")


if __name__ == "__main__":
    main()
