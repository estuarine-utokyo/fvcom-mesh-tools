"""PoC #28: channel-width / h ratio detector for under-resolved channels.

The existing :func:`fvcom_mesh_tools.diagnostics.thin_chain_elements_flag`
catches only 1-cell-wide channels (every triangle has 3 boundary
vertices). FVCOM flow representation in narrow water bodies typically
needs at least 3 cells across, so a metric that flags 2-cell-wide
channels too is useful. PoC #28 explores a medial-axis-style
channel-width / h ratio detector.

Approach
--------

1. Project mesh from (lon, lat) to local metric coordinates via a
   flat-earth approximation around the mesh's mid-latitude.
2. Sample each boundary polyline (open + land) at fine spacing.
3. For each element, query the distance from its centroid to each
   polyline (one ``cKDTree`` per polyline; ``k=1``).
4. Channel width estimate at the centroid: sum of distances to the two
   nearest *distinct* polylines.
5. Local h: median edge length of the element's 3 edges, in metres.
6. ``w_h_ratio = channel_width / h_local``.
7. Flag elements with ``w_h_ratio < min_w_h`` (default 3 = "channel
   should be at least 3 cells across").

Compares the new flag with ``thin_chain_elements_flag`` to see whether
the medial-axis detector catches additional under-resolved channels
that are 2- or 3-cell wide rather than only 1-cell.

Inputs:
    outputs/19_tokyo_bay_oceanmesh.14            (raw, 144 components)
    outputs/19_tokyo_bay_oceanmesh_cleaned.14    (after A+B+C+D)

Outputs:
    outputs/28_channel_width_summary.txt
    outputs/28_channel_width_<name>_map.png
"""

from __future__ import annotations

from pathlib import Path

import matplotlib

matplotlib.use("Agg")

import matplotlib.pyplot as plt  # noqa: E402
import numpy as np  # noqa: E402
from scipy.spatial import cKDTree  # noqa: E402

from fvcom_mesh_tools.diagnostics import (  # noqa: E402
    face_face_adjacency,
    thin_chain_elements_flag,
    thin_elements_flag,
)
from fvcom_mesh_tools.io import Fort14Mesh, read_fort14  # noqa: E402
from fvcom_mesh_tools.plotting import MESH_PNG_DPI  # noqa: E402

REPO_ROOT = Path(__file__).resolve().parent.parent
OUT_DIR = REPO_ROOT / "outputs"
SUMMARY_TXT = OUT_DIR / "28_channel_width_summary.txt"

CASES: list[tuple[str, Path]] = [
    ("19_raw",     OUT_DIR / "19_tokyo_bay_oceanmesh.14"),
    ("19_cleaned", OUT_DIR / "19_tokyo_bay_oceanmesh_cleaned.14"),
]

EARTH_R_M = 6_371_000.0
SAMPLE_DS_M = 50.0     # boundary-sample spacing in metres
THRESHOLDS = (2, 3, 5)
DEFAULT_THRESHOLD = 3


# ---------------------------------------------------------------------------
# Geometry helpers
# ---------------------------------------------------------------------------


def to_metric(nodes_lonlat: np.ndarray, *, lat0: float, lon0: float) -> np.ndarray:
    """Flat-earth projection (lon, lat) -> (x_m, y_m) about (lon0, lat0)."""
    cos_lat0 = np.cos(np.deg2rad(lat0))
    x = (nodes_lonlat[:, 0] - lon0) * np.deg2rad(1.0) * EARTH_R_M * cos_lat0
    y = (nodes_lonlat[:, 1] - lat0) * np.deg2rad(1.0) * EARTH_R_M
    return np.column_stack([x, y])


def edge_lengths_m(nodes_m: np.ndarray, elements: np.ndarray) -> np.ndarray:
    """Per-element ``(NE, 3)`` edge lengths in metres (already metric input)."""
    p0 = nodes_m[elements[:, 0]]
    p1 = nodes_m[elements[:, 1]]
    p2 = nodes_m[elements[:, 2]]
    return np.column_stack([
        np.linalg.norm(p1 - p0, axis=1),
        np.linalg.norm(p2 - p1, axis=1),
        np.linalg.norm(p0 - p2, axis=1),
    ])


def sample_polyline(
    nodes_m: np.ndarray, seg_ids: np.ndarray, ds: float,
) -> np.ndarray:
    """Sample a polyline (sequence of node ids) at ~``ds`` spacing in m.

    Returns ``(N, 2)`` metric coords; always at least 2 samples per edge.
    """
    pts: list[np.ndarray] = []
    seg_ids = np.asarray(seg_ids, dtype=np.int64)
    for i in range(len(seg_ids) - 1):
        p0 = nodes_m[seg_ids[i]]
        p1 = nodes_m[seg_ids[i + 1]]
        L = float(np.linalg.norm(p1 - p0))
        n = max(2, int(np.ceil(L / ds)) + 1)
        t = np.linspace(0.0, 1.0, n)
        pts.append(p0 + t[:, None] * (p1 - p0))
    if not pts:
        return np.empty((0, 2), dtype=np.float64)
    return np.vstack(pts)


# ---------------------------------------------------------------------------
# Channel-width detector
# ---------------------------------------------------------------------------


def detect_channel_width(
    mesh: Fort14Mesh, *,
    sample_ds_m: float = SAMPLE_DS_M,
) -> dict[str, np.ndarray]:
    """Compute per-element channel width, local h, and the w/h ratio.

    Returns a dict with keys ``"channel_width_m"``, ``"h_local_m"``,
    ``"w_h_ratio"``, and ``"d_to_polylines"`` (distance to each polyline,
    shape ``(NE, n_polylines)``).
    """
    if not mesh.open_boundaries and not mesh.land_boundaries:
        raise ValueError("mesh has no boundaries; cannot compute channel width")

    lat0 = float(mesh.nodes[:, 1].mean())
    lon0 = float(mesh.nodes[:, 0].mean())
    nodes_m = to_metric(mesh.nodes, lat0=lat0, lon0=lon0)

    polylines: list[np.ndarray] = list(mesh.open_boundaries)
    polylines += [seg for _ib, seg in mesh.land_boundaries]
    n_polylines = len(polylines)

    trees: list[cKDTree] = []
    for seg in polylines:
        samples = sample_polyline(nodes_m, seg, sample_ds_m)
        trees.append(cKDTree(samples))

    centroids = nodes_m[mesh.elements].mean(axis=1)
    dists = np.empty((mesh.n_elements, n_polylines), dtype=np.float64)
    for k, tree in enumerate(trees):
        d, _ = tree.query(centroids, k=1)
        dists[:, k] = d

    # Channel width = sum of distances to 2 nearest distinct polylines.
    if n_polylines >= 2:
        sorted_dists = np.sort(dists, axis=1)
        channel_width = sorted_dists[:, 0] + sorted_dists[:, 1]
    else:
        channel_width = 2.0 * dists[:, 0]   # only one polyline; use 2x as proxy

    edge_lens = edge_lengths_m(nodes_m, mesh.elements)
    h_local = np.median(edge_lens, axis=1)

    ratio = channel_width / np.where(h_local > 0, h_local, 1.0)
    return {
        "channel_width_m": channel_width,
        "h_local_m": h_local,
        "w_h_ratio": ratio,
        "d_to_polylines": dists,
    }


# ---------------------------------------------------------------------------
# Driver
# ---------------------------------------------------------------------------


def _summary_lines(name: str, mesh: Fort14Mesh, det: dict[str, np.ndarray]) -> list[str]:
    ratio = det["w_h_ratio"]
    cw = det["channel_width_m"]
    h = det["h_local_m"]

    thin = thin_elements_flag(mesh)
    adj = face_face_adjacency(mesh.elements)
    chain = thin_chain_elements_flag(adj, thin, min_chain_length=3)

    lines = [
        f"=== {name} ===",
        f"NP={mesh.n_nodes:,}  NE={mesh.n_elements:,}",
        f"channel_width_m  p10={np.percentile(cw, 10):8.0f}  "
        f"p50={np.percentile(cw, 50):8.0f}  "
        f"p90={np.percentile(cw, 90):8.0f}",
        f"h_local_m        p10={np.percentile(h, 10):8.0f}  "
        f"p50={np.percentile(h, 50):8.0f}  "
        f"p90={np.percentile(h, 90):8.0f}",
        f"w_h_ratio        p10={np.percentile(ratio, 10):8.2f}  "
        f"p50={np.percentile(ratio, 50):8.2f}  "
        f"p90={np.percentile(ratio, 90):8.2f}",
        "",
        "[Element counts at thresholds: w_h_ratio < N]",
    ]
    for thr in THRESHOLDS:
        n_below = int((ratio < thr).sum())
        pct = n_below / mesh.n_elements * 100
        lines.append(f"  ratio < {thr}: {n_below:6,d} elements  ({pct:5.2f} %)")
    lines += [
        "",
        "[Comparison with existing detectors]",
        f"  thin elements (any 3-bdy):   {int(thin.sum()):,}",
        f"  thin chain (>= 3):           {int(chain.sum()):,}",
        f"  ratio < 3 (this PoC):        {int((ratio < 3).sum()):,}",
        "",
        "[Co-occurrence (ratio < 3 AND thin chain): "
        f"{int(((ratio < 3) & chain).sum()):,} elements]",
        f"  ratio < 3 only (not thin chain): "
        f"{int(((ratio < 3) & ~chain).sum()):,}",
        f"  thin chain only (not ratio < 3): "
        f"{int((chain & ~(ratio < 3)).sum()):,}",
    ]
    return lines


def _plot_map(
    name: str, mesh: Fort14Mesh, det: dict[str, np.ndarray], png: Path,
) -> None:
    ratio = det["w_h_ratio"]
    fig, ax = plt.subplots(figsize=(10, 8), dpi=120)
    # Cap at p99 for readable colour scale
    vmax = float(np.percentile(ratio, 99))
    tpc = ax.tripcolor(
        mesh.nodes[:, 0], mesh.nodes[:, 1], mesh.elements,
        facecolors=np.minimum(ratio, vmax),
        cmap="magma", edgecolors="none",
        vmin=0, vmax=vmax,
    )
    cb = fig.colorbar(tpc, ax=ax, shrink=0.7)
    cb.set_label(f"w/h ratio (capped at p99 = {vmax:.1f})")
    # Highlight under-resolved (ratio < 3) elements with a thick edge.
    bad = ratio < DEFAULT_THRESHOLD
    if bad.any():
        # Plot only flagged elements as outlined triangles
        from matplotlib.collections import PolyCollection
        tris = mesh.nodes[mesh.elements[bad]]
        pc = PolyCollection(
            tris, facecolor="none", edgecolor="tab:cyan",
            linewidths=0.4, alpha=0.6,
        )
        ax.add_collection(pc)

    for seg in mesh.open_boundaries:
        s = np.asarray(seg)
        ax.plot(mesh.nodes[s, 0], mesh.nodes[s, 1], "-", color="tab:red", lw=0.6)
    for ib, seg in mesh.land_boundaries:
        s = np.asarray(seg)
        col = "tab:blue" if int(ib) == 21 else "0.4"
        lw = 1.0 if int(ib) == 21 else 0.3
        ax.plot(mesh.nodes[s, 0], mesh.nodes[s, 1], "-", color=col, lw=lw)

    ax.set_aspect("equal")
    ax.set_xlabel("lon (deg)")
    ax.set_ylabel("lat (deg)")
    ax.set_title(
        f"PoC #28 channel-width / h ratio  {name}  "
        f"NP={mesh.n_nodes:,} NE={mesh.n_elements:,}"
        f"\nratio<{DEFAULT_THRESHOLD}: "
        f"{int((ratio < DEFAULT_THRESHOLD).sum()):,} elements (cyan outline)"
    )
    fig.tight_layout()
    fig.savefig(png, bbox_inches="tight", dpi=MESH_PNG_DPI)
    plt.close(fig)


def main() -> None:
    OUT_DIR.mkdir(parents=True, exist_ok=True)
    sections: list[str] = []
    for name, path in CASES:
        if not path.exists():
            sections.append(f"=== {name} ===\nSKIP (missing {path})")
            continue
        print(f"[28] === case={name}  reading {path.name} ===")
        mesh = read_fort14(path)
        det = detect_channel_width(mesh)
        section = "\n".join(_summary_lines(name, mesh, det))
        print(section)
        print()
        sections.append(section)
        png = OUT_DIR / f"28_channel_width_{name}_map.png"
        _plot_map(name, mesh, det, png)
        print(f"[28] wrote {png}")
    summary = "\n\n".join(sections) + "\n"
    SUMMARY_TXT.write_text(summary, encoding="utf-8")
    print(f"\n[28] wrote {SUMMARY_TXT}")


if __name__ == "__main__":
    main()
