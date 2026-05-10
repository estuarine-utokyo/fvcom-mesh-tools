"""PoC #25: characterise the 440 over-connected nodes in PoC #16.

Background: PoC #24's mesh-check found 440 nodes with valence > 8 in
``16_tokyo_bay_with_rivers.14`` (max valence = 26). FVCOM's
``MAX_NBR_ELEM`` cap is typically 8-12 depending on the build, so this
mesh is borderline-to-fatal. PoC #16 was built with the OCSMesh+gmsh
engine plus ``--quality-pass 6 --refine-min-angle 20`` and 5 river
inflows; PoC #19 (oceanmesh engine, otherwise similar inputs) had only
3 over-connected nodes (max valence 9). The cause is therefore in the
OCSMesh+post-processing path, not in the inputs themselves.

This PoC characterises *where* the high-valence nodes sit and *what*
their local triangulation looks like, so that an upstream mitigation
can be proposed without re-running multiple mesh-build ablations.

Approach:
    1. Classify each over-connected node by boundary type (open / land
       coast / river ibtype=21 / interior).
    2. Distance to the nearest river-inflow input point (CSV).
    3. Valence histogram.
    4. Bay-wide map: all OC nodes coloured by valence, river points
       overlaid.
    5. Zoom panels: 1-ring patch around the top-N worst-valence nodes,
       so the local triangulation pattern is visible.

Inputs:
    outputs/16_tokyo_bay_with_rivers.14
    outputs/24_mesh_check_16_tokyo_bay_with_rivers_diag.json
    data/rivers/tokyo_bay/tokyo_bay_rivers.csv

Outputs:
    outputs/25_overconnected_summary.txt
    outputs/25_overconnected_map.png
    outputs/25_overconnected_zoom.png
"""

from __future__ import annotations

import json
from pathlib import Path

import matplotlib

matplotlib.use("Agg")

import matplotlib.pyplot as plt  # noqa: E402
import numpy as np  # noqa: E402
from matplotlib.collections import PolyCollection  # noqa: E402

from fvcom_mesh_tools.io import load_river_points, read_fort14  # noqa: E402
from fvcom_mesh_tools.plotting import MESH_PNG_DPI  # noqa: E402

REPO_ROOT = Path(__file__).resolve().parent.parent
MESH_PATH = REPO_ROOT / "outputs" / "16_tokyo_bay_with_rivers.14"
DIAG_JSON = (
    REPO_ROOT / "outputs" / "24_mesh_check_16_tokyo_bay_with_rivers_diag.json"
)
RIVERS_CSV = REPO_ROOT / "data" / "rivers" / "tokyo_bay" / "tokyo_bay_rivers.csv"

OUT_DIR = REPO_ROOT / "outputs"
SUMMARY_TXT = OUT_DIR / "25_overconnected_summary.txt"
MAP_PNG = OUT_DIR / "25_overconnected_map.png"
ZOOM_PNG = OUT_DIR / "25_overconnected_zoom.png"

EARTH_R_M = 6_371_000.0
N_ZOOM_PATCHES = 6   # how many worst-valence nodes to plot zoom-in for
ZOOM_RADIUS = 0.005  # degrees padding around centre


def haversine_m(p0: np.ndarray, p1: np.ndarray) -> np.ndarray:
    """Pairwise great-circle distance in metres for matched (lon, lat) arrays."""
    lon0 = np.deg2rad(p0[..., 0])
    lat0 = np.deg2rad(p0[..., 1])
    lon1 = np.deg2rad(p1[..., 0])
    lat1 = np.deg2rad(p1[..., 1])
    dlat = lat1 - lat0
    dlon = lon1 - lon0
    a = np.sin(dlat / 2) ** 2 + np.cos(lat0) * np.cos(lat1) * np.sin(dlon / 2) ** 2
    return 2 * EARTH_R_M * np.arcsin(np.sqrt(a))


def min_distance_to_points_m(query: np.ndarray, anchors: np.ndarray) -> np.ndarray:
    """For each lon/lat in ``query``, return the minimum great-circle
    distance (m) to the lon/lat in ``anchors``.
    """
    if anchors.size == 0:
        return np.full(len(query), np.inf)
    out = np.empty(len(query), dtype=np.float64)
    for i, q in enumerate(query):
        out[i] = haversine_m(np.broadcast_to(q, anchors.shape), anchors).min()
    return out


def classify_boundary_type(mesh, node_ids: np.ndarray) -> dict[str, np.ndarray]:
    """Classify each node id as 'open', 'river' (ibtype==21), 'coast' (other
    land), or 'interior'. Returns a dict mapping the 4 categories to bool
    masks aligned with ``node_ids``.
    """
    open_set = set()
    for seg in mesh.open_boundaries:
        open_set.update(int(n) for n in np.asarray(seg).ravel())
    river_set = set()
    coast_set = set()
    for ib, seg in mesh.land_boundaries:
        target = river_set if int(ib) == 21 else coast_set
        target.update(int(n) for n in np.asarray(seg).ravel())

    is_open = np.fromiter(
        (int(i) in open_set for i in node_ids), dtype=bool, count=len(node_ids),
    )
    # An OB node may also appear in a land segment at the OB-coast junction;
    # priority is open > river > coast > interior so the classification is
    # one-hot.
    is_river = np.fromiter(
        (int(i) in river_set for i in node_ids), dtype=bool, count=len(node_ids),
    ) & ~is_open
    is_coast = np.fromiter(
        (int(i) in coast_set for i in node_ids), dtype=bool, count=len(node_ids),
    ) & ~is_open & ~is_river
    is_interior = ~(is_open | is_river | is_coast)
    return {
        "open": is_open, "river": is_river,
        "coast": is_coast, "interior": is_interior,
    }


def one_ring_elements(mesh, node_id: int) -> np.ndarray:
    """Element ids that contain ``node_id``."""
    return np.where((mesh.elements == node_id).any(axis=1))[0]


def _draw_patch(ax, mesh, node_id: int, valence: int) -> None:
    elems = one_ring_elements(mesh, node_id)
    coords = mesh.nodes[mesh.elements[elems]]
    pc = PolyCollection(
        coords, facecolor="white", edgecolor="0.3", linewidths=0.6,
    )
    ax.add_collection(pc)
    centre = mesh.nodes[node_id]
    ax.scatter(*centre, color="tab:red", s=40, zorder=4)
    # mark each one-ring neighbour
    nbr_ids = np.unique(mesh.elements[elems])
    nbr_ids = nbr_ids[nbr_ids != node_id]
    ax.scatter(
        mesh.nodes[nbr_ids, 0], mesh.nodes[nbr_ids, 1],
        color="tab:blue", s=14, zorder=3,
    )
    # boundary segments locally
    for seg in mesh.open_boundaries:
        s = np.asarray(seg)
        ax.plot(mesh.nodes[s, 0], mesh.nodes[s, 1], "-", color="tab:red", lw=1.0)
    for ib, seg in mesh.land_boundaries:
        s = np.asarray(seg)
        col = "tab:blue" if int(ib) == 21 else "0.5"
        lw = 1.4 if int(ib) == 21 else 0.7
        ax.plot(mesh.nodes[s, 0], mesh.nodes[s, 1], "-", color=col, lw=lw)
    pad = ZOOM_RADIUS
    ax.set_xlim(centre[0] - pad, centre[0] + pad)
    ax.set_ylim(centre[1] - pad, centre[1] + pad)
    ax.set_aspect("equal")
    ax.set_title(
        f"node {node_id}  valence={valence}  "
        f"({centre[0]:.4f}, {centre[1]:.4f})",
        fontsize=9,
    )
    ax.tick_params(labelsize=7)


def main() -> None:
    for p in (MESH_PATH, DIAG_JSON, RIVERS_CSV):
        if not p.exists():
            raise SystemExit(f"required input missing: {p}")
    OUT_DIR.mkdir(parents=True, exist_ok=True)

    mesh = read_fort14(MESH_PATH)
    diag = json.loads(DIAG_JSON.read_text(encoding="utf-8"))
    oc_records = diag["detectors"]["overconnected_nodes"]["nodes"]
    if not oc_records:
        SUMMARY_TXT.write_text("no over-connected nodes flagged\n", encoding="utf-8")
        print("no over-connected nodes flagged; nothing to do")
        return
    node_ids = np.asarray([r["id"] for r in oc_records], dtype=np.int64)
    valences = np.asarray([r["valence"] for r in oc_records], dtype=np.int64)
    coords = np.asarray(
        [[r["lon"], r["lat"]] for r in oc_records], dtype=np.float64,
    )
    rivers = load_river_points([RIVERS_CSV])
    print(f"[25] {len(node_ids)} over-connected nodes; {len(rivers)} river inputs")

    cats = classify_boundary_type(mesh, node_ids)
    dist_to_river = min_distance_to_points_m(coords, rivers)

    # Valence histogram.
    bins = np.arange(valences.min(), valences.max() + 2)

    # --- summary text ---
    lines = [
        f"mesh:   {MESH_PATH}",
        f"diag:   {DIAG_JSON}",
        f"rivers: {RIVERS_CSV}",
        f"n_over_connected = {len(node_ids)}",
        f"max valence      = {int(valences.max())}",
        "",
        "[Boundary classification]",
    ]
    n = len(node_ids)
    for label in ("open", "river", "coast", "interior"):
        m = int(cats[label].sum())
        lines.append(f"  {label:<9s}: {m:4d}  ({m / n * 100:5.1f} %)")
    lines += [
        "",
        "[Valence histogram]",
    ]
    counts, edges = np.histogram(valences, bins=bins)
    for v, c in zip(edges[:-1], counts):
        lines.append(f"  valence={int(v):2d}: {int(c):4d}")
    lines += [
        "",
        "[Distance to nearest river-inflow point (m)]",
        f"  min/p25/p50/p75/max : "
        f"{dist_to_river.min():.0f} / "
        f"{np.percentile(dist_to_river, 25):.0f} / "
        f"{np.percentile(dist_to_river, 50):.0f} / "
        f"{np.percentile(dist_to_river, 75):.0f} / "
        f"{dist_to_river.max():.0f}",
        f"  fraction within 1 km : "
        f"{(dist_to_river < 1000).mean() * 100:.1f} %",
        f"  fraction within 5 km : "
        f"{(dist_to_river < 5000).mean() * 100:.1f} %",
        "",
        "[Top-10 worst-valence nodes (id, valence, lon, lat, dist_to_river_m, class)]",
    ]
    order = np.argsort(-valences)[:10]
    for k in order:
        cls = (
            "open" if cats["open"][k]
            else "river" if cats["river"][k]
            else "coast" if cats["coast"][k]
            else "interior"
        )
        lines.append(
            f"  {int(node_ids[k]):6d}  v={int(valences[k]):2d}  "
            f"({coords[k, 0]:.4f}, {coords[k, 1]:.4f})  "
            f"d_river={dist_to_river[k]:7.0f} m  {cls}"
        )
    summary = "\n".join(lines)
    print(summary)
    SUMMARY_TXT.write_text(summary + "\n", encoding="utf-8")
    print(f"[25] wrote {SUMMARY_TXT}")

    # --- bay-wide map ---
    fig, ax = plt.subplots(figsize=(11, 9), dpi=120)
    ax.triplot(
        mesh.nodes[:, 0], mesh.nodes[:, 1], mesh.elements,
        color="0.85", lw=0.1,
    )
    for seg in mesh.open_boundaries:
        s = np.asarray(seg)
        ax.plot(mesh.nodes[s, 0], mesh.nodes[s, 1], "-", color="tab:red", lw=0.9)
    for ib, seg in mesh.land_boundaries:
        s = np.asarray(seg)
        col = "tab:blue" if int(ib) == 21 else "0.4"
        lw = 1.4 if int(ib) == 21 else 0.4
        ax.plot(mesh.nodes[s, 0], mesh.nodes[s, 1], "-", color=col, lw=lw)
    sc = ax.scatter(
        coords[:, 0], coords[:, 1], c=valences, cmap="magma",
        s=18, edgecolors="black", linewidths=0.3, zorder=4,
    )
    cb = fig.colorbar(sc, ax=ax, shrink=0.7)
    cb.set_label("node valence")
    ax.scatter(
        rivers[:, 0], rivers[:, 1], marker="x", color="tab:purple",
        s=80, linewidths=2.0, label="river input pt", zorder=5,
    )
    ax.set_aspect("equal")
    ax.set_xlabel("lon (deg)")
    ax.set_ylabel("lat (deg)")
    ax.set_title(
        f"PoC #25 over-connected nodes  "
        f"n={len(node_ids)}  max_valence={int(valences.max())}"
    )
    ax.legend(loc="lower right", fontsize=8)
    fig.tight_layout()
    fig.savefig(MAP_PNG, bbox_inches="tight", dpi=MESH_PNG_DPI)
    plt.close(fig)
    print(f"[25] wrote {MAP_PNG}")

    # --- zoom panels for the worst N ---
    n_zoom = min(N_ZOOM_PATCHES, len(node_ids))
    cols = 3
    rows = (n_zoom + cols - 1) // cols
    fig, axes = plt.subplots(
        rows, cols, figsize=(4.5 * cols, 4.5 * rows), dpi=120,
    )
    axes_arr = np.atleast_1d(axes).ravel()
    for ax in axes_arr:
        ax.set_visible(False)
    worst = np.argsort(-valences)[:n_zoom]
    for ax, k in zip(axes_arr, worst):
        ax.set_visible(True)
        _draw_patch(ax, mesh, int(node_ids[k]), int(valences[k]))
    fig.suptitle(
        "Local 1-ring around top "
        f"{n_zoom} highest-valence nodes (red=node, blue=neighbour)",
        y=1.0, fontsize=11,
    )
    fig.tight_layout()
    fig.savefig(ZOOM_PNG, bbox_inches="tight", dpi=MESH_PNG_DPI)
    plt.close(fig)
    print(f"[25] wrote {ZOOM_PNG}")


if __name__ == "__main__":
    main()
