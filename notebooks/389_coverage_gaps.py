"""Measure resolvable omitted water; run via jobs/octopus/389_coverage_gaps.sh.

The domain follows 325_sample_repro, clipped again at the certified mesh's
actual open boundary. Original OSM water, policy fill and defects stay distinct.
"""
from __future__ import annotations

import argparse
import hashlib
import json
import os
import sys
from pathlib import Path

import geopandas as gpd
import matplotlib

matplotlib.use("Agg")
import matplotlib.pyplot as plt  # noqa: E402
import numpy as np  # noqa: E402
from matplotlib.patches import Rectangle  # noqa: E402
from pyproj import Transformer  # noqa: E402
from shapely.geometry import LineString, Polygon, box, shape  # noqa: E402
from shapely.ops import split, transform, unary_union  # noqa: E402

ROOT = Path(__file__).resolve().parents[1]
sys.path.insert(0, str(ROOT / "src"))
from fvcom_mesh_tools.coverage import coverage_gaps  # noqa: E402
from fvcom_mesh_tools.gridref import TOKYO_BAY_GRID  # noqa: E402
from fvcom_mesh_tools.io import read_fort14  # noqa: E402
from fvcom_mesh_tools.plotting import add_atlas_grid, use_readable_style  # noqa: E402

OUT = ROOT / "outputs/sample_repro"
FIG = ROOT / "outputs/figures/389_coverage_gaps.png"
UTM = Transformer.from_crs(4326, 32654, always_xy=True)
LL = Transformer.from_crs(32654, 4326, always_xy=True)


def domain_of_interest(mesh):
    """325's domain, with any seaward sliver at the actual OBC removed."""
    arc = [[139.6713, 35.1396], [139.6737, 35.1288], [139.6772, 35.1168],
           [139.6816, 35.1031], [139.6871, 35.0877], [139.6946, 35.0705],
           [139.7000, 35.0576], [139.7069, 35.0445], [139.7134, 35.0327],
           [139.7216, 35.0184], [139.7289, 35.0047], [139.7373, 34.9916],
           [139.7497, 34.9750]]
    domain = transform(UTM.transform, Polygon(
        [[139.83, 34.973], [140.12, 34.973], [140.12, 35.75],
         [139.60, 35.75], [139.60, 35.20], [139.6642, 35.1546]] + arc))
    if len(mesh.open_boundaries) != 1:
        raise ValueError("certified runner expects exactly one open-boundary polyline")
    points = mesh.nodes[mesh.open_boundaries[0]]
    if len(points) < 2:
        raise ValueError("open boundary needs at least two nodes")
    span = 4 * np.hypot(domain.bounds[2] - domain.bounds[0],
                        domain.bounds[3] - domain.bounds[1])
    start, end = points[0] - points[1], points[-1] - points[-2]
    line = LineString(np.vstack([points[0] + span * start / np.linalg.norm(start),
                                points, points[-1] + span * end / np.linalg.norm(end)]))
    parts = list(split(domain, line).geoms)
    covered = unary_union([Polygon(p) for p in mesh.nodes[mesh.elements]])
    # The bay side contains essentially all mesh area. Do not use a hull:
    # original bays beyond the meshed shoreline must remain in the domain.
    return max(parts, key=lambda p: p.intersection(covered).area)


def locate(rec):
    lon, lat = LL.transform(*rec["xy_m"])
    rec["lon_lat"] = [lon, lat]
    try:
        rec["gridref"] = TOKYO_BAY_GRID.point_to_subcell(lon, lat)
    except ValueError:
        rec["gridref"] = "outside-atlas"


def draw(mesh, land, report):
    use_readable_style()
    worst = report["regions"][:9]
    fig = plt.figure(figsize=(23, 19))
    grid = fig.add_gridspec(3, 4)
    overview = fig.add_subplot(grid[:, 0])
    xy = np.column_stack(LL.transform(mesh.nodes[:, 0], mesh.nodes[:, 1]))
    coast = gpd.GeoSeries([transform(LL.transform, land)], crs=4326)
    colours = {"clear-defect": "#d7191c", "marginal": "#e66101",
               "subtarget": "#2c7bb6", "unmeasured": "#7b3294"}

    def base(ax, bounds, lw):
        coast.clip(box(*bounds)).plot(ax=ax, color="0.88", edgecolor="0.55", linewidth=0.4)
        ax.triplot(xy[:, 0], xy[:, 1], mesh.elements, color="0.4", linewidth=lw, zorder=4)
        ax.set_xlim(bounds[0], bounds[2])
        ax.set_ylim(bounds[1], bounds[3])
        ax.set_aspect(1 / np.cos(np.deg2rad(35.35)))
        add_atlas_grid(ax, crs="EPSG:4326")

    bounds = transform(LL.transform, shape(report["domain_geometry"])).bounds
    base(overview, bounds, 0.25)
    overview.set_title(f"Omitted original water\n{report['clear_defects']} clear defects")
    for i, rec in enumerate(worst):
        poly = transform(LL.transform, shape(rec["geometry"]))
        x0, y0, x1, y1 = poly.bounds
        pad = max(0.002, (x1 - x0) * 0.1, (y1 - y0) * 0.1)
        bounds = (x0 - pad, y0 - pad, x1 + pad, y1 + pad)
        colour = colours[rec["classification"]]
        overview.add_patch(Rectangle(bounds[:2], bounds[2] - bounds[0], bounds[3] - bounds[1],
                                     fill=False, edgecolor=colour, linewidth=2, zorder=6))
        overview.text(bounds[0], bounds[3], str(rec["id"]), color=colour,
                      fontweight="bold", bbox=dict(fc="white", ec=colour, alpha=0.9), zorder=7)
        ax = fig.add_subplot(grid[i // 3, 1 + i % 3])
        base(ax, bounds, 0.65)
        gpd.GeoSeries([poly], crs=4326).plot(ax=ax, color=colour, alpha=0.6, zorder=3)
        policy = [transform(LL.transform, shape(p["geometry"])) for p in report["policy_filled"]
                  if shape(p["geometry"]).intersects(shape(rec["geometry"]).buffer(1000))]
        if policy:
            clipped = gpd.GeoSeries(policy, crs=4326).clip(box(*bounds))
            clipped = clipped[~clipped.is_empty]
            if len(clipped):   # geopandas cannot set an aspect from nothing
                clipped.plot(ax=ax, color="#888888", alpha=0.5, hatch="///", zorder=2)
        ratio = rec.get("w_h_p50")
        detail = f"w/h={ratio:.2f}" if ratio is not None else "unmeasured"
        ax.set_title(f"[{rec['id']}] {rec['gridref']} {rec['classification']}\n{detail}",
                     color=colour, fontsize=12)
    fig.suptitle("Coverage gaps: coloured original water, mesh edges overlaid; "
                 "grey hatching = policy-filled", fontsize=16)
    fig.tight_layout(rect=(0, 0, 1, 0.97))
    FIG.parent.mkdir(parents=True, exist_ok=True)
    fig.savefig(FIG, dpi=140)
    plt.close(fig)


def main():
    parser = argparse.ArgumentParser(description=__doc__)
    parser.add_argument("--pixel-size-m", type=float, default=None)
    parser.add_argument("--no-figure", action="store_true")
    args = parser.parse_args()
    mesh_path = OUT / "sample_repro_final.14"
    original_path = ROOT / "outputs/tb_varres_3r/land_osm_wide.shp"
    policy_path = OUT / "land_channel_adj.shp"
    normalize_path = OUT / "normalize.json"
    destination = OUT / "coverage_gaps.json"
    for path in [destination] + ([] if args.no_figure else [FIG]):
        if path.exists():
            raise FileExistsError(f"preserve existing output; move it before rerunning: {path}")
    mesh = read_fort14(mesh_path)
    land = unary_union(gpd.read_file(original_path).to_crs(32654).geometry)
    policy = unary_union(gpd.read_file(policy_path).to_crs(32654).geometry)
    normalization = json.loads(normalize_path.read_text())
    domain = domain_of_interest(mesh)
    report = coverage_gaps(mesh, land, domain, policy_land=policy,
                           pixel_size_m=args.pixel_size_m,
                           h_target_m=float(os.environ.get("SR_H0", 290.0)))
    from shapely.geometry import mapping

    report["domain_geometry"] = mapping(domain)
    report["crs"] = "EPSG:32654"
    report["mesh_counts"] = {"nodes": mesh.n_nodes, "elements": mesh.n_elements}
    report["inputs"] = {str(p): {"sha256": hashlib.sha256(p.read_bytes()).hexdigest(),
                                 "mtime_ns": p.stat().st_mtime_ns}
                        for p in (mesh_path, original_path, policy_path, normalize_path)}
    report["normalization_log"] = normalization
    report["policy_note"] = (
        "Footprints use adjusted land minus original land, restricted to unmeshed water. "
        "normalize.json centres/areas are retained as context, not exact footprints. "
        "Input hashes record provenance but cannot prove these files came from the same run.")
    for rec in report["regions"] + report["policy_filled"]:
        locate(rec)
    print("ID atlas       lon       lat       class          area_m2 h_m W50_m w/h rows "
          "C1 C2 connected length_m inland_m", flush=True)
    for r in report["regions"][:15]:
        def number(key):
            value = r.get(key)
            if value is None:
                return "NA"
            if key == "rows_p50":
                return str(value)
            return f"{value:.3f}" if key == "w_h_p50" else f"{value:.1f}"

        corridor = r.get("corridors", {})
        flags = [str(int(corridor[str(n)]["feasible"])) if str(n) in corridor else "NA"
                 for n in (1, 2)]
        print(f"{r['id']:2} {r['gridref']:10} {r['lon_lat'][0]:.5f} {r['lon_lat'][1]:.5f} "
              f"{r['classification']:14} {r['area_m2']:.0f} {number('h_m')} "
              f"{number('width_p50_m')} {number('w_h_p50')} {number('rows_p50')} "
              f"{' '.join(flags)} {r['touches_mesh']} {number('length_m')} "
              f"{number('inland_reach_m')}", flush=True)
    print(f"clear defects={report['clear_defects']}; unmeasured={report['unmeasured']}; "
          f"below threshold={report['below_area_threshold']}", flush=True)
    for r in report["policy_filled"]:
        print(f"policy-filled {r['id']} {r['gridref']} {r['lon_lat']} "
              f"area_m2={r['area_m2']:.1f}", flush=True)
    destination.write_text(json.dumps(report, indent=2, allow_nan=False) + "\n")
    print(f"saved {destination}", flush=True)
    if not args.no_figure:
        draw(mesh, land, report)
        print(f"saved {FIG}", flush=True)


if __name__ == "__main__":
    main()
