"""PoC #67 — smooth-section boundary-conformity metric.

The user's acceptance criterion is NOT the aggregate boundary-to-
coastline distance (curved harbour geometry makes a chord error of
~h²·κ/8 unavoidable at any element size), but exact conformity on the
SMOOTH sections that a coarse mesh can represent:

* the reference polylines are resampled at ``DS``;
* a sample is "smooth at mesh scale" when the line deviates from the
  chord spanning ±``HMIN/2`` of arc by less than ``SMOOTH_TOL``
  (i.e. an hmin-long element edge can lie on it within that tol);
* every mesh land-boundary node whose nearest reference sample is
  smooth contributes to the metric.

Reports p50/p90/max and a pass count against ``CONFORM_TOL`` for each
mesh given on the command line (defaults compare 60c / 64 / 66).
"""

from __future__ import annotations

import json
import sys
from pathlib import Path

import numpy as np

REPO = Path(__file__).resolve().parents[1]
COAST = REPO / "outputs" / "osm_shoreline" / "osm_true_land_tokyo_bay.shp"
OUT_JSON = REPO / "outputs" / "67_smooth_fidelity.json"

HMIN = 300.0        # m — mesh scale the smoothness test is relative to
DS = 50.0           # m — reference sampling
SMOOTH_TOL = 15.0   # m — max line-vs-chord deviation to call it smooth
CONFORM_TOL = 20.0  # m — smooth-section nodes should sit within this
BAND_DEG = 0.01     # artificial-cut exclusion at domain edges

DEFAULT_MESHES = [
    "outputs/60c_v4_raw_300m_osm.14",
    "outputs/64_v4_raw_opened_h0100.14",
    "outputs/66_v4_raw_opened_mult36.14",
]


def _sample_lines(lines, ds):
    pts, smooth = [], []
    half = HMIN / 2.0
    for ln in lines:
        L = ln.length
        if L < HMIN:
            continue
        s = np.arange(0.0, L + ds / 2, ds)
        P = np.asarray([ln.interpolate(float(t)).coords[0] for t in s])
        A = np.asarray([ln.interpolate(max(0.0, float(t) - half)).coords[0]
                        for t in s])
        B = np.asarray([ln.interpolate(min(L, float(t) + half)).coords[0]
                        for t in s])
        ab = B - A
        nrm = np.hypot(ab[:, 0], ab[:, 1])
        nrm = np.where(nrm > 0, nrm, 1.0)
        dev = np.abs(
            (P[:, 0] - A[:, 0]) * ab[:, 1] - (P[:, 1] - A[:, 1]) * ab[:, 0]
        ) / nrm
        pts.append(P)
        smooth.append(dev <= SMOOTH_TOL)
    return np.vstack(pts), np.concatenate(smooth)


def main() -> int:
    import geopandas as gpd
    import shapely
    from pyproj import Transformer
    from scipy.spatial import cKDTree
    from shapely.strtree import STRtree

    from fvcom_mesh_tools.io import read_fort14

    gdf = gpd.read_file(COAST).to_crs(32654)
    lines = []
    for g in gdf.geometry:
        if g is None or g.is_empty:
            continue
        b = g.boundary
        lines.extend(list(b.geoms) if hasattr(b, "geoms") else [b])
    samples, smooth_mask = _sample_lines(lines, DS)
    print(f"[67] reference samples: {len(samples):,} "
          f"(smooth at {HMIN:.0f} m scale: {smooth_mask.mean() * 100:.1f}%)",
          flush=True)
    sample_tree = cKDTree(samples)
    line_tree = STRtree(lines)
    line_arr = np.array(lines, dtype=object)
    tr = Transformer.from_crs("EPSG:4326", "EPSG:32654", always_xy=True)

    results = {}
    mesh_paths = sys.argv[1:] or DEFAULT_MESHES
    for mp in mesh_paths:
        mp = Path(mp)
        if not mp.exists():
            print(f"[67] skip missing {mp}", flush=True)
            continue
        m = read_fort14(mp)
        lon, lat = m.nodes[:, 0], m.nodes[:, 1]
        land_nodes = np.unique(np.concatenate(
            [np.asarray(s) for _ib, s in m.land_boundaries]
        ))
        keep = (
            (lat[land_nodes] > lat[land_nodes].min() + BAND_DEG)
            & (lon[land_nodes] > lon[land_nodes].min() + BAND_DEG)
            & (lon[land_nodes] < lon[land_nodes].max() - BAND_DEG)
        )
        nodes = land_nodes[keep]
        px, py = tr.transform(lon[nodes], lat[nodes])
        P = np.column_stack([px, py])
        _d_s, idx_s = sample_tree.query(P, workers=-1)
        on_smooth = smooth_mask[idx_s]
        pts = shapely.points(P[on_smooth, 0], P[on_smooth, 1])
        d = shapely.distance(pts, line_arr[line_tree.nearest(pts)])
        res = {
            "n_land_nodes": int(nodes.size),
            "n_smooth_nodes": int(on_smooth.sum()),
            "p50_m": float(np.percentile(d, 50)),
            "p90_m": float(np.percentile(d, 90)),
            "max_m": float(d.max()),
            "n_beyond_conform_tol": int((d > CONFORM_TOL).sum()),
            "frac_within_tol": float((d <= CONFORM_TOL).mean()),
        }
        results[mp.stem] = res
        print(f"[67] {mp.stem}: smooth-nodes {res['n_smooth_nodes']:,}/"
              f"{res['n_land_nodes']:,}  p50={res['p50_m']:.1f} m  "
              f"p90={res['p90_m']:.1f} m  within {CONFORM_TOL:.0f} m: "
              f"{res['frac_within_tol'] * 100:.1f}%", flush=True)

    OUT_JSON.write_text(json.dumps({
        "params": {"hmin_m": HMIN, "ds_m": DS, "smooth_tol_m": SMOOTH_TOL,
                   "conform_tol_m": CONFORM_TOL},
        "results": results,
    }, indent=2), encoding="utf-8")
    print(f"[67] wrote {OUT_JSON}", flush=True)
    return 0


if __name__ == "__main__":
    raise SystemExit(main())
