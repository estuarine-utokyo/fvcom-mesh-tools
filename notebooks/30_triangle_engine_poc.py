"""PoC #30: Tokyo-Bay mesh build via ocsmesh's Triangle engine.

Why
---
``--engine ocsmesh`` currently calls ``MeshDriver(engine_name='gmsh')``
under the hood. PoC #25/26 showed gmsh itself produces ~380
over-connected nodes (max valence 26) on Tokyo Bay before any
post-processing — the dominant source of quality degradation in the
ocsmesh path. ``ocsmesh`` also ships a ``triangle`` engine wrapping
Shewchuk's Triangle (constrained Delaunay + refinement), which is
gmsh-independent and has a different quality profile. If Triangle's
output approaches the oceanmesh path's quality, it would let us drop
the gmsh dependency entirely (license simplification + quality win).

This PoC re-runs PoC #16's exact ocsmesh build configuration on the
same Tokyo Bay inputs, but invokes the Triangle engine instead of
gmsh. It writes a fort.14, computes the same quality metrics
(``alpha mean``, ``frac<20°``, ``min angle p50``, max valence,
over-connected count, flipped triangles), and prints a side-by-side
comparison against PoC #16 (gmsh) and PoC #19 (oceanmesh) reference
numbers.

Inputs (same as PoC #16):
    data/bathymetry/tokyo_bay/dem_00_01_change.nc
    data/coastline/tokyo_bay/MLIT_C23/C23-06_TOKYOBAY.shp
    data/rivers/tokyo_bay/tokyo_bay_rivers.csv

Build parameters (same as PoC #16): hmin=200 m, hmax=5000 m,
coast-target-size=200, coast-expansion-rate=0.005,
min-polygon-area=1e6 m², min-island-area=1e5 m².

Outputs:
    outputs/30_tokyo_bay_triangle.14
    outputs/30_triangle_engine_summary.json
    outputs/30_triangle_engine_summary.txt

The engine is invoked by editing one keyword (``engine_name``) in an
ad-hoc copy of the ocsmesh build path; the production
``mesh_engine/ocsmesh.py`` is left untouched. Promoting Triangle to
the production path is deferred until the comparison numbers are in.
"""
from __future__ import annotations

import json
import time
from pathlib import Path

import matplotlib

matplotlib.use("Agg")

import numpy as np  # noqa: E402

from fvcom_mesh_tools.algorithms import (  # noqa: E402
    alpha_quality,
    min_interior_angle,
    signed_areas,
)
from fvcom_mesh_tools.diagnostics import (  # noqa: E402
    node_valence,
    overconnected_nodes_flag,
)
from fvcom_mesh_tools.io import Fort14Mesh, write_fort14  # noqa: E402

REPO_ROOT = Path(__file__).resolve().parent.parent
DEM = REPO_ROOT / "data" / "bathymetry" / "tokyo_bay" / "dem_00_01_change.nc"
COASTLINE = (
    REPO_ROOT / "data" / "coastline" / "tokyo_bay"
    / "MLIT_C23" / "C23-06_TOKYOBAY.shp"
)
OUT_DIR = REPO_ROOT / "outputs"
F14_OUT = OUT_DIR / "30_tokyo_bay_triangle.14"
SUMMARY_JSON = OUT_DIR / "30_triangle_engine_summary.json"
SUMMARY_TXT = OUT_DIR / "30_triangle_engine_summary.txt"

# PoC #16 build configuration (Tokyo Bay), reused verbatim.
HMIN_M = 200.0
HMAX_M = 5000.0
COAST_TARGET_SIZE_M = 200.0
COAST_EXPANSION_RATE = 0.005
MIN_POLYGON_AREA_M2 = 1_000_000.0
MIN_ISLAND_AREA_M2 = 100_000.0
ZMAX = 0.0


def build_with_triangle() -> tuple[np.ndarray, np.ndarray, float]:
    """Run ocsmesh's Geom + Hfun pipeline and drive Triangle.

    Returns ``(coords_lonlat, cells, wall_seconds)``. ``coords`` is
    ``(NP, 2)`` in EPSG:4326 even when ocsmesh internally meshes in a
    metric CRS; ``cells`` is ``(NE, 3)`` 0-based.
    """
    from ocsmesh import Geom, Hfun, MeshDriver, Raster
    from pyproj import CRS, Transformer

    from fvcom_mesh_tools.io import (
        filter_multipolygon_by_area,
        load_coastline_as_lines,
    )

    raster = Raster(str(DEM))
    geom = Geom(raster, zmax=ZMAX)

    # Match PoC #16 area filtering.
    mp = geom.get_multipolygon()
    mp_filt = filter_multipolygon_by_area(
        mp,
        src_crs=raster.crs,
        min_polygon_area_m2=MIN_POLYGON_AREA_M2,
        min_island_area_m2=MIN_ISLAND_AREA_M2,
    )
    print(
        f"[30] geom: polygons {len(mp.geoms)} -> {len(mp_filt.geoms)} "
        f"after area filter"
    )
    geom = Geom(mp_filt, crs=raster.crs)

    hfun = Hfun(raster, hmin=HMIN_M, hmax=HMAX_M)
    coast = load_coastline_as_lines([COASTLINE], bbox=None)
    print(
        f"[30] coastline: {len(coast.geoms)} line strings; "
        f"target={COAST_TARGET_SIZE_M:g} expansion={COAST_EXPANSION_RATE:g}"
    )
    hfun.add_feature(
        feature=coast,
        expansion_rate=COAST_EXPANSION_RATE,
        target_size=COAST_TARGET_SIZE_M,
    )

    print("[30] generate_mesh (engine_name='triangle') ...")
    t0 = time.perf_counter()
    driver = MeshDriver(geom, hfun=hfun, engine_name="triangle")
    mesh = driver.run()
    wall = time.perf_counter() - t0
    print(f"[30] triangle wall: {wall:.2f} s")

    coords = np.asarray(mesh.coord, dtype=np.float64)
    cells = np.asarray(mesh.triangles, dtype=np.int64)
    if cells.size == 0:
        raise RuntimeError("Triangle produced zero triangles.")

    src_crs = mesh.crs
    if src_crs is not None and not CRS(src_crs).equals(CRS.from_epsg(4326)):
        print(f"[30] projecting coords {src_crs} -> EPSG:4326")
        transformer = Transformer.from_crs(
            src_crs, CRS.from_epsg(4326), always_xy=True,
        )
        lon, lat = transformer.transform(coords[:, 0], coords[:, 1])
        coords = np.column_stack([lon, lat])
    return coords, cells, wall


def quality_metrics(mesh: Fort14Mesh) -> dict:
    """Compute the same metrics PoC #16 / #19 reported."""
    nodes = mesh.nodes
    elems = mesh.elements
    NP, NE = mesh.n_nodes, mesh.n_elements
    alpha = alpha_quality(nodes, elems)
    min_ang = min_interior_angle(nodes, elems)
    sa = signed_areas(nodes, elems)
    frac_lt20 = float((min_ang < np.deg2rad(20.0)).sum() / NE)

    valence = node_valence(elems, n_nodes=NP)
    overconn_flag, _ = overconnected_nodes_flag(elems, n_nodes=NP, max_nbr=8)
    return {
        "NP": int(NP),
        "NE": int(NE),
        "alpha_mean": float(alpha.mean()),
        "alpha_p50": float(np.median(alpha)),
        "min_angle_p50_deg": float(np.degrees(np.median(min_ang))),
        "frac_lt_20deg": frac_lt20,
        "max_valence": int(valence.max()),
        "n_overconnected": int(overconn_flag.sum()),
        "n_flipped": int((sa < 0).sum()),
    }


def comparison_table(triangle_metrics: dict) -> str:
    """Side-by-side text comparison against PoC #16 (gmsh) and #19 (oceanmesh).

    Reference numbers are pinned from docs/architecture.md §2.97-110
    so the table can be read without external lookup.
    """
    refs = {
        "ocsmesh (gmsh, PoC #16)": {
            "NP": 19362, "NE": 27609,
            "alpha_mean": 0.847, "min_angle_p50_deg": 40.7,
            "frac_lt_20deg": 0.0113, "max_valence": 26,
            "n_overconnected": 440, "n_flipped": 0,
        },
        "oceanmesh (PoC #19)": {
            "NP": 31771, "NE": 53203,
            "alpha_mean": 0.959, "min_angle_p50_deg": 51.05,
            "frac_lt_20deg": 0.0010, "max_valence": 9,
            "n_overconnected": 3, "n_flipped": 0,
        },
        "ocsmesh (triangle, PoC #30)": triangle_metrics,
    }
    keys = ["NP", "NE", "alpha_mean", "min_angle_p50_deg",
            "frac_lt_20deg", "max_valence", "n_overconnected", "n_flipped"]
    fmt = {
        "NP": "{:>10,d}", "NE": "{:>10,d}",
        "alpha_mean": "{:>10.3f}", "min_angle_p50_deg": "{:>10.2f}",
        "frac_lt_20deg": "{:>10.4f}", "max_valence": "{:>10d}",
        "n_overconnected": "{:>10d}", "n_flipped": "{:>10d}",
    }
    col_w = 30
    out = []
    header = "metric".ljust(20) + "".join(name.rjust(col_w) for name in refs)
    out.append(header)
    out.append("-" * len(header))
    for k in keys:
        row = k.ljust(20)
        for name in refs:
            v = refs[name].get(k)
            row += fmt[k].format(v).rjust(col_w) if v is not None else "n/a".rjust(col_w)
        out.append(row)
    return "\n".join(out)


def main() -> int:
    OUT_DIR.mkdir(parents=True, exist_ok=True)
    if not DEM.exists() or not COASTLINE.exists():
        raise SystemExit(
            f"required inputs missing: {DEM=}, {COASTLINE=}"
        )

    coords, cells, wall = build_with_triangle()

    # Construct a minimal Fort14Mesh (no boundaries set; this PoC
    # measures element-shape quality, not the boundary classification).
    mesh = Fort14Mesh(
        title="poc30_triangle",
        nodes=coords,
        depths=np.zeros(coords.shape[0], dtype=np.float64),
        elements=cells,
        open_boundaries=[],
        land_boundaries=[],
    )
    write_fort14(mesh, F14_OUT)
    print(f"[30] wrote {F14_OUT}")

    metrics = quality_metrics(mesh)
    metrics["wall_seconds"] = wall
    print(f"[30] metrics: {metrics}")

    table = comparison_table({k: v for k, v in metrics.items()
                              if k != "wall_seconds"})
    print()
    print("=== Quality comparison (Triangle vs gmsh vs oceanmesh) ===")
    print(table)

    SUMMARY_TXT.write_text(table + "\n", encoding="utf-8")
    SUMMARY_JSON.write_text(
        json.dumps(
            {
                "input_dem": str(DEM.resolve()),
                "input_coastline": str(COASTLINE.resolve()),
                "build_config": {
                    "hmin_m": HMIN_M, "hmax_m": HMAX_M,
                    "coast_target_size_m": COAST_TARGET_SIZE_M,
                    "coast_expansion_rate": COAST_EXPANSION_RATE,
                    "min_polygon_area_m2": MIN_POLYGON_AREA_M2,
                    "min_island_area_m2": MIN_ISLAND_AREA_M2,
                    "zmax": ZMAX,
                    "engine_name": "triangle",
                },
                "metrics": metrics,
            },
            indent=2,
        ),
        encoding="utf-8",
    )
    print(f"[30] wrote {SUMMARY_TXT}")
    print(f"[30] wrote {SUMMARY_JSON}")
    return 0


if __name__ == "__main__":
    raise SystemExit(main())
