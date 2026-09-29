"""Build a whole-bay base mesh from a base recipe, and record what it used.

    python notebooks/440_base_mesh.py recipes/base/tokyo_bay_tool.yaml [OUTDIR]

On OCTOPUS: ``qsub -v FMESH_RECIPE=... jobs/octopus/440_base_mesh.sh``.

The steps, all from the recipe and the raw data under $DATA_DIR:

1. land: OSM land polygons minus sea-connected inland water, for the
   recipe's window (``prep.fetch_true_land``), written to ``land/``;
2. generation: ``notebooks/325_sample_repro.py`` (oceanmesh DistMesh with the
   open boundary constrained), into ``generate/``;
3. finishing: ``notebooks/331_finish2.py`` (channel policy, open-boundary
   finishing, coastline fit);
4. depths and the FVCOM case: ``notebooks/422_tool_base.py``.

``report.json`` holds the provenance (commits, input hashes, the effective
settings, seeds), the hashes of the products, and whether the fort.14 is the
recipe's reference byte for byte. Nothing is read from an earlier run: a
source that is missing stops the build rather than being replaced.
"""

import importlib.util
import json
import os
import subprocess
import sys
import time
from pathlib import Path

REPO = Path(__file__).resolve().parents[1]
sys.path.insert(0, str(REPO / "src"))

import numpy as np  # noqa: E402

from fvcom_mesh_tools.base_recipe import compare_to_reference, load_base_recipe  # noqa: E402
from fvcom_mesh_tools.provenance import collect, dataset_files, file_sha256  # noqa: E402

T0 = time.time()


def say(msg):
    print(f"[base] {msg} +{time.time() - T0:.0f}s", flush=True)


recipe = load_base_recipe(sys.argv[1] if len(sys.argv) > 1 else
                          os.environ.get("FMESH_RECIPE", "recipes/base/tokyo_bay_tool.yaml"))
OUT = Path(sys.argv[2] if len(sys.argv) > 2 else REPO / "outputs" / f"base_{recipe['name']}")
OUT = OUT.resolve()
if OUT.exists() and any(OUT.iterdir()):
    raise SystemExit(f"{OUT} is not empty; move it first -- a build never mixes with another")
OUT.mkdir(parents=True, exist_ok=True)
say(f"recipe {recipe['recipe_path']} -> {OUT}")

# ------------------------------------------------------------ the raw data
DATA = os.environ.get("DATA_DIR")
if not DATA:
    raise SystemExit("DATA_DIR is not set")
DATA = Path(DATA)
SOURCES = {
    "osm_land": DATA / "geodata/OSM/land-polygons-split-4326/land_polygons.shp",
    "osm_water": DATA / "geodata/OSM/geofabrik_kanto/gis_osm_water_a_free_1.shp",
    "osm_waterways": DATA / "geodata/OSM/geofabrik_kanto/gis_osm_waterways_free_1.shp",
    "srtm15": DATA / "geodata/bathymetry/tokyo_bay/SRTM15_kanto_15s.nc",
    "m7001_fine": DATA / "geodata/bathymetry/M7001/TP/M7001_dem_tokyobay.nc",
    "m7001_wide": DATA / "geodata/bathymetry/tokyo_bay/kanto_M7001_srtm_15s.nc",
}
absent = [f"{k}: {p}" for k, p in SOURCES.items() if not p.exists()]
if absent:
    # fetch_true_land falls back to xcoast's own download and skips the water
    # subtraction when a source is missing; a base mesh must not
    raise SystemExit("missing source data:\n  " + "\n  ".join(absent))

# ------------------------------------------------------------------ 1. land
from fvcom_mesh_tools.prep import fetch_true_land  # noqa: E402

land_dir = OUT / "land"
land_dir.mkdir()
say(f"land: OSM for {recipe['land']['bbox']}")
land = fetch_true_land(tuple(recipe["land"]["bbox"]),
                       land_shp_path=SOURCES["osm_land"],
                       min_water_area_deg2=float(recipe["land"]["min_water_area_deg2"]),
                       cache_dir=land_dir / "xcoast_cache", force=True)
LAND_SHP = land_dir / "land_osm.shp"
land.to_file(LAND_SHP)
say(f"land: {len(land)} polygon(s)")

# ------------------------------------------------- 2-3. generate and finish
gen = OUT / "generate"
env = {k: v for k, v in os.environ.items()
       if not k.startswith("SR_")}                  # nothing leaks in from the shell
env.update(recipe["settings"])
env.update(SR_OUT=str(gen), SR_LAND=str(LAND_SHP), SR_OBC_FILE=recipe["open_boundary"],
           SR_DOMAIN_FILE=recipe["domain"], SR_EDITS_DIR=recipe["edits"],
           FMESH_OVERWRITE="1", PYTHONPATH=str(REPO / "src"))
for script in ("325_sample_repro.py", "331_finish2.py"):
    say(f"run {script}")
    subprocess.run([sys.executable, str(REPO / "notebooks" / script)], cwd=REPO, env=env,
                   check=True)

# ------------------------------------------------------ 4. depths, the case
say("run 422_tool_base.py")
subprocess.run([sys.executable, str(REPO / "notebooks" / "422_tool_base.py"),
                str(gen / "sample_repro_final.14"), str(OUT), recipe["case"]],
               cwd=REPO, env=env, check=True)

# ------------------------------------------------------------- the record
report_path = OUT / "report.json"
report = json.loads(report_path.read_text())
case = recipe["case"]
products = {p.name: file_sha256(p) for p in sorted(OUT.glob(f"{case}*"))}
fort14 = products.get(f"{case}.14")
ref = recipe.get("reference")
qa = json.loads((OUT / f"{case}_qa.json").read_text())


def wet_area_km2(grd):
    lines = Path(grd).read_text().splitlines()
    nn, ne = (int(lines[k].split("=")[1]) for k in (0, 1))
    tri = np.array([ln.split()[1:4] for ln in lines[2:2 + ne]], dtype=np.int64) - 1
    xy = np.array([ln.split()[1:3] for ln in lines[2 + ne:2 + ne + nn]], dtype=float)
    a, b = xy[tri[:, 1]] - xy[tri[:, 0]], xy[tri[:, 2]] - xy[tri[:, 0]]
    return float(np.abs(0.5 * (a[:, 0] * b[:, 1] - a[:, 1] * b[:, 0])).sum() / 1e6)


summary = {"fort14_sha256": fort14, "n_nodes": qa["mesh"]["n_nodes"],
           "n_elements": qa["mesh"]["n_elements"], "n_obc_nodes": qa["mesh"]["n_obc_nodes"],
           "wet_area_km2": wet_area_km2(OUT / f"{case}_grd.dat")}
verdict = None if ref is None else compare_to_reference(summary, ref, qa["passed"])


def spec_path(name):
    spec = importlib.util.find_spec(name)       # located, not imported
    return spec.origin if spec is not None else None


edits = sorted(Path(recipe["edits"]).glob("*.json"))
report.update(
    recipe=recipe["recipe_path"],
    products_sha256=products,
    reference=ref,
    reproduction=verdict,
    settings=recipe["settings"],
    seeds={"generate": int(recipe["settings"]["SR_GEN_SEED"]),
           "finish": int(recipe["settings"]["SR_FIN_SEED"])},
    threads={k: os.environ.get(k) for k in ("OMP_NUM_THREADS", "NUMBA_NUM_THREADS",
                                            "OPENBLAS_NUM_THREADS", "MKL_NUM_THREADS")},
    land_polygons=int(len(land)),
    provenance=collect(
        code={"fvcom_mesh_tools": str(REPO / "src" / "fvcom_mesh_tools" / "__init__.py"),
              "driver": __file__,
              "oceanmesh": spec_path("oceanmesh"),
              "xcoast": spec_path("xcoast")},
        files={"recipe": recipe["recipe_path"],
               "open_boundary": recipe["open_boundary"],
               "domain": recipe["domain"],
               "edits": [str(p) for p in edits],
               **{k: dataset_files(p) for k, p in SOURCES.items()}},
    ),
)
report_path.write_text(json.dumps(report, indent=1, default=str))
if verdict is None:
    say(f"fort.14 sha256 {fort14}")
else:
    rel = ", ".join(f"{k} {v:+.2%}" for k, v in verdict["relative_difference"].items())
    say(("REPRODUCES" if verdict["reproduces"] else "DOES NOT REPRODUCE")
        + f" the recipe's reference ({rel}; QA "
        + ("passed" if verdict["qa_passed"] else "FAILED") + "; "
        + ("byte-identical" if verdict["byte_identical"] else "not byte-identical") + ")")
say("done")
if verdict is not None and not verdict["reproduces"]:
    sys.exit(3)
