"""Build a wide mesh: a base mesh as it is, plus the sea out to a new boundary.

    python notebooks/445_extend_mesh.py recipes/extend/tokyo_bay_enshu.yaml [OUTDIR]

On OCTOPUS: ``qsub -v FMESH_RECIPE=... jobs/octopus/445_extend_mesh.sh``.

1. generation of the outer sea: ``notebooks/446_extend_generate.py``
   (oceanmesh DistMesh; the base's open boundary and the new one fixed);
2. finishing, merge onto the base, depths of the new nodes, the FVCOM case
   and QA: ``notebooks/447_extend_merge.py``.

``report.json`` holds the provenance (commits, input hashes, the effective
settings), the hashes of the products and the stage reports.
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

from fvcom_mesh_tools.dem.sources import source_files  # noqa: E402
from fvcom_mesh_tools.extend_recipe import load_extend_recipe  # noqa: E402
from fvcom_mesh_tools.outdir import reserve  # noqa: E402
from fvcom_mesh_tools.provenance import collect, dataset_files, file_sha256  # noqa: E402

T0 = time.time()


def say(msg):
    print(f"[extend] {msg} +{time.time() - T0:.0f}s", flush=True)


recipe = load_extend_recipe(sys.argv[1] if len(sys.argv) > 1 else os.environ["FMESH_RECIPE"])
OUT = Path(sys.argv[2] if len(sys.argv) > 2 else REPO / "outputs" / f"extend_{recipe['name']}")
OUT = OUT.resolve()
# reserve the output atomically: checking that it is empty and then creating
# it let two builds into the same directory (review F20)
OUT = reserve(OUT)
say(f"recipe {recipe['recipe_path']} -> {OUT}")
DATA = Path(os.environ.get("DATA_DIR") or sys.exit("DATA_DIR is not set"))
osm_land = DATA / "geodata/OSM/land-polygons-split-4326/land_polygons.shp"
names = sorted(set(recipe["bathymetry"]["sizing"]) | set(recipe["bathymetry"]["depths"]))
bathy = source_files(names)
missing = [str(p) for p in [osm_land, *[f for v in bathy.values() for f in v]] if not p.exists()]
if missing:
    raise SystemExit("missing source data:\n  " + "\n  ".join(missing))

gen = OUT / "generate"
env = dict(os.environ, PYTHONPATH=str(REPO / "src"))
# a failed stage still leaves a report with the provenance and the stage that
# failed (review round 2 F18): the stages are run first, the report written
# either way, and the exit code is the first failure's
failed = None
for script, args in (("446_extend_generate.py", [recipe["recipe_path"], str(gen)]),
                     ("447_extend_merge.py", [recipe["recipe_path"], str(gen), str(OUT)])):
    say(f"run {script}")
    rc = subprocess.run([sys.executable, str(REPO / "notebooks" / script), *args], cwd=REPO,
                        env=env).returncode
    if rc != 0:
        failed = {"stage": script, "returncode": rc}
        break


def spec_path(name):
    spec = importlib.util.find_spec(name)
    return spec.origin if spec is not None else None


case = recipe["case"]
b = Path(recipe["base"]) / recipe["base_case"]
report = {
    "recipe": recipe["recipe_path"],
    "products_sha256": {p.name: file_sha256(p) for p in sorted(OUT.glob(f"{case}*"))},
    "base_sha256": {k: file_sha256(Path(f"{b}_{k}.dat")) for k in ("grd", "dep", "obc")},
    "settings": recipe["settings"], "depths": recipe["depths"],
    "bathymetry": recipe["bathymetry"],
    "status": "failed" if failed else "ok", "failure": failed,
    "generate": (json.loads((gen / "generate.json").read_text())
                 if (gen / "generate.json").exists() else None),
    "merge": (json.loads((OUT / "merge.json").read_text())
              if (OUT / "merge.json").exists() else None),
    "threads": {k: os.environ.get(k) for k in ("OMP_NUM_THREADS", "NUMBA_NUM_THREADS",
                                               "OPENBLAS_NUM_THREADS", "MKL_NUM_THREADS")},
    "provenance": collect(
        code={"fvcom_mesh_tools": str(REPO / "src" / "fvcom_mesh_tools" / "__init__.py"),
              "driver": __file__, "oceanmesh": spec_path("oceanmesh")},
        files={"recipe": recipe["recipe_path"], "open_boundary": recipe["open_boundary"],
               "osm_land": dataset_files(osm_land),
               **{f"bathymetry_{k}": [str(p) for p in v] for k, v in bathy.items()}},
    ),
}
(OUT / "report.json").write_text(json.dumps(report, indent=1, default=str))
if failed:
    raise SystemExit(f"{failed['stage']} failed (exit {failed['returncode']}); "
                     f"report in {OUT / 'report.json'}")
qa = report["merge"]["qa"]
say(f"QA {qa['n_gate_total'] - qa['n_gate_failed']}/{qa['n_gate_total']}; "
    f"NP={report['merge']['n_nodes']:,} NE={report['merge']['n_elements']:,}")
say("done")
