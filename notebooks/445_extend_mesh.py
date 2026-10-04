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

import atexit
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
from fvcom_mesh_tools.extend_recipe import EXPECT_ENV, load_extend_recipe  # noqa: E402
from fvcom_mesh_tools.outdir import TOKEN_ENV, new_token, reserve  # noqa: E402
from fvcom_mesh_tools.provenance import (  # noqa: E402
    changed_files,
    changed_inventory,
    code_identity,
    collect,
    dataset_files,
    file_sha256,
)

T0 = time.time()


def say(msg):
    print(f"[extend] {msg} +{time.time() - T0:.0f}s", flush=True)


recipe = load_extend_recipe(sys.argv[1] if len(sys.argv) > 1 else os.environ["FMESH_RECIPE"])
OUT = Path(sys.argv[2] if len(sys.argv) > 2 else REPO / "outputs" / f"extend_{recipe['name']}")
OUT = OUT.resolve()
if len(sys.argv) <= 2 and (REPO / "outputs").resolve() not in OUT.parents:
    raise SystemExit(f"the default output {OUT} is not under {REPO / 'outputs'}")
# A report is written whatever happens, from the moment the output is
# reserved (review round 3 F6, round 4 F7): an exit handler writes a failure
# report unless the final one was written. A failure before the provenance
# is taken says so ("provenance": null). The handler is registered before
# anything else that can fail, a log line included (round 5 F8).
STATE = {"done": False, "stage": "inputs", "provenance": None, "extra": {}}
# reserve the output atomically: checking that it is empty and then creating
# it let two builds into the same directory (review F20)
TOKEN = new_token()           # proves to the stages which run they belong to
OUT = reserve(OUT, token=TOKEN)


def _on_exit():
    if not STATE["done"]:
        (OUT / "report.json").write_text(json.dumps(
            {"recipe": recipe["recipe_path"], "status": "failed", "stage": STATE["stage"],
             **STATE["extra"], "provenance": STATE["provenance"],
             "provenance_complete": STATE["provenance"] is not None}, indent=1, default=str))


atexit.register(_on_exit)
say(f"recipe {recipe['recipe_path']} -> {OUT}")
if not os.environ.get("DATA_DIR"):
    raise SystemExit("DATA_DIR is not set")
DATA = Path(os.environ["DATA_DIR"])
osm_land = DATA / "geodata/OSM/land-polygons-split-4326/land_polygons.shp"
names = sorted(set(recipe["bathymetry"]["sizing"]) | set(recipe["bathymetry"]["depths"]))
bathy = source_files(names)
missing = [str(p) for p in [osm_land, *[f for v in bathy.values() for f in v]] if not p.exists()]
if missing:
    STATE["extra"]["missing"] = missing
    raise SystemExit("missing source data:\n  " + "\n  ".join(missing))


def spec_path(name):
    spec = importlib.util.find_spec(name)
    return spec.origin if spec is not None else None


# what makes this build, captured BEFORE it runs. The stages read the
# recipe, the open boundary and the base themselves; these are hashed again
# at the end and a change fails the build (round 4 F8).
STATE["stage"] = "provenance"
code = {"fvcom_mesh_tools": str(REPO / "src" / "fvcom_mesh_tools" / "__init__.py"),
        "driver": __file__}
if spec_path("oceanmesh"):
    code["oceanmesh"] = spec_path("oceanmesh")
b = Path(recipe["base"]) / recipe["base_case"]
INPUTS = {"recipe": recipe["recipe_path"], "open_boundary": recipe["open_boundary"],
          "base": [f"{b}_{k}.dat" for k in ("grd", "dep", "obc")]}
if Path(recipe["open_boundary"]).with_suffix(".json").exists():     # 444's report
    INPUTS["open_boundary_report"] = str(Path(recipe["open_boundary"]).with_suffix(".json"))
STATE["provenance"] = PROV = collect(
    code=code,
    files={**INPUTS, "osm_land": dataset_files(osm_land),
           **{f"bathymetry_{k}": [str(p) for p in v] for k, v in bathy.items()}})
# the recipe and boundary hashed are the ones parsed above, and the stages
# must read the same (review round 5 F3)
for key, name in (("recipe_sha256", "recipe"), ("open_boundary_sha256", "open_boundary")):
    if PROV["files"][name]["sha256"] != recipe[key]:
        raise SystemExit(f"{name} changed between reading and recording it")

gen = OUT / "generate"
env = dict(os.environ, PYTHONPATH=str(REPO / "src"), **{TOKEN_ENV: TOKEN},
           **{var: recipe[key] for key, var in EXPECT_ENV.items()})
# a failed stage still leaves a report with the provenance and the stage that
# failed (review round 2 F18): the stages are run first, the report written
# either way, and the exit code is the first failure's
failed = None
for script, args in (("446_extend_generate.py", [recipe["recipe_path"], str(gen)]),
                     ("447_extend_merge.py", [recipe["recipe_path"], str(gen), str(OUT)])):
    say(f"run {script}")
    STATE["stage"] = script
    rc = subprocess.run([sys.executable, str(REPO / "notebooks" / script), *args], cwd=REPO,
                        env=env).returncode
    if rc != 0:
        failed = {"stage": script, "returncode": rc}
        break
STATE["stage"] = "report"
# both stages must have left their reports before the build can be "ok"
# (review round 27 F7)
if not failed:
    for need in (gen / "generate.json", OUT / "merge.json"):
        if not need.exists():
            failed = {"stage": need.name, "missing": str(need)}
            break


def _count(value, what):
    """A whole number from a report: not a float, a bool or a string (round 29 F2)."""
    if isinstance(value, bool) or not isinstance(value, int):
        raise ValueError(f"{what} is not an integer: {value!r}")
    return value


def _acceptance_problem():
    """Why a build whose stages exited 0 is still not "ok", or None (review rounds 28 F1, 29 F2)."""
    try:
        merged = json.loads((OUT / "merge.json").read_text())
        generated = json.loads((gen / "generate.json").read_text())
        for name, rep in (("merge.json", merged), ("generate.json", generated)):
            if not isinstance(rep, dict):
                raise ValueError(f"{name} is not an object")
        if not isinstance(generated.get("settings"), dict) or not generated["settings"]:
            raise ValueError("generate.json lacks its settings")
        # the generation must name the very recipe and boundary this build
        # parsed, and carry every identity the merge accepted against
        # (review round 30 F1)
        ident = generated.get("inputs")
        if not isinstance(ident, dict) or not all(
                ident.get(k) for k in ("recipe_sha256", "open_boundary_sha256", "base_sha256",
                                       "land_sha256", "outer_utm14_sha256")):
            raise ValueError("generate.json lacks its input identities")
        for k in ("recipe_sha256", "open_boundary_sha256"):
            if ident[k] != recipe[k]:
                raise ValueError(f"generate.json's {k} is not this build's")
        for k in ("n_nodes", "n_elements"):
            if _count(generated.get(k), f"generate {k}") <= 0 or _count(
                    merged.get(k), f"merge {k}") <= 0:
                raise ValueError(f"{k} is not positive")
        qa = merged["qa"]
        total, nfail = _count(qa["n_gate_total"], "n_gate_total"), _count(
            qa["n_gate_failed"], "n_gate_failed")
        problems = merged["problems"]
        if not isinstance(problems, list):
            raise ValueError("problems is not a list")
    except (OSError, ValueError, KeyError, TypeError) as exc:
        return f"unreadable stage report ({type(exc).__name__}: {exc})"
    if total <= 0 or nfail != 0:
        return f"QA {total - nfail}/{total}"
    if problems:
        return f"merge problems: {problems}"
    # the products parse, and agree with what the merge reported
    from fvcom_mesh_tools.io.fvcom_native import read_dep, read_grd, read_obc

    try:
        case = recipe["case"]
        nodes, els = read_grd(OUT / f"{case}_grd.dat")
        read_dep(OUT / f"{case}_dep.dat")
        read_obc(OUT / f"{case}_obc.dat")
    except (OSError, ValueError, IndexError) as exc:
        return f"unreadable product ({type(exc).__name__}: {exc})"
    if (len(nodes), len(els)) != (merged["n_nodes"], merged["n_elements"]):
        return (f"products hold NP={len(nodes)} NE={len(els)}, the merge reported "
                f"NP={merged['n_nodes']} NE={merged['n_elements']}")
    return None


if not failed and (why := _acceptance_problem()):
    failed = {"stage": "acceptance", "problem": why}
# every input the provenance names, the bathymetry and land data included
# (review round 6 F8), and the datasets listed again: a file that appeared
# since is a change too (round 7 F9)
changed = sorted(set(changed_files(PROV, list(PROV["files"]))) | set(changed_inventory(
    PROV, {"osm_land": dataset_files(osm_land),
           **{f"bathymetry_{k}": v for k, v in source_files(names).items()}})))
# the stages ran the live code: it must still be the code recorded (the
# package, this driver, oceanmesh; round 7 F10)
_now = collect(code=code, libraries=())["code"]
if any(code_identity(_now[k]) != code_identity(PROV["code"][k]) for k in PROV["code"]):
    changed.append("code")
if changed and not failed:
    failed = {"stage": "inputs", "changed_during_build": changed}

case = recipe["case"]
report = {
    "recipe": recipe["recipe_path"],
    "products_sha256": {p.name: file_sha256(p) for p in sorted(OUT.glob(f"{case}*"))},
    "settings": recipe["settings"], "depths": recipe["depths"],
    "bathymetry": recipe["bathymetry"],
    "status": "failed" if failed else "ok", "failure": failed,
    "inputs_changed_during_build": changed,
    "generate": (json.loads((gen / "generate.json").read_text())
                 if (gen / "generate.json").exists() else None),
    "merge": (json.loads((OUT / "merge.json").read_text())
              if (OUT / "merge.json").exists() else None),
    "threads": {k: os.environ.get(k) for k in ("OMP_NUM_THREADS", "NUMBA_NUM_THREADS",
                                               "OPENBLAS_NUM_THREADS", "MKL_NUM_THREADS")},
    "provenance": PROV,
}
(OUT / "report.json").write_text(json.dumps(report, indent=1, default=str))
STATE["done"] = True
if failed:
    why = failed.get("returncode", failed.get("missing", failed.get("problem", "inputs changed")))
    raise SystemExit(f"{failed['stage']} failed ({why}); "
                     f"report in {OUT / 'report.json'}")
qa = report["merge"]["qa"]
say(f"QA {qa['n_gate_total'] - qa['n_gate_failed']}/{qa['n_gate_total']}; "
    f"NP={report['merge']['n_nodes']:,} NE={report['merge']['n_elements']:,}")
say("done")
