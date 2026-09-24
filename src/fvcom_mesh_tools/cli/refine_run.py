"""``fmesh-refine``: refine a region of an existing FVCOM mesh from a recipe.

The generator itself is ``notebooks/420_local_refine.py``. It fills the hole
with DistMesh from the GPL-3.0 ``oceanmesh`` fork, and the licence policy of
this Apache-2.0 package forbids importing that, so this command runs it as a
**subprocess** -- which is also how a batch job runs it
(``jobs/octopus/417_hires_refine.sh``).

    fmesh-refine recipes/refine/kimitsu_port_hires.yaml
    fmesh-refine my_port.yaml --out outputs/my_port --seeds 0,1,2 --land land.shp

It is heavy (minutes to tens of minutes, several GB): on OCTOPUS run it inside
a batch job, never on a login node. It refuses to write over an existing
output, because an old ``report.json`` makes a failed rerun look like success.
"""

from __future__ import annotations

import argparse
import os
import subprocess
import sys
from pathlib import Path

DEFAULT_LAND = ("{DATA_DIR}/geodata/OSM/coastmask_cache/"
                "custom_139.55_34.9_140.3_35.75_minarea1e-05/land.shp")


def repo_root() -> Path:
    """The repository this package is installed from (an editable install)."""
    root = Path(__file__).resolve().parents[3]
    if not (root / "notebooks" / "420_local_refine.py").exists():
        raise FileNotFoundError(
            "notebooks/420_local_refine.py not found next to the package; "
            "fmesh-refine needs an editable install of the repository "
            "(pip install -e .)")
    return root


def build_parser() -> argparse.ArgumentParser:
    p = argparse.ArgumentParser(
        prog="fmesh-refine",
        description="Refine a region of an existing FVCOM mesh from a recipe "
                    "(runs notebooks/420_local_refine.py as a subprocess).")
    p.add_argument("recipe", type=Path, help="refinement recipe (YAML)")
    p.add_argument("--out", type=Path, default=None,
                   help="output directory (default outputs/refine_<recipe name>)")
    p.add_argument("--seeds", default="0,1,2,3,4",
                   help="DistMesh seeds to try, best kept (default 0,1,2,3,4)")
    p.add_argument("--land", type=Path, default=None,
                   help="OSM land polygons (default: $FMESH_LAND, else the "
                        "Tokyo Bay coastmask under $DATA_DIR)")
    p.add_argument("--dry-run", action="store_true",
                   help="print the command and environment, run nothing")
    return p


def main(argv: list[str] | None = None) -> int:
    args = build_parser().parse_args(argv)
    try:
        root = repo_root()
    except FileNotFoundError as exc:
        print(f"fmesh-refine: {exc}", file=sys.stderr)
        return 2
    recipe = args.recipe.expanduser().resolve()
    if not recipe.exists():
        print(f"fmesh-refine: no such recipe: {recipe}", file=sys.stderr)
        return 2
    out = (args.out or root / "outputs" / f"refine_{recipe.stem}").expanduser().resolve()
    # ANY content counts: a run that failed before its first report leaves
    # the rim, the constraints and the filtered shoreline behind, and a rerun
    # into that directory would mix two runs (review F10).
    if out.exists() and any(out.iterdir()):
        print(f"fmesh-refine: {out} is not empty; move it first", file=sys.stderr)
        return 2
    land = args.land or os.environ.get("FMESH_LAND") or DEFAULT_LAND.format(
        DATA_DIR=os.environ.get("DATA_DIR", "/octfs/work/G16445/share/Data"))
    land = Path(land).expanduser().resolve()
    if not land.exists():
        print(f"fmesh-refine: no land polygons at {land} (--land)", file=sys.stderr)
        return 2
    seeds = ",".join(s.strip() for s in args.seeds.replace(":", ",").split(",") if s.strip())
    env = {**os.environ, "LR_OUT": str(out), "LR_SEEDS": seeds, "FMESH_LAND": str(land),
           "MPLBACKEND": "Agg"}
    cmd = [sys.executable, str(root / "notebooks" / "420_local_refine.py"), str(recipe)]
    print(f"fmesh-refine: {' '.join(cmd)}\n  LR_OUT={out}\n  LR_SEEDS={seeds}\n"
          f"  FMESH_LAND={land}", flush=True)
    if args.dry_run:
        return 0
    return subprocess.call(cmd, env=env, cwd=root)


if __name__ == "__main__":  # pragma: no cover
    raise SystemExit(main())
