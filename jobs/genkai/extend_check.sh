#!/bin/bash
#PJM -L rscgrp=a-batch
#PJM -L vnode=1
#PJM -L vnode-core=32
#PJM -L elapse=03:00:00
#PJM -j
#PJM -X
#PJM -N fmesh_extcheck
#============================================================================
# GENKAI: verify the extension tools end to end on the real Tokyo Bay /
# Enshu-nada case -- the unit tests, the open-boundary design (444), the
# build (445) and the re-depth (453) -- into a fresh directory, and print the
# fingerprints to compare with an earlier run (grd sha256, NP, NE, QA).
#
#   pjsub -x FMESH_CHECK=extcheck28 jobs/genkai/extend_check.sh
# Output: $WORK_DIR/scratch/$FMESH_CHECK (must not exist yet)
# Log:    logs/$FMESH_CHECK.<jobid>.log
#============================================================================
set -euo pipefail    # the stage failures that are collected carry an explicit `|| rc=1`
cd "${PJM_O_WORKDIR:?Submit from the repository root}"
NAME=${FMESH_CHECK:?pjsub -x FMESH_CHECK=<name>}
# one plain name: the run directory and the log must stay where they are meant to
[[ $NAME =~ ^[A-Za-z0-9][A-Za-z0-9._-]*$ ]] || { echo "FMESH_CHECK '$NAME' is not a plain name"; exit 2; }
. jobs/genkai/common.sh "$NAME" 32
R=$WORK_DIR/scratch/$NAME
mkdir -p "$WORK_DIR/scratch"
mkdir -- "$R" || { echo "$R exists: pick a new FMESH_CHECK"; exit 2; }
echo "commit=$(git rev-parse HEAD) dirty=$(git status --porcelain --untracked-files=all | wc -l)"
rc=0
echo "== pytest $(date -Is)"
python -m pytest -q -p no:cacheprovider tests > "$R/pytest.log" 2>&1 || rc=1
tail -3 "$R/pytest.log" || true
echo "== 444 $(date -Is)"
# the later stages build on the boundary it designs, so they cannot go on without it
python notebooks/444_design_obc.py recipes/extend/tokyo_bay_enshu_obc_design.yaml \
    "$R/obc/obc.csv" > "$R/444.log" 2>&1 || { tail -5 "$R/444.log"; echo "444 failed"; exit 1; }
tail -3 "$R/444.log"
# a run-local recipe: the repository recipe with the base and the boundary
# just designed named absolutely, so that 445 and 453 consume that boundary
# and not the checked-in CSV (review round 28 F6)
python - "$R" recipes/extend/tokyo_bay_enshu.yaml outputs/base_tokyo_bay_tool <<'PY'
import sys
from pathlib import Path

import yaml

from fvcom_mesh_tools.yaml_strict import load_unique

run, recipe, base = (Path(a).resolve() for a in sys.argv[1:4])
doc = load_unique(recipe.read_text())      # the strict loader: a repeated key is refused
# the paths are assigned as values and emitted by YAML, never spliced into
# text: a path may hold '#' or ': ' (review round 29 F6)
doc["base"] = str(base)
doc["open_boundary"] = str(run / "obc" / "obc.csv")
(run / "recipe.yaml").write_text(yaml.safe_dump(doc, sort_keys=False, allow_unicode=True))
PY
RECIPE=$R/recipe.yaml
echo "== 445 $(date -Is)"
python notebooks/445_extend_mesh.py "$RECIPE" "$R/build" \
    > "$R/445.log" 2>&1 || rc=1
grep -E "ladders|QA|WARNING|NP=|Error|failed" "$R/445.log" | tail -15 || true
echo "== boundary consumed"
python - "$R" <<'PY' || rc=1
import hashlib, json, sys
from pathlib import Path
run = Path(sys.argv[1])
want = hashlib.sha256((run / "obc" / "obc.csv").read_bytes()).hexdigest()
got = json.loads((run / "build" / "report.json").read_text())["provenance"]["files"]["open_boundary"]["sha256"]
print("obc.csv sha256", want[:16], "build used", got[:16])
sys.exit(0 if got == want else "the build did not use the boundary stage 444 designed")
PY
echo "== 453 $(date -Is)"
python notebooks/453_redepth_extended.py "$RECIPE" "$R/build" \
    "$R/redepth" --sources cao_shutochokka_2025,m7001,srtm15plus > "$R/453.log" 2>&1 || rc=1
tail -4 "$R/453.log" | cut -c1-400 || true
echo "== fingerprints"
find "$R/build" -name '*_grd.dat' -exec sha256sum {} + | sed "s|$R/||"
echo "== end rc=$rc $(date -Is)"
exit $rc
