#!/bin/bash
#PJM -L rscgrp=a-batch
#PJM -L vnode=1
#PJM -L vnode-core=32
#PJM -L elapse=03:00:00
#PJM -j
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
set -uo pipefail
cd "${PJM_O_WORKDIR:?Submit from the repository root}"
NAME=${FMESH_CHECK:?pjsub -x FMESH_CHECK=<name>}
. jobs/genkai/common.sh "$NAME" 32
R=$WORK_DIR/scratch/$NAME
mkdir -p "$WORK_DIR/scratch"
mkdir -- "$R" || { echo "$R exists: pick a new FMESH_CHECK"; exit 2; }
echo "commit=$(git rev-parse HEAD) dirty=$(git status --porcelain --untracked-files=all | wc -l)"
rc=0
echo "== pytest $(date -Is)"
python -m pytest -q -p no:cacheprovider tests > "$R/pytest.log" 2>&1 || rc=1
tail -3 "$R/pytest.log"
echo "== 444 $(date -Is)"
python notebooks/444_design_obc.py recipes/extend/tokyo_bay_enshu_obc_design.yaml \
    "$R/obc/obc.csv" > "$R/444.log" 2>&1 || rc=1
tail -3 "$R/444.log"
echo "== 445 $(date -Is)"
python notebooks/445_extend_mesh.py recipes/extend/tokyo_bay_enshu.yaml "$R/build" \
    > "$R/445.log" 2>&1 || rc=1
grep -E "ladders|QA|WARNING|NP=|Error|failed" "$R/445.log" | tail -15
echo "== 453 $(date -Is)"
python notebooks/453_redepth_extended.py recipes/extend/tokyo_bay_enshu.yaml "$R/build" \
    "$R/redepth" --sources cao_shutochokka_2025,m7001,srtm15plus > "$R/453.log" 2>&1 || rc=1
tail -4 "$R/453.log" | cut -c1-400
echo "== fingerprints"
find "$R/build" -name '*_grd.dat' -exec sha256sum {} + | sed "s|$R/||"
echo "== end rc=$rc $(date -Is)"
exit $rc
