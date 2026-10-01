#!/bin/bash
#PBS -q OCT-S
#PBS --group=G16445
#PBS -l cpunum_job=8
#PBS -l memsz_job=32GB
#PBS -l elapstim_req=00:30:00
#PBS -N fmesh_414prep
#PBS -j o
#PBS -o logs/414_refine_m2_prep.pbs.log
#PBS -r n
# Stage the two M2 cases of a local-refinement test: the base, and the base
# refined. Eight cores is all staging needs; the integrations ask for 64 each
# and start in their own shared slots.
#
# Required: FMESH_RUN_ROOT, FMESH_BASE, FMESH_REFINED (case prefixes, without
# _grd.dat).  Optional: FMESH_DTE.
set -euo pipefail
cd "${PBS_O_WORKDIR:?Submit from the repository root}"
. jobs/octopus/common.sh 414_refine_m2_prep 8
RUN_ROOT=${FMESH_RUN_ROOT:?set FMESH_RUN_ROOT}
BASE=${FMESH_BASE:?set FMESH_BASE}
REFINED=${FMESH_REFINED:?set FMESH_REFINED}
fmesh_stage_lock "$RUN_ROOT" || exit 2
# the cases are about to be replaced: every acceptance of the old ones goes
# first, before anything is written (review round 14 F3)
rm -f "$RUN_ROOT/STAGED" "$RUN_ROOT/SMOKE_OK" "$RUN_ROOT"/*/RUN_OK
for case in base refined; do
    if compgen -G "$RUN_ROOT/$case/output/m2_*.nc" >/dev/null; then
        echo "Existing model output: $RUN_ROOT/$case/output; archive before rerunning"
        exit 2
    fi
done
echo "base    = $BASE"
echo "refined = $REFINED"
python notebooks/414_refine_m2_prep.py --root "$RUN_ROOT" \
    --base "$BASE" --refined "$REFINED" ${FMESH_DTE:+--dte "$FMESH_DTE"}
echo "end=$(date -Is)"
