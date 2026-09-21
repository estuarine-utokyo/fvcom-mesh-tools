#!/bin/bash
#PBS -q OCT-S
#PBS --group=G16445
#PBS -l cpunum_job=8
#PBS -l memsz_job=32GB
#PBS -l elapstim_req=00:30:00
#PBS -N fmesh_411prep
#PBS -j o
#PBS -o logs/411_m2_prep.pbs.log
#PBS -r n
# Stage the three M2 cases.  Split out of the single 192-core job so that each
# stage asks the scheduler only for what it needs: three 64-core run jobs find
# three shared slots far sooner than one 192-core job finds room on one node
# (job 115307 sat in O1 queued), and nothing holds 192 cores during staging or
# analysis.  Chained with qsub --after; see jobs/octopus/410_m2_chain.sh.
#
# Required: FMESH_RUN_ROOT.  Optional: FMESH_B_MESH.
set -euo pipefail
cd "${PBS_O_WORKDIR:?Submit from the repository root}"
. jobs/octopus/common.sh 411_m2_prep 8
REPO=$(pwd)
RUN_ROOT=${FMESH_RUN_ROOT:?set FMESH_RUN_ROOT}
B_MESH=${FMESH_B_MESH:-$REPO/outputs/verify_409.115302/fit/sample_repro_final.14}
[ -f "$B_MESH" ] || { echo "mesh not found: $B_MESH"; exit 2; }
for case in A B_own B_m7001; do
    if compgen -G "$RUN_ROOT/$case/output/m2_*.nc" >/dev/null; then
        echo "Existing model output: $RUN_ROOT/$case/output; archive before rerunning"
        exit 2
    fi
done
echo "B mesh = $B_MESH"
head -2 "$B_MESH"
python notebooks/383_m2_case_prep.py --root "$RUN_ROOT" --mesh "$B_MESH"
echo "end=$(date -Is)"
