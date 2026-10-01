#!/bin/bash
#PBS -q OCT-S
#PBS --group=G16445
#PBS -l cpunum_job=8
#PBS -l memsz_job=64GB
#PBS -l elapstim_req=01:00:00
#PBS -N fmesh_413ana
#PBS -j o
#PBS -o logs/413_m2_analysis.pbs.log
#PBS -r n
# Analyse the three finished integrations.
# Required: FMESH_RUN_ROOT.  Optional: FMESH_OUT.
set -euo pipefail
cd "${PBS_O_WORKDIR:?Submit from the repository root}"
. jobs/octopus/common.sh 413_m2_analysis 8
REPO=$(pwd)
RUN_ROOT=${FMESH_RUN_ROOT:?set FMESH_RUN_ROOT}
# `--after` fires when the predecessor terminates, not when it succeeds, so
# a failed prep or run would otherwise reach this as a puzzling traceback.
[ -f "$RUN_ROOT/manifest.json" ] || { echo "not staged: $RUN_ROOT"; exit 2; }
# the cases are the manifest's: 410 runs A, B_own, B_m7001; 415 base and
# refined (review round 8 F13)
CASES=$(python -c 'import json, sys; print(" ".join(json.load(open(sys.argv[1]))["runs"]))' \
    "$RUN_ROOT/manifest.json") || { echo "manifest without runs: $RUN_ROOT"; exit 2; }
[ -n "$CASES" ] || { echo "manifest lists no runs: $RUN_ROOT"; exit 2; }
for c in $CASES; do
    [ -f "$RUN_ROOT/$c/RUN_OK" ] || { echo "run $c did not pass: no $RUN_ROOT/$c/RUN_OK"; exit 2; }
done
OUT=${FMESH_OUT:-$REPO/outputs/m2_$(basename "$RUN_ROOT")}
python notebooks/384_m2_analysis.py --root "$RUN_ROOT" --output "$OUT"
echo "wrote $OUT"
echo "end=$(date -Is)"
