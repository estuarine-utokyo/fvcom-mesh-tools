#!/bin/bash
# Submit the M2 validation as a dependency chain of small requests:
#
#   411 prep (8 cores)  ->  412 run x3 (64 cores each, concurrent)  ->  413 analysis
#
# Each stage asks only for what it needs, so the three integrations take three
# shared 64-core slots instead of one 192-core slot on a single node (job
# 115307 sat queued in O1).  NQSV chains them with `qsub --after`.
#
# Run from the repository root (a login-node shell, not a batch job):
#   bash jobs/octopus/410_m2_chain.sh [path/to/candidate.14]
set -euo pipefail
cd "$(dirname "$0")/../.."
REPO=$(pwd)
B_MESH=${1:-$REPO/outputs/verify_409.115302/fit/sample_repro_final.14}
[ -f "$B_MESH" ] || { echo "mesh not found: $B_MESH"; exit 2; }
STAMP=$(date +%Y%m%d_%H%M%S)
RUN_ROOT=/octfs/work/G16445/v61021/scratch/m2_$STAMP
mkdir -p "$RUN_ROOT"
echo "run root: $RUN_ROOT"
echo "B mesh  : $B_MESH"

id() { sed 's/^\([0-9:]*\.[a-z]*\).*/\1/;s/^0://' <<<"$1" | tr -d ' '; }

prep=$(qsub -v "FMESH_RUN_ROOT=$RUN_ROOT,FMESH_B_MESH=$B_MESH" jobs/octopus/411_m2_prep.sh)
prep=$(grep -oE '[0-9]+\.[a-z]+' <<<"$prep" | head -1)
echo "prep     = $prep"

runs=()
for case in A B_own B_m7001; do
    r=$(qsub --after "$prep" -N "m2_$case" \
        -v "FMESH_RUN_ROOT=$RUN_ROOT,FMESH_CASE=$case" jobs/octopus/412_m2_run.sh)
    r=$(grep -oE '[0-9]+\.[a-z]+' <<<"$r" | head -1)
    runs+=("$r")
    echo "run $case = $r"
done

ana=$(qsub --after "$(IFS=,; echo "${runs[*]}")" \
    -v "FMESH_RUN_ROOT=$RUN_ROOT" jobs/octopus/413_m2_analysis.sh)
ana=$(grep -oE '[0-9]+\.[a-z]+' <<<"$ana" | head -1)
echo "analysis = $ana"
echo "outputs  -> $REPO/outputs/m2_$STAMP"
