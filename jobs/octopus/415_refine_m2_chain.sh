#!/bin/bash
# Submit the M2 test of a local refinement as a chain of small requests:
#
#   414 prep (8 cores) -> 412 run x2 (64 cores each, concurrent) -> 413 analysis
#
# The two integrations differ in ONE thing -- the patch -- so they share the
# open boundary, the forcing, the sponge, the sigma levels and the external
# step. The step is the refined mesh's: comparing each mesh at its own step
# would compare two time steps as well as two meshes.
#
# Run from the repository root on a LOGIN node:
#   bash jobs/octopus/415_refine_m2_chain.sh \
#       outputs/base_tool/TokyoBayTool \
#       outputs/refine_futtsu_nori_tool/fvcom/futtsu_nori_tool [DTE]
set -euo pipefail
cd "$(dirname "$0")/../.."
REPO=$(pwd)
BASE=${1:-$REPO/outputs/base_tool/TokyoBayTool}
REFINED=${2:-$REPO/outputs/refine_futtsu_nori_tool/fvcom/futtsu_nori_tool}
# Empty by default: 414 works the step out from the finer mesh and rounds it
# DOWN. A hard-coded 1.5 s belonged to one experiment and is wrong for the
# next mesh in either direction.
DTE=${3:-}
for p in "${BASE}_grd.dat" "${BASE}_dep.dat" "${BASE}_obc.dat" \
         "${REFINED}_grd.dat" "${REFINED}_dep.dat" "${REFINED}_obc.dat"; do
    [ -f "$p" ] || { echo "not found: $p"; exit 2; }
done
STAMP=$(date +%Y%m%d_%H%M%S)
RUN_ROOT=/octfs/work/G16445/v61021/scratch/m2r_$STAMP
mkdir -p "$RUN_ROOT"
echo "run root: $RUN_ROOT"
echo "base    : $BASE"
echo "refined : $REFINED"
echo "DTE     : $DTE s"

# `--after` starts a dependent when its predecessor TERMINATES, not when it
# succeeds -- this NQSV has no afterok -- so every dependent checks its own
# input and exits 2 rather than burning a 64-core slot on a failed prep.
prep=$(qsub -v "FMESH_RUN_ROOT=$RUN_ROOT,FMESH_BASE=$BASE,FMESH_REFINED=$REFINED,FMESH_DTE=$DTE" \
    jobs/octopus/414_refine_m2_prep.sh)
prep=$(grep -oE '[0-9]+\.[a-z]+' <<<"$prep" | head -1)
echo "prep     = $prep"

runs=()
for case in base refined; do
    r=$(qsub --after "$prep" -N "m2r_$case" \
        -v "FMESH_RUN_ROOT=$RUN_ROOT,FMESH_CASE=$case" jobs/octopus/412_m2_run.sh)
    r=$(grep -oE '[0-9]+\.[a-z]+' <<<"$r" | head -1)
    runs+=("$r")
    echo "run $case = $r"
done

ana=$(qsub --after "$(IFS=,; echo "${runs[*]}")" \
    -v "FMESH_RUN_ROOT=$RUN_ROOT,FMESH_OUT=$REPO/outputs/m2r_$STAMP" \
    jobs/octopus/413_m2_analysis.sh)
ana=$(grep -oE '[0-9]+\.[a-z]+' <<<"$ana" | head -1)
echo "analysis = $ana"
echo "outputs  -> $REPO/outputs/m2r_$STAMP"
echo "jobs     : $prep ${runs[*]} $ana"
