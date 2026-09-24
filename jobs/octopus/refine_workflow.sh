#!/bin/bash
# The whole local-refinement workflow for one recipe, as a chain of batch jobs.
#
# Run on an OCTOPUS LOGIN node from the repository root; it only submits.
#
#   bash jobs/octopus/refine_workflow.sh recipes/refine/my_port.yaml
#
#   1. refine    417_hires_refine.sh       the mesh, QA, report     (1 core)
#   2. figures   418 + notebooks/434        final_mesh_*.png          (after 1)
#   3. depths    421_finish_and_run.sh      fmesh-finish-depths, and  (after 1)
#                                           an M2 pair staged against the base
#   4. smoke     423_m2_smoke.sh            2-day FVCOM run, both      (after 3)
#   5. M2        412 x 2 + 413              20-day M2 and the gauge    (after 4,
#                                           comparison                  FMESH_M2=1)
#
# Options (environment):
#   FMESH_HMIN (3)  FMESH_HMAX (300)  FMESH_RFACTOR (0.2)   the depth product
#   FMESH_VIEWS     close-ups for the figures, "name:x0:x1:y0:y1+name2:..."
#   FMESH_M2        1 to go on to the 20-day M2 comparison (default 0)
#   LR_SEEDS        DistMesh seeds, colon-separated (default 0:1:2:3:4)
#
# NQSV has no afterok: a stage starts when the one before it ENDS, whatever
# the outcome. So each stage writes a marker only when it passed, and the
# next requires it: ACCEPTED (refine) -> STAGED (depths) -> SMOKE_OK (smoke)
# -> RUN_OK per case (M2). Read the logs in order -- and always start a
# monitor (the command is printed at the end).
set -euo pipefail
cd "$(dirname "$0")/../.."
RECIPE=${1:?usage: bash jobs/octopus/refine_workflow.sh RECIPE.yaml}
[ -f "$RECIPE" ] || { echo "no such recipe: $RECIPE"; exit 2; }
NAME=$(basename "${RECIPE%.yaml}")
OUT=$(pwd)/outputs/refine_$NAME
if [ -d "$OUT" ] && [ -n "$(ls -A "$OUT")" ]; then
    echo "$OUT is not empty; move it first (mv $OUT $OUT.old)"
    exit 2
fi
STAMP=$(date +%Y%m%d_%H%M%S)
RUN_ROOT=/octfs/work/G16445/v61021/scratch/m2_${NAME}_$STAMP
# FVCOM truncates INPUT_DIR at 80 characters (seen on this machine)
if [ $(( ${#RUN_ROOT} + 23 )) -gt 80 ]; then   # + "/smoke/refined/output/"
    RUN_ROOT=/octfs/work/G16445/v61021/scratch/m2_$STAMP
fi
# qsub -v separates variables with commas, so a comma inside a value would
# silently become another variable (review F14)
for v in "$RECIPE" "$OUT" "$RUN_ROOT" "${FMESH_VIEWS:-}"; do
    case "$v" in *,*) echo "a comma in '$v' cannot pass through qsub -v; rename it"; exit 2 ;; esac
done
id() { grep -oE '[0-9]+\.[a-z]+' | head -1; }

r=$(qsub -N "fm_$NAME" -v "FMESH_RECIPE=$RECIPE,LR_SEEDS=${LR_SEEDS:-0:1:2:3:4}" \
    jobs/octopus/417_hires_refine.sh | id)
f=$(qsub --after "$r" -v "FMESH_SCRIPT=434_final_mesh.py,FMESH_OUT=$OUT,FMESH_VIEWS=${FMESH_VIEWS:-}" \
    jobs/octopus/418_hires_coastline_check.sh | id)
d=$(qsub --after "$r" -v "FMESH_OUT=$OUT,FMESH_RUN_ROOT=$RUN_ROOT,FMESH_HMIN=${FMESH_HMIN:-3},FMESH_HMAX=${FMESH_HMAX:-300},FMESH_RFACTOR=${FMESH_RFACTOR:-0.2}" \
    jobs/octopus/421_finish_and_run.sh | id)
s=$(qsub --after "$d" -v "FMESH_RUN_ROOT=$RUN_ROOT" jobs/octopus/423_m2_smoke.sh | id)
echo "recipe   $RECIPE"
echo "output   $OUT"
echo "refine   $r    logs/417_hires_refine.${r%%.*}.log"
echo "figures  $f    logs/418_hires_coastline_check.${f%%.*}.log"
echo "depths   $d    logs/421_finish_and_run.${d%%.*}.log"
echo "smoke    $s    logs/423_m2_smoke.${s%%.*}.log"
last=$s
if [ "${FMESH_M2:-0}" = 1 ]; then
    runs=()
    for c in base refined; do
        runs+=("$(qsub --after "$s" -N "m2_$c" \
            -v "FMESH_RUN_ROOT=$RUN_ROOT,FMESH_CASE=$c,FMESH_REQUIRE_SMOKE=1" \
            jobs/octopus/412_m2_run.sh | id)")
    done
    a=$(qsub --after "$(IFS=,; echo "${runs[*]}")" \
        -v "FMESH_RUN_ROOT=$RUN_ROOT,FMESH_OUT=$(pwd)/outputs/m2_${NAME}_$STAMP" \
        jobs/octopus/413_m2_analysis.sh | id)
    echo "M2 runs  ${runs[*]}"
    echo "M2 cmp   $a    -> outputs/m2_${NAME}_$STAMP"
    last=$a
fi
echo "run root $RUN_ROOT"
echo
echo "monitor:  until ! qstat | grep -q ${last%%.*}; do sleep 60; done; tail -5 logs/*.${r%%.*}.log"
