#!/bin/bash
#PBS -q OCT-S
#PBS --group=G16445
#PBS -l cpunum_job=8
#PBS -l memsz_job=64GB
#PBS -l elapstim_req=01:00:00
#PBS -N fmesh_453redepth
#PBS -j o
#PBS -o logs/453_redepth_extended.pbs.log
#PBS -r n
# New depths for the new nodes of an extended case from another bathymetry
# stack (notebooks/453_redepth_extended.py).
#
# Required: FMESH_RECIPE, FMESH_BUILT (the built case dir), FMESH_OUT,
#           FMESH_SOURCES (joined by "+": qsub -v cannot pass commas),
#           DATA_DIR (qsub -v).
# Optional: FMESH_CASE_NAME.
set -euo pipefail
cd "${PBS_O_WORKDIR:?Submit from the repository root}"
. jobs/octopus/common.sh 453_redepth_extended 8
case $(hostname -s) in oct-cpu*) ;; *) echo 'Compute nodes only'; exit 1 ;; esac
: "${FMESH_SOURCES:?set FMESH_SOURCES}"
name=()
[ -n "${FMESH_CASE_NAME:-}" ] && name=(--case-name "$FMESH_CASE_NAME")
python notebooks/453_redepth_extended.py "${FMESH_RECIPE:?set FMESH_RECIPE}" \
    "${FMESH_BUILT:?set FMESH_BUILT}" "${FMESH_OUT:?set FMESH_OUT}" \
    --sources "${FMESH_SOURCES//+/,}" "${name[@]}"
echo "end=$(date -Is)"
