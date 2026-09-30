#!/bin/bash
#PBS -q OCT-S
#PBS --group=G16445
#PBS -l cpunum_job=8
#PBS -l memsz_job=32GB
#PBS -l elapstim_req=01:00:00
#PBS -N fmesh_451map
#PBS -j o
#PBS -o logs/451_tide_model_map.pbs.log
#PBS -r n
# Maps of a 449 tide run against the tide model that forces it
# (notebooks/451_tide_model_map.py).
#
# Required: FMESH_RUN (a 449 run dir), FMESH_OUT, DATA_DIR (qsub -v).
# Optional: FMESH_LINES ("label=obc.dat+label=obc.dat"), FMESH_REFERENCE
#           (forcing | nao99jb | tpxo10 | fes2022).
set -euo pipefail
cd "${PBS_O_WORKDIR:?Submit from the repository root}"
. jobs/octopus/common.sh 451_tide_model_map 4
case $(hostname -s) in oct-cpu*) ;; *) echo 'Compute nodes only'; exit 1 ;; esac
args=()
if [ -n "${FMESH_LINES:-}" ]; then
    IFS=+ read -ra lines <<< "$FMESH_LINES"
    for l in "${lines[@]}"; do args+=(--line "$l"); done
fi
python notebooks/451_tide_model_map.py --run "${FMESH_RUN:?set FMESH_RUN}" \
    --out "${FMESH_OUT:?set FMESH_OUT}" --reference "${FMESH_REFERENCE:-forcing}" "${args[@]}"
echo "end=$(date -Is)"
