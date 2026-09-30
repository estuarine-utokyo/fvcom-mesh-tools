#!/bin/bash
#PBS -q OCT-S
#PBS --group=G16445
#PBS -l cpunum_job=8
#PBS -l memsz_job=32GB
#PBS -l elapstim_req=01:00:00
#PBS -N fmesh_452ts
#PBS -j o
#PBS -o logs/452_tide_timeseries.pbs.log
#PBS -r n
# Time series of 449 tide runs against tide gauges
# (notebooks/452_tide_timeseries.py). Optional: FMESH_EXTRA (more arguments,
# e.g. '--start 2021-03-01 --days 14 --stations tokyo,yokosuka,mera').
#
# Required: FMESH_RUNS ("label=dir+label=dir+..."), FMESH_OUT, DATA_DIR
# (pass them with qsub -v).
set -euo pipefail
cd "${PBS_O_WORKDIR:?Submit from the repository root}"
. jobs/octopus/common.sh 452_tide_timeseries 4
case $(hostname -s) in oct-cpu*) ;; *) echo 'Compute nodes only'; exit 1 ;; esac
args=()
IFS=+ read -ra runs <<< "${FMESH_RUNS:?set FMESH_RUNS}"
for r in "${runs[@]}"; do args+=(--run "$r"); done
python notebooks/452_tide_timeseries.py "${args[@]}" --out "${FMESH_OUT:?set FMESH_OUT}" ${FMESH_EXTRA:-}
echo "end=$(date -Is)"
