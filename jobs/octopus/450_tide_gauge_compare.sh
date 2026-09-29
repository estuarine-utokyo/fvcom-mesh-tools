#!/bin/bash
#PBS -q OCT-S
#PBS --group=G16445
#PBS -l cpunum_job=8
#PBS -l memsz_job=32GB
#PBS -l elapstim_req=01:00:00
#PBS -N fmesh_450cmp
#PBS -j o
#PBS -o logs/450_tide_gauge_compare.pbs.log
#PBS -r n
# Harmonic constants of 449 tide runs against tide gauges
# (notebooks/450_tide_gauge_compare.py).
#
# Required: FMESH_RUNS ("label=dir+label=dir+..."), FMESH_OUT, DATA_DIR
# (pass them with qsub -v).
set -euo pipefail
cd "${PBS_O_WORKDIR:?Submit from the repository root}"
. jobs/octopus/common.sh 450_tide_gauge_compare 4
case $(hostname -s) in oct-cpu*) ;; *) echo 'Compute nodes only'; exit 1 ;; esac
args=()
IFS=+ read -ra runs <<< "${FMESH_RUNS:?set FMESH_RUNS}"
for r in "${runs[@]}"; do args+=(--run "$r"); done
python notebooks/450_tide_gauge_compare.py "${args[@]}" --out "${FMESH_OUT:?set FMESH_OUT}"
echo "end=$(date -Is)"
