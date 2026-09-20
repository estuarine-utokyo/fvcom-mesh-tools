#!/bin/bash
#PBS -q OCT-S
#PBS --group=G16445
#PBS -l cpunum_job=1
#PBS -l memsz_job=16GB
#PBS -l elapstim_req=00:30:00
#PBS -N fmesh_391
#PBS -j o
#PBS -o logs/391_sizing_report.pbs.log
#PBS -r n
# Submit from repository root: qsub jobs/octopus/391_sizing_report.sh
set -euo pipefail
cd "${PBS_O_WORKDIR:?Submit from the repository root}"
. jobs/octopus/common.sh 391_sizing_report 1
python notebooks/391_sizing_report.py \
    "${SIZING_MESH:-outputs/sample_repro/sample_repro.14}" \
    "${SR_SIZING:-recipes/sizing/tokyo_bay.yaml}" \
    > "logs/391_sizing_report.${JOBID}.csv"
