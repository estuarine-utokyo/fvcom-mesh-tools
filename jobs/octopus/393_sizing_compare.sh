#!/bin/bash
#PBS -q OCT-S
#PBS --group=G16445
#PBS -l cpunum_job=8
#PBS -l memsz_job=32GB
#PBS -l elapstim_req=00:40:00
#PBS -N fmesh_393
#PBS -j o
#PBS -o logs/393_sizing_compare.pbs.log
#PBS -r n
# Side-by-side meshes of the sizing sweep (owner: show the mesh every time).
set -euo pipefail
cd "${PBS_O_WORKDIR:?Submit from the repository root}"
. jobs/octopus/common.sh 393_sizing_compare 8
python notebooks/393_sizing_compare.py "${1:-outputs/sizing_392}"
echo "end=$(date -Is)"
