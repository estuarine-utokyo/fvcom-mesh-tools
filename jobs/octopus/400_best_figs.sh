#!/bin/bash
#PBS -q OCT-S
#PBS --group=G16445
#PBS -l cpunum_job=8
#PBS -l memsz_job=32GB
#PBS -l elapstim_req=00:40:00
#PBS -N fmesh_400
#PBS -j o
#PBS -o logs/400_best_figs.pbs.log
#PBS -r n
# Figures for the optimised mesh (achieved targets + SR_OBC_SKIP=2).
set -euo pipefail
cd "${PBS_O_WORKDIR:?Submit from the repository root}"
. jobs/octopus/common.sh 400_best_figs 8
mkdir -p outputs/best_pair/optimised
cp -p outputs/best_399/sample_repro_final.14 outputs/best_pair/optimised/
python notebooks/393_sizing_compare.py outputs/best_pair
python notebooks/394_spec_match.py outputs/best_399/sample_repro_final.14 optimised
echo "end=$(date -Is)"
