#!/bin/bash
#PBS -q OCT-S
#PBS --group=G16445
#PBS -l cpunum_job=8
#PBS -l memsz_job=32GB
#PBS -l elapstim_req=00:40:00
#PBS -N fmesh_394
#PBS -j o
#PBS -o logs/394_spec_match.pbs.log
#PBS -r n
# How closely each mesh follows its sizing spec (owner 2026-09-20).
set -euo pipefail
cd "${PBS_O_WORKDIR:?Submit from the repository root}"
. jobs/octopus/common.sh 394_spec_match 8
python notebooks/394_spec_match.py outputs/sample_repro/sample_repro_final.14 certified
SR_GRADE=0.12 python notebooks/394_spec_match.py outputs/sizing_392/A/sample_repro_final.14 A
SR_H0=250 python notebooks/394_spec_match.py outputs/sizing_392/B/sample_repro_final.14 B
echo "end=$(date -Is)"
