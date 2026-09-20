#!/bin/bash
#PBS -q OCT-S
#PBS --group=G16445
#PBS -l cpunum_job=4
#PBS -l memsz_job=32GB
#PBS -l elapstim_req=01:00:00
#PBS -N fmesh_389
#PBS -j o
#PBS -o logs/389_coverage_gaps.pbs.log
#PBS -r n
# Measure certified-mesh omitted water and draw the overview plus zoom panels.
# Usage from repository root: qsub jobs/octopus/389_coverage_gaps.sh
set -euo pipefail
cd "${PBS_O_WORKDIR:?Submit from the repository root}"
. jobs/octopus/common.sh 389_coverage_gaps 4
python notebooks/389_coverage_gaps.py
