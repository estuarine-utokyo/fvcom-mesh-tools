#!/bin/bash
#PBS -q OCT-S
#PBS --group=G16445
#PBS -l cpunum_job=4
#PBS -l memsz_job=32GB
#PBS -l elapstim_req=01:00:00
#PBS -N fmesh_444obc
#PBS -j o
#PBS -o logs/444_design_obc.pbs.log
#PBS -r n
# Design an open boundary from a design recipe (notebooks/444_design_obc.py)
# and write its nodes as a CSV, with a report (.json) and a figure (.png).
#
# Required: FMESH_DESIGN (the design YAML), FMESH_OBC_CSV (the CSV to write).
#   qsub -v FMESH_DESIGN=recipes/extend/tokyo_bay_enshu_obc_design.yaml,FMESH_OBC_CSV=recipes/extend/tokyo_bay_enshu_obc.csv jobs/octopus/444_design_obc.sh
set -euo pipefail
cd "${PBS_O_WORKDIR:?Submit from the repository root}"
case $(hostname -s) in oct-cpu*) ;; *) echo "compute nodes only"; exit 1 ;; esac
. jobs/octopus/common.sh 444_design_obc 4
echo "commit=$(git rev-parse HEAD) dirty=$(git status --porcelain --untracked-files=all | wc -l)"
python notebooks/444_design_obc.py "${FMESH_DESIGN:?set FMESH_DESIGN}" "${FMESH_OBC_CSV:?set FMESH_OBC_CSV}"
echo "end=$(date -Is)"
