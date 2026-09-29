#!/bin/bash
#PBS -q OCT-S
#PBS --group=G16445
#PBS -l cpunum_job=16
#PBS -l memsz_job=96GB
#PBS -l elapstim_req=04:00:00
#PBS -N fmesh_445ext
#PBS -j o
#PBS -o logs/445_extend_mesh.pbs.log
#PBS -r n
# Build a wide mesh from an extension recipe (notebooks/445_extend_mesh.py):
# the base mesh as it is, plus the sea out to the recipe's open boundary.
#
# Required: FMESH_RECIPE (e.g. recipes/extend/tokyo_bay_enshu.yaml).
# Optional: FMESH_OUT (default outputs/extend_<name>; must be empty or absent).
#   qsub -v FMESH_RECIPE=recipes/extend/tokyo_bay_enshu.yaml jobs/octopus/445_extend_mesh.sh
set -euo pipefail
cd "${PBS_O_WORKDIR:?Submit from the repository root}"
case $(hostname -s) in oct-cpu*) ;; *) echo "compute nodes only"; exit 1 ;; esac
. jobs/octopus/common.sh 445_extend_mesh 16
RECIPE=${FMESH_RECIPE:?set FMESH_RECIPE}
echo "commit=$(git rev-parse HEAD) dirty=$(git status --porcelain --untracked-files=all | wc -l)"
python notebooks/445_extend_mesh.py "$RECIPE" ${FMESH_OUT:+"$FMESH_OUT"}
echo "end=$(date -Is)"
