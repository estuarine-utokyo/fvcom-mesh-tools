#!/bin/bash
#PBS -q OCT-S
#PBS --group=G16445
#PBS -l cpunum_job=16
#PBS -l memsz_job=96GB
#PBS -l elapstim_req=04:00:00
#PBS -N fmesh_440base
#PBS -j o
#PBS -o logs/440_base_mesh.pbs.log
#PBS -r n
# Build a whole-bay base mesh from a base recipe (notebooks/440_base_mesh.py):
# OSM land from $DATA_DIR, generation, finishing, M7001 depths, the FVCOM
# case, and report.json with the provenance and product hashes.
#
# Required: FMESH_RECIPE (e.g. recipes/base/tokyo_bay_tool.yaml).
# Optional: FMESH_OUT (default outputs/base_<name>; must be empty or absent).
#
# Submit from the repository root on a login node:
#   qsub -v FMESH_RECIPE=recipes/base/tokyo_bay_tool.yaml jobs/octopus/440_base_mesh.sh
# Threads are fixed at 16 (the count the reference mesh was built with); the
# report records them.
set -euo pipefail
cd "${PBS_O_WORKDIR:?Submit from the repository root}"
case $(hostname -s) in oct-cpu*) ;; *) echo "compute nodes only"; exit 1 ;; esac
. jobs/octopus/common.sh 440_base_mesh 16
RECIPE=${FMESH_RECIPE:?set FMESH_RECIPE}
echo "commit=$(git rev-parse HEAD) dirty=$(git status --porcelain --untracked-files=all | wc -l)"
python notebooks/440_base_mesh.py "$RECIPE" ${FMESH_OUT:+"$FMESH_OUT"}
echo "end=$(date -Is)"
