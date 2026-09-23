#!/bin/bash
#PBS -q OCT-S
#PBS --group=G16445
#PBS -l cpunum_job=4
#PBS -l memsz_job=16GB
#PBS -l elapstim_req=00:20:00
#PBS -N fmesh_419valid
#PBS -j o
#PBS -o logs/419_mesh_validity.pbs.log
#PBS -r n
# Why a mesh that verify_patch accepts is one matplotlib's TriFinder refuses.
# Required: FMESH_MESH (a .14).
set -euo pipefail
cd "${PBS_O_WORKDIR:?Submit from the repository root}"
. jobs/octopus/common.sh 419_mesh_validity 4
case $(hostname -s) in oct-cpu*) ;; *) echo 'Compute nodes only'; exit 1 ;; esac
python notebooks/422_mesh_validity.py "${FMESH_MESH:?set FMESH_MESH}"
echo "end=$(date -Is)"
