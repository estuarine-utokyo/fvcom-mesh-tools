#!/bin/bash
#PBS -q OCT-S
#PBS --group=G16445
#PBS -l cpunum_job=4
#PBS -l memsz_job=16GB
#PBS -l elapstim_req=00:20:00
#PBS -N fmesh_424struct
#PBS -j o
#PBS -o logs/424_structure_survey.pbs.log
#PBS -r n
# What the width filter removed from a harbour.  Required: FMESH_SITE
# ("lon:lat:radius:h0") and FMESH_MESH_DIR (a refinement output dir).
set -euo pipefail
cd "${PBS_O_WORKDIR:?Submit from the repository root}"
. jobs/octopus/common.sh 424_structure_survey 4
case $(hostname -s) in oct-cpu*) ;; *) echo 'Compute nodes only'; exit 1 ;; esac
export FMESH_LAND=${FMESH_LAND:-"$DATA_DIR/geodata/OSM/coastmask_cache/custom_139.55_34.9_140.3_35.75_minarea1e-05/land.shp"}
IFS=: read -r a b c d <<< "${FMESH_SITE:?set FMESH_SITE}"
python notebooks/425_structure_survey.py "$a" "$b" "$c" "$d"
echo "end=$(date -Is)"
