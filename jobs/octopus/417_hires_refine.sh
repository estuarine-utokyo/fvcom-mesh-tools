#!/bin/bash
#PBS -q OCT-S
#PBS --group=G16445
#PBS -l cpunum_job=8
#PBS -l memsz_job=64GB
#PBS -l elapstim_req=02:00:00
#PBS -N fmesh_417hires
#PBS -j o
#PBS -o logs/417_hires_refine.pbs.log
#PBS -r n
# The hires branch of the local refinement: coastline resolved against OSM,
# depths from the Tokyo Bay ladder, nothing floored or smoothed.
#
# Memory, not cores, is what this needs: the ladder's 30 m grid is 2466 x 1982
# and the M7001 sounding file is 3.95 million rows.  DistMesh is serial.
#
# Required: FMESH_RECIPE.  Optional: LR_OUT, LR_SEEDS, FMESH_LAND.
#
# FMESH_LAND defaults to the xcoast "true land" product -- OSM land polygons
# minus inland water, which DATA_INVENTORY.md records as coastline precedence
# #1 and which the owner named on 2026-09-23.  A `resolve` run that silently
# fell back to some other shoreline would be a run whose coastline claim is
# about the wrong line, so the job prints what it used.
set -euo pipefail
cd "${PBS_O_WORKDIR:?Submit from the repository root}"
. jobs/octopus/common.sh 417_hires_refine 8
case $(hostname -s) in oct-cpu*) ;; *) echo 'Compute nodes only'; exit 1 ;; esac
RECIPE=${FMESH_RECIPE:?set FMESH_RECIPE}
[ -f "$RECIPE" ] || { echo "no such recipe: $RECIPE"; exit 2; }
export FMESH_LAND=${FMESH_LAND:-"$DATA_DIR/geodata/OSM/coastmask_cache/custom_139.55_34.9_140.3_35.75_minarea1e-05/land.shp"}
[ -f "$FMESH_LAND" ] || { echo "no shoreline: $FMESH_LAND"; exit 2; }
export LR_OUT=${LR_OUT:-"outputs/refine_$(basename "${RECIPE%.yaml}")"}
export LR_SEEDS=${LR_SEEDS:-0,1,2,3,4}
if [ -e "$LR_OUT/report.json" ]; then
    # Output that already exists makes a rerun look like a success.
    echo "Existing output: $LR_OUT/report.json; move it before rerunning"
    exit 2
fi
echo "recipe=$RECIPE land=$FMESH_LAND out=$LR_OUT seeds=$LR_SEEDS"
python notebooks/420_local_refine.py "$RECIPE"
echo "end=$(date -Is)"
