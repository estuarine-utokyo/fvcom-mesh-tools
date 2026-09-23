#!/bin/bash
#PBS -q OCT-S
#PBS --group=G16445
#PBS -l cpunum_job=4
#PBS -l memsz_job=32GB
#PBS -l elapstim_req=00:30:00
#PBS -N fmesh_418coast
#PBS -j o
#PBS -o logs/418_hires_coastline_check.pbs.log
#PBS -r n
# Where a `resolve` coastline departs from OSM, and whether that is the
# shoreline or the resampler.  Required: FMESH_OUT (a refinement output dir).
set -euo pipefail
cd "${PBS_O_WORKDIR:?Submit from the repository root}"
. jobs/octopus/common.sh 418_hires_coastline_check 4
case $(hostname -s) in oct-cpu*) ;; *) echo 'Compute nodes only'; exit 1 ;; esac
OUTDIR=${FMESH_OUT:?set FMESH_OUT}
[ -d "$OUTDIR" ] || { echo "no such directory: $OUTDIR"; exit 2; }
export FMESH_LAND=${FMESH_LAND:-"$DATA_DIR/geodata/OSM/coastmask_cache/custom_139.55_34.9_140.3_35.75_minarea1e-05/land.shp"}
python notebooks/${FMESH_SCRIPT:-421_hires_coastline_check.py} "$OUTDIR"
echo "end=$(date -Is)"
