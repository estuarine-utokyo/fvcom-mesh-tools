#!/bin/bash
#PBS -q OCT-S
#PBS --group=G16445
#PBS -l cpunum_job=4
#PBS -l memsz_job=32GB
#PBS -l elapstim_req=00:30:00
#PBS -N fmesh_420site
#PBS -j o
#PBS -o logs/420_site_survey.pbs.log
#PBS -r n
# Size a refinement site before declaring it.
# Required: FMESH_SITE -- one or more "x:y:radius[:crs]", separated by "+".
# Not spaces, because `qsub -v` splits its argument on whitespace and treats
# the rest as a script name; not commas, because that is how it separates
# VARIABLES.  Both mistakes have now been made in this repository.
set -euo pipefail
cd "${PBS_O_WORKDIR:?Submit from the repository root}"
. jobs/octopus/common.sh 420_site_survey 4
case $(hostname -s) in oct-cpu*) ;; *) echo 'Compute nodes only'; exit 1 ;; esac
export FMESH_LAND=${FMESH_LAND:-"$DATA_DIR/geodata/OSM/coastmask_cache/custom_139.55_34.9_140.3_35.75_minarea1e-05/land.shp"}
SITES=${FMESH_SITE:?set FMESH_SITE}
for site in ${SITES//+/ }; do
    IFS=: read -r sx sy sr scrs <<< "$site"
    echo "=================== $sx $sy $sr ${scrs:-utm}"
    python notebooks/423_site_survey.py "$sx" "$sy" "$sr" "${scrs:-utm}"
done
echo "end=$(date -Is)"
