#!/bin/bash
#PBS -q OCT-S
#PBS --group=G16445
#PBS -l cpunum_job=2
#PBS -l memsz_job=8GB
#PBS -l elapstim_req=00:15:00
#PBS -N fmesh_zoom
#PBS -j o
#PBS -o logs/zoom.pbs.log
#PBS -r n
# Zoom on the constraints and the mesh at several points (notebook 430).
# Required: FMESH_OUT, FMESH_POINTS ("x:y:half+x:y:half...").
set -euo pipefail
cd "${PBS_O_WORKDIR:?}"
. jobs/octopus/common.sh zoom 2
P=${FMESH_POINTS:?}
for pt in ${P//+/ }; do
    IFS=: read -r x y h <<< "$pt"
    echo "=== $x $y"
    python notebooks/430_constraint_zoom.py "${FMESH_OUT:?}" "$x" "$y" "$h" | tail -25
done
