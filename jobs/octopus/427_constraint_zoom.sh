#!/bin/bash
#PBS -q OCT-S
#PBS --group=G16445
#PBS -l cpunum_job=2
#PBS -l memsz_job=8GB
#PBS -l elapstim_req=00:10:00
#PBS -N fmesh_zoom
#PBS -j o
#PBS -o logs/zoom.pbs.log
#PBS -r n
set -euo pipefail
cd "${PBS_O_WORKDIR:?}"
. jobs/octopus/common.sh zoom 2
python notebooks/430_constraint_zoom.py outputs/refine_kimitsu_port_hires 392622 3909130 80
