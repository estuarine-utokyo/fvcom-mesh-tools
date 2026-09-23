#!/bin/bash
#PBS -q OCT-S
#PBS --group=G16445
#PBS -l cpunum_job=4
#PBS -l memsz_job=16GB
#PBS -l elapstim_req=00:15:00
#PBS -N fmesh_422rfac
#PBS -j o
#PBS -o logs/422_base_rfactor.pbs.log
#PBS -r n
set -euo pipefail
cd "${PBS_O_WORKDIR:?Submit from the repository root}"
. jobs/octopus/common.sh 422_base_rfactor 4
case $(hostname -s) in oct-cpu*) ;; *) echo 'Compute nodes only'; exit 1 ;; esac
G=$HOME/Github/TB-FVCOM/input/goto2023/grid
python notebooks/424_base_rfactor.py "$G/TokyoBay_grd.dat" \
    "$G/TokyoBay_dep.dat" "$G/TokyoBay_dep_m7001tp_rfac0p2_cap300.dat"
echo "--- the finished patch itself"
F=outputs/refine_kimitsu_port_hires/fvcom_finished
python notebooks/424_base_rfactor.py "$F"/*_grd.dat "$F"/*_dep.dat
echo "end=$(date -Is)"
