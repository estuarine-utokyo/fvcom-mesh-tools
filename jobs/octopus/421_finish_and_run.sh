#!/bin/bash
#PBS -q OCT-S
#PBS --group=G16445
#PBS -l cpunum_job=8
#PBS -l memsz_job=32GB
#PBS -l elapstim_req=00:40:00
#PBS -N fmesh_421finish
#PBS -j o
#PBS -o logs/421_finish_and_run.pbs.log
#PBS -r n
# Finish a hires refinement's depths and stage an M2 pair against the base.
# Required: FMESH_OUT (a refinement output dir), FMESH_RUN_ROOT.
# Optional: FMESH_HMIN (3), FMESH_RFACTOR (0.2).
set -euo pipefail
cd "${PBS_O_WORKDIR:?Submit from the repository root}"
. jobs/octopus/common.sh 421_finish_and_run 8
case $(hostname -s) in oct-cpu*) ;; *) echo 'Compute nodes only'; exit 1 ;; esac
OUTDIR=${FMESH_OUT:?set FMESH_OUT}
RUN_ROOT=${FMESH_RUN_ROOT:?set FMESH_RUN_ROOT}
python -m fvcom_mesh_tools.cli.refine_depths "$OUTDIR" \
    --hmin "${FMESH_HMIN:-3}" --rfactor "${FMESH_RFACTOR:-0.2}"
FIN=$(ls "$OUTDIR"/fvcom_finished/*_grd.dat)
FIN=${FIN%_grd.dat}
echo "finished case = $FIN"
python notebooks/414_refine_m2_prep.py --root "$RUN_ROOT" \
    --base "$HOME/Github/TB-FVCOM/input/goto2023/grid/TokyoBay" \
    --refined "$FIN"
echo "end=$(date -Is)"
