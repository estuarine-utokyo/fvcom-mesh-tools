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
# The base has to carry the SAME depth product the refinement inherited.
# goto2023's grid directory holds two: TokyoBay_dep.dat, which the grid was
# written with (max r 0.8841 over its edges, 647 m deep), and the baseline's
# own TokyoBay_dep_m7001tp_rfac0p2_cap300.dat (max r 0.2000, capped at 300).
# The refinement inherited the second; staging the base with the first would
# make the comparison a comparison of bathymetry products as well as meshes.
BASEDIR=$RUN_ROOT/base_case
mkdir -p "$BASEDIR"
G=$HOME/Github/TB-FVCOM/input/goto2023/grid
cp "$G/TokyoBay_grd.dat" "$BASEDIR/TokyoBayB_grd.dat"
cp "$G/TokyoBay_obc.dat" "$BASEDIR/TokyoBayB_obc.dat"
[ -f "$G/TokyoBay_cor.dat" ] && cp "$G/TokyoBay_cor.dat" "$BASEDIR/TokyoBayB_cor.dat"
cp "$G/TokyoBay_dep_m7001tp_rfac0p2_cap300.dat" "$BASEDIR/TokyoBayB_dep.dat"
python notebooks/414_refine_m2_prep.py --root "$RUN_ROOT" \
    --base "$BASEDIR/TokyoBayB" --refined "$FIN"
echo "end=$(date -Is)"
