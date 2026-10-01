#!/bin/bash
#PBS -q OCT-S
#PBS --group=G16445
#PBS -l cpunum_job=8
#PBS -l memsz_job=16GB
#PBS -l elapstim_req=01:00:00
#PBS -N fmesh_454chan
#PBS -j o
#PBS -o logs/454_flather_channel.pbs.log
#PBS -r n
# The idealised channel of notebooks/454_flather_channel.py: stage one mode
# (clamped | flather | flather0), run FVCOM, check the run.
#
# Required: FMESH_MODE, FMESH_RUN_ROOT (short: FVCOM paths <= 80 chars),
#           FMESH_FVCOM (an FVCOM with the Flather boundary), WORK_DIR.
set -euo pipefail
cd "${PBS_O_WORKDIR:?Submit from the repository root}"
rm -f "${FMESH_RUN_ROOT:?set FMESH_RUN_ROOT}/RUN_OK"
. jobs/octopus/common.sh 454_flather_channel 1
case $(hostname -s) in oct-cpu*) ;; *) echo 'Compute nodes only'; exit 1 ;; esac
RUN_ROOT=$FMESH_RUN_ROOT
RANKS=${FMESH_RANKS:-8}
FVCOM=${FMESH_FVCOM:-${WORK_DIR:?set WORK_DIR}/Github/FVCOM/src/fvcom}
python notebooks/454_flather_channel.py stage --root "$RUN_ROOT" --mode "${FMESH_MODE:?set FMESH_MODE}"
CASE_DIR=$RUN_ROOT
set +u; conda deactivate; set -u
if ! type module >/dev/null 2>&1; then
    for init in /etc/profile.d/modules.sh /usr/share/Modules/init/bash /usr/share/lmod/lmod/init/bash; do
        if [[ -r $init ]]; then source "$init"; break; fi
    done
fi
module purge
module load BaseCPU/2026
module load hdf5/1.14.6 netcdf-c/4.9.3 netcdf-fortran/4.6.2
INSTALLDIR=$FVCOM_LIBS
export INSTALLDIR OMP_NUM_THREADS=1 PROJ_DATA=/usr/share/proj
export LD_LIBRARY_PATH="$INSTALLDIR/lib:$INSTALLDIR/lib64:${LD_LIBRARY_PATH:-}"
ulimit -s unlimited
fail=0
t0=$(date +%s)
if ( cd "$CASE_DIR" && mpiexec -np "$RANKS" "$FVCOM" --casename=m2 > fvcom.log 2>&1 ); then
    rc=0; else rc=$?; fi
echo "[454] exit=$rc seconds=$(( $(date +%s) - t0 ))"
if grep -Ei 'fatal|non[ -]?finite|floating exception|segmentation|nan detected' \
        "$CASE_DIR/fvcom.log" >/dev/null; then
    echo "[454] UNHEALTHY log"; tail -25 "$CASE_DIR/fvcom.log"; fail=1
fi
[ "$rc" -eq 0 ] || { tail -25 "$CASE_DIR/fvcom.log"; fail=1; }
module purge
unset LD_LIBRARY_PATH
set +u; conda activate "${FMESH_ENV:-fvcom-mesh-tools}"; set -u
python -m fvcom_mesh_tools.cli.check_run "$CASE_DIR" || fail=1
[ "$fail" -eq 0 ] && date -Is > "$RUN_ROOT/RUN_OK"
echo "end=$(date -Is) fail=$fail"
exit $fail
