#!/bin/bash
#PBS -q OCT-S
#PBS --group=G16445
#PBS -l cpunum_job=64
#PBS -l memsz_job=120GB
#PBS -l elapstim_req=10:00:00
#PBS -N fmesh_449tide
#PBS -j o
#PBS -o logs/449_tide_nao_run.pbs.log
#PBS -r n
# An astronomical-tide hindcast of one case: stage it
# (notebooks/449_tide_nao_run.py: NAO.99Jb, 8 constituents, on the open
# boundary), integrate it, and judge the run as 448 does: exit code, the log
# grepped for fatal signatures, and fmesh-check-run.
#
# Required: FMESH_CASE (a case prefix), FMESH_RUN_ROOT (under $WORK_DIR),
#           WORK_DIR, DATA_DIR (pass them with qsub -v).
# Optional: FMESH_DAYS (200), FMESH_START (2021-01-01), FMESH_RANKS (64),
#           FMESH_EQUI=1 (add the tidal potential; then FMESH_FVCOM must be an
#           FVCOM built with -DEQUI_TIDE), FMESH_FVCOM (the FVCOM binary).
set -euo pipefail
cd "${PBS_O_WORKDIR:?Submit from the repository root}"
rm -f "${FMESH_RUN_ROOT:?set FMESH_RUN_ROOT}/RUN_OK"
. jobs/octopus/common.sh 449_tide_nao_run 1
case $(hostname -s) in oct-cpu*) ;; *) echo 'Compute nodes only'; exit 1 ;; esac
RUN_ROOT=$FMESH_RUN_ROOT
RANKS=${FMESH_RANKS:-64}
FVCOM=${FMESH_FVCOM:-${WORK_DIR:?set WORK_DIR}/Github/FVCOM/src/fvcom}
equi=()
[ "${FMESH_EQUI:-0}" = 1 ] && equi=(--equilibrium)
python notebooks/449_tide_nao_run.py --case "${FMESH_CASE:?set FMESH_CASE}" --root "$RUN_ROOT" \
    --days "${FMESH_DAYS:-200}" --start "${FMESH_START:-2021-01-01}" "${equi[@]}"
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
INSTALLDIR=/octfs/work/G16445/share/local/fvcom/libs/install-oneapi-2025.3.1
export INSTALLDIR OMP_NUM_THREADS=1 PROJ_DATA=/usr/share/proj
export LD_LIBRARY_PATH="$INSTALLDIR/lib:$INSTALLDIR/lib64:${LD_LIBRARY_PATH:-}"
ulimit -s unlimited
fail=0
t0=$(date +%s)
if ( cd "$CASE_DIR" && mpiexec -np "$RANKS" "$FVCOM" --casename=m2 > fvcom.log 2>&1 ); then
    rc=0; else rc=$?; fi
echo "[449] exit=$rc seconds=$(( $(date +%s) - t0 ))"
if grep -Ei 'fatal|non[ -]?finite|floating exception|segmentation|nan detected' \
        "$CASE_DIR/fvcom.log" >/dev/null; then
    echo "[449] UNHEALTHY log"; tail -25 "$CASE_DIR/fvcom.log"; fail=1
fi
[ "$rc" -eq 0 ] || { tail -25 "$CASE_DIR/fvcom.log"; fail=1; }
module purge
unset LD_LIBRARY_PATH
set +u; conda activate "${FMESH_ENV:-fvcom-mesh-tools}"; set -u
python -m fvcom_mesh_tools.cli.check_run "$CASE_DIR" || fail=1
[ "$fail" -eq 0 ] && date -Is > "$RUN_ROOT/RUN_OK"
echo "end=$(date -Is) fail=$fail"
exit $fail
