#!/bin/bash
#PBS -q OCT-S
#PBS --group=G16445
#PBS -l cpunum_job=64
#PBS -l memsz_job=120GB
#PBS -l elapstim_req=03:00:00
#PBS -N fmesh_412run
#PBS -j o
#PBS -o logs/412_m2_run.pbs.log
#PBS -r n
# One FVCOM integration.  64 ranks = half a 128-core CPU, the owner's
# production choice; OCT-S (O1SS) takes exactly 64 on a shared node.
# Required: FMESH_RUN_ROOT, FMESH_CASE.
set -euo pipefail
cd "${PBS_O_WORKDIR:?Submit from the repository root}"
# INVALIDATE first, before anything that can fail (review 2, R1)
rm -f "${FMESH_RUN_ROOT:?set FMESH_RUN_ROOT}/${FMESH_CASE:?set FMESH_CASE}/RUN_OK"
. jobs/octopus/common.sh "412_m2_run" 1
case $(hostname -s) in oct-cpu*) ;; *) echo 'Compute nodes only'; exit 1 ;; esac
RUN_ROOT=${FMESH_RUN_ROOT:?set FMESH_RUN_ROOT}
CASE=${FMESH_CASE:?set FMESH_CASE}
RANKS=${FMESH_RANKS:-64}
FVCOM=/octfs/work/G16445/v61021/Github/FVCOM/src/fvcom
[ -f "$RUN_ROOT/$CASE/m2_run.nml" ] || { echo "not staged: $RUN_ROOT/$CASE"; exit 2; }
# In a chain with a smoke test, a failed smoke must stop the long runs: the
# namelists were staged before it ran, so their existence proves nothing.
if [ "${FMESH_REQUIRE_SMOKE:-0}" = 1 ] && [ ! -f "$RUN_ROOT/SMOKE_OK" ]; then
    echo "smoke test did not pass: no $RUN_ROOT/SMOKE_OK"; exit 2
fi
set +u
conda deactivate
set -u
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
echo "case=$CASE ranks=$RANKS cores=$(nproc) start=$(date -Is)"
status=0
# an earlier attempt's output must not be judged as this attempt's (review 2, R2)
rm -f "$RUN_ROOT/$CASE/output/"*.nc
( cd "$RUN_ROOT/$CASE" && mpiexec -np "$RANKS" "$FVCOM" --casename=m2 > fvcom.log 2>&1 ) || status=1
tail -20 "$RUN_ROOT/$CASE/fvcom.log"
# Some Fortran STOP paths return zero: a clean exit is not enough.
if grep -Ei 'fatal|non[ -]?finite|floating.*exception|segmentation|nan detected' \
        "$RUN_ROOT/$CASE/fvcom.log"; then
    status=1
fi
# the run is judged on its output as well as its log (review F6)
module purge
unset LD_LIBRARY_PATH
set +u; conda activate "${FMESH_ENV:-oceanmesh-bench}"; set -u
# RUN_OK only when the solver's own exit AND the output check both pass
# (review 2, R2); the check always runs, for its diagnosis
if [ "$status" -eq 0 ]; then
    python -m fvcom_mesh_tools.cli.check_run "$RUN_ROOT/$CASE" \
        --marker "$RUN_ROOT/$CASE/RUN_OK" || status=1
else
    python -m fvcom_mesh_tools.cli.check_run "$RUN_ROOT/$CASE" || true
fi
echo "case=$CASE end=$(date -Is) status=$status"
exit "$status"
