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
. jobs/octopus/common.sh "412_m2_run" 1
case $(hostname -s) in oct-cpu*) ;; *) echo 'Compute nodes only'; exit 1 ;; esac
RUN_ROOT=${FMESH_RUN_ROOT:?set FMESH_RUN_ROOT}
CASE=${FMESH_CASE:?set FMESH_CASE}
RANKS=${FMESH_RANKS:-64}
FVCOM=/octfs/work/G16445/v61021/Github/FVCOM/src/fvcom
[ -f "$RUN_ROOT/$CASE/m2_run.nml" ] || { echo "not staged: $RUN_ROOT/$CASE"; exit 2; }
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
( cd "$RUN_ROOT/$CASE" && mpiexec -np "$RANKS" "$FVCOM" --casename=m2 > fvcom.log 2>&1 ) || status=1
tail -20 "$RUN_ROOT/$CASE/fvcom.log"
# Some Fortran STOP paths return zero: a clean exit is not enough.
if grep -Ei 'fatal|non[ -]?finite|floating.*exception|segmentation|nan detected' \
        "$RUN_ROOT/$CASE/fvcom.log"; then
    status=1
fi
echo "case=$CASE end=$(date -Is) status=$status"
exit "$status"
