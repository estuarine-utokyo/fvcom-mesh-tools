#!/bin/bash
#PBS -q OCT-S
#PBS --group=G16445
#PBS -l cpunum_job=8
#PBS -l memsz_job=32GB
#PBS -l elapstim_req=02:00:00
#PBS -N fmesh_383
#PBS -j o
#PBS -o logs/383_m2.pbs.log
#PBS -r n
# Submit from the fvcom-mesh-tools repository root: qsub jobs/octopus/383_m2.sh
set -euo pipefail
cd "${PBS_O_WORKDIR:?Submit from repository root}"
. jobs/octopus/common.sh 383_m2 1
trap 'echo "ERROR line $LINENO: $BASH_COMMAND" >&2' ERR
case $(hostname -s) in oct-cpu*) ;; *) echo 'Compute nodes only'; exit 1 ;; esac
REPO=$(pwd)
PYTHON=$(command -v python)
RUN_ROOT=/octfs/work/G16445/v61021/scratch/m2_383
FVCOM=/octfs/work/G16445/v61021/Github/FVCOM/src/fvcom
export PYTHONDONTWRITEBYTECODE=1
# Avoid stale output being interpreted as a successful rerun. Move previous
# scratch run directories aside manually before resubmitting a completed case.
for case in A B_own B_Adepth; do
    if compgen -G "$RUN_ROOT/$case/output/m2_*.nc" >/dev/null; then
        echo "Existing model output: $RUN_ROOT/$case/output; archive before rerunning"
        exit 2
    fi
done
# Capture Python/conda first; keep its shared libraries out of the MPI runtime.
"$PYTHON" notebooks/383_m2_case_prep.py --root "$RUN_ROOT"
set +u
conda deactivate
set -u
# Exact module/library setup from FVCOM/octopus/common.sh. Do not source that
# file: its logging/build paths would write to the read-only FVCOM repository.
if ! type module >/dev/null 2>&1; then
    for init in /etc/profile.d/modules.sh /usr/share/Modules/init/bash /usr/share/lmod/lmod/init/bash; do
        if [[ -r $init ]]; then source "$init"; break; fi
    done
fi
module purge
module load BaseCPU/2026
module load hdf5/1.14.6 netcdf-c/4.9.3 netcdf-fortran/4.6.2
module list
INSTALLDIR=/octfs/work/G16445/share/local/fvcom/libs/install-oneapi-2025.3.1
export INSTALLDIR OMP_NUM_THREADS=1 PROJ_DATA=/usr/share/proj
export LD_LIBRARY_PATH="$INSTALLDIR/lib:$INSTALLDIR/lib64:${LD_LIBRARY_PATH:-}"
ulimit -s unlimited
mpiexec --version
ldd "$FVCOM"
sha256sum "$FVCOM"
status=0
for case in A B_own B_Adepth; do
    echo "START $case $(date -Is)"
    (
        cd "$RUN_ROOT/$case"
        mpiexec -np 8 "$FVCOM" --casename=m2 > fvcom.log 2>&1
    ) || status=1
    tail -30 "$RUN_ROOT/$case/fvcom.log"
    # Some Fortran STOP paths return zero: analysis also requires complete output.
    if grep -Ei 'fatal|non[ -]?finite|floating.*exception|segmentation|nan detected' "$RUN_ROOT/$case/fvcom.log"; then
        status=1
    fi
    echo "END $case $(date -Is)"
done
# Restore the Python environment without allowing module libraries to override conda.
module purge
unset LD_LIBRARY_PATH
set +u
conda activate "${FMESH_ENV:-oceanmesh-bench}"
set -u
cd "$REPO"
"$PYTHON" notebooks/384_m2_analysis.py --root "$RUN_ROOT" || status=1
echo "end=$(date -Is) status=$status"
exit "$status"
