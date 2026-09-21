#!/bin/bash
#PBS -q OCT-S
#PBS --group=G16445
#PBS -l cpunum_job=8
#PBS -l memsz_job=32GB
#PBS -l elapstim_req=02:00:00
#PBS -N fmesh_410
#PBS -j o
#PBS -o logs/410_m2_optimised.pbs.log
#PBS -r n
# M2 tidal validation of the OPTIMISED mesh (achieved-size targets, OBC skip 2,
# coastline fit), against the goto2023 production mesh.
#   A         goto2023 production mesh and depths
#   B_own     the mesh under test with its own depths
#   B_Adepth  the mesh under test carrying A's depths (isolates geometry)
# The previous round (job 6203xxx, certified 3,393-node mesh) gave
# B_own - A = +0.006 to +0.009 m in M2 amplitude and -0.3 to -0.4 deg in phase.
# Submit from the fvcom-mesh-tools repository root:
#   qsub jobs/octopus/410_m2_optimised.sh
set -euo pipefail
cd "${PBS_O_WORKDIR:?Submit from repository root}"
. jobs/octopus/common.sh 410_m2_optimised 1
trap 'echo "ERROR line $LINENO: $BASH_COMMAND" >&2' ERR
case $(hostname -s) in oct-cpu*) ;; *) echo 'Compute nodes only'; exit 1 ;; esac
REPO=$(pwd)
PYTHON=$(command -v python)
RUN_ROOT=/octfs/work/G16445/v61021/scratch/m2_410.${JOBID}
# The mesh under test: the verified optimised + coastline-fitted mesh from
# job 115302.  Override with FMESH_B_MESH to validate a different candidate.
B_MESH=${FMESH_B_MESH:-$REPO/outputs/verify_409.115302/fit/sample_repro_final.14}
FVCOM=/octfs/work/G16445/v61021/Github/FVCOM/src/fvcom
export PYTHONDONTWRITEBYTECODE=1
# Avoid stale output being interpreted as a successful rerun. Move previous
# scratch run directories aside manually before resubmitting a completed case.
# The run root carries the job id, so a resubmission never lands on an
# existing result; this stays as a guard against a hand-edited RUN_ROOT.
for case in A B_own B_Adepth; do
    if compgen -G "$RUN_ROOT/$case/output/m2_*.nc" >/dev/null; then
        echo "Existing model output: $RUN_ROOT/$case/output; archive before rerunning"
        exit 2
    fi
done
[ -f "$B_MESH" ] || { echo "mesh not found: $B_MESH"; exit 2; }
echo "B mesh = $B_MESH"
head -2 "$B_MESH"
# Capture Python/conda first; keep its shared libraries out of the MPI runtime.
"$PYTHON" notebooks/383_m2_case_prep.py --root "$RUN_ROOT" --mesh "$B_MESH"
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
"$PYTHON" notebooks/384_m2_analysis.py --root "$RUN_ROOT" \
    --output "$REPO/outputs/m2_410.${JOBID}" || status=1
echo "end=$(date -Is) status=$status"
exit "$status"
