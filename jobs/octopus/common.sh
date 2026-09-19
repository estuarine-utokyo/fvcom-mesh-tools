# Sourced by the OCTOPUS (NQSV) job scripts in this directory.
#
# Sets up logging to logs/<name>.<jobid>.log in the repository root,
# the conda env, DATA_DIR and thread counts. Callers must `cd` to the
# repository root first (PBS_O_WORKDIR = submission directory).
#
# Usage inside a job script:
#   cd "${PBS_O_WORKDIR:-$PWD}"
#   . jobs/octopus/common.sh <log-name> [threads]

_name="${1:?usage: . jobs/octopus/common.sh NAME [THREADS]}"
_threads="${2:-4}"

[ -f pyproject.toml ] && [ -d src/fvcom_mesh_tools ] || {
    echo "ERROR: submit from the fvcom-mesh-tools repository root" >&2
    exit 2
}

# NQSV job ids look like "0:114895.oct"; keep the numeric part.
JOBID="${PBS_JOBID:-interactive.$$}"
JOBID="${JOBID##*:}"
JOBID="${JOBID%%.*}"
mkdir -p logs
LOG="$(pwd)/logs/${_name}.${JOBID}.log"
exec > "${LOG}" 2>&1
echo "job=${JOBID} host=$(hostname) start=$(date -Is)"

export OMP_NUM_THREADS="${_threads}" OPENBLAS_NUM_THREADS="${_threads}"
export MKL_NUM_THREADS="${_threads}" NUMBA_NUM_THREADS="${_threads}"
export PYTHONUNBUFFERED=1 MPLBACKEND=Agg
export DATA_DIR="${DATA_DIR:-/octfs/work/G16445/share/Data}"

# conda's activate scripts are not `set -u` clean.
set +u
. /octfs/work/G16445/v61021/miniforge3/etc/profile.d/conda.sh
conda activate "${FMESH_ENV:-oceanmesh-bench}"
set -u
echo "python=$(command -v python) DATA_DIR=${DATA_DIR}"
