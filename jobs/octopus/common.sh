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
# Paths come from the login profile only (owner rule, 2026-09-29): stop when
# they are missing rather than fall back to a path right on one machine.
: "${DATA_DIR:?DATA_DIR is not set (login profile)}"
: "${WORK_DIR:?WORK_DIR is not set (login profile)}"
export DATA_DIR WORK_DIR
# OCTOPUS site layer: the shared FVCOM library install, beside $DATA_DIR in
# the group area (override with FVCOM_LIBS).
# Made absolute here, before any job changes directory: a relative override
# would otherwise name another place inside the case (review round 9 F14).
export FVCOM_LIBS
FVCOM_LIBS=$(realpath -m -- "${FVCOM_LIBS:-$(dirname "$DATA_DIR")/local/fvcom/libs/install-oneapi-2025.3.1}")

# conda's activate scripts are not `set -u` clean.
set +u
. "$WORK_DIR/miniforge3/etc/profile.d/conda.sh"
conda activate "${FMESH_ENV:-fvcom-mesh-tools}"
set -u
echo "python=$(command -v python) DATA_DIR=${DATA_DIR} WORK_DIR=${WORK_DIR}"

# The FVCOM executable: FMESH_FVCOM, or the repository build under WORK_DIR,
# resolved to an absolute path here -- before a job changes into its run
# directory, where a relative override would name another file (review of
# the extend tools, round 4 F12) -- and required to be executable.
fmesh_fvcom() {
    local want exe
    want=${FMESH_FVCOM:-$WORK_DIR/Github/FVCOM/src/fvcom}
    exe=$(realpath -e -- "$want") || { echo "FVCOM executable not found: $want" >&2; return 1; }
    [[ -f $exe && -x $exe ]] || { echo "not an executable file: $exe" >&2; return 1; }
    printf '%s\n' "$exe"
}

# The run-root staging lock (421 and 423 take it inline, before sourcing
# this file): one stage at a time may write or read a whole run root, and
# never while a 412 case run holds <case>/.running (review rounds 10 F3,
# 13 F2). Released on exit, or earlier with fmesh_stage_unlock.
fmesh_stage_lock() {
    mkdir -p -- "$1"
    FMESH_STAGE_LOCK=$1/.staging
    mkdir -- "$FMESH_STAGE_LOCK" 2>/dev/null || {
        echo "another stage holds $FMESH_STAGE_LOCK (remove it if that job is gone)"
        FMESH_STAGE_LOCK=; return 1; }
    trap '[ -n "${FMESH_STAGE_LOCK:-}" ] && rmdir -- "$FMESH_STAGE_LOCK"' EXIT
    local c
    for c in "$1"/*/.running; do
        if [ -e "$c" ]; then echo "a run is active: $c"; return 1; fi
    done
}
fmesh_stage_unlock() {
    rmdir -- "$FMESH_STAGE_LOCK"
    FMESH_STAGE_LOCK=
}
