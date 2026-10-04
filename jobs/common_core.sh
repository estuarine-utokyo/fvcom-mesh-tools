# The machine-independent part of the job set-up, sourced by each site's
# jobs/<site>/common.sh after it has set JOBID (owner's portability rule,
# 2026-09-29: the computation is one code path; only the batch system and
# the install locations differ per machine).
#
# Sets up logging to logs/<name>.<jobid>.log in the repository root, the
# conda env, thread counts, and checks DATA_DIR / WORK_DIR. Defines
# fmesh_fvcom, fmesh_stage_lock and fmesh_stage_unlock.

_name="${1:?usage: . jobs/<site>/common.sh NAME [THREADS]}"
_threads="${2:-4}"
: "${JOBID:?the site layer sets JOBID before sourcing jobs/common_core.sh}"

[ -f pyproject.toml ] && [ -d src/fvcom_mesh_tools ] || {
    echo "ERROR: submit from the fvcom-mesh-tools repository root" >&2
    exit 2
}

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
# conda's activate scripts are not `set -u` clean. The conda install is
# $WORK_DIR/miniforge3 unless the site layer sets FMESH_CONDA_ROOT.
set +u
# A failed set-up stops the job: a script that goes on with the system python
# checks the wrong software (review round 28 F5).
_env=${FMESH_ENV:-fvcom-mesh-tools}
. "${FMESH_CONDA_ROOT:-$WORK_DIR/miniforge3}/etc/profile.d/conda.sh" \
    || { echo "ERROR: cannot source conda from ${FMESH_CONDA_ROOT:-$WORK_DIR/miniforge3}"; exit 2; }
conda activate "$_env" || { echo "ERROR: cannot activate conda env $_env"; exit 2; }
# the active env is the requested one: by name, or by prefix when FMESH_ENV is
# a path (conda then reports the env's name, review round 29 F5)
if [[ $_env == */* ]]; then
    # each side resolved on its own: two failed resolutions are two empty
    # strings, which compare equal (review round 30 F3)
    _want=$(realpath -e -- "$_env") || { echo "ERROR: no such conda env path: $_env"; exit 2; }
    _have=$(realpath -e -- "${CONDA_PREFIX:-}") \
        || { echo "ERROR: CONDA_PREFIX '${CONDA_PREFIX:-}' does not exist"; exit 2; }
    [ -n "$_want" ] && [ "$_want" = "$_have" ] \
        || { echo "ERROR: conda prefix is '$_have', not $_want"; exit 2; }
else
    [ "${CONDA_DEFAULT_ENV:-}" = "$_env" ] \
        || { echo "ERROR: conda env is '${CONDA_DEFAULT_ENV:-}', not $_env"; exit 2; }
fi
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
