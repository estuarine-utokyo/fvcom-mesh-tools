# Sourced by the OCTOPUS (NQSV) job scripts in this directory.
#
# Sets up logging to logs/<name>.<jobid>.log in the repository root,
# the conda env, DATA_DIR and thread counts. Callers must `cd` to the
# repository root first (PBS_O_WORKDIR = submission directory).
#
# Usage inside a job script:
#   cd "${PBS_O_WORKDIR:-$PWD}"
#   . jobs/octopus/common.sh <log-name> [threads]

# NQSV job ids look like "0:114895.oct"; keep the numeric part.
# Only a scheduler id is parsed; the interactive fallback keeps its process id,
# so two shells do not share a log (review round 49 F3).
if [ -n "${PBS_JOBID:-}" ]; then
    JOBID="${PBS_JOBID##*:}"
    JOBID="${JOBID%%.*}"
else
    JOBID="interactive.$$"
fi
. jobs/common_core.sh "$@"

# OCTOPUS site layer: the shared FVCOM library install, beside $DATA_DIR in
# the group area (override with FVCOM_LIBS).
# Made absolute here, before any job changes directory: a relative override
# would otherwise name another place inside the case (review round 9 F14).
export FVCOM_LIBS
FVCOM_LIBS=$(realpath -m -- "${FVCOM_LIBS:-$(dirname "$DATA_DIR")/local/fvcom/libs/install-oneapi-2025.3.1}")

