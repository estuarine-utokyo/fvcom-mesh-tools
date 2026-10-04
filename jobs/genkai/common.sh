# Sourced by the GENKAI (Fujitsu TCS, pjsub) job scripts in this directory:
# the GENKAI site layer over jobs/common_core.sh. Callers must `cd` to the
# repository root first (PJM_O_WORKDIR = submission directory).
#
# Usage inside a job script:
#   cd "${PJM_O_WORKDIR:?Submit from the repository root}"
#   . jobs/genkai/common.sh <log-name> [threads]

JOBID="${PJM_JOBID:-interactive.$$}"
# GENKAI's conda: $WORK_DIR/miniforge3, else the older $WORK_DIR/mambaforge.
if [ -z "${FMESH_CONDA_ROOT:-}" ] && [ ! -d "${WORK_DIR:?WORK_DIR is not set (login profile)}/miniforge3" ] \
        && [ -d "$WORK_DIR/mambaforge" ]; then
    FMESH_CONDA_ROOT=$WORK_DIR/mambaforge
fi
. jobs/common_core.sh "$@"

# The FVCOM libraries have no agreed place on GENKAI yet: a run that needs
# them passes FVCOM_LIBS, made absolute here before any job changes directory.
if [ -n "${FVCOM_LIBS:-}" ]; then
    export FVCOM_LIBS
    FVCOM_LIBS=$(realpath -m -- "$FVCOM_LIBS")
fi
