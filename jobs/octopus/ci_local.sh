#!/bin/bash
#PBS -q OCT-S
#PBS --group=G16445
#PBS -l cpunum_job=8
#PBS -l memsz_job=32GB
#PBS -l elapstim_req=00:45:00
#PBS -N fmesh_ci_local
#PBS -j o
#PBS -o logs/ci_local.pbs.log
#PBS -r n
#============================================================================
# OCTOPUS header for tools/ci_local/run_ci_local.sh (the GitHub CI job reproduced locally: ruff,
# pytest without oceanmesh and WORK_DIR/DATA_DIR, pytest with the full environment, pre-commit,
# CLI smoke); the GENKAI header is jobs/genkai/ci_local.sh. Required before every commit that
# touches src/, tests/ or notebooks/ (CLAUDE.md).
#
#   qsub jobs/octopus/ci_local.sh        # from the repository root
#
# pre-commit keeps its hook environments in PRE_COMMIT_HOME, by default ~/.cache/pre-commit;
# the home directory is not written on OCTOPUS, so it is pointed under $WORK_DIR here.
#============================================================================
set -euo pipefail
cd "${PBS_O_WORKDIR:?Submit from the repository root}"
. jobs/octopus/common.sh ci_local 8
case $(hostname -s) in oct-cpu*) ;; *) echo 'Compute nodes only'; exit 1 ;; esac
export PRE_COMMIT_HOME=${PRE_COMMIT_HOME:-"${WORK_DIR:?}/.cache/pre-commit"}
echo "commit: $(git rev-parse HEAD)  dirty files: $(git status --porcelain --untracked-files=no | wc -l)"
echo "PRE_COMMIT_HOME=$PRE_COMMIT_HOME"
tools/ci_local/run_ci_local.sh
