#!/bin/bash
# Reproduce the GitHub CI job (.github/workflows/ci.yml) locally, BEFORE every commit.
# GitHub's CI has no oceanmesh fork and no WORK_DIR/DATA_DIR; tests that need oceanmesh must carry
# @pytest.mark.needs_oceanmesh or they fail there. This script runs the suite twice:
#   1. CI mode: oceanmesh hidden (tools/ci_local/hide_oceanmesh/ shadows it), WORK_DIR/DATA_DIR unset;
#   2. lab mode: the full environment, oceanmesh present.
# plus ruff (repository and maintained notebooks), pre-commit on all files and the CLI smoke checks.
# Usage (GENKAI: in a batch job, see jobs/ci_local.pjsub): tools/ci_local/run_ci_local.sh
# Exit status 0 only if every step passed.
set -uo pipefail
cd "$(dirname "$0")/../.."
FAIL=0
step() { echo "=== $1"; shift; "$@"; rc=$?; echo "--- rc=$rc"; [[ $rc -eq 0 ]] || FAIL=1; }
step "ruff (repository)" ruff check .
step "ruff (notebooks 42*, 43*)" ruff check notebooks/42*.py notebooks/43*.py
step "pytest, CI mode (no oceanmesh, no WORK_DIR/DATA_DIR)" \
  env -u WORK_DIR -u DATA_DIR PYTHONPATH="$PWD/tools/ci_local/hide_oceanmesh${PYTHONPATH:+:$PYTHONPATH}" \
  pytest -q -p no:cacheprovider
step "pytest, lab mode (full environment)" pytest -q -p no:cacheprovider
step "pre-commit (all files)" pre-commit run --all-files --show-diff-on-failure
smoke() {
  local rc=0 c
  for c in $(sed -n '/CLI smoke checks/,$p' .github/workflows/ci.yml | grep -oE '^\s+fmesh-[a-z0-9-]+' | tr -d ' '); do
    "$c" --help > /dev/null || { echo "FAIL: $c --help"; rc=1; }
  done
  return $rc
}
step "CLI smoke checks (--help of the entry points in ci.yml)" smoke
echo "=== ci_local: $([[ $FAIL -eq 0 ]] && echo PASS || echo FAIL)"
exit $FAIL
