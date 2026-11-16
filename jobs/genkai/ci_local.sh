#!/bin/bash
#PJM -L rscgrp=a-batch
#PJM -L vnode=1
#PJM -L vnode-core=8
#PJM -L elapse=00:30:00
#PJM -j
#PJM -X
#PJM -N fmesh_ci_local
#============================================================================
# Run tools/ci_local/run_ci_local.sh (the GitHub CI job reproduced locally: ruff, pytest without
# oceanmesh and WORK_DIR/DATA_DIR, pytest with the full environment, pre-commit, CLI smoke) on a
# compute node. Required before every commit that touches src/, tests/ or notebooks/ (CLAUDE.md).
#
#   pjsub jobs/genkai/ci_local.sh        # from the repository root
#============================================================================
set -euo pipefail
cd "${PJM_O_WORKDIR:?Submit from the repository root}"
. jobs/genkai/common.sh ci_local 1
echo "commit: $(git rev-parse HEAD)  dirty files: $(git status --porcelain --untracked-files=no | wc -l)"
tools/ci_local/run_ci_local.sh
