#!/bin/bash
#PJM -L rscgrp=a-batch
#PJM -L vnode=1
#PJM -L vnode-core=4
#PJM -L elapse=00:15:00
#PJM -j
#PJM -X
#PJM -N fmesh_455
#============================================================================
# GENKAI: read-only analysis of how the extended mesh follows the coastline.
#   pjsub -x BUILD=<445 build dir>,OUT=<new analysis dir> jobs/genkai/455_coast_analysis.sh
#============================================================================
set -euo pipefail
cd "${PJM_O_WORKDIR:?Submit from the repository root}"
: "${BUILD:?set BUILD}" "${OUT:?set OUT}"
. jobs/genkai/common.sh 455_coast_analysis 4
python notebooks/455_coast_fit_analysis.py "$BUILD" "$OUT"
