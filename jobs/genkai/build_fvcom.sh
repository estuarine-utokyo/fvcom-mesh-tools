#!/bin/bash
#PJM -L rscgrp=a-batch
#PJM -L vnode=1
#PJM -L vnode-core=8
#PJM -L elapse=00:30:00
#PJM -j
#PJM -X
#PJM -N fmesh_fvcom_build
#============================================================================
# GENKAI: rebuild the FVCOM executable that the mesh-test runs use
# ($WORK_DIR/Github/FVCOM/src/fvcom, the default of fmesh_fvcom), from the
# source and the local src/make.inc as they are. Only step 2 of
# FVCOM/genkai_rebuild_all.sh: no FABM rebuild (SKIP_FABM), and the binary is
# not copied into TB-FVCOM run directories.
#
#   pjsub jobs/genkai/build_fvcom.sh
# Log: logs/build_fvcom.<jobid>.log
# The previous binary and library are kept as fvcom.<date>.prev / libfvcom.a.<date>.prev.
#============================================================================
set -euo pipefail
cd "${PJM_O_WORKDIR:?Submit from the repository root}"
: "${WORK_DIR:?WORK_DIR is not set (login profile)}"
mkdir -p logs
exec > "logs/build_fvcom.${PJM_JOBID:-interactive.$$}.log" 2>&1

module load intel/2025.1.3 impi/2021.15 netcdf/4.9.2 netcdf-fortran/4.6.1 hdf5/1.14.4

SRC="$WORK_DIR/Github/FVCOM/src"
cd "$SRC"
echo "FVCOM $(git -C .. rev-parse --short HEAD) branch $(git -C .. branch --show-current); make.inc sha $(sha256sum make.inc | cut -c1-12)"
stamp=$(date +%Y%m%d)
[ -f fvcom ] && cp -p fvcom "fvcom.${stamp}.prev"
[ -f libfvcom.a ] && cp -p libfvcom.a "libfvcom.a.${stamp}.prev"
make clean
make -j 8 SYSTEM_TYPE=GENKAI
ls -la fvcom libfvcom.a
echo "built $(date -Is)"
