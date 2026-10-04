#!/bin/bash
#PJM -L rscgrp=a-batch
#PJM -L vnode=1
#PJM -L vnode-core=8
#PJM -L elapse=01:00:00
#PJM -j
#PJM -X
#PJM -N fmesh_build
#============================================================================
# GENKAI: compile and install the local repositories into the conda env
# (the GENKAI twin of jobs/octopus/build_env.sh).
#
# Create the env first on the LOGIN node, from conda-forge only:
#   mamba env create -n fvcom-mesh-tools -f environment.yml
# Then, from the repository root:
#   pjsub jobs/genkai/build_env.sh
# Log: logs/build_env.<jobid>.log
#============================================================================
set -euo pipefail
cd "${PJM_O_WORKDIR:?Submit from the repository root}"
. jobs/genkai/common.sh build_env 8

REPO="$(pwd)"
GH="$(dirname "${REPO}")"

for d in "${GH}/oceanmesh" "${REPO}" "${GH}/xcoast"; do
    echo "=== pip install -e ${d} ==="
    (cd "${d}" && python -m pip install -e . --no-deps --no-build-isolation -v 2>&1 \
        | grep -vE '^\s*(Created temporary|Removed build tracker)' | tail -40)
done

echo "=== import smoke test ==="
python - <<'PY'
import importlib
for m in ["oceanmesh", "fvcom_mesh_tools", "xcoast",
          "geopandas", "rasterio", "netCDF4", "numba", "skfmm", "yaml"]:
    mod = importlib.import_module(m)
    print(f"{m:18s} {getattr(mod, '__version__', '-'):12s} {mod.__file__}")
PY

echo "=== pytest ==="
python -m pytest -q -x -p no:cacheprovider tests | tail -15
echo "end=$(date -Is)"
