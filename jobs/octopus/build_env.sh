#!/bin/bash
#PBS -q OCT-S
#PBS --group=G16445
#PBS -l cpunum_job=8
#PBS -l memsz_job=32GB
#PBS -l elapstim_req=01:00:00
#PBS -N fmesh_build
#PBS -j o
#PBS -o logs/build_env.pbs.log
#PBS -r n
#============================================================================
# OCTOPUS: compile and install the local repositories into the conda env.
#
# The env itself is created on the LOGIN node (compute nodes have no
# network), from conda-forge only:
#   mamba env create -n oceanmesh-bench -f environment.yml
# This job then builds, without network and without PyPI:
#   ../oceanmesh   (our fork; C++/CGAL extensions compiled in place)
#   .              (fvcom-mesh-tools)
#   ../xcoast      (coastline plotting helper)
# with `pip install -e . --no-deps --no-build-isolation`, so pip only
# runs the local build; every dependency comes from conda-forge.
# Finally it imports everything and runs the unit tests.
#
# Usage (from the repository root):
#   qsub jobs/octopus/build_env.sh
# Log: logs/build_env.<jobid>.log
#============================================================================
set -euo pipefail
cd "${PBS_O_WORKDIR:-$PWD}"
. jobs/octopus/common.sh build_env 8

REPO="$(pwd)"
GH="$(dirname "${REPO}")"

for d in "${GH}/oceanmesh" "${REPO}" "${GH}/xcoast"; do
    echo "=== pip install -e ${d} ==="
    (cd "${d}" && python -m pip install -e . --no-deps --no-build-isolation -v 2>&1 \
        | grep -vE '^\s*(Created temporary|Removed build tracker)' | tail -40)
done

echo "=== extension linkage (expect the env's libgmp/libmpfr) ==="
for so in "${GH}"/oceanmesh/_*.so; do
    echo "--- $(basename "${so}")"
    ldd "${so}" | grep -E 'gmp|mpfr|stdc\+\+' || true
done

echo "=== import smoke test ==="
python - <<'PY'
import importlib
for m in ["oceanmesh", "fvcom_mesh_tools", "xcoast", "ocsmesh", "jigsawpy",
          "geopandas", "rasterio", "netCDF4", "numba", "skfmm", "yaml"]:
    mod = importlib.import_module(m)
    print(f"{m:18s} {getattr(mod, '__version__', '-'):12s} {mod.__file__}")
PY

echo "=== CLI entry points ==="
fmesh-mesh-qa --help | head -3

echo "=== pytest ==="
python -m pytest -q -x -p no:cacheprovider tests | tail -15

echo "end=$(date -Is)"
