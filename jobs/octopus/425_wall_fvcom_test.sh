#!/bin/bash
#PBS -q OCT-S
#PBS --group=G16445
#PBS -l cpunum_job=8
#PBS -l memsz_job=32GB
#PBS -l elapstim_req=01:00:00
#PBS -N fmesh_425wall
#PBS -j o
#PBS -o logs/425_wall_fvcom_test.pbs.log
#PBS -r n
# A split wall against an unsplit one, in FVCOM (notebook 426).
# Required: FMESH_RUN_ROOT (a fresh directory).
set -euo pipefail
cd "${PBS_O_WORKDIR:?Submit from the repository root}"
. jobs/octopus/common.sh 425_wall_fvcom_test 8
case $(hostname -s) in oct-cpu*) ;; *) echo 'Compute nodes only'; exit 1 ;; esac
ROOT=${FMESH_RUN_ROOT:?set FMESH_RUN_ROOT}
[ -e "$ROOT" ] && { echo "exists: $ROOT"; exit 2; }
python notebooks/426_wall_fvcom_test.py prep "$ROOT"
FVCOM=/octfs/work/G16445/v61021/Github/FVCOM/src/fvcom
(
set +u; conda deactivate; set -u
if ! type module >/dev/null 2>&1; then
    for init in /etc/profile.d/modules.sh /usr/share/Modules/init/bash /usr/share/lmod/lmod/init/bash; do
        if [[ -r $init ]]; then source "$init"; break; fi
    done
fi
module purge
module load BaseCPU/2026
module load hdf5/1.14.6 netcdf-c/4.9.3 netcdf-fortran/4.6.2
INSTALLDIR=/octfs/work/G16445/share/local/fvcom/libs/install-oneapi-2025.3.1
export INSTALLDIR OMP_NUM_THREADS=1 PROJ_DATA=/usr/share/proj
export LD_LIBRARY_PATH="$INSTALLDIR/lib:$INSTALLDIR/lib64:${LD_LIBRARY_PATH:-}"
ulimit -s unlimited
for case in split edges; do
    if ( cd "$ROOT/$case" && mpiexec -np 8 "$FVCOM" --casename=m2 > fvcom.log 2>&1 ); then
        rc=0; else rc=$?; fi
    echo "[425] $case exit=$rc"
    grep -Eci 'fatal|non[ -]?finite|floating exception|segmentation|nan detected' \
        "$ROOT/$case/fvcom.log" || true
    tail -3 "$ROOT/$case/fvcom.log"
done
)
python notebooks/426_wall_fvcom_test.py analyze "$ROOT"
echo "end=$(date -Is)"
