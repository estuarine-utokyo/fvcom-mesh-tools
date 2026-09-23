#!/bin/bash
#PBS -q OCT-S
#PBS --group=G16445
#PBS -l cpunum_job=64
#PBS -l memsz_job=120GB
#PBS -l elapstim_req=01:00:00
#PBS -N fmesh_426wallp
#PBS -j o
#PBS -o logs/426_wall_port_test.pbs.log
#PBS -r n
# Split walls in the finished port case, which is known to run (notebook 427).
# Required: FMESH_CASE (a staged, already-run 2-day case), FMESH_OUT (new dir).
set -euo pipefail
cd "${PBS_O_WORKDIR:?Submit from the repository root}"
. jobs/octopus/common.sh 426_wall_port_test 1
case $(hostname -s) in oct-cpu*) ;; *) echo 'Compute nodes only'; exit 1 ;; esac
CASE=${FMESH_CASE:?}; OUT=${FMESH_OUT:?}
python notebooks/427_wall_port_test.py prep "$CASE" "$OUT"
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
if ( cd "$OUT" && mpiexec -np 64 "$FVCOM" --casename=m2 > fvcom.log 2>&1 ); then rc=0; else rc=$?; fi
echo "[426] exit=$rc"; tail -3 "$OUT/fvcom.log"
[ "$rc" -eq 0 ] && grep -q TADA "$OUT/fvcom.log"
) || { echo "[426] the integration did not finish; no verdict"; exit 1; }
python notebooks/427_wall_port_test.py analyze "$OUT" "$CASE"
echo "end=$(date -Is)"
