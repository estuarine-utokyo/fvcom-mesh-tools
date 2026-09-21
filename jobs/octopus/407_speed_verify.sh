#!/bin/bash
#PBS -q OCT-S
#PBS --group=G16445
#PBS -l cpunum_job=16
#PBS -l memsz_job=64GB
#PBS -l elapstim_req=02:00:00
#PBS -N fmesh_407
#PBS -j o
#PBS -o logs/407_speed_verify.pbs.log
#PBS -r n
# Verify the detect_waterways speed-up: the eps land buffer is loop-invariant
# and used to be rebuilt once per network (1.9 s each on the 168k-vertex
# coastline).  Job 115299 measured the input stage at 1,427 s of 1,487 s wall,
# with 1,185 s of it inside shapely buffer.  This run must produce a mesh
# IDENTICAL to the reference and report a much shorter input stage.
# Submit from the repository root:
#   qsub jobs/octopus/407_speed_verify.sh
set -euo pipefail
cd "${PBS_O_WORKDIR:?Submit from the repository root}"
. jobs/octopus/common.sh 407_speed_verify 16

REPO=$(pwd)
REF=${FMESH_REF:?Set FMESH_REF to the reference sample_repro_final.14}
WORK=/octfs/work/G16445/v61021/scratch/speed_407.${JOBID}
mkdir -p "$WORK"
for d in notebooks recipes; do rsync -a --delete "$REPO/$d/" "$WORK/$d/"; done
mkdir -p "$WORK/outputs/figures"
rsync -a "$REPO/outputs/tb_varres_3r/" "$WORK/outputs/tb_varres_3r/"
rsync -a "$REPO/outputs/sample_repro/" "$WORK/outputs/sample_repro/"
rm -f "$WORK/outputs/sample_repro/sample_repro_final.14"

cd "$WORK"
export SR_NORMALIZE=on SR_OBC_H1=1680 FMESH_OVERWRITE=1
export SR_H_TARGET=achieved SR_OBC_SKIP=2 SR_ONE_WIDE=forbid SR_COAST_FIT=on
export OMP_NUM_THREADS=16 OPENBLAS_NUM_THREADS=16 MKL_NUM_THREADS=16 NUMBA_NUM_THREADS=16
t0=$(date +%s)
python notebooks/325_sample_repro.py > "$WORK/325.log" 2>&1
t1=$(date +%s)
python notebooks/331_finish2.py > "$WORK/331.log" 2>&1
t2=$(date +%s)
grep -E "inputs \+|sizing done" "$WORK/325.log" || true
grep -E "coast fit" "$WORK/331.log" || true
echo "325 wall = $((t1 - t0)) s   331 wall = $((t2 - t1)) s"
fmesh-mesh-qa outputs/sample_repro/sample_repro_final.14 | tail -4 || true
echo "=== identity check against $REF ==="
if cmp -s outputs/sample_repro/sample_repro_final.14 "$REF"; then
    echo "IDENTICAL to the reference"
else
    echo "DIFFERS from the reference:"
    python - "$REF" outputs/sample_repro/sample_repro_final.14 <<'PY'
import sys
import numpy as np
def load(p):
    L = p.read_text().splitlines() if hasattr(p, "read_text") else open(p).read().splitlines()
    ne, nn = map(int, L[1].split()[:2])
    xy = np.array([L[2 + i].split()[1:3] for i in range(nn)], float)
    tri = np.array([L[2 + nn + i].split()[2:5] for i in range(ne)], int)
    return xy, tri
a, ta = load(sys.argv[1])
b, tb = load(sys.argv[2])
print(f"  nodes {len(a)} vs {len(b)}, elements {len(ta)} vs {len(tb)}")
if a.shape == b.shape:
    d = np.linalg.norm(a - b, axis=1)
    print(f"  node displacement: max {d.max():.6f} m, >1 m: {(d > 1).sum()}")
PY
fi
echo "end=$(date -Is)"
