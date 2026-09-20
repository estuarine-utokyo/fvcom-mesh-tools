#!/bin/bash
#PBS -q OCT-S
#PBS --group=G16445
#PBS -l cpunum_job=64
#PBS -l memsz_job=192GB
#PBS -l elapstim_req=03:00:00
#PBS -N fmesh_392
#PBS -j o
#PBS -o logs/392_sizing_sweep.pbs.log
#PBS -r n
#============================================================================
# Coastal-resolution sweep at fixed inputs (owner 2026-09-20: optimise the
# present configuration before adding region targets).
#
# Diagnostic 391 showed the global time step is set by the DEEP MOUTH
# (depth 727 m, edge 1436 m) while shallow coastal elements allow 31-44 s,
# so refining the coast is free in dt terms. What actually holds the coast
# at ~400 m is the distance sizing: h = SR_H0 + SR_GRADE * distance.
#
# Variants, all else identical to the certified chain, run in parallel in
# private copies (325 writes to fixed repository-relative paths):
#   A  SR_GRADE=0.12                 gentler growth away from the coast
#   B  SR_H0=250                     finer base size
#   C  SR_H0=250 SR_GRADE=0.12       both
#
# Usage (from the repository root): qsub jobs/octopus/392_sizing_sweep.sh
#============================================================================
set -euo pipefail
cd "${PBS_O_WORKDIR:?Submit from the repository root}"
. jobs/octopus/common.sh 392_sizing_sweep 16

REPO=$(pwd)
WORK=/octfs/work/G16445/v61021/scratch/sizing_392.${JOBID}
COLLECT="$REPO/outputs/sizing_392"
mkdir -p "$WORK" "$COLLECT"

run_variant() {   # run_variant NAME "ENV=VAL ..."
    local name="$1"
    local envs="$2"
    local dir="$WORK/$name"
    mkdir -p "$dir"
    for d in notebooks recipes; do
        rsync -a --delete "$REPO/$d/" "$dir/$d/"
    done
    mkdir -p "$dir/outputs/figures"
    rsync -a "$REPO/outputs/tb_varres_3r/" "$dir/outputs/tb_varres_3r/"
    rsync -a "$REPO/outputs/sample_repro/" "$dir/outputs/sample_repro/"
    (
        cd "$dir"
        export SR_NORMALIZE=on SR_OBC_H1=1680 FMESH_OVERWRITE=1
        export OMP_NUM_THREADS=16 OPENBLAS_NUM_THREADS=16 MKL_NUM_THREADS=16 NUMBA_NUM_THREADS=16
        # shellcheck disable=SC2163
        for kv in $envs; do export "$kv"; done
        echo "=== variant $name ($envs) ==="
        python notebooks/325_sample_repro.py
        python notebooks/331_finish2.py
        fmesh-mesh-qa outputs/sample_repro/sample_repro_final.14 || true
        python notebooks/342_connectivity_check.py || true
        python notebooks/346_one_wide_flags.py || true
        python notebooks/364_over_resolution.py || true
        python notebooks/389_coverage_gaps.py || true
    ) > "$dir/run.log" 2>&1
    echo "variant $name finished rc=$? $(date -Is)"
}

echo "=== launching A B C in parallel ==="
pids=()
run_variant A "SR_GRADE=0.12" & pids+=($!)
run_variant B "SR_H0=250" & pids+=($!)
run_variant C "SR_H0=250 SR_GRADE=0.12" & pids+=($!)
status=0
for p in "${pids[@]}"; do wait "$p" || status=1; done
echo "=== variants done (status=$status) $(date -Is) ==="

for v in A B C; do
    mkdir -p "$COLLECT/$v"
    cp -p "$WORK/$v/run.log" "$COLLECT/$v/" 2>/dev/null || true
    for f in sample_repro_final.14 sample_repro_final_qa.json one_wide_cells.json \
             over_resolution.json coverage_gaps.json; do
        [ -f "$WORK/$v/outputs/sample_repro/$f" ] && cp -p "$WORK/$v/outputs/sample_repro/$f" "$COLLECT/$v/"
    done
done

echo "=== summary ==="
python notebooks/392_sizing_summary.py "$COLLECT" | tee "$COLLECT/summary.txt"
echo "end=$(date -Is) status=$status"
exit "$status"
