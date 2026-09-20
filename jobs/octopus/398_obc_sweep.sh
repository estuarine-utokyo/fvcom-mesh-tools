#!/bin/bash
#PBS -q OCT-S
#PBS --group=G16445
#PBS -l cpunum_job=48
#PBS -l memsz_job=144GB
#PBS -l elapstim_req=03:00:00
#PBS -N fmesh_398
#PBS -j o
#PBS -o logs/398_obc_sweep.pbs.log
#PBS -r n
#============================================================================
# The size targets as ACHIEVED edge lengths (owner 2026-09-20, default from
# now on) versus as SIZING-FIELD values (the certified chain), and with
# one-element-wide channels allowed.
#
#   field     SR_H_TARGET=field                     = certified behaviour
#   achieved  SR_H_TARGET=achieved                  new default
#   ach1w     SR_H_TARGET=achieved SR_ONE_WIDE=allow
#
# Usage (from the repository root): qsub jobs/octopus/395_target_mode.sh
#============================================================================
set -euo pipefail
cd "${PBS_O_WORKDIR:?Submit from the repository root}"
. jobs/octopus/common.sh 398_obc_sweep 16

REPO=$(pwd)
WORK=/octfs/work/G16445/v61021/scratch/target_398.${JOBID}
COLLECT="$REPO/outputs/target_398"
mkdir -p "$WORK" "$COLLECT"

run_variant() {
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
        python notebooks/394_spec_match.py outputs/sample_repro/sample_repro_final.14 "$name" || true
    ) > "$dir/run.log" 2>&1
    echo "variant $name finished rc=$? $(date -Is)"
}

# Close the remaining QA failures of the achieved-target mesh
# (OBC orthogonality 45.9 deg at one node, C1 x3, C4 x1).
pids=()
# Every QA failure of the achieved mesh sits at the NW OBC junction,
# so vary the band geometry there.
run_variant k100   "SR_H_TARGET=achieved SR_OBC_K=1.0" & pids+=($!)
run_variant k150   "SR_H_TARGET=achieved SR_OBC_K=1.5" & pids+=($!)
run_variant hs120  "SR_H_TARGET=achieved SR_OBC_HSCALE=1.2" & pids+=($!)
run_variant h0_500 "SR_H_TARGET=achieved SR_OBC_H0=500" & pids+=($!)
status=0
for p in "${pids[@]}"; do wait "$p" || status=1; done
echo "=== variants done (status=$status) $(date -Is) ==="

for v in k100 k150 hs120 h0_500; do
    mkdir -p "$COLLECT/$v"
    cp -p "$WORK/$v/run.log" "$COLLECT/$v/" 2>/dev/null || true
    for f in sample_repro_final.14 sample_repro_final_qa.json one_wide_cells.json \
             over_resolution.json coverage_gaps.json; do
        [ -f "$WORK/$v/outputs/sample_repro/$f" ] && cp -p "$WORK/$v/outputs/sample_repro/$f" "$COLLECT/$v/"
    done
done
echo "=== summary ==="
python notebooks/392_sizing_summary.py "$COLLECT" | tee "$COLLECT/summary.txt"
python notebooks/393_sizing_compare.py "$COLLECT" || true
echo "end=$(date -Is) status=$status"
exit "$status"
