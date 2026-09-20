#!/bin/bash
#PBS -q OCT-S
#PBS --group=G16445
#PBS -l cpunum_job=32
#PBS -l memsz_job=96GB
#PBS -l elapstim_req=02:00:00
#PBS -N fmesh_388
#PBS -j o
#PBS -o logs/388_cert_regen.pbs.log
#PBS -r n
#============================================================================
# Regenerate the CERTIFIED configuration (all four edits ON) in a private
# copy, re-measure it and redraw the issue map from a self-consistent set.
#
# Needed because job 115210 (all-edits-off ablation) left the auxiliary
# registries -- land_channel_adj.shp, widen_ops.json, one_wide_cells.json --
# describing the edit-free mesh while the .14 files were restored, which
# mis-classifies intended vs unintended land elements and wall crossings.
#
# The repository working tree is NOT modified except for the collected
# results under outputs/cert_regen/ and the two issue figures.
#
# Usage (from the repository root): qsub jobs/octopus/388_cert_regen.sh
#============================================================================
set -euo pipefail
cd "${PBS_O_WORKDIR:?Submit from the repository root}"
. jobs/octopus/common.sh 388_cert_regen 32

REPO=$(pwd)
WORK=/octfs/work/G16445/v61021/scratch/cert_regen.${JOBID}
COLLECT="$REPO/outputs/cert_regen"
mkdir -p "$WORK" "$COLLECT"

for d in notebooks recipes; do
    rsync -a --delete "$REPO/$d/" "$WORK/$d/"
done
mkdir -p "$WORK/outputs/figures"
rsync -a "$REPO/outputs/tb_varres_3r/" "$WORK/outputs/tb_varres_3r/"
rsync -a "$REPO/outputs/sample_repro/" "$WORK/outputs/sample_repro/"

cd "$WORK"
export SR_NORMALIZE=on SR_OBC_H1=1680
echo "=== generate (certified configuration: all four edits ON) ==="
python notebooks/325_sample_repro.py
echo "=== finish ==="
python notebooks/331_finish2.py
echo "=== QA ==="
fmesh-mesh-qa outputs/sample_repro/sample_repro_final.14 || true
echo "=== connectivity comparator ==="
python notebooks/342_connectivity_check.py || true
echo "=== one-wide ledger ==="
python notebooks/346_one_wide_flags.py || true
echo "=== over-resolution / wall integrity ==="
python notebooks/364_over_resolution.py || true
echo "=== OW05 axis check ==="
python notebooks/372_axis_check.py || true
echo "=== issue map (self-consistent) ==="
python notebooks/386_issue_map.py

echo "=== collect ==="
rsync -a "$WORK/outputs/sample_repro/" "$COLLECT/sample_repro/"
cp -p "$WORK/outputs/figures/386_issue_map.png" "$REPO/outputs/figures/386_issue_map.png"
cp -p "$WORK/outputs/figures/386_issue_zooms.png" "$REPO/outputs/figures/386_issue_zooms.png"
echo "=== quality: certified 115144 vs this regeneration ==="
fmesh-mesh-quality "$REPO/outputs/sample_repro/cert_octopus_115144/sample_repro_final.14" \
    "$WORK/outputs/sample_repro/sample_repro_final.14" || true
echo "end=$(date -Is)"
