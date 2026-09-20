#!/bin/bash
#PBS -q OCT-S
#PBS --group=G16445
#PBS -l cpunum_job=8
#PBS -l memsz_job=64GB
#PBS -l elapstim_req=02:00:00
#PBS -N fmesh_385
#PBS -j o
#PBS -o logs/385_no_edits.pbs.log
#PBS -r n
#============================================================================
# How far does the chain get WITHOUT the four human-judgment edits?
#
# Regenerates the sample-repro mesh with SR_EDITS_EXCLUDE set to all of
# recipes/edits/sample_repro/, then runs the same acceptance measurements as
# the certified run (QA gate, connectivity comparator, one-wide ledger,
# over-resolution/wall integrity) and compares the two meshes.
#
# The certified outputs are moved aside and restored afterwards; the
# edit-free results are kept in outputs/sample_repro_noedits/.
#
# Usage (from the repository root): qsub jobs/octopus/385_no_edits.sh
#============================================================================
set -euo pipefail
cd "${PBS_O_WORKDIR:?Submit from the repository root}"
. jobs/octopus/common.sh 385_no_edits 4

SR=outputs/sample_repro
CERT="$SR/cert_octopus_115144"
NOED=outputs/sample_repro_noedits
CERT_FILES="sample_repro.14 sample_repro_utm.14 sample_repro_final.14
            sample_repro_qa.json sample_repro_utm_qa.json sample_repro_final_qa.json"

mkdir -p "$CERT" "$NOED"
for f in $CERT_FILES; do
    [ -f "$SR/$f" ] && cp -p "$SR/$f" "$CERT/$f"
done
restore() {
    for f in $CERT_FILES; do
        [ -f "$CERT/$f" ] && cp -p "$CERT/$f" "$SR/$f"
    done
    echo "certified outputs restored into $SR"
}
trap restore EXIT

export SR_NORMALIZE=on
export SR_OBC_H1=1680
export SR_EDITS_EXCLUDE=edit_001_haneda_d_runway,edit_003_west_edge_crack,edit_004_ow05_harbor,edit_005_ow05_urayasu

echo "=== generate WITHOUT the four manual edits ==="
python notebooks/325_sample_repro.py
echo "=== finish ==="
python notebooks/331_finish2.py
echo "=== QA (edit-free mesh) ==="
fmesh-mesh-qa "$SR/sample_repro_final.14" || true
echo "=== connectivity comparator ==="
python notebooks/342_connectivity_check.py || true
echo "=== one-wide ledger ==="
python notebooks/346_one_wide_flags.py || true
echo "=== over-resolution / wall integrity ==="
python notebooks/364_over_resolution.py || true

echo "=== keep the edit-free results ==="
for f in $CERT_FILES; do
    [ -f "$SR/$f" ] && cp -p "$SR/$f" "$NOED/$f"
done
for f in issue_registry.json ow_registry.json one_wide_cells.json normalize.json axis_check.json; do
    [ -f "$SR/$f" ] && cp -p "$SR/$f" "$NOED/$f"
done

echo "=== certified vs edit-free ==="
fmesh-mesh-quality "$CERT/sample_repro_final.14" "$NOED/sample_repro_final.14" || true
echo "=== QA of the certified mesh, for the same report ==="
fmesh-mesh-qa "$CERT/sample_repro_final.14" || true
echo "end=$(date -Is)"
