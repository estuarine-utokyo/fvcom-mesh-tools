#!/bin/bash
# Phase 3 of the OceanMesh2D -> $DATA_DIR migration: DESTRUCTIVE.
#
# Replaces source paths under $HOME/Github/OceanMesh2D with symlinks that
# point into $DATA_DIR. After this runs, the OceanMesh2D MATLAB scripts
# continue to work transparently because every original path resolves via
# a symlink to the canonical copy in $DATA_DIR.
#
# PREREQUISITES:
#   - tools/migrate_oceanmesh2d_data.pjsub has completed successfully and
#     all "[verify] OK" lines are present in its log.
#
# Invocation:
#   bash tools/replace_oceanmesh2d_with_symlinks.sh           # dry-run (default)
#   bash tools/replace_oceanmesh2d_with_symlinks.sh --apply   # actually do it
#
# This is fast (just rm + ln) so it does not need to be a batch job. Run on
# the login node.

set -euo pipefail

: "${HOME:?HOME not set}"
: "${DATA_DIR:?DATA_DIR not set}"

OM2D="${HOME}/Github/OceanMesh2D"
[[ -d "${OM2D}" ]] || { echo "ERROR: ${OM2D} not found" >&2; exit 1; }

APPLY=0
case "${1-}" in
  --apply) APPLY=1 ;;
  ""|--dry-run) APPLY=0 ;;
  -h|--help)
    sed -n '2,18p' "$0"; exit 0 ;;
  *) echo "unknown arg: $1" >&2; exit 2 ;;
esac

if (( APPLY )); then
  echo "[symlink] APPLY mode — will rm and ln -s the originals."
else
  echo "[symlink] DRY RUN — pass --apply to perform the changes."
fi
echo

# Each entry: <original path under OM2D>  <target under DATA_DIR>
mapping=(
  "${OM2D}/datasets/GEBCO_2024.nc                    ${DATA_DIR}/bathymetry/GEBCO/GEBCO_2024.nc"
  "${OM2D}/datasets/SRTM15+.nc                       ${DATA_DIR}/bathymetry/SRTM15plus/SRTM15+.nc"
  "${OM2D}/datasets/CUDEMS                           ${DATA_DIR}/bathymetry/CUDEM"
  "${OM2D}/datasets/TokyoBay/dem                     ${DATA_DIR}/bathymetry/tokyo_bay"
  "${OM2D}/datasets/GSHHS_shp                        ${DATA_DIR}/coastline/GSHHS"
  "${OM2D}/datasets/TokyoBay/shp/Futtsu_coastline    ${DATA_DIR}/coastline/tokyo_bay/Futtsu"
  "${OM2D}/datasets/TokyoBay/shp/SMS                 ${DATA_DIR}/coastline/tokyo_bay/SMS"
  "${OM2D}/datasets/TokyoBay/shp/Wang_coastline      ${DATA_DIR}/coastline/tokyo_bay/Wang"
  "${OM2D}/datasets/TokyoBay/shp/tokyo.nc            ${DATA_DIR}/coastline/tokyo_bay/misc/tokyo.nc"
  "${OM2D}/datasets/TokyoBay/shp/tokyo_0001.nc       ${DATA_DIR}/coastline/tokyo_bay/misc/tokyo_0001.nc"
  "${OM2D}/Tokyo_Bay/data/Futtsu_coastline           ${DATA_DIR}/coastline/tokyo_bay/Futtsu"
  "${OM2D}/tb_mesh/tb_futtsu.14                      ${DATA_DIR}/mesh/reference/tokyo_bay/tb_futtsu.14"
  "${OM2D}/tb_mesh/tb_futtsu.mat                     ${DATA_DIR}/mesh/reference/tokyo_bay/tb_futtsu.mat"
  "${OM2D}/tb_mesh/tb_futtsu20220311.14              ${DATA_DIR}/mesh/reference/tokyo_bay/tb_futtsu20220311.14"
  "${OM2D}/tb_mesh/tokyobay_futtsu5.0                ${DATA_DIR}/mesh/reference/tokyo_bay/tokyobay_futtsu5.0"
  "${OM2D}/tb_mesh/tokyobay_futtsu5.0.14             ${DATA_DIR}/mesh/reference/tokyo_bay/tokyobay_futtsu5.0.14"
)

# MLIT C23 shapefiles: replace each component file individually because
# the source layout is flat (multiple shapefile stems share Tokyo_Bay/data/).
declare -a mlit_files=()
for f in "${OM2D}/Tokyo_Bay/data/"C23-06_TOKYOBAY*; do
  [[ -e "${f}" ]] || continue
  base="${f##*/}"
  mlit_files+=("${f}    ${DATA_DIR}/coastline/tokyo_bay/MLIT_C23/${base}")
done

run_one () {
  local original="$1" target="$2"
  if [[ ! -e "${target}" ]]; then
    echo "[symlink] SKIP (target missing): ${target}"
    return
  fi
  if [[ -L "${original}" ]]; then
    echo "[symlink] already symlink: ${original} -> $(readlink "${original}")"
    return
  fi
  if [[ ! -e "${original}" ]]; then
    echo "[symlink] SKIP (original missing): ${original}"
    return
  fi
  echo "[symlink] ${original}"
  echo "       -> ${target}"
  if (( APPLY )); then
    rm -rf "${original}"
    ln -s "${target}" "${original}"
  fi
}

for line in "${mapping[@]}"; do
  read -r original target <<< "${line}"
  run_one "${original}" "${target}"
done

echo
echo "[symlink] === MLIT C23 (per-file) ==="
for line in "${mlit_files[@]}"; do
  read -r original target <<< "${line}"
  run_one "${original}" "${target}"
done

echo
if (( APPLY )); then
  echo "[symlink] DONE. OceanMesh2D paths now resolve via $DATA_DIR."
  echo "          Verify with: ls -la ${OM2D}/datasets ${OM2D}/Tokyo_Bay/data ${OM2D}/tb_mesh"
else
  echo "[symlink] DRY RUN done. Re-run with --apply to perform the replacements."
fi
