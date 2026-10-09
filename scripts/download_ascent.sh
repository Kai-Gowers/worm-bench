#!/usr/bin/env bash
# Download the ASCENT in-house recording and its ground-truth CSVs from Zenodo
# (record 17561700) into a data directory outside the Git repo, and verify
# the MD5 checksums published by Zenodo.
#
# Usage:  bash scripts/download_ascent.sh [DEST_DIR]
# Default DEST_DIR is data/ascent (a symlink to scratch in this repo).
set -euo pipefail

DEST="${1:-data/ascent}"
REC="https://zenodo.org/api/records/17561700/files"
mkdir -p "$DEST"

# name  md5 (from https://zenodo.org/api/records/17561700)
FILES=(
  "InHouse_Opterra_coords_GT.csv 37bd1234033d360dfd888b46daa369b1"
  "InHouse_Opterra_tracks_GT.csv 5a7c414b95365bb13cc2674df03d5277"
  "InHouse_Opterra_CH_RG_Day1_OP50Stim_1100frames.h5 9581892c969bb0e7f3d84eba84988133"
)

for entry in "${FILES[@]}"; do
  name="${entry%% *}"; md5="${entry##* }"
  out="$DEST/$name"
  if [[ -f "$out" ]] && echo "$md5  $out" | md5sum -c --quiet 2>/dev/null; then
    echo "ok (cached): $name"; continue
  fi
  echo "downloading: $name"
  # -f fails on HTTP errors (Zenodo sometimes returns 504 pages), -C - resumes,
  # --retry-all-errors retries on those failures.
  if [[ "$name" == *.h5 ]]; then
    # the 6 GB movie: Zenodo gives ~0.3 MB/s per connection, so fetch in parallel ranges
    bash "$(dirname "$0")/parallel_fetch.sh" "$REC/$name/content" "$out" 16 128
  else
    curl -sS -f -L -C - --retry 20 --retry-delay 15 --retry-all-errors \
         -o "$out" "$REC/$name/content"
  fi
  echo "$md5  $out" | md5sum -c
done
echo "done: $(du -sh "$DEST" | cut -f1) in $DEST"
