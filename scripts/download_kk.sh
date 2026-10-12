#!/usr/bin/env bash
# Download one NWB recording from the KK dataset (DANDI 000692, Kimura lab)
# into a data directory outside the Git repo and verify the SHA-256 digest
# published in the DANDI asset metadata.
#
# Usage:  bash scripts/download_kk.sh [SUBJECT] [DEST_DIR]
#   SUBJECT   e.g. sub-231013-02 (default)
#   DEST_DIR  default data/kk (a symlink to scratch in this repo)
#
# Files are ~9 GB each.  DANDI serves ~50 MB/s on a single connection, so a plain resumable curl is
# enough (no parallel range fetch needed, unlike Zenodo).
set -euo pipefail

SUBJECT="${1:-sub-231013-02}"
DEST="${2:-data/kk}"
DANDISET=000692
VERSION=0.240402.2118
API="https://api.dandiarchive.org/api"
mkdir -p "$DEST"

# look up the asset by path prefix; the API returns asset_id, size, path
listing=$(curl -sS -f "$API/dandisets/$DANDISET/versions/$VERSION/assets/?path=$SUBJECT/&page_size=10")
asset_id=$(echo "$listing" | python3 -P -c 'import json,sys; r=json.load(sys.stdin)["results"]; assert len(r)==1, r; print(r[0]["asset_id"])')
path=$(echo "$listing" | python3 -P -c 'import json,sys; print(json.load(sys.stdin)["results"][0]["path"])')
name=$(basename "$path")
meta=$(curl -sS -f "$API/dandisets/$DANDISET/versions/$VERSION/assets/$asset_id/")
sha=$(echo "$meta" | python3 -P -c 'import json,sys; print(json.load(sys.stdin)["digest"]["dandi:sha2-256"])')
size=$(echo "$meta" | python3 -P -c 'import json,sys; print(json.load(sys.stdin)["contentSize"])')
echo "asset $asset_id  $path  $size bytes  sha256 $sha"
echo "$meta" > "$DEST/$name.asset.json"

out="$DEST/$name"
if [[ -f "$out" ]] && echo "$sha  $out" | sha256sum -c --quiet 2>/dev/null; then
  echo "ok (cached): $name"; exit 0
fi
curl -sS -f -L -C - --retry 10 --retry-delay 10 --retry-all-errors \
     -o "$out" "$API/assets/$asset_id/download/"
echo "$sha  $out" | sha256sum -c
echo "done: $(du -sh "$out" | cut -f1) at $out"
