#!/usr/bin/env bash
# Download one NWB recording from the Kato-lab dandisets used by WormID-Bench
# (SK1 = DANDI 000565, SK2 = DANDI 000472) into a data directory outside the Git
# repo and verify the SHA-256 digest published in the DANDI asset metadata.
#
# Usage:  bash scripts/download_sk.sh DANDISET SUBJECT [DEST_DIR]
#   DANDISET  000565 (SK1) or 000472 (SK2)
#   SUBJECT   e.g. sub-20230322-18-57-04 or sub-2022-04-26-w00-NP1
#   DEST_DIR  default data/sk1 or data/sk2 (symlinks to scratch in this repo)
#
# Only 4 of the 21 SK1 assets and 2 of the 10 SK2 assets contain a calcium
# movie (see scripts/survey_dandi_sk.py); the rest are NeuroPAL-only stacks.
set -euo pipefail

DANDISET="${1:?dandiset id (000565 or 000472)}"
SUBJECT="${2:?subject, e.g. sub-20230322-18-57-04}"
case "$DANDISET" in
  000565) VERSION=0.241009.1504; DEFAULT_DEST=data/sk1 ;;
  000472) VERSION=0.241009.1502; DEFAULT_DEST=data/sk2 ;;
  *) echo "unknown dandiset $DANDISET" >&2; exit 1 ;;
esac
DEST="${3:-$DEFAULT_DEST}"
API="https://api.dandiarchive.org/api"
mkdir -p "$DEST"

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
