#!/usr/bin/env bash
# Download one large file using many parallel HTTP range requests.
# Zenodo serves ~0.3 MB/s per connection but honours Range headers, so
# 16 connections give roughly 16x the throughput.
#
# Usage: bash scripts/parallel_fetch.sh URL OUT [N_CONN] [CHUNK_MB]
set -euo pipefail
URL="$1"; OUT="$2"; NCONN="${3:-16}"; CHUNK_MB="${4:-128}"

SIZE=$(curl -sS -I -L "$URL" | tr -d '\r' | awk 'tolower($1)=="content-length:"{s=$2} END{print s}')
[[ -n "$SIZE" && "$SIZE" -gt 0 ]] || { echo "could not get content-length" >&2; exit 1; }
CHUNK=$((CHUNK_MB * 1024 * 1024))
NPARTS=$(( (SIZE + CHUNK - 1) / CHUNK ))
PARTDIR="$OUT.parts"; mkdir -p "$PARTDIR"
echo "size=$SIZE bytes, $NPARTS parts of ${CHUNK_MB}MB, $NCONN connections"

fetch_part() {  # idx
  local i=$1 start end have part want
  start=$(( i * CHUNK )); end=$(( start + CHUNK - 1 )); (( end >= SIZE )) && end=$(( SIZE - 1 ))
  part="$PARTDIR/part_$(printf %05d "$i")"
  want=$(( end - start + 1 ))
  for attempt in $(seq 1 30); do
    have=0; [[ -f "$part" ]] && have=$(stat -c %s "$part")
    (( have >= want )) && return 0
    # append the remaining byte range; -f fails on 5xx error pages so we never append HTML
    curl -sS -f -L -r "$(( start + have ))-$end" "$URL" >> "$part" || sleep $(( attempt * 2 ))
  done
  echo "part $i failed after 30 attempts" >&2; return 1
}
export -f fetch_part; export URL SIZE CHUNK PARTDIR

seq 0 $(( NPARTS - 1 )) | xargs -P "$NCONN" -I{} bash -c 'fetch_part {}'

# verify every part has exactly the expected size, then concatenate in order
for i in $(seq 0 $(( NPARTS - 1 ))); do
  start=$(( i * CHUNK )); end=$(( start + CHUNK - 1 )); (( end >= SIZE )) && end=$(( SIZE - 1 ))
  want=$(( end - start + 1 )); part="$PARTDIR/part_$(printf %05d "$i")"
  [[ $(stat -c %s "$part") -eq $want ]] || { echo "size mismatch in $part" >&2; exit 1; }
done
cat "$PARTDIR"/part_* > "$OUT"
[[ $(stat -c %s "$OUT") -eq $SIZE ]] && rm -rf "$PARTDIR"
echo "assembled $OUT ($SIZE bytes)"
