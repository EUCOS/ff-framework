#!/bin/bash
set -euo pipefail

ROOT="$(cd "$(dirname "${BASH_SOURCE[0]}")/../.." && pwd)"
OUTDIR="${ROOT}/AI_GENERATED/ai_fast_gauge_searcher"
WORKDIR="${OUTDIR}/comparison_output"
ORDER="${ORDER:-2}"

"${OUTDIR}/build.sh" >/dev/null

rm -rf "${WORKDIR}"
mkdir -p "${WORKDIR}"

"${OUTDIR}/HumanGaugeBaseline" \
  --dimension 4 \
  --output "${WORKDIR}/human_order${ORDER}.txt" \
  "${ORDER}"

"${OUTDIR}/AIFastGaugeSearcher" \
  --dimension 4 \
  --orders "${ORDER}" \
  --output "${WORKDIR}/ai_order${ORDER}.txt" \
  --shard-index 0 \
  --shard-count 1 \
  --threads "${THREADS:-1}" \
  --checkpoint "${WORKDIR}/ai_order${ORDER}.checkpoint"

normalize() {
  sed '/^Total time:/d;/^Models built per sec:/d' "$1"
}

normalize "${WORKDIR}/human_order${ORDER}.txt" > "${WORKDIR}/human.normalized"
normalize "${WORKDIR}/ai_order${ORDER}.txt" > "${WORKDIR}/ai.normalized"

if ! diff -u "${WORKDIR}/human.normalized" "${WORKDIR}/ai.normalized"; then
  echo "Comparison failed: AI output differs from human baseline for order ${ORDER}." >&2
  exit 1
fi

echo "Comparison passed for D=4 order ${ORDER}."
