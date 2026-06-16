#!/bin/bash
set -euo pipefail

ROOT="$(cd "$(dirname "${BASH_SOURCE[0]}")/../.." && pwd)"
OUTDIR="${ROOT}/AI_GENERATED/ai_fast_gauge_searcher"
CXX="${CXX:-g++}"
CXXFLAGS="${CXXFLAGS:--std=c++11 -O3 -DNDEBUG -Wall -pedantic}"

mkdir -p "${OUTDIR}/build"

needs_build() {
  local output="$1"
  shift
  [[ ! -x "${output}" ]] && return 0
  local input
  for input in "$@"; do
    [[ "${input}" -nt "${output}" ]] && return 0
  done
  return 1
}

FRAMEWORK_SOURCES=("${ROOT}"/src/*.cpp)

if needs_build "${OUTDIR}/AIFastGaugeSearcher" \
    "${OUTDIR}/AIFastGaugeSearcher.cpp" "${FRAMEWORK_SOURCES[@]}"; then
  "${CXX}" ${CXXFLAGS} -I"${ROOT}/include" -pthread \
    -o "${OUTDIR}/AIFastGaugeSearcher" \
    "${OUTDIR}/AIFastGaugeSearcher.cpp" \
    "${FRAMEWORK_SOURCES[@]}"
fi

if needs_build "${OUTDIR}/HumanGaugeBaseline" \
    "${OUTDIR}/HumanGaugeBaseline.cpp" "${FRAMEWORK_SOURCES[@]}"; then
  "${CXX}" ${CXXFLAGS} -I"${ROOT}/include" -pthread \
    -o "${OUTDIR}/HumanGaugeBaseline" \
    "${OUTDIR}/HumanGaugeBaseline.cpp" \
    "${FRAMEWORK_SOURCES[@]}"
fi

echo "${OUTDIR}/AIFastGaugeSearcher"
