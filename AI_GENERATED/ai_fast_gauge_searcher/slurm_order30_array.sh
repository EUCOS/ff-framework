#!/bin/bash
#SBATCH --job-name=ai-fast-gauge-o30
#SBATCH --output=AI_GENERATED/ai_fast_gauge_searcher/logs/%x-%A_%a.out
#SBATCH --error=AI_GENERATED/ai_fast_gauge_searcher/logs/%x-%A_%a.err
#SBATCH --array=0-127
#SBATCH --ntasks=1
#SBATCH --cpus-per-task=8
#SBATCH --time=24:00:00
#SBATCH --mem=8G

set -euo pipefail

ROOT="${SLURM_SUBMIT_DIR:-$(pwd)}"
APP="${ROOT}/AI_GENERATED/ai_fast_gauge_searcher/AIFastGaugeSearcher"
ORDERS="${ORDERS:-30}"
THREADS="${SLURM_CPUS_PER_TASK:-8}"
SHARD_COUNT="${SLURM_ARRAY_TASK_COUNT:-128}"
SHARD_INDEX="${SLURM_ARRAY_TASK_ID:-0}"
OUTDIR="${OUTDIR:-${ROOT}/AI_GENERATED/ai_fast_gauge_searcher/results/order${ORDERS}}"

mkdir -p "${OUTDIR}" "${ROOT}/AI_GENERATED/ai_fast_gauge_searcher/logs"

"${APP}" \
  --dimension 4 \
  --orders "${ORDERS}" \
  --output "${OUTDIR}/shard-${SHARD_INDEX}-of-${SHARD_COUNT}.txt" \
  --shard-index "${SHARD_INDEX}" \
  --shard-count "${SHARD_COUNT}" \
  --threads "${THREADS}" \
  --checkpoint "${OUTDIR}/shard-${SHARD_INDEX}-of-${SHARD_COUNT}.checkpoint" \
  --resume
