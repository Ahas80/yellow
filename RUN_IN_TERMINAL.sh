#!/bin/bash

# Terminal wrapper for the only supported complete analysis entrypoint.
set -euo pipefail

cd "$(dirname "$0")"
mkdir -p logs

run_stamp="$(date '+%Y%m%dT%H%M%S%z')"
export RUN_ID="${RUN_ID:-${run_stamp}}"
export PIPELINE_LOG_PATH="${PWD}/logs/complete_analysis_${RUN_ID}.log"

echo "Writing the complete run log to ${PIPELINE_LOG_PATH}"
# Reuse a core-SNP result only when its module confirms that every required
# output and the exact SHA-256 input manifest still match.
FORCE_RERUN_CORE_SNP="${FORCE_RERUN_CORE_SNP:-0}" \
RUN_LEGACY_PUBLICATION_FIGURES="${RUN_LEGACY_PUBLICATION_FIGURES:-0}" \
    bash RUN_COMPLETE_ANALYSIS.sh 2>&1 | tee "${PIPELINE_LOG_PATH}"
