#!/bin/bash

set -uo pipefail

cd "$(dirname "$0")"
mkdir -p logs results/pipeline

stamp="$(date '+%Y%m%dT%H%M%S%z')"
log_file="$(pwd)/logs/complete_analysis_${stamp}.log"
pid_file="results/pipeline/runner.pid"
log_pointer="results/pipeline/runner.logpath"
exit_file="results/pipeline/runner.exit"

echo "$$" > "${pid_file}"
echo "${log_file}" > "${log_pointer}"
rm -f "${exit_file}"

# A detached continuation must reuse a manifest-matching core-SNP result when
# available.  The core-SNP module still reruns automatically if its required
# outputs or exact input hash are absent/stale.
FORCE_RERUN_CORE_SNP="${FORCE_RERUN_CORE_SNP:-0}" \
    bash RUN_COMPLETE_ANALYSIS.sh >> "${log_file}" 2>&1
status=$?
{
    echo "exit_status=${status}"
    echo "ended=$(date '+%Y-%m-%d %H:%M:%S %z')"
} > "${exit_file}"
exit "${status}"
