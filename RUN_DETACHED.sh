#!/bin/bash

# Submit the complete analysis to the user's launchd session so the long-running
# worker survives terminal and Codex-thread lifecycle changes.

set -euo pipefail

cd "$(dirname "$0")"

label="com.openai.ruti.longcycler"
domain="gui/$(id -u)"
worker="$(pwd)/RUN_LAUNCHD_WORKER.sh"

if launchctl print "${domain}/${label}" >/dev/null 2>&1; then
    echo "A launchd job named ${label} already exists. Inspect it before relaunching."
    exit 1
fi

launchctl submit -l "${label}" -- /bin/bash "${worker}"
echo "Submitted complete analysis to launchd."
echo "Job: ${domain}/${label}"

