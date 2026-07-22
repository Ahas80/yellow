#!/usr/bin/env bash
set -euo pipefail

# Complete, resumable CGE VirulenceFinder sensitivity workflow.
#
# ALL scientific, cohort, software, database, threshold, bootstrap and normal
# execution settings live in the commented file below:
#   config/virulencefinder_sensitivity.toml
# Each user-changeable value is labelled USER-TUNABLE and explained there.
#
# Optional launch-only override:
#   VIRULENCEFINDER_CONFIG=/path/to/another.toml bash RUN_VIRULENCEFINDER_SENSITIVITY.sh
# This is useful for a deliberately separate sensitivity configuration. The
# default file preserves the approved scientific contract.

PROJECT_ROOT="$(cd "$(dirname "${BASH_SOURCE[0]}")" && pwd)"
cd "$PROJECT_ROOT"

CONFIG="${VIRULENCEFINDER_CONFIG:-config/virulencefinder_sensitivity.toml}"
RUNNER="scripts/run_virulencefinder_batch.py"
ANALYSIS="scripts/analyze_virulencefinder_sensitivity.R"

echo "[1/7] Validating the pinned environment and all 532 selected FASTAs"
python3 "$RUNNER" --config "$CONFIG" --preflight-only

echo "[2/7] Running runner regression tests"
python3 -m unittest -v tests/test_virulencefinder_batch.py

echo "[3/7] Running the deterministic six-job pilot and cache-reuse gate"
python3 "$RUNNER" --config "$CONFIG" --pilot-only

echo "[4/7] Running/resuming all 1,064 single-assembly jobs"
python3 "$RUNNER" --config "$CONFIG" --resume

echo "[5/7] Revalidating and reconciling every cached JSON and combined table"
python3 "$RUNNER" --config "$CONFIG" --verify-only

echo "[6/7] Running the separate RQ06-RQ08 sensitivity analysis"
Rscript "$ANALYSIS"

echo "[7/7] Checking final release markers"
test -f results/virulencefinder_cge_3_2_1/BATCH_COMPLETE.txt
test -f results/virulencefinder_cge_3_2_1/RUN_COMPLETE.txt
grep -q "All acceptance checks: PASS" results/virulencefinder_cge_3_2_1/RUN_COMPLETE.txt

echo "VirulenceFinder sensitivity workflow complete: PASS"
