#!/usr/bin/env bash
set -euo pipefail

cd "$(dirname "$0")/.."

marker="results/amr/RUN_COMPLETE.txt"
if [[ ! -f "${marker}" ]] || ! grep -q '^status=complete$' "${marker}"; then
  echo "Script 29 has not published a successful genomic-AMR completion marker." >&2
  exit 1
fi

echo "[1/9] Proving Script 29 cache-safe idempotence..."
AMR_CALLER_WORKERS="${AMR_CALLER_WORKERS:-10}" \
AMRFINDER_THREADS_PER_CALL="${AMRFINDER_THREADS_PER_CALL:-1}" \
AMR_BOOTSTRAP_REPS="${AMR_BOOTSTRAP_REPS:-10000}" \
  Rscript 29_vf_amr_combined_profile.R

echo "[2/9] Refreshing mechanism-first results from Script 29 profiles..."
Rscript 33_mechanism_first_addon.R

echo "[3/9] Rebuilding project summary tables and the AMR results narrative..."
Rscript 30_vf_project_summary_tables.R

echo "[4/9] Rebuilding the canonical/supplementary figure pack, including FigS07..."
Rscript -e 'options(warn = 2); source("35_final_figure_pack.R", chdir = FALSE)'

echo "[5/9] Validating final figures..."
Rscript scripts/validate_final_figures.R

echo "[6/9] Rebuilding figure QA derivatives..."
Rscript scripts/visual_qa_final_figures.R

echo "[7/9] Refreshing the figure audit..."
Rscript scripts/build_figure_audit.R

echo "[8/9] Rebuilding the audited claim registry..."
Rscript scripts/build_longcycler_release_claim_registry.R

echo "[9/9] Running the final Longcycler-only release gate..."
Rscript scripts/verify_longcycler_only_pipeline.R --stage final

echo "Genomic-AMR results integration and final release validation completed."
