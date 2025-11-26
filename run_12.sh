#!/usr/bin/env bash
set -euo pipefail

# ---- config you can tweak ----
THREADS="${THREADS:-8}"
PIDS="${PIDS:-ALL}"              # e.g. "P001,P007" or "ALL"
MODE="${MODE:-auto}"             # auto|reads|assemblies
ENV_NAME="${ENV_NAME:-yellow-wgs}"

# ---- activate conda env with pinned tools ----
# (create once with the YAML below)
if command -v conda >/dev/null 2>&1; then
  eval "$(conda shell.bash hook)"
  safe_conda() { set +u; conda "$@"; local rc=$?; set -u; return $rc; }
  safe_conda activate "$ENV_NAME"
fi

# ---- harden environment ----
export OMP_NUM_THREADS=1
export MKL_NUM_THREADS=1
export WGS_PARALLEL=0
export LC_ALL=C
export LANG=C

# ---- preflight binaries ----
need=(minimap2 samtools bcftools mash snp-dists nucmer show-diff mosdepth)
miss=()
for b in "${need[@]}"; do
  command -v "$b" >/dev/null 2>&1 || miss+=("$b")
done
# Pan-genome: at least one of these
if ! command -v panaroo >/dev/null 2>&1 && ! command -v roary >/dev/null 2>&1; then
  echo "Note: panaroo/roary not found - stage C (pangenome) will be skipped."
fi
# Optional but nice:
if ! command -v gubbins >/dev/null 2>&1 && ! command -v run_gubbins.py >/dev/null 2>&1; then
  echo "Note: gubbins not found — recombination masking will be skipped."
fi

if ((${#miss[@]})); then
  echo "Missing tools: ${miss[*]}"
  echo "Activate the right conda env or install them, then rerun."
  exit 1
fi

# ---- directories & logging ----
mkdir -p results/wgs/logs results/plots
ts="$(date +%Y%m%d_%H%M%S)"
LOG="results/wgs/logs/run12_${ts}.log"

# ---- run the R script ----
echo ">>> Starting 12_wgs_exact_compare.R (mode=${MODE}, threads=${THREADS}, pids=${PIDS})" | tee -a "$LOG"
Rscript 12_wgs_exact_compare.R \
  --mode "${MODE}" \
  --threads "${THREADS}" \
  --pids "${PIDS}" 2>&1 | tee -a "$LOG"

echo ">>> Done. See $LOG and results/wgs/reports/ for outputs."
