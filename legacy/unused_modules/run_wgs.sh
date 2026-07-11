#!/usr/bin/env bash
set -Eeo pipefail

# --------- configurable defaults (you can also pass CLI args) -----------------
PIDS_DEFAULT="20002"
MODE_DEFAULT="assemblies"   # assemblies|reads|auto
THREADS_DEFAULT="8"
AUTO_PATCH_HAPLOID="1"      # 1 = add '--ploidy 1' to bcftools call if missing

# --------- usage ---------------------------------------------------------------
if [[ "${1:-}" == "-h" || "${1:-}" == "--help" ]]; then
  cat <<'USAGE'
Usage:
  ./run_wgs.sh [PIDS [MODE [THREADS]]]

Examples:
  ./run_wgs.sh
  ./run_wgs.sh 20002 assemblies 8
  ./run_wgs.sh ALL auto 16

Run this from the rUTIs folder (where 12_wgs_exact_compare.R lives).
USAGE
  exit 0
fi

# --------- locate project & args ----------------------------------------------
SCRIPT_DIR="$( cd "$( dirname "${BASH_SOURCE[0]}" )" && pwd )"
cd "$SCRIPT_DIR"

PIDS="${1:-$PIDS_DEFAULT}"
MODE="${2:-$MODE_DEFAULT}"
THREADS="${3:-$THREADS_DEFAULT}"

RFILE="12_wgs_exact_compare.R"
[[ -f "$RFILE" ]] || { echo "❌ Can't find $RFILE in: $SCRIPT_DIR"; exit 1; }

# --------- conda env (create if missing, then activate) -----------------------
if ! command -v conda >/dev/null 2>&1; then
  echo "❌ 'conda' not found. Install Miniconda/Conda first, then re-run."
  exit 1
fi
eval "$(conda shell.bash hook)"

if ! conda env list | awk '{print $1}' | grep -qx "wgs"; then
  echo "🧪 Creating Conda env 'wgs' (first-time only)…"
  conda create -y -n wgs -c conda-forge -c bioconda \
    minimap2 samtools bcftools mash snp-dists mummer4 mosdepth bedtools gubbins \
    panaroo roary r-base r-essentials
fi

echo "📦 Activating env 'wgs'…"
conda activate wgs

# --------- versions (nice for troubleshooting) --------------------------------
echo "🔎 Tool versions:"
for t in minimap2 samtools bcftools mash snp-dists nucmer show-diff show-coords mosdepth bedtools gubbins panaroo roary Rscript; do
  if command -v "$t" >/dev/null 2>&1; then
    v="$("$t" --version 2>/dev/null | head -n1 || true)"
    printf "  - %-12s %s\n" "$t" "$v"
  else
    printf "  - %-12s (not found)\n" "$t"
  fi
done

# --------- optional one-time auto-patch: haploid bcftools call ----------------
if [[ "$AUTO_PATCH_HAPLOID" == "1" && -w "$RFILE" ]]; then
  if ! grep -q -- '--ploidy[[:space:]]*1' "$RFILE"; then
    echo "🩹 Patching $RFILE: adding '--ploidy 1' to bcftools call…"
    cp -n "$RFILE" "${RFILE}.bak.$(date +%Y%m%d-%H%M%S)"
    perl -0777 -pe 's/\bcall\s+-mv\b/call -m --ploidy 1/g' -i "$RFILE"
  fi
fi

# --------- quick syntax check & self-test -------------------------------------
echo "🧪 Parsing $RFILE…"
Rscript -e "parse('$RFILE')" >/dev/null

echo "🧪 R self-test…"
Rscript "$RFILE" --selftest

# --------- run the pipeline ---------------------------------------------------
echo "🚀 Running WGS pipeline…"
echo "    PIDs=$PIDS  MODE=$MODE  THREADS=$THREADS"
Rscript "$RFILE" --pids "$PIDS" --mode "$MODE" --threads "$THREADS"

echo "✅ Done. Check results/wgs/logs and results/plots/_progress.png"
