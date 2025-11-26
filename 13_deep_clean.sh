#!/usr/bin/env bash
# ==============================================================================
# 13_deep_clean.sh
# ------------------------------------------------------------------------------
# Purpose: Aggressive cleanup of root directory clutter.
# Moves artifacts to results/cleanup_<date>/ and organizes inputs/outputs.
# ==============================================================================

set -euo pipefail
shopt -s nullglob

DATE=$(date +%Y%m%d)
CLEANUP_DIR="results/cleanup_${DATE}"
INPUT_DIR="data/inputs"
STATS_DIR="results/stats"
PLOTS_DIR="results/plots"
LEGACY_DIR="scripts/legacy"

mkdir -p "$CLEANUP_DIR"
mkdir -p "$INPUT_DIR"
mkdir -p "$STATS_DIR"
mkdir -p "$PLOTS_DIR"
mkdir -p "$LEGACY_DIR"

echo "Starting deep cleanup..."

# 1. Move Temporary Directories (tmp*)
# These are likely leftover from parallel processes or crashes.
echo "Moving temporary directories..."
mv tmp* "$CLEANUP_DIR/" 2>/dev/null || true

# 2. Move Logs and Intermediate Files
echo "Moving logs and intermediate files..."
mv *.log *.err *.out "$CLEANUP_DIR/" 2>/dev/null || true
mv *.sam *.bam "$CLEANUP_DIR/" 2>/dev/null || true
mv *.list *.sorted "$CLEANUP_DIR/" 2>/dev/null || true
mv triangle.ids triangle.sorted "$CLEANUP_DIR/" 2>/dev/null || true

# 3. Move Legacy/Backup Scripts
echo "Moving legacy scripts and backups..."
mv *.bak *.R.bak.* "$LEGACY_DIR/" 2>/dev/null || true
mv 12 allin1.r InitialScript.R InitialChatGPTupdated.R Plotting.R analyze_assemblies*.R "$LEGACY_DIR/" 2>/dev/null || true
mv r_patch.R "$LEGACY_DIR/" 2>/dev/null || true

# 4. Move Plots from Root
echo "Moving plots..."
mv *.png "$PLOTS_DIR/" 2>/dev/null || true

# 5. Move Stats/Summary CSVs
echo "Moving stats and summaries..."
mv stats_*.csv "$STATS_DIR/" 2>/dev/null || true
mv summary_*.csv "$STATS_DIR/" 2>/dev/null || true
mv persistence_by_participant.csv "$STATS_DIR/" 2>/dev/null || true
mv audit_manifest.csv "$STATS_DIR/" 2>/dev/null || true

# 6. Move Inputs (Excel and Batch CSVs)
# Note: We will update 00_process_clinical_data.R to look here.
echo "Moving input data..."
mv *.xlsx "$INPUT_DIR/" 2>/dev/null || true
mv batch*.csv "$INPUT_DIR/" 2>/dev/null || true

# 7. Move other miscellaneous files
mv .DS_Store "$CLEANUP_DIR/" 2>/dev/null || true
mv .RData .Rhistory "$CLEANUP_DIR/" 2>/dev/null || true

echo "✓ Deep cleanup complete. Check $CLEANUP_DIR for moved artifacts."
