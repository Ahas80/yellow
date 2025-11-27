#!/bin/bash
# WGS Pipeline Execution Guide - Optimized for Terminal
# Copy-paste these commands into your terminal at ~/Desktop/rUTIs

# ============================================================================
# SETUP
# ============================================================================

cd ~/Desktop/rUTIs

# Activate conda environment (adjust name if different)
conda activate asm-snp-x86

# Create logs directory
mkdir -p logs

# ============================================================================
# STEP 1: Check if 12b is still running
# ============================================================================

# Check Parsnp status
ps aux | grep parsnp_core | grep -v grep

# If it's still running and you want to let it finish:
# - It's been running ~70+ minutes
# - Should complete eventually (large dataset = slow)
# - Monitor with: tail -f logs/12b_core_snp_pass2.log

# If you want to kill it and skip to critical steps:
# pkill -f parsnp_core
# rm -rf results/wgs/core/parsnp_out results/wgs/core/temp_fastas

# ============================================================================
# STEP 2: Run 12b in background (if not already done or if restarting)
# ============================================================================

# Run with output to log file, in background
nohup Rscript 12b_core_snp.R > logs/12b_core_snp_pass2.log 2>&1 &

# Note the PID (process ID) it prints, or find it:
ps aux | grep "12b_core_snp.R" | grep -v grep

# Monitor progress:
tail -f logs/12b_core_snp_pass2.log
# Press Ctrl+C to stop monitoring (process keeps running)

# ============================================================================
# STEP 3: Run 12c Panaroo (AFTER 12b completes)
# ============================================================================

# Check 12b completion first:
ls -lh results/wgs/core/strain_pairs.csv
# If this file exists, 12b is done!

# Run Panaroo
nohup Rscript 12c_panaroo.R > logs/12c_panaroo_pass2.log 2>&1 &

# Monitor:
tail -f logs/12c_panaroo_pass2.log

# Panaroo takes ~10-30 mins for 361 genomes

# ============================================================================
# STEP 4: Run 13 Selection Visualization
# ============================================================================

nohup Rscript 13_visualise_panaroo_selection.R > logs/13_visualise_pass2.log 2>&1 &

# This is quick (~1-2 mins)

# ============================================================================
# STEP 5: Run 11 Strain Comparison (CRITICAL - Can run without 12b!)
# ============================================================================

# This script does NOT depend on 12b core SNPs!
# It can run with just VF/MLST/plasmid data

nohup Rscript 11_compare_strains.R --participants ALL > logs/11_compare_strains_pass2.log 2>&1 &

# Monitor:
tail -f logs/11_compare_strains_pass2.log

# Takes ~5-15 mins

# CRITICAL VERIFICATION after completion:
head -20 results/strain_compare/pairwise_metrics.csv
# Check that Participant_id_A and Participant_id_B have TRUE IDs (e.g., 20002, 100009)
# NOT batch IDs (PR0010)

# ============================================================================
# STEP 6: Run 14 GWAS (CRITICAL - Also independent of 12b!)
# ============================================================================

nohup Rscript 14_genotype_phenotype_model.R > logs/14_gwas_pass2.log 2>&1 &

# Monitor:
tail -f logs/14_gwas_pass2.log

# Takes ~5-10 mins

# CRITICAL VERIFICATION after completion:
head -20 results/models/gwas_univariable_stats.csv
# Should have populated results

# ============================================================================
# MONITORING COMMANDS
# ============================================================================

# See all running R processes:
ps aux | grep Rscript | grep -v grep

# Check system resources:
top -o cpu

# See recent log output for all scripts:
tail -20 logs/*pass2.log

# ============================================================================
# QUICK STATUS CHECK
# ============================================================================

# Run this to see what's complete:
echo "=== Pipeline Status ==="
echo "12a QC: $([ -f results/wgs/qc_summary.csv ] && echo '✅' || echo '❌')"
echo "12b Core SNPs: $([ -f results/wgs/core/strain_pairs.csv ] && echo '✅' || echo '❌')"
echo "12c Panaroo: $([ -f results/wgs/pan/gene_presence_absence.csv ] && echo '✅' || echo '❌')"
echo "11 Strain Compare: $([ -f results/strain_compare/pairwise_metrics.csv ] && echo '✅' || echo '❌')"
echo "14 GWAS: $([ -f results/models/gwas_univariable_stats.csv ] && echo '✅' || echo '❌')"

# ============================================================================
# RECOMMENDED APPROACH (Most Efficient)
# ============================================================================

# Since 12b is very slow and 11 & 14 don't need it:

# 1. Let 12b keep running in background (or skip it for now)
# 2. Run 11 (strain comparison) NOW - independently
# 3. Run 14 (GWAS) NOW - independently
# 4. These two are CRITICAL for verifying the Participant ID fix
# 5. Come back to 12b/12c/13 later if needed

# Parallel execution (run both critical scripts at once):
nohup Rscript 11_compare_strains.R --participants ALL > logs/11_compare_strains_pass2.log 2>&1 &
nohup Rscript 14_genotype_phenotype_model.R > logs/14_gwas_pass2.log 2>&1 &

# Monitor both:
tail -f logs/11_compare_strains_pass2.log logs/14_gwas_pass2.log

# ============================================================================
# VERIFICATION (After scripts complete)
# ============================================================================

# 1. Check Participant IDs are correct (not batch IDs):
cut -d',' -f1,7 results/strain_compare/pairwise_metrics.csv | head -10
# Should show real participant IDs like: 20002, 100009, etc.

# 2. Ensure no batch IDs in outputs:
grep -r "PR0010" results/strain_compare/*.csv results/models/*.csv
# Should return nothing (empty)

# 3. Check within-host pairs are truly within-host:
awk -F',' 'NR>1 && $18=="TRUE" {if ($1 != $7) print "MISMATCH:", $1, $7}' \
  results/strain_compare/pairwise_metrics.csv
# Should return nothing (empty) = all good!

# ============================================================================
# IF ANYTHING FAILS
# ============================================================================

# Check the log file for errors:
grep -i "error\|failed\|stop" logs/[script_name]_pass2.log

# Common fixes:
# - Missing tools: conda activate asm-snp-x86
# - Missing packages: install.packages(c("tidyverse", "lme4", "broom.mixed"))
# - Permission errors: chmod +x *.R

echo "✅ Commands ready! Start with the RECOMMENDED APPROACH section."
