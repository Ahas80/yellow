#!/bin/bash
# ==============================================================================
# RUN_COMPLETE_ANALYSIS.sh
# ------------------------------------------------------------------------------
# Complete pipeline execution script for rUTIs project
# Runs all 4 phases in order
# ==============================================================================

set -e  # Exit on error
set -u  # Exit on undefined variable

echo "=========================================="
echo "rUTIs Complete Analysis Pipeline"
echo "Started: $(date)"
echo "=========================================="

# Check conda environment
if [ -z "${CONDA_DEFAULT_ENV:-}" ]; then
    echo "⚠️  WARNING: No conda environment active. Activate with:"
    echo "   conda activate asm-snp-x86"
    exit 1
fi

echo "✓ Conda environment: $CONDA_DEFAULT_ENV"

# ============================================================================
# Phase 0: Clinical Data Foundation
# ============================================================================
echo ""
echo "━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━"
echo "PHASE 0: Clinical Data Foundation [Est: 2 min]"
echo "━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━"

echo "[1/3] Loading clinical data..."
Rscript 00a_load_clean_clinical.R
echo "[2/3] Classifying episodes..."
Rscript 00b_classify_episodes.R
echo "[3/3] Generating clinical plots..."
Rscript 00c_plot_clinical_summary.R

echo "✓ Phase 0 complete: status_map.csv created"

# ============================================================================
# Phase 1: WGS Processing (LONGEST PHASE)
# ============================================================================
echo ""
echo "━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━"
echo "PHASE 1: WGS Processing [Est: 1-2 hours total]"
echo "━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━"

echo "[1/6] Assembly QC [~2 min]..."
Rscript 12a_wgs_qc.R
echo "✓ QC complete"

echo "[2/6] Core SNP calling with Parsnp [~30-60 min - SLOWEST STEP]..."
Rscript 12b_core_snp.R
echo "✓ Core SNPs complete"

echo "[3/6] Pangenome analysis with Panaroo [~15 min]..."
Rscript 12c_panaroo.R
echo "✓ Pangenome complete"

echo "[4/6] Selection visualization [~2 min]..."
Rscript 13_visualise_panaroo_selection.R
echo "✓ Visualization complete"

echo "[5/6] Gene presence/absence matrix [~5 min]..."
Rscript 02_gene_presence_analysis.R
echo "✓ Gene presence matrix created"

echo "[6/6] MLST typing [~3 min]..."
Rscript 06_MLST.R
echo "✓ MLST typing complete"

echo "✓ Phase 1 complete: WGS data processed"

# ============================================================================
# Phase 1b: Additional Plots & Exploration (OPTIONAL)
# ============================================================================
echo ""
echo "━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━"
echo "PHASE 1b: Additional Plots [Est: 10 min total]"
echo "━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━"

Rscript 03_plotting.R
echo "✓ VF heatmaps and summary plots created"

Rscript 04_gene_breakdown.R
echo "✓ Focused gene analysis complete"

Rscript 05_gene_overview_plots.R
echo "✓ Gene distribution plots created"

Rscript 07_explore_MLST.R
echo "✓ MLST exploration plots created"

Rscript 08_core_vs_plasmid.R
echo "✓ Core vs plasmid comparison complete"

Rscript 09_inc_plasmid_network.R
echo "✓ Plasmid network visualization created"

Rscript 10_replicon_heatmap.R
echo "✓ Plasmid replicon heatmap created"

echo "✓ Phase 1b complete: Additional plots generated"

# ============================================================================
# Phase 2: Comparative Genomics
# ============================================================================
echo ""
echo "━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━"
echo "PHASE 2: Comparative Genomics [Est: 15 min]"
echo "━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━"

echo "[1/3] Within-host strain comparison [~5 min]..."
Rscript 11_compare_strains.R --participants ALL
echo "✓ Strain comparison complete"

echo "[2/3] GWAS for UTI-associated genes [~8 min]..."
Rscript 14_genotype_phenotype_model.R
echo "✓ GWAS complete"

echo "[3/3] Lineage risk analysis [~2 min]..."
Rscript 17_lineage_analysis.R
echo "✓ Lineage analysis complete"

echo "✓ Phase 2 complete: Comparative analysis done"

# ============================================================================
# Phase 3: Longitudinal & Mechanism
# ============================================================================
echo ""
echo "━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━"
echo "PHASE 3: Longitudinal & Mechanism [Est: 5 min]"
echo "━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━"

echo "[1/6] Reconstructing patient timelines [~1 min]..."
Rscript 15_longitudinal_patterns.R
echo "✓ Timelines reconstructed"

echo "[2/6] Within-host evolution analysis [~1 min]..."
Rscript 16_within_host_evolution.R
echo "✓ Evolution analysis complete"

echo "[3/6] Basic variant annotation [~1 min]..."
Rscript 18_annotate_variants.R
echo "✓ Variant annotation complete"

echo "[4/6] Deep variant annotation with GFF [~1 min]..."
Rscript 20_variant_annotation_deep.R
echo "✓ Deep annotation complete"

echo "[5/6] Host context analysis [~1 min]..."
Rscript 19_host_context.R
echo "✓ Host context analysis complete"

echo "[6/6] Publication figures [~1 min]..."
Rscript 21_publication_figures.R
echo "✓ Publication figures generated"

echo "✓ Phase 3 complete: Mechanism identified"

# ============================================================================
# Phase 4: Virulence Factor (VF) Deep Analysis
# ============================================================================
echo ""
echo "━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━"
echo "PHASE 4: VF Deep Analysis [Est: 1 min]"
echo "━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━"

echo "[1/4] Building canonical VF dataset..."
Rscript 22_vf_build_analysis_dataset.R
echo "✓ VF dataset built"

echo "[2/4] Cross-sectional analysis..."
Rscript 23_vf_cross_sectional.R
echo "✓ VF cross-sectional done"

echo "[3/4] Longitudinal dynamics..."
Rscript 24_vf_longitudinal_dynamics.R
echo "✓ VF longitudinal done"

echo "[4/4] Lineage confounding check..."
Rscript 25_vf_lineage_vf_interaction.R
echo "✓ Lineage interaction checked"

echo "✓ Phase 4 complete: VF Deep Analysis"

# ============================================================================
# Done
# ============================================================================
echo ""
echo "=========================================="
echo "✅ COMPLETE ANALYSIS FINISHED"
echo "Ended: $(date)"
echo "Total Runtime: ~1.5-2.5 hours"
echo "=========================================="
echo ""
echo "📊 Key Outputs:"
echo "   - results/clinical/status_map.csv"
echo "   - results/vf/vf_pa_all.csv"
echo "   - results/strain_compare/pairwise_metrics.csv"
echo "   - results/longitudinal/variant_annotation_detailed.csv"
echo "   - plots/publication/Fig*.png"
echo "   - plots/vf/ (VF heatmaps)"
echo "   - plots/mlst/ (MLST plots)"
echo "   - plots/plasmids/ (Plasmid networks)"
echo ""
echo "📖 Next steps:"
echo "   1. Review REPO_CLEANUP_MANIFEST.md"
echo "   2. Check results/ and plots/ after run"
echo "   3. See README.md for workflow details"
