#!/bin/bash
# ==============================================================================
# RUN_COMPLETE_ANALYSIS.sh
# ------------------------------------------------------------------------------
# Complete Longcycler-linked pipeline execution script for the rUTIs project.
# The completion marker is published only after every numbered phase, RQ01-RQ10,
# and the final provenance/content gates pass.
# ==============================================================================

set -euo pipefail

# Use deterministic local tool locations in GUI terminals as well as Codex.
export PATH="/Applications/ChatGPT.app/Contents/Resources:/usr/local/bin:/opt/homebrew/bin:${PATH}"

RUN_ID="${RUN_ID:-$(date '+%Y%m%dT%H%M%S%z')}"
PIPELINE_DIR="results/pipeline"
COMPLETE_MARKER="${PIPELINE_DIR}/RUN_COMPLETE.txt"
FAILED_MARKER="${PIPELINE_DIR}/RUN_FAILED.txt"
mkdir -p "${PIPELINE_DIR}"
rm -f "${COMPLETE_MARKER}" "${FAILED_MARKER}"

on_error() {
    status=$?
    rm -f "${COMPLETE_MARKER}"
    {
        echo "Complete analysis failed"
        echo "Run ID: ${RUN_ID}"
        echo "Ended: $(date)"
        echo "Exit status: ${status}"
    } > "${FAILED_MARKER}"
    exit "${status}"
}
trap on_error ERR

publish_complete_marker() {
    require_complete_release_prerequisites
    {
        echo "Complete analysis: PASS"
        echo "Run ID: ${RUN_ID}"
        echo "Ended: $(date)"
        echo "RQ release: results/research_questions/RUN_COMPLETE.txt"
        echo "Final figures: results/figure_audit/final_figure_manifest.csv"
        echo "Figure validation: results/figure_audit/automated_validation_results.txt"
        echo "Claim registry: results/pipeline/longcycler_release_claim_registry.json"
        echo "Verification: results/qc/longcycler_only_pipeline_verification.csv"
        echo "Predicted plasmids: results/plasmids/mob_suite/RUN_COMPLETE.txt"
    } > "${COMPLETE_MARKER}"
    rm -f "${FAILED_MARKER}"
}

require_completed_genomic_amr() {
    local marker="results/amr/RUN_COMPLETE.txt"
    local required
    if [[ ! -f "${marker}" ]]; then
        echo "Genomic-AMR gate failed: ${marker} is missing." >&2
        echo "Complete/resume Script 29 with AMR_ONLY=1 before starting the integrated cohort run." >&2
        return 1
    fi
    for required in \
        "status=complete" \
        "episodes=532" \
        "residents=161" \
        "adjacent_pairs=371" \
        "focused_transitions=9"
    do
        if ! grep -Fqx "${required}" "${marker}"; then
            echo "Genomic-AMR gate failed: ${marker} lacks '${required}'." >&2
            return 1
        fi
    done
}

require_complete_release_prerequisites() {
    local required_file
    local pf_marker="results/plasmids/PLASMIDFINDER_RUN_COMPLETE.txt"
    local mob_marker="results/plasmids/mob_suite/RUN_COMPLETE.txt"
    local rq_marker="results/research_questions/RUN_COMPLETE.txt"
    require_completed_genomic_amr
    for required_file in \
        "${pf_marker}" \
        "${mob_marker}" \
        "results/plasmids/mob_suite/episode_mechanism_profiles.csv" \
        "results/plasmids/mob_suite/plasmid_gene_locations_long.csv" \
        "results/plasmids/mob_suite/plasmid_gene_location_validation.csv" \
        "${rq_marker}" \
        "results/figure_audit/final_figure_manifest.csv" \
        "results/figure_audit/automated_validation_results.txt" \
        "results/pipeline/longcycler_release_claim_registry.json" \
        "results/qc/longcycler_only_pipeline_verification.csv"
    do
        if [[ ! -s "${required_file}" ]]; then
            echo "Release gate failed: required output is missing or empty: ${required_file}" >&2
            return 1
        fi
    done
    for required_file in \
        "episodes=532" \
        "primary_hits=1257" \
        "primary_gene_labels=42" \
        "positive_episodes=422" \
        "valid_no_hit_episodes=110"
    do
        if ! grep -Fqx "${required_file}" "${pf_marker}"; then
            echo "Release gate failed: PlasmidFinder marker lacks '${required_file}'." >&2
            return 1
        fi
    done
    if ! grep -Fqx "episodes=532" "${mob_marker}"; then
        echo "Release gate failed: MOB-suite marker does not cover 532 episodes." >&2
        return 1
    fi
    if ! grep -Fq "Research-question analysis runner: PASS" "${rq_marker}"; then
        echo "Release gate failed: RQ01-RQ10 marker is not PASS." >&2
        return 1
    fi
}

if [[ "${FINALIZE_ONLY:-0}" == "1" ]]; then
    require_completed_genomic_amr
    echo "Finalizing an already completed analysis after refreshing figure audit and final gates..."
    Rscript scripts/validate_final_figures.R
    if [[ "${FINALIZE_REUSE_REVIEWED_FIGURE_QA:-0}" == "1" ]]; then
        echo "↪ Reusing the current reviewed visual-QA derivatives."
        echo "  The consolidated validation gate will still reject missing, stale, or incomplete derivatives."
    else
        Rscript scripts/visual_qa_final_figures.R
    fi
    Rscript scripts/build_figure_audit.R
    echo "✓ Final figure validation, visual-QA derivatives, and inventory refreshed"
    Rscript scripts/verify_uti_not_uti_alignment.R
    echo "✓ Operational-status alignment verified"
    Rscript scripts/build_longcycler_release_claim_registry.R
    echo "✓ Longcycler release claim registry written"
    Rscript scripts/verify_longcycler_only_pipeline.R --stage final
    echo "✓ Complete Longcycler-only release verified"
    publish_complete_marker
    Rscript scripts/write_figure_audit_validation_results.R
    echo "✓ Consolidated figure and pipeline validation report written"
    trap - ERR
    echo "✅ COMPLETE ANALYSIS FINALIZATION FINISHED"
    exit 0
fi

echo "=========================================="
echo "rUTIs Complete Analysis Pipeline"
echo "Run ID: $RUN_ID"
echo "Started: $(date)"
echo "=========================================="

echo "Validating the prerequisite 532-episode genomic-AMR completion marker..."
require_completed_genomic_amr
echo "✓ Genomic-AMR prerequisite is complete"

echo "Preparing the active generated-output roots..."
Rscript scripts/prepare_longcycler_release.R --apply

echo "Preparing the pinned MOB-suite 3.1.9 runtime and database snapshot..."
bash scripts/setup_mob_suite_runtime.sh
echo "✓ MOB-suite runtime complete"

echo "Running dependency, database, policy, and disk preflight..."
REQUIRE_COMPLETED_AMR=1 Rscript scripts/preflight_complete_analysis.R
echo "✓ Preflight complete"

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
echo "[3/3] Deriving display-only poster timepoints..."
Rscript 00d_derive_plot_timepoints.R

echo "✓ Phase 0 complete: primary status maps created"

# ============================================================================
# Phase 1: WGS Processing (LONGEST PHASE)
# ============================================================================
echo ""
echo "━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━"
echo "PHASE 1: WGS Processing [Est: 18-24 hours on this workstation]"
echo "━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━"

echo "[1/8] Refreshing Longcycler assembly metadata..."
Rscript 00_make_assembly_metadata.r
echo "✓ Assembly metadata refreshed"

echo "[2/8] Assembly QC and Longcycler-only canonical selection [~2 min]..."
Rscript 12a_wgs_qc.R
echo "✓ QC and Longcycler-linked clinical cohort complete"

echo "[3/8] Generating clinical plots from selected Longcycler assemblies..."
Rscript 00c_plot_clinical_summary.R
echo "✓ Clinical and assembly-QC plots complete"

echo "[4/8] Core SNP calling with Parsnp..."
FORCE_RERUN_CORE_SNP="${FORCE_RERUN_CORE_SNP:-1}" Rscript 12b_core_snp.R
echo "✓ Core SNPs complete"

echo "[5/8] Pangenome analysis with Panaroo [historically the slowest step]..."
Rscript 12c_panaroo.R
echo "✓ Pangenome complete"

echo "[5b/8] Re-sweeping generated roots before any downstream analysis..."
Rscript scripts/prepare_longcycler_release.R --apply
echo "✓ Pre-downstream Longcycler-only generated-output sweep complete"

echo "[6/8] Selection visualization [~2 min]..."
Rscript 13_visualise_panaroo_selection.R
echo "✓ Visualization complete"

echo "[7/8] Gene presence/absence matrix [~5 min]..."
Rscript 02_gene_presence_analysis.R
echo "✓ Gene presence matrix created"

echo "[8/8] Active Longcycler-only provider/local MLST integration [~3 min]..."
Rscript 06_MLST.R
echo "✓ Active RIVM/provider MLST complete"

echo "✓ Phase 1 complete: WGS data processed"

# ============================================================================
# Phase 1b: Plasmid reconstruction and additional plots
# ============================================================================
echo ""
echo "━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━"
echo "PHASE 1b: Plasmid reconstruction and additional plots [Est: 8-12 hours]"
echo "━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━"

if [ "${RUN_LEGACY_EXPLORATORY_PLOTS:-0}" = "1" ]; then
    Rscript 03_plotting.R
    echo "✓ Legacy exploratory plots created"
else
    echo "↪ Skipping 03_plotting.R legacy exploratory plots."
    echo "  Canonical UTI-vs-Not_UTI VF figures are generated by scripts 23-30."
fi

Rscript 04_gene_breakdown.R
echo "✓ Focused gene analysis complete"

Rscript 05_gene_overview_plots.R
echo "✓ Gene distribution plots created"

Rscript 07_explore_MLST.R
echo "✓ MLST exploration plots created"

Rscript 09_inc_plasmid_network.R
echo "✓ Canonical gene-level PlasmidFinder layer created"

Rscript 09b_mob_plasmid_reconstruction.R
echo "✓ Assembly-first predicted plasmid reconstruction complete"

Rscript 08_core_vs_plasmid.R
echo "✓ Descriptive ST-replicon context complete"

Rscript 10_replicon_heatmap.R
echo "✓ Gene-level replicon-marker heatmap created"

echo "Verifying all selected upstream WGS, VF, MLST, plasmid, core, and GFF inputs..."
Rscript scripts/verify_longcycler_only_pipeline.R --stage upstream
echo "✓ Upstream Longcycler-only release gate passed"

echo "✓ Phase 1b complete: Additional plots generated"

# ============================================================================
# Phase 2: Comparative Genomics
# ============================================================================
echo ""
echo "━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━"
echo "PHASE 2: Comparative Genomics [Est: 15 min]"
echo "━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━"

echo "[1/5] Within-host strain comparison [~5 min]..."
Rscript 11_compare_strains.R --participants ALL
echo "✓ Strain comparison complete"

echo "[2/5] Building canonical VF dataset for UTI vs Not_UTI modelling..."
Rscript 22_vf_build_analysis_dataset.R
echo "✓ VF dataset built"

echo "[3/5] Publishing canonical Longcycler-linked adjacent transitions..."
Rscript scripts/rebuild_longcycler_sensitivity.R
echo "✓ Canonical transition export written"

echo "[4/5] Genotype-phenotype model: UTI vs Not_UTI [~8 min]..."
Rscript 14_genotype_phenotype_model.R
echo "✓ UTI vs Not_UTI model complete"

echo "[5/5] Lineage risk analysis [~2 min]..."
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

echo "[1/4] Reconstructing patient timelines [~1 min]..."
Rscript 15_longitudinal_patterns.R
echo "✓ Timelines reconstructed"

echo "[2/4] Within-host evolution analysis [~1 min]..."
Rscript 16_within_host_evolution.R
echo "✓ Evolution analysis complete"

echo "[3/4] Basic variant annotation [~1 min]..."
Rscript 18_annotate_variants.R
echo "✓ Variant annotation complete"

echo "[4/4] Host context analysis [~1 min]..."
Rscript 19_host_context.R
echo "✓ Host context analysis complete"

if [ "${RUN_LEGACY_PUBLICATION_FIGURES:-0}" = "1" ]; then
    echo "[legacy optional] Generating obsolete script-21 publication figures by explicit request..."
    Rscript 21_publication_figures.R
    echo "✓ Legacy script-21 figures generated for provenance only"
else
    echo "↪ Skipping obsolete 21_publication_figures.R output (default)."
    echo "  The canonical numbered thesis figures are generated only after RQ01-RQ10."
fi

echo "✓ Phase 3 complete: longitudinal candidate events and host context generated"

# ============================================================================
# Phase 4: Virulence Factor (VF) Deep Analysis
# ============================================================================
echo ""
echo "━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━"
echo "PHASE 4: VF Deep Analysis and canonical transitions [Est: 30-90 min]"
echo "━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━"

echo "[1/14] Cross-sectional analysis..."
Rscript 23_vf_cross_sectional.R
echo "✓ VF cross-sectional done"

echo "[2/14] Longitudinal dynamics..."
Rscript 24_vf_longitudinal_dynamics.R
echo "✓ VF longitudinal done"

echo "[3/14] Lineage confounding check..."
Rscript 25_vf_lineage_vf_interaction.R
echo "✓ Lineage interaction checked"

echo "[4/14] VF gene modules..."
Rscript 26_vf_define_gene_modules.R
echo "✓ VF modules defined"

echo "[5/14] VF score framework..."
Rscript 27_vf_score_framework.R
echo "✓ VF scores analysed"

echo "[6/14] Transition case studies..."
Rscript 28_vf_transition_case_studies.R
echo "✓ Transition case studies done"

echo "[7/14] Deep variant annotation against the canonical 9-transition case index..."
Rscript 20_variant_annotation_deep.R
echo "✓ Deep annotation complete"

echo "[8/14] Genomic AMR and VF/plasmid integration..."
Rscript 29_vf_amr_combined_profile.R
echo "✓ Genomic AMR and VF/plasmid integration done"

echo "[9/14] UTI/Not_UTI diagnostic statistics..."
Rscript 32_uti_not_uti_diagnostic_stats.R
echo "✓ UTI/Not_UTI diagnostics written"

echo "[10/14] Rebuilding the selected-cohort UTI/Not_UTI count audit..."
Rscript scripts/audit_uti_status_count_explanation.R
echo "✓ Selected Longcycler count audit written"

echo "[11/14] Mechanism-first add-on..."
Rscript 33_mechanism_first_addon.R
echo "✓ Mechanism-first add-on written"

echo "[12/14] Robustness-first add-on..."
Rscript 34_robustness_first_addon.R
echo "✓ Robustness-first add-on written"

echo "[13/14] Targeted statistical sensitivity add-on..."
Rscript 36_statistical_sensitivity_addon.R
echo "✓ Targeted statistical sensitivity add-on written"

echo "[14/14] Project summary tables..."
Rscript 30_vf_project_summary_tables.R
echo "✓ VF project summary tables written"

echo "Archiving stale ASB-vs-UTI generated outputs..."
Rscript scripts/archive_legacy_asb_uti_outputs.R
echo "✓ Legacy ASB-vs-UTI generated outputs archived"

echo "✓ Phase 4 complete: VF Deep Analysis"

# ============================================================================
# Phase 5: Prespecified Research Questions and Canonical Figure Release
# ============================================================================
echo ""
echo "━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━"
echo "PHASE 5: Research Questions and canonical thesis figures [Est: 1-3 hours]"
echo "━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━"

echo "[1/6] Running prespecified research questions RQ01-RQ10..."
Rscript scripts/research_questions/run_all.R
echo "✓ RQ01-RQ10 release complete"

echo "[2/6] Preparing exact-reference-aware variant coordinates for Fig08..."
Rscript scripts/prepare_reference_aware_variants.R
echo "✓ Reference-aware variant validation and plotting table written"

echo "[3/6] Generating the numbered canonical thesis figure pack..."
Rscript -e 'options(warn = 2); source("35_final_figure_pack.R", chdir = FALSE)'
echo "✓ Canonical thesis figure pack written"

echo "[4/6] Running automated final-figure acceptance checks..."
Rscript scripts/validate_final_figures.R
echo "✓ Final figure validation passed"

echo "[5/6] Rendering PDF, thesis-size, greyscale, and colour-vision QA derivatives..."
Rscript scripts/visual_qa_final_figures.R
echo "✓ Reproducible visual-QA derivatives written"

echo "[6/6] Rebuilding the repository-wide figure audit and inventory..."
Rscript scripts/build_figure_audit.R
echo "✓ Figure audit and artifact census refreshed"

echo "✓ Phase 5 complete: RQ01-RQ10 and canonical numbered figures released"

# ============================================================================
# Phase 6: Final Release Gates
# ============================================================================
echo ""
echo "━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━"
echo "PHASE 6: Final Longcycler-only release gates [Est: 5-10 min]"
echo "━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━"

echo "[1/3] Verifying operational-status and Longcycler-linked output alignment..."
Rscript scripts/verify_uti_not_uti_alignment.R
echo "✓ Operational-status alignment verified"

echo "[2/3] Publishing the single audited claim registry for decks and handouts..."
Rscript scripts/build_longcycler_release_claim_registry.R
echo "✓ Longcycler release claim registry written"

echo "[3/3] Verifying the complete Longcycler-only release..."
Rscript scripts/verify_longcycler_only_pipeline.R --stage final
echo "✓ Complete Longcycler-only release verified"

echo "✓ Phase 6 complete: final release gates passed"

# ============================================================================
# Done
# ============================================================================
echo ""
echo "=========================================="
publish_complete_marker
Rscript scripts/write_figure_audit_validation_results.R
echo "✓ Consolidated figure and pipeline validation report written"
trap - ERR

echo "✅ COMPLETE ANALYSIS FINISHED"
echo "Ended: $(date)"
echo "Expected workstation runtime: ~20-30 hours"
echo "=========================================="
echo ""
echo "📊 Key Outputs:"
echo "   - results/clinical/status_map.csv"
echo "   - results/clinical/analysis_cohort_longcycler.csv"
echo "   - results/vf/vf_pa_all.csv"
echo "   - results/vf/vf_analysis_ready.csv"
echo "   - results/vf/vf_binary_uti_ready.csv"
echo "   - results/strain_compare/pairwise_metrics.csv"
echo "   - results/plasmids/mob_suite/RUN_COMPLETE.txt"
echo "   - results/plasmids/mob_suite/plasmids_long.csv"
echo "   - results/plasmids/mob_suite/plasmid_gene_locations_long.csv"
echo "   - results/models/model_dataset_denominator.csv"
echo "   - results/longitudinal/variant_annotation_detailed.csv"
echo "   - results/research_questions/RUN_COMPLETE.txt"
echo "   - results/figure_audit/final_figure_manifest.csv"
echo "   - results/figure_audit/figure_inventory.csv"
echo "   - results/figure_audit/automated_validation_results.txt"
echo "   - results/figure_audit/validation_results.txt"
echo "   - results/pipeline/RUN_COMPLETE.txt"
echo "   - plots/final/Fig01-Fig08.{png,pdf}"
echo "   - plots/final/supplementary/FigS01-FigS10.{png,pdf}"
echo "   - plots/vf/ (VF heatmaps)"
echo "   - plots/mlst/ (MLST plots)"
echo "   - plots/plasmids/ (Plasmid networks)"
echo ""
echo "📖 Next steps:"
echo "   1. Review results/summary/final_key_results_summary.md"
echo "   2. Review results/models/model_interpretation_warnings.txt"
echo "   3. Review plots/final and plots/final/supplementary for the canonical thesis figures"
