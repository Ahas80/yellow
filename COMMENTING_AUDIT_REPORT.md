# COMMENTING AUDIT REPORT

## 1) Summary

- **Canonical anchor reviewed:** `RUN_COMPLETE_ANALYSIS.sh`.
- **Total scripts reviewed:** 35 (30 active pipeline scripts + 5 active helper/config scripts sourced by the pipeline).
- **Total scripts edited:** 19.
- **Total scripts left unchanged:** 16.

### Unchanged scripts and rationale
These scripts already contained clear scientific headers and adequate block-level comments for goals, IO, and key analysis context; no behavior-preserving comment edits were needed:

- `00a_load_clean_clinical.R`
- `00b_classify_episodes.R`
- `00c_plot_clinical_summary.R`
- `02_gene_presence_analysis.R`
- `04_gene_breakdown.R`
- `05_gene_overview_plots.R`
- `06_MLST.R`
- `11_compare_strains.R`
- `14_genotype_phenotype_model.R`
- `22_vf_build_analysis_dataset.R`
- `23_vf_cross_sectional.R`
- `24_vf_longitudinal_dynamics.R`
- `25_vf_lineage_vf_interaction.R`
- `00_config.R`
- `R/clinical_helpers.R`
- `R/plot_helpers.R`

## 2) Script-by-script table

| Script | Pipeline phase | Commenting status before | What was improved | Logic changed? | Validation status |
|---|---|---|---|---|---|
| RUN_COMPLETE_ANALYSIS.sh | Orchestration | Good but brief | Expanded canonical-entry header incl. conda rationale, phase purpose, reproducibility notes | No | `bash -n` pass |
| 00a_load_clean_clinical.R | Phase 0 | Already strong | No change required | No | Not parsed (Rscript unavailable) |
| 00b_classify_episodes.R | Phase 0 | Already strong | No change required | No | Not parsed (Rscript unavailable) |
| 00c_plot_clinical_summary.R | Phase 0 | Already strong | No change required | No | Not parsed (Rscript unavailable) |
| 12a_wgs_qc.R | Phase 1 | Good | Added explicit key decisions/pipeline position/limitations language | No | Not parsed (Rscript unavailable) |
| 12b_core_snp.R | Phase 1 | Good | Added rationale for QC-gated SNP distances and caveats | No | Not parsed (Rscript unavailable) |
| 12c_panaroo.R | Phase 1 | Good | Added pangenome interpretation caveats and position notes | No | Not parsed (Rscript unavailable) |
| 13_visualise_panaroo_selection.R | Phase 1 | Good | Added bias-diagnostic framing and limitations | No | Not parsed (Rscript unavailable) |
| 02_gene_presence_analysis.R | Phase 1 | Already strong | No change required | No | Not parsed (Rscript unavailable) |
| 06_MLST.R | Phase 1 | Already strong | No change required | No | Not parsed (Rscript unavailable) |
| 03_plotting.R | Phase 1b | Moderate | Added exploratory status, pipeline placement, repeated-measure caveat | No | Not parsed (Rscript unavailable) |
| 04_gene_breakdown.R | Phase 1b | Already strong | No change required | No | Not parsed (Rscript unavailable) |
| 05_gene_overview_plots.R | Phase 1b | Already strong | No change required | No | Not parsed (Rscript unavailable) |
| 07_explore_MLST.R | Phase 1b | Moderate | Added canonical-join rationale and MLST limitations | No | Not parsed (Rscript unavailable) |
| 08_core_vs_plasmid.R | Phase 1b | Moderate | Added key decisions, placement, and observational caveats | No | Not parsed (Rscript unavailable) |
| 09_inc_plasmid_network.R | Phase 1b | Moderate | Added cache rationale and replicon interpretation caveats | No | Not parsed (Rscript unavailable) |
| 10_replicon_heatmap.R | Phase 1b | Moderate | Added descriptive-only interpretation and clustering caveats | No | Not parsed (Rscript unavailable) |
| 11_compare_strains.R | Phase 2 | Already strong | No change required | No | Not parsed (Rscript unavailable) |
| 14_genotype_phenotype_model.R | Phase 2 | Already strong | No change required | No | Not parsed (Rscript unavailable) |
| 17_lineage_analysis.R | Phase 2 | Moderate | Added ST-risk stability constraints and repeated-measures caveat | No | Not parsed (Rscript unavailable) |
| 15_longitudinal_patterns.R | Phase 3 | Good | Added graph-clustering design rationale and timeline caveats | No | Not parsed (Rscript unavailable) |
| 16_within_host_evolution.R | Phase 3 | Good | Added explicit hypothesis-generating framing | No | Not parsed (Rscript unavailable) |
| 18_annotate_variants.R | Phase 3 | Good | Added provenance and annotation-dependence caveats | No | Not parsed (Rscript unavailable) |
| 20_variant_annotation_deep.R | Phase 3 | Good | Added scoped-candidate rationale and generalizability caveat | No | Not parsed (Rscript unavailable) |
| 19_host_context.R | Phase 3 | Good | Added host-confounding interpretation notes | No | Not parsed (Rscript unavailable) |
| 21_publication_figures.R | Phase 3 | Moderate | Replaced duplicate “Purpose” with design/position/limitations sections | No | Not parsed (Rscript unavailable) |
| 22_vf_build_analysis_dataset.R | Phase 4 | Already strong | No change required | No | Not parsed (Rscript unavailable) |
| 23_vf_cross_sectional.R | Phase 4 | Already strong | No change required | No | Not parsed (Rscript unavailable) |
| 24_vf_longitudinal_dynamics.R | Phase 4 | Already strong | No change required | No | Not parsed (Rscript unavailable) |
| 25_vf_lineage_vf_interaction.R | Phase 4 | Already strong | No change required | No | Not parsed (Rscript unavailable) |
| 00_config.R | Helper/config | Already strong | No change required | No | Not parsed (Rscript unavailable) |
| R/clinical_helpers.R | Helper/config | Already strong | No change required | No | Not parsed (Rscript unavailable) |
| R/plot_helpers.R | Helper/config | Already strong | No change required | No | Not parsed (Rscript unavailable) |
| R/wgs_helpers.R | Helper/config | Sparse header | Added full structured helper header for workflow role | No | Not parsed (Rscript unavailable) |
| 11_compare_strains_helpers.R | Helper/config | Basic header | Added full structured helper header and pipeline context | No | Not parsed (Rscript unavailable) |

## 3) Important scientific explanations added

- Clarified that multiple Phase 1b plots are **exploratory/descriptive** and are not independent-sample inference.
- Added explicit warnings that **presence/absence does not imply expression or causal mechanism** (pangenome/plasmid sections).
- Added lineage and repeated-measures caveats in ST-risk and plasmid-lineage scripts.
- Added longitudinal interpretation framing: graph-connected “Same strain” clusters are used to track within-host continuity.
- Added explicit canonical-entry rationale in `RUN_COMPLETE_ANALYSIS.sh` to support reproducible, auditable reruns.

## 4) Remaining concerns

- Environment currently lacks `Rscript`, so parse checks could not be executed in this runtime.
- Because parse checks were blocked, syntax validation of edited R files remains pending until run in the intended conda environment.
- Existing workflow still depends on availability of external tools/binaries (Panaroo, Parsnp, Abricate, etc.) as before.

## 5) Validation commands run

- `Rscript --version` → **failed**: `/bin/bash: Rscript: command not found`
- `bash -n RUN_COMPLETE_ANALYSIS.sh` → **passed**
- Intended but blocked due missing Rscript:
  - `Rscript -e "parse(file='03_plotting.R')"`
  - `Rscript -e "parse(file='07_explore_MLST.R')"`
  - `Rscript -e "parse(file='08_core_vs_plasmid.R')"`
  - `Rscript -e "parse(file='09_inc_plasmid_network.R')"`
  - `Rscript -e "parse(file='10_replicon_heatmap.R')"`
  - `Rscript -e "parse(file='11_compare_strains_helpers.R')"`
  - `Rscript -e "parse(file='12a_wgs_qc.R')"`
  - `Rscript -e "parse(file='12b_core_snp.R')"`
  - `Rscript -e "parse(file='12c_panaroo.R')"`
  - `Rscript -e "parse(file='13_visualise_panaroo_selection.R')"`
  - `Rscript -e "parse(file='15_longitudinal_patterns.R')"`
  - `Rscript -e "parse(file='16_within_host_evolution.R')"`
  - `Rscript -e "parse(file='17_lineage_analysis.R')"`
  - `Rscript -e "parse(file='18_annotate_variants.R')"`
  - `Rscript -e "parse(file='19_host_context.R')"`
  - `Rscript -e "parse(file='20_variant_annotation_deep.R')"`
  - `Rscript -e "parse(file='21_publication_figures.R')"`
  - `Rscript -e "parse(file='R/wgs_helpers.R')"`
