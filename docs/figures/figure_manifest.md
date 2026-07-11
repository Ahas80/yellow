# Current Figure Manifest

All listed figures are current UTI-vs-Not_UTI outputs unless explicitly labelled legacy.

| Key Figure | Figure ID | Filename | Script | Description |
| :---: | :--- | :--- | :--- | :--- |
| ✅ | F001 | `plots/clinical/trajectories_heatmap.png` | `00c_plot_clinical_summary.R` | Longitudinal primary UTI status heatmap per participant. |
| ✅ | F002 | `plots/clinical/transitions_alluvial_or_heatmap.png` | `00c_plot_clinical_summary.R` | Transitions between primary `Not_UTI` and `UTI` states. |
| ✅ | F003 | `plots/clinical/waterfall_counts.png` | `00c_plot_clinical_summary.R` | Catheter-aware UTI definition waterfall. |
| ✅ | F004 | `plots/clinical/uti_reclassification_heatmap.png` | `00c_plot_clinical_summary.R` | Legacy ASB / UTI / Negative rows reclassified into primary UTI or Not_UTI subgroups. |
| ✅ | F005 | `plots/clinical/not_uti_subgroup_by_batch_event.png` | `00c_plot_clinical_summary.R` | Composition of heterogeneous `Not_UTI` comparator by batch and sampling context. |
| ✅ | F006 | `plots/clinical/uti_symptom_rule_provenance.png` | `00c_plot_clinical_summary.R` | Catheter-rule and symptom-rule provenance for primary UTI classification. |
| ✅ | F007 | `plots/clinical/uti_cfu_threshold_provenance.png` | `00c_plot_clinical_summary.R` | CFU threshold/source provenance for primary UTI classification. |
| ✅ | F008 | `plots/vf/vf_burden_by_status.png` | `23_vf_cross_sectional.R` | Exploratory VF burden by primary UTI status. |
| ✅ | F009 | `plots/vf/vf_gene_prevalence_difference_uti_not_uti.png` | `23_vf_cross_sectional.R` | Gene prevalence differences for UTI vs Not_UTI. |
| ✅ | F010 | `results/models/plots/volcano_plot_UTI_vs_Not_UTI.png` | `14_genotype_phenotype_model.R` | Exploratory genotype-phenotype association plot for primary contrast. |
| ✅ | F011 | `plots/publication/Fig1_Swimmer_Plot.png` | `21_publication_figures.R` | Publication swimmer plot using primary UTI status. |
| ✅ | F012 | `plots/publication/Fig2_Mutation_Map.png` | `21_publication_figures.R` | Mutation map for primary-status switch candidates, if variants are available. |

Legacy ASB-vs-UTI files, if present, belong under `plots/legacy/old_asb_uti_outputs/` or `results/legacy/old_asb_uti_outputs/` and should not be used as current figures.
