# Publication-Grade Code and Research-Logic Audit

**Project:** YELLOW RoUTIne / rUTIs analysis  
**Audit date:** 2026-07-11 (Europe/Amsterdam)  
**Review type:** findings-only, read-only review of code, data flow, generated artifacts, Git history, and research logic  
**Publication verdict:** **NOT READY for primary-result publication in the present state**  
**Permitted use now:** the explicitly supported, claim-limited results in the stop/go table below; all other genomic headline claims remain exploratory, superseded, withdrawn, or unverifiable  
**Files changed by this audit:** this README only

## 1. Executive verdict

The clinical counts are reproducible, and the current Longcycler selection contract is internally coherent. The repository as a whole is not publication-ready because its active source now requires a **532-genome Longcycler-only** analysis, while most core-genome, pangenome, VF, model, pairwise, transition, summary, and figure artifacts were generated from the earlier **556-genome mixed Longcycler/Flye** analysis. The repository therefore contains two incompatible definitions of “current.”

The scientific primary question in the published protocol is also not the cross-sectional `UTI` versus heterogeneous `Not_UTI` comparison currently implemented. The protocol asks whether time-varying ASB—especially low-virulence *E. coli* ASB—predicts subsequent UTI during follow-up. The available isolate-linked dataset can support descriptive and exploratory *E. coli* analyses, but it cannot by itself estimate the protocol's full-cohort ASB-to-subsequent-UTI estimand.

“Perfect sense” is therefore interpreted here as internally consistent, protocol-aware, source-supported, reproducible, and appropriately limited. It is not a guarantee of scientific certainty.

### Claim-specific stop/go table

| Claim or artifact | Decision | What is supportable now | Publication gate |
|---|---|---|---|
| Operational clinical denominator: 583 included episodes, 166 participants, 18 UTI and 565 Not_UTI | **GO with qualification** | These counts are independently reproduced under `catheter_adjusted_sns_cfu1e3_v1`. Call this the **implemented operational phenotype**, not the protocol phenotype. | Resolve the phenotype departures in P0-2 before calling it protocol-defined UTI. |
| Current analysis assembly manifest: 532 Longcycler genomes, 161 participants, 16 UTI and 516 Not_UTI | **GO** | The manifest is internally valid, unique by participant-timepoint, QC-passing, file-backed, and Longcycler-only. | Record content hashes and regenerate every downstream primary artifact from this exact manifest. |
| Current Longcycler pair universe: 893 all within-participant pairs | **GO as a reconciled denominator only** | The number follows exactly from `sum(n_i choose 2)` and the 893 Longcycler-only rows inside the older pair table match that universe. | Write a fresh official 893-row pair table and pass the final manifest verifier. |
| 371 Longcycler adjacent transitions; 140 at operational `dnadiff` <=25 SNPs; 5/9 Not_UTI-to-UTI at <=25 | **GO as exploratory operational results** | Independently reproduced with SHA-bound `dnadiff` provenance. “<=25” is an operational assembly-to-assembly boundary, not method-equivalent validation of same strain. | Fix missing-WGS/direct-pair logic, relabel the former “sensitivity” artifact as current primary if that remains the policy, and retain the method caveat. |
| Earlier 556 genomes, 963 all-pairs, 394 mixed transitions, 138 transitions at <=25 | **SUPERSEDED / sensitivity only** | Reproducible as the prior mixed-assembler analysis. It is not current-canonical after the 532-row manifest was written. | Preserve under a clearly labelled historical/sensitivity namespace; do not merge with current-primary outputs. |
| Eight current phenotype-switch candidates | **GO for candidate selection only** | All eight have direct Longcycler-to-Longcycler comparisons, `Classification == "Same"`, and <=25 SNPs (maximum 20). | Do not attach the current mutation/gene annotations; regenerate variants by current hash-bound files and contig-aware coordinates. |
| Current mutation loci or mechanisms, including named `rpoD`/`lpxL` stories | **STOP / withdraw** | No locus-specific mechanism is currently supportable. | Fix P0-5, rerun annotation, validate manually, and use non-causal wording. |
| Core-SNP tree/distances and Panaroo outputs as current | **STOP** | Historical mixed-556 outputs only. Existing GREEN reports refer to the old manifests. | Regenerate from the declared policy, require complete outputs and content hashes, and connect the branch to a stated estimand. |
| VF prevalence, burden, module, diagnostic, or genotype-phenotype estimates from 556 rows | **STOP as current primary; exploratory historical sensitivity only** | The old tables can document what the former analysis produced. They are stale under the new manifest. | Reprofile/rebuild the 532 rows and rerun every downstream table/model/figure. |
| “No FDR-significant VF association” | **EXPLORATORY, prior mixed analysis** | The old 556-row output has 0/119 univariable and 0/50 selected GLMM results at FDR <0.05. This is not evidence of equivalence or absence of association. | Rerun the 532-row analysis with a prespecified sparse-outcome/population-structure strategy and correct testing family. |
| Complete VF systems/UPEC modules | **STOP as complete-system claims** | The current module fields are “at least one component detected” summaries. | Define component-completeness rules or rename fields as partial module evidence. |
| AMR findings | **STOP / not performed** | Script 29 correctly records that dedicated AMR screening was not integrated. | Add a versioned AMRFinderPlus/ResFinder workflow and distinguish genotype from phenotype. |
| Plasmid carriage/location | **DESCRIPTIVE replicon evidence only** | PlasmidFinder-style replicon presence can be described. It does not prove a complete plasmid or locate VF/AMR genes on it. | Add contig/read-context validation if plasmid localization or transmission is claimed. |
| ASB protects against subsequent UTI, bacterial changes cause a phenotype switch, or a gene is a treatment target | **STOP / causal claim unsupported** | Hypothesis-generating language only. | Estimate the protocol-compatible longitudinal estimand with the required cohort data and use an explicit causal framework if a causal claim is intended. |

## 2. Frozen reviewed state and provenance

### 2.1 Authoritative freeze

- Branch: `codex-source-only-github-update`
- HEAD: `69c4eb2812d69071e3bfd955ecdf0209754eca75`
- HEAD subject/date: `Refactor application logic and update related tests`, 2026-07-11 04:29:28 +02:00
- Pre-review dirty inventory: 30 entries—28 modified and 2 untracked source/test files
- Pre-review status fingerprint: `03c87805a7e646dbbd1f486207b466102482afde493acde942a7dfe0608a95e2`
- Dirty diff at freeze: 28 tracked files, 1,419 insertions and 407 deletions
- `git diff --check`: PASS
- Source fingerprint: `fefb583a603238ad23af0efdd7250d3849defcb7ead2eba995f4ef1a363b8a5f`
- Fingerprint scope: 83 files—46 root numbered R/r scripts, 6 `R/` helpers, 26 top-level `scripts/` R files, 3 tests, and 2 shell runners

The planning snapshot had been based on `d14d0e5`; before the execution freeze, the branch advanced to `69c4eb2` and the Longcycler-only manifest was generated. This review uses the later state only. The state change is material reproducibility evidence: a manuscript analysis should be run from an immutable commit plus run-specific output directory, not a mutable shared results tree.

### 2.2 Pre-existing dirty-tree inventory

These changes pre-dated the audit README and were preserved without modification:

```text
 M 00_config.R
 M 00c_plot_clinical_summary.R
 M 02_gene_presence_analysis.R
 M 06_MLST.R
 M 07_explore_MLST.R
 M 11_compare_strains.R
 M 11_compare_strains_helpers.R
 M 12a_wgs_qc.R
 M 12b_core_snp.R
 M 22_vf_build_analysis_dataset.R
 M 29_vf_amr_combined_profile.R
 M 30_vf_project_summary_tables.R
 M 32_uti_not_uti_diagnostic_stats.R
 M 34_robustness_first_addon.R
 M 35_final_figure_pack.R
 M 36_statistical_sensitivity_addon.R
 M R/pipeline_qc_helpers.R
 M R/provider_mlst_integration.R
 M R/wgs_helpers.R
 M RUN_COMPLETE_ANALYSIS.sh
 M scripts/audit_uti_status_count_explanation.R
 M scripts/compare_mlst_sources.R
 M scripts/create_workflow_case_count_flowchart.R
 M scripts/rebuild_longcycler_sensitivity.R
 M scripts/run_local_mlst_deprecated.R
 M scripts/verify_mlst_source_usage.R
 M scripts/verify_uti_not_uti_alignment.R
 M tests/test_pairwise_cache_provenance.R
?? scripts/verify_longcycler_only_pipeline.R
?? tests/test_longcycler_only_selection.R
```

### 2.3 Review scope

Primary review scope comprised active numbered scripts `00` through `36`, shared `R/` helpers, `RUN_COMPLETE_ANALYSIS.sh`, `RUN_IN_TERMINAL.sh`, all top-level invoked audit/reconciliation R scripts, three tests, current results, and the source-facing documentation. Legacy/generated material was read only to trace superseded logic and withdrawn claims. No raw row-level participant data are reproduced in this README.

### 2.4 Key artifact timestamps and hashes

The timestamps show the incompatible generations directly.

| Artifact | Modified | SHA-256 | Interpretation at freeze |
|---|---|---|---|
| `results/clinical/status_map.csv` | 2026-05-29 23:54:58 +02:00 | `e7fd159d5ce3f413fd835c63745d38a3b73e289ef03dd17b644e96466c2b8455` | Current operational clinical map |
| `results/qc/canonical_assembly_selection.csv` | 2026-07-11 04:50:20 +02:00 | `1da1c6a7d0fbd2573b2b7fa01cf39b58e546a579cbf34099671903a0c52eb050` | Current 532-row Longcycler selection |
| `results/qc/analysis_assembly_manifest.csv` | 2026-07-11 04:50:20 +02:00 | `fb35354cb5253a21853fbf1801d5d2151f4ea68470a9f6445ccb36eb5ebe9970` | Current analysis contract |
| `results/wgs/core/core_snp_input_manifest.csv` | 2026-05-29 23:55:52 +02:00 | not current-manifest bound | Stale mixed 556 input |
| `results/wgs/pan/panaroo_input_manifest.csv` | 2026-05-29 23:55:55 +02:00 | not current-manifest bound | Stale mixed 556 input |
| `results/vf/vf_analysis_ready.csv` | 2026-05-30 00:12:12 +02:00 | `21a049a4661f26d5453d9aeb4b7ac63ac25a1fa32d6512c8d087a7363135a47d` | Stale mixed 556 rows |
| `results/strain_compare/pairwise_metrics.csv` | 2026-07-10 02:20:58 +02:00 | `3efe7de83e0ac697c2439afe9fc10c5db7b641d0eefec612745ab73ef5346b19` | SHA-audited but mixed 963 rows |
| `results/vf/vf_longitudinal_transitions.csv` | 2026-07-10 02:25 +02:00 | `5caf69379373fd14c09a07423ef09ac02bb04fc9865e7fa8757e0c85ba97d321` | Prior mixed 394 adjacencies |
| `results/sensitivity/longcycler_only/longcycler_transitions.csv` | 2026-07-10 02:26:26 +02:00 | `d0bddb265af2796aaea22df865a7226ad8287ded27e7a93e1e60ee76f0a623e6` | Numerically current with 532-row policy; hierarchy label stale |
| `results/longitudinal/annotated_snps.csv` | 2026-05-30 00:13:30 +02:00 | `5e28e6ada6352acec51d47b309252b84399d71051b6d71e90c0e901dfe614d80` | Withdrawn/stale variant source |
| `results/longitudinal/variant_annotation_detailed.csv` | 2026-05-30 00:13:40 +02:00 | `42f8bc5b3290b73a73051eab93bad092618d42d0c56c404d5e9fb9ff2393da54` | Withdrawn/stale and contig-unsafe |
| `results/summary/final_key_results_summary.md` | 2026-07-10 02:54:20 +02:00 | `37b0e8d67c738168d049c88c4ad14671a9175d1d665a57edb7b143c851e4c8ec` | Clinical lines usable; WGS/VF/core/Panaroo lines stale |

### 2.5 Available environment provenance

Current shell provenance is not proof of the environment that generated the May/July artifacts, because output-level environment snapshots were not stored.

- Active Conda environment: `base`
- R: 4.3.3
- Key R packages: `dplyr` 1.1.4, `readr` 2.1.5, `testthat` 3.2.3, `lintr` 3.2.0, `lme4` 1.1-37, `broom.mixed` 0.2.9.6, `future` 1.40.0, `furrr` 0.3.1, `igraph` 2.1.4, `digest` 0.6.37, `jsonlite` 2.0.0
- Available commands: ABRicate 1.0.1, DNAdiff 1.3, Mash 2.2.2, MLST 2.19.0, snp-dists 0.7.0, Nucmer/MUMmer
- Missing from the current PATH: `parsnp`, `panaroo`, `prokka`, `plasmidfinder`
- Current ABRicate database listing: VFDB 2,597 sequences, ResFinder 3,077, PlasmidFinder 460, NCBI 5,386; database date 2025-10-28
- Missing reproducibility files: `install_r_packages.R`—despite `README.md:19-26` requiring it—and `renv.lock`
- Environment mismatch: `README.md:9-17,44-51` names `yellow-wgs-x86`; `RUN_COMPLETE_ANALYSIS.sh:17-24` merely accepts any active environment and instructs `asm-snp-x86`; `env-wgs.yml:5-51` omits several tools the runner invokes

## 3. Reconciled source of truth

Units are deliberately separated; assembly candidates, episodes, all pairs, and adjacent transitions are not interchangeable.

| Layer | Unit | Rows | Participants | Primary-status detail | Status |
|---|---|---:|---:|---|---|
| Classified clinical map | Clinical episode | 585 | 167 | 18 UTI, 567 Not_UTI before primary exclusions | Audit universe only |
| Included clinical analysis | Clinical episode | 583 | 166 | 18 UTI, 565 Not_UTI | Current operational clinical denominator |
| WGS candidate table | Assembly candidate | 1,291 | 166 | Assembler alternatives are separate rows | Current QC audit universe |
| Current selected WGS | Selected participant-timepoint genome | 532 | 161 | 16 UTI, 516 Not_UTI after joining status | Current Longcycler-only contract |
| Prior mixed WGS/VF | Selected participant-timepoint genome/VF row | 556 | 162 | 17 UTI, 539 Not_UTI | Superseded primary; retain as mixed sensitivity/history |
| Current all-pair universe | Unordered within-participant genome pair | 893 | 139 with repeated WGS | 277 have <=25 `dnadiff` SNPs in the current-manifest subset of the older table | Denominator reconciled; fresh official table missing |
| Prior mixed all-pair table | Unordered within-participant genome pair | 963 | 144 | 70 rows involve a now-nonmanifest Flye endpoint | Superseded mixed sensitivity |
| Current adjacent Longcycler transitions | Adjacent retained participant-timepoint pair | 371 | 139 | 140 at <=25; 9 Not_UTI-to-UTI, of which 5 at <=25 | Supported exploratory operational result |
| Prior mixed adjacent transitions | Adjacent selected participant-timepoint pair | 394 | 144 | 138 at <=25 | Superseded mixed sensitivity |

The 371 current adjacencies are not obtained by simply filtering 394 rows. Of the old mixed adjacencies, 362 retain two current Longcycler endpoints; removing intervening Flye-only timepoints creates 9 new Longcycler-to-Longcycler adjacencies, yielding 371. This is recorded at `results/sensitivity/longcycler_only/denominator_counts.csv:18-24` and `longcycler_sensitivity_summary.txt:9-13`.

## 4. Prioritized findings

## P0 — blocks or invalidates affected results

### P0-1. Active source and generated genomic outputs use incompatible assembly policies

**Evidence.** The active policy is `ANALYSIS_ASSEMBLER <- "longcycler"` with no fallback (`00_config.R:166-173`). Selection eligibility and canonical choice enforce that contract (`R/pipeline_qc_helpers.R:1544-1585`), and manifest validation requires selected, QC-passing, existing, unique Longcycler rows (`R/pipeline_qc_helpers.R:434-513`). The current denominator log records 532 rows at `results/qc/pipeline_denominator_summary.csv:28-29`. The same file records core SNP, Panaroo, VF, model, module, and score branches at 556 at lines 16-25. The old core report still claims GREEN for 556 rows (`results/wgs/core/core_snp_staleness_report.txt:1-9`), while the prospective current core-manifest hash is `357d19743d568346f63058e54c9618db`, not stored `f1425fc48328a50024a4c551846e77b2`.

**Reproduced consequence.** The 556-row artifacts contain all 532 current rows plus 24 former Flye fallback rows. The 963-row pair table contains exactly 893 current-manifest pairs plus 70 pairs involving 23 Flye endpoints. VF-ready has 24 rows outside the current manifest. The strict final verifier has no generated report and would fail its 532-row equality checks (`scripts/verify_longcycler_only_pipeline.R:48-142`).

**Affected claims/outputs.** Current core tree/distances, Panaroo, VF/MLST integrations, genotype-phenotype models, mixed pairwise totals, mixed transition totals, summary tables, workbooks, and figures. The independent `dnadiff` values for the 893 retained Longcycler pairs are not invalidated by the 70 extra rows, but the official file is not a current-primary artifact.

**Scientific rationale.** Assembler selection changes eligibility, missingness, allele/gene detection, pair adjacency, and outcome composition. The Landman benchmark supports Longcycler for molecular typing but recommends Flye for AMR gene and plasmid-replicon detection, so a universal assembler policy is a scientific decision, not merely file housekeeping.

**Remediation.** Freeze the intended policy in a versioned analysis plan. If the present source is authoritative, regenerate every active downstream branch from the 532-row manifest in a new run directory; preserve the mixed 556 analysis as a labelled sensitivity. Do not overwrite historical artifacts in place.

**Verification test.** Require exact equality of episode keys and normalized FASTA paths between every active genomic artifact and `analysis_assembly_manifest.csv`; require pair count `sum(n_i * (n_i - 1) / 2) == 893`; require adjacency count `sum(n_i - 1) == 371`; then run the final Longcycler verifier and fail on any mismatch.

### P0-2. The implemented UTI phenotype is not the published protocol phenotype

**Evidence.** Suprapubic pain is detected but explicitly marked descriptive and excluded from the rule (`R/clinical_helpers.R:276-290`). The non-catheter local rule uses dysuria, urgency, frequency, incontinence, and pus, plus flank pain with systemic symptoms (`R/clinical_helpers.R:424-465`). No recent-onset or alternative-infectious-focus field exists in `status_map.csv`. Any episode not satisfying both culture and symptom rules becomes binary `Not_UTI`, including unknown/indeterminate inputs by construction (`R/clinical_helpers.R:470-508`; `00b_classify_episodes.R:9-18`). Episode aggregation uses `any()` across source rows (`00b_classify_episodes.R:62-66,255-339`), allowing culture and symptoms from different source rows to be synthesized into one episode.

The published [YELLOW RoUTIne protocol](https://pmc.ncbi.nlm.nih.gov/articles/PMC11363575/) is the prespecified benchmark. Its clinical framework includes suprapubic pain, recent-onset symptoms, exclusion of another infectious focus, catheter-aware assessment, and specimen/collection handling. IDSA, revised McGeer, and Loeb are useful sensitivity/context definitions, not replacements for that protocol.

**Reproduced consequence.** Adding suprapubic pain to the existing non-catheter local-symptom rule—without changing any other condition—identifies two included, culture-supported, suprapubic-only episodes currently labelled Not_UTI. The conditional counterfactual would change clinical UTI from 18 to 20 and current Longcycler-linked UTI from 16 to 18. This is not an adjudicated reclassification; it quantifies the effect of one omitted protocol element. All 18 current operational UTIs meet `>=10^4` as well as `>=10^3`, so the universal lower threshold alone does not explain the present UTI count. The included frozen dataset happens to have no `unknown_or_indeterminate` rows, but the code would force future unknowns into controls.

**Affected claims/outputs.** Any statement that the 18/565 phenotype is protocol-defined; every status-stratified clinical, VF, model, lineage, and transition output; sensitivity/specificity language about `Not_UTI`.

**Scientific rationale.** In a population with prevalent bacteriuria, small symptom-definition differences materially alter case/control classification. `Not_UTI` is not synonymous with ASB, and absence of documented evidence is not evidence of absence.

**Remediation.** Obtain clinical/protocol-owner adjudication and implement an explicit versioned decision table containing: symptom recency, suprapubic pain, alternative focus, collection method, specimen-specific culture interpretation, catheter status, and a retained indeterminate state. Keep protocol, IDSA, McGeer, and Loeb classifications in separate named columns.

**Verification test.** Add unit fixtures for suprapubic-only culture-supported episodes, old/non-recent symptoms, another infectious focus, catheter cases, each specimen threshold, conflicting source rows, and missing inputs. Assert that unknown stays indeterminate and that protocol/sensitivity definitions never silently overwrite one another.

### P0-3. The available isolate dataset cannot estimate the protocol's primary longitudinal estimand

**Evidence.** The protocol's main hypothesis asks whether low-virulence *E. coli* ASB reduces subsequent UTI risk during follow-up and prespecifies time-varying analyses with repeated visits nested in residents and nursing homes. The repository's principal model is a contemporaneous `UTI` versus `Not_UTI` isolate-profile comparison (`14_genotype_phenotype_model.R:510-540`). The available analysis table has no nursing-home/facility identifier, only 13 participants with both genomic primary states, and strong sampling-occasion structure: in the older 556-row VF table, 16/17 UTI rows are UTI-event samples, while 521/539 Not_UTI rows are routine samples.

**Reproduced consequence.** The 9 adjacent Longcycler Not_UTI-to-UTI transitions are not the protocol estimand: `Not_UTI` includes more than ASB; adjacency among available isolates is not a prespecified three-month risk window; residents without an *E. coli* isolate are absent from the genomic exposure set; and facility-level clustering cannot be fitted. UTI incidence, ASB exposure prevalence, no-ASB comparisons, and protective effects are therefore not identifiable from this dataset alone.

**Affected claims/outputs.** “ASB protects,” “low virulence reduces UTI risk,” protocol objective 2C, full-cohort UTI incidence, and causal wording. The current cross-sectional and transition analyses remain possible as explicitly secondary/exploratory analyses.

**Scientific rationale.** Conditioning on successful culture, *E. coli* selection, sequencing, assembly QC, and available adjacent isolates changes the target population and can create selection bias. A contemporaneous bacterial comparison is not a prospective risk model.

**Remediation.** Link the full resident-visit cohort, including visits without ASB, without *E. coli*, and without WGS; define ASB and low-virulence exposure before outcome; construct the prespecified follow-up window; include death/discharge/loss-to-follow-up handling; and model visits within residents within nursing homes. If those data are unavailable, state that protocol objective 2C is not answered.

**Verification test.** Produce a target-trial/estimand table specifying population, time zero, exposure, comparator, outcome window, censoring, clustering, and missingness. Require every eligible resident-visit to appear exactly once in an at-risk dataset before fitting the model.

### P0-4. Longitudinal graph logic turns missing WGS and transitive links into false strain states

**Evidence.** Script 15 builds graph edges from the composite `Classification == "Same"` rule (`15_longitudinal_patterns.R:68-86`), makes every clinical episode a graph vertex (`:88-123`), and equates component equality with adjacent same strain (`:181-199`). It does not call the primary-analysis filter after reading the status map (`:54-60`). The composite classification itself treats missing accessory Jaccards as passing and can call “Same” without SNP/identity evidence when ST/accessory evidence permits (`11_compare_strains_helpers.R:569-605`).

**Reproduced consequence.** The frozen timeline reports zero missing `Strain_ID`, although 53/585 timeline episodes are absent from the current Longcycler manifest. Among 418 adjacent timeline rows, 63 have at least one current-manifest-missing endpoint, and all 63 are labelled `Strain_Replacement`. The graph contains 22 directly compared pairs above 25 SNPs—maximum 139—inside the same connected component; five adjacent rows are called same strain despite a direct distance above 25. It also includes two rows excluded from the 583-row primary clinical denominator.

**Affected claims/outputs.** Persistence, replacement, and clearance counts; `participant_timelines.csv`; `transitions.csv`; swimmer/timeline summaries; any claim that every episode has strain evidence. The eight current phenotype-switch candidates are protected from this particular defect because each has a direct Longcycler comparison at <=25 SNPs.

**Scientific rationale.** Pairwise relatedness is not necessarily transitive at a fixed distance threshold. A graph component means “connected by a chain,” not “every pair satisfies same-strain criteria.” Missing sequence cannot demonstrate replacement.

**Remediation.** Filter primary clinical rows explicitly. Represent missing sequence as `Missing_WGS`, never a singleton biological strain. For adjacent transitions, use the direct pair's SNP evidence first; use graph components only as separately labelled connected-lineage context, or impose a cluster-diameter/complete-linkage constraint.

**Verification test.** Assert that an absent endpoint yields `Missing_WGS`; a directly compared pair above 25 never yields strict same strain; excluded rows never enter primary timelines; and every component's maximum direct distance is reported. Add an A-B<=25, B-C<=25, A-C>25 counterexample.

### P0-5. Variant and locus annotations are stale, hash-disconnected, and unsafe for multi-contig references

**Evidence.** SHA-bound `dnadiff` prefixes and sidecars are constructed in `11_compare_strains_helpers.R:375-403,463-499`. Script 16 calls `run_dnadiff()` (`16_within_host_evolution.R:91-96`) but its result table drops the generated report/SNP paths (`:142-160`). Script 18 reconstructs a legacy unsuffixed `key.snps` path (`18_annotate_variants.R:101-115`) and keeps only the first four SNP columns, discarding reference/query contig tags (`:76-96`). Script 20 reads GFF `seqid` (`20_variant_annotation_deep.R:207-223`) but matches only numeric coordinate (`:249-254`). The runner executes script 20 at `RUN_COMPLETE_ANALYSIS.sh:169-171`, before preferred case/summary sources from scripts 28 and 30 at lines 211-213 and 235-237.

**Reproduced consequence.** Eight current hashed candidate SNP files and eight older unsuffixed counterparts exist. One counterpart differs from the current file by nine substitutions. Three of eight current comparisons have multiple reference contigs, up to four. `annotated_snps.csv` and `variant_annotation_detailed.csv` are May artifacts, while the current candidate set and SHA caches are July artifacts. Numeric positions can therefore be assigned to the wrong contig/gene.

**Affected claims/outputs.** All named mutations, gene/locus annotations, “mechanism” stories, mutation-map figures, and claims that a phenotype switch was caused by a particular genomic change.

**Scientific rationale.** A coordinate is unique only within a sequence/contig and a specific reference build. Provenance must bind the variant to both input assemblies, comparison direction, contig, and annotation version.

**Remediation.** Persist the exact `.snps` path and both FASTA hashes from script 16; parse and retain reference/query contig identifiers; join GFF by `(reference_fasta_hash, seqid, coordinate)`; run annotation only after candidate sources exist; archive legacy unsuffixed caches.

**Verification test.** Use a fixture with two contigs sharing the same coordinate and assert the correct gene is selected. Assert every annotated variant traces to a current sidecar and FASTA hashes, and that no legacy unsuffixed file is read.

### P0-6. `vf_jaccard` is silently erased in the transition case summary

**Evidence.** `28_vf_transition_case_studies.R:507-528` initializes a `vf_jaccard` column to `NA`. A scalar with the same name is computed at line 542, then `mutate(vf_jaccard = vf_jaccard)` at lines 576-586 resolves the data-mask column to itself instead of the environment scalar.

**Reproduced consequence.** All 417 rows in `results/vf/vf_transition_case_summary.csv` have missing `vf_jaccard`; 386 rows have non-missing `module_jaccard`. This is a silent wrong-result bug, not merely missing documentation.

**Affected claims/outputs.** Transition case tables, casebooks, and any case-specific VF-stability statement sourced from this field. The separate Jaccard calculation in script 24 is not invalidated by this assignment bug, but remains stale under P0-1 until rerun.

**Remediation.** Use an unambiguous scalar name or `.env$vf_jaccard_value`, add a non-missing assertion when both VF endpoints exist, and fail rather than publish an all-NA metric.

**Verification test.** A two-vector fixture with known intersection/union must produce the expected scalar in the written summary; assert `has_vf_pair => !is.na(vf_jaccard)`.

### P0-7. “Current/final” documentation mixes stale, withdrawn, and active claims

**Evidence.** `README.md:137-148` still reports a highly significant long-polar-fimbriae hit, two exact phenotype-switch mechanisms, and named mutation stories. `docs/FINAL_SUMMARY.md:16-55,59-99` reports 361 genomes, 522 pairs, significant genes, causal/treatment implications, and publication readiness. `docs/DOCUMENTATION_INDEX.md:73-87` reports 274 episodes and 42 UTI. `results/KEY_FINDINGS.md:14-23` and `results/ANALYSIS_README.md:15-19` report the obsolete 579/12/557 denominators. Even `results/summary/final_key_results_summary.md:10-21,36-40,54-64` calls the 556-row layer current and says the old core/Panaroo reports are GREEN.

**Reproduced consequence.** A reader can select mutually incompatible “current” denominators and conclusions from repository-root, docs, results-root, and results-summary files. The old model has no FDR-significant feature, so the root significant-lpf claim is directly contradicted. Old Longcycler handouts report 116/371 and 7/9 (`outputs/longcycler_only_methods_summary/Longcycler_only_methods_handout.md:20-30`), versus current verified 140/371 and 5/9.

**Affected claims/outputs.** Manuscript text, presentations, supervisor handouts, thesis briefs, publication templates, and automated summaries.

**Scientific rationale.** Claim provenance is part of reproducibility. A correct pipeline cannot rescue a manuscript assembled from superseded artifacts.

**Remediation.** After the full rerun, generate one run-scoped claim registry with denominator, unit, manifest hash, source file, status, and approved wording. Move obsolete narratives to a dated archive and remove “current/final” labels from historical files.

**Verification test.** Search all non-archive documentation for forbidden old anchors/claims and fail CI unless each occurrence is explicitly labelled historical/withdrawn. Require every headline number to resolve to one current run manifest.

## P1 — material bias or reproducibility risk

### P1-1. The <=25 SNP threshold is borrowed from a non-equivalent method

**Evidence.** `SAME_STRAIN_SNP_THRESHOLD <- 25` is set at `00_config.R:206-208`. The nursing-home study used Illumina reads, Snippy, isolates of the same ST, ST-specific references, and episodes at least 30 days apart; it considered <=25 SNPs same strain. The current project applies the number to assembly-to-assembly MUMmer/DNAdiff SNP totals across adjacent retained observations.

**Consequence.** The 140/371 and 5/9 figures are reproducible, but their biological calibration is not validated by the cited study. Parsnp core-SNP distances and `dnadiff` totals must not be presented as interchangeable.

**Affected claim.** “Same strain” as a definitive biological label. Use “meets the operational <=25 DNAdiff-SNP rule” unless validated.

**Source.** [Recurrent *E. coli* UTIs in nursing homes](https://pmc.ncbi.nlm.nih.gov/articles/PMC9686610/).

**Remediation/test.** Validate the threshold against a read-based, recombination-aware or ST-stratified reference workflow on a representative subset; report sensitivity across thresholds and methods; assert method name beside every threshold result.

### P1-2. WGS QC does not implement all configured quality dimensions, and assembler choice is not task-neutral

**Evidence.** The configured QC includes completeness >=95% and contamination <=5% (`R/wgs_helpers.R:367-375`), but script 12a applies only contig count, N50, total size, and read errors (`12a_wgs_qc.R:134-160`). Raw-read coverage, basecalling model/version, polishing, taxonomic contamination, completeness, and provider laboratory provenance are not available. The current selection excludes every QC-passing Flye fallback (`R/pipeline_qc_helpers.R:1552-1585`).

**Consequence.** Genome inclusion is conditional on a limited assembly-QC screen. The current QC report shows 49 Not_UTI and 2 UTI episodes without a passing selected Longcycler genome (`results/qc/qc_selection_bias_report.txt:4-11`). The Longcycler policy is plausible for molecular typing but may be suboptimal for AMR/plasmid tasks.

**Source.** The [Landman long-read benchmark](https://pmc.ncbi.nlm.nih.gov/articles/PMC11587635/) recommends at least 40x coverage and Miniasm/Longcycler for molecular typing, but Flye for AMR gene and plasmid-replicon detection.

**Remediation/test.** Add read/assembly provenance and implemented completeness/contamination metrics; justify a task-specific assembly policy; compare key VF/AMR/plasmid endpoints across assemblers. Assert that every configured QC field is calculated or removed from configuration with rationale.

### P1-3. Sparse-outcome modelling, outcome-driven selection, and population/sampling structure limit inference

**Evidence.** The old model has 17 UTI among 556 profiles and warns that fewer than 20 events are exploratory (`14_genotype_phenotype_model.R:518-533`). Fisher screens ignore participant clustering (`:538-540`). GLMM features are chosen using outcome-dependent `p <0.1` or top 50 (`:618-623`), then BH correction is applied only to the selected GLMM results (`:786-803`). Models adjust Timepoint and Batch but deliberately omit ST/population structure (`:719-729`). Singular fits are accepted, and failed GLMMs fall back to unclustered GLM (`:733-763`).

**Reproduced consequence.** Existing mixed results: 119 univariable tests, minimum FDR 0.152, zero significant; 50 selected GLMM results, minimum FDR 0.112, zero significant; 6 singular fits and 37 sparse/separation flags. There were no GLM fallbacks in this particular output. Only 13 VF participants contribute both statuses, and sampling event strongly overlaps outcome.

**Scientific rationale/sources.** Bacterial GWAS must address clonal population structure; [pyseer](https://pmc.ncbi.nlm.nih.gov/articles/PMC6289128/) was designed for that problem. Sparse/separated logistic estimates can require penalized methods; [Heinze and Schemper](https://pubmed.ncbi.nlm.nih.gov/12210625/) describe Firth bias reduction for separation.

**Remediation.** Prespecify a small primary endpoint family; treat genome-wide screens as exploratory; correct across the full tested family rather than a data-selected subset; model sampling occasion/event, resident clustering, facility, and population structure; use a prespecified sparse-outcome method. Do not mix GLMM and GLM fallback estimates under one inferential label.

**Verification test.** Simulate null phenotypes on the observed lineage/event structure and assess calibration; require the exact multiplicity family and fallback policy before fitting; report events per parameter, convergence, singularity, separation, and effective within-resident information.

### P1-4. VF calling can silently lose zero-hit or failed assemblies, and its cache is not input-bound

**Evidence.** The cache name includes only FASTA basename and identity/coverage thresholds (`02_gene_presence_analysis.R:233-253`), not FASTA hash, ABRicate version, or VFDB version. `safely()` converts errors to `NULL`, and rows with `NULL` or zero hits are filtered before the P/A matrix is built (`:256-266,283-295`).

**Reproduced consequence.** The older output happens to include all 556 selected profiles and the 532 current keys, so no observed row is currently lost for zero hits. The code nevertheless cannot distinguish a true all-zero genome from tool failure and can silently change denominators after a database/tool/input update.

**Source.** [VFDB 2022](https://academic.oup.com/nar/article/50/D1/D912/6446532) is a curated virulence-factor knowledgebase; detection remains version- and threshold-dependent and does not prove expression or activity.

**Remediation/test.** Scaffold the P/A matrix from the manifest, left-join hits, encode verified zero-hit rows as all zero, and stop or quarantine tool failures. Bind caches to FASTA SHA-256, command, tool version, database name/version/hash, and thresholds. Test one zero-hit FASTA and one forced ABRicate failure.

### P1-5. Core-SNP cache, command, completeness checks, and downstream role need repair

**Evidence.** Script 12b fingerprints path, size, and mtime (`12b_core_snp.R:74-89`); the helper MD5s that metadata table rather than FASTA contents (`R/pipeline_qc_helpers.R:247-258`). A matching hash skips even when `FORCE_RERUN_CORE_SNP=1` (`12b_core_snp.R:117-120`). “Outputs exist” checks only tree and distance, not alignment or `strain_pairs.csv` (`:89`). Missing alignment or snp-dists can warn and finish successfully (`:216-282`). The command comments say `-n` is passed and `-x` filters duplicates, but the command omits `-n`; official Parsnp documentation defines `-x` as recombination detection/filtering (`:163-180`). No active primary R/shell consumer reads the generated `snp_dists.tsv`, `strain_pairs.csv`, or `core_genome.tree`; script 30 reads only the status text.

**Source.** [Parsnp/Harvest paper](https://pmc.ncbi.nlm.nih.gov/articles/PMC4262987/) and [official Harvest option documentation](https://harvest.readthedocs.io/en/v1.0/content/parsnp/tutorial.html).

**Remediation/test.** Use content hashes; honor true force; require alignment, distance, tree, pair table, command log, and manifest equality; invoke commands with argument vectors rather than shell strings; correct option documentation; declare the branch optional/report-only or connect it to a stated analysis. A test must mutate FASTA content without changing its path and prove cache invalidation.

### P1-6. Runner order and environment setup are history-dependent

**Evidence.** Script 20 runs before scripts 28 and 30 even though it preferentially consumes their outputs (`RUN_COMPLETE_ANALYSIS.sh:169-171,211-213,235-237`). The runner has `set -e` and `set -u` (`:9-10`) but only checks that some Conda environment is active (`:17-24`), not tool versions or the expected environment. The current shell lacks four required commands. The README refers to a missing installer and a different environment name.

**Consequence.** A clean run and an incremental run can use different variant candidate sources; the pipeline can start in an incomplete environment and fail late; published output cannot be reconstructed from a lockfile.

**Remediation/test.** Topologically order stages, add an exact preflight, pin Conda and R dependencies, record session/tool/database versions per run, and use run-scoped output directories. Execute a clean-room dry run that fails before any output write when a dependency is absent.

### P1-7. Module, AMR, and plasmid labels overstate what is measured

**Evidence.** A module is marked present when any one assigned component gene is detected (`26_vf_define_gene_modules.R:394-421`). Script 29 searches for dedicated AMR inputs and explicitly reports none (`29_vf_amr_combined_profile.R:94-126`); its output sets `true_amr_integration_performed = FALSE` (`:333-352`). Replicon detection is not plasmid reconstruction or gene localization.

**Sources.** Operational ExPEC/UPEC marker definitions require specified combinations, not arbitrary one-gene modules ([Tetzschner et al.](https://journals.asm.org/doi/10.1128/jcm.01269-20)). Dedicated AMR methods include [AMRFinderPlus](https://pubmed.ncbi.nlm.nih.gov/34135355/) and [ResFinder 4.0](https://pmc.ncbi.nlm.nih.gov/articles/PMC7662176/). [PlasmidFinder](https://pubmed.ncbi.nlm.nih.gov/24777092/) detects/characterizes replicon sequences.

**Remediation/test.** Rename partial modules or implement explicit component-completeness rules; preserve established ExPEC/UPEC definitions separately; add versioned AMR calling; describe plasmid outputs as replicon evidence unless contig/read context proves more. Unit-test full, partial, and absent operons.

### P1-8. CFU parsing does not support scientific notation despite its stated purpose

**Evidence.** The parser's output fields claim threshold-aware numeric bounds (`R/clinical_helpers.R:269-274`). Synthetic checks show `10^5`, `1 x 10^5`, and `1×10^5` are parsed with lower bound `1`, while `100000` is parsed correctly.

**Consequence.** No frozen status-map row uses caret or multiplication-sign scientific notation, so this does not change the current 18/565 result. It is a material future-data defect and can silently turn a high-count culture into below-threshold.

**Remediation/test.** Normalize multiplication/caret notation and add tests for exponents, inequalities, ranges, locale separators, text categories, and ambiguous inputs; ambiguous values must remain unknown rather than false.

## P2 — maintainability, documentation, or test-coverage concern

### P2-1. Output publication is not consistently atomic or failure-safe

`safe_write_csv()` writes a fixed `.tmp`, deletes the previous file, and ignores rename failure (`11_compare_strains_helpers.R:42-48`). Many scripts write final CSVs directly. A crash can leave a missing/partial current artifact, and concurrent runs can collide.

**Remediation/test.** Use unique same-directory temporary files, flush/close, check rename, preserve the old file until replacement succeeds, and write a completion manifest last. Force a rename failure and assert the prior artifact remains intact.

### P2-2. Tests cover important helpers but not the result-invalidating paths

The current suite passes 105 expectations across 20 test blocks and three files. It now covers clinical helpers, Longcycler selection, and pairwise cache provenance. Missing tests include: protocol phenotype fixtures, graph transitivity, missing WGS, primary filtering, output-manifest equality, zero-hit/error VF rows, variant hash/contig joins, runner dependency order, Jaccard assignment, true core force behavior, complete-output checks, sparse-model family/FDR behavior, and stale-claim scanning.

**Remediation/test.** Add small synthetic end-to-end fixtures that execute each critical branch and validate the written artifact, not only helper return values.

### P2-3. Lint debt obscures higher-value diagnostics

All 81 reviewed R files parse, but default `lintr` reports 8,285 diagnostics: 4,795 line-length, 1,778 indentation, 1,050 object-usage, 281 object-name, 229 commented-code, and smaller categories. Most are style, but the object-usage volume makes genuine name/data-mask mistakes—such as P0-6—harder to see.

**Remediation/test.** Adopt a staged lint baseline: first make object-usage, vector-logic, and parse-risk linters blocking; then mechanically format style debt. Do not mix mass formatting with scientific-logic changes.

### P2-4. Dead/disconnected and legacy stages need an explicit lifecycle

Core-SNP outputs are largely report-only, optional plotting expects a different path, and numerous “current/final” documents are historical. The runner still executes broad exploratory plotting and reporting stages without a machine-readable declaration of primary versus optional outputs.

**Remediation/test.** Maintain a stage registry with input contract, output contract, primary/secondary/legacy status, consumers, and deprecation date. Fail documentation generation when a listed primary artifact has no current producer or hash.

## 5. Verified logic and positive findings

The following logic was independently confirmed and should be retained while remediating the blockers:

- **Clinical denominator:** 585 classified rows/167 participants; after the two declared exclusions, 583/166 with 18 UTI and 565 Not_UTI.
- **Collection strata:** spontaneous 426 episodes (15 UTI), incontinence-material 153 (1 UTI), catheter 4 (2 UTI). These reproduce exactly, but the catheter stratum is too small for stable separate inference.
- **Current Longcycler selection:** 532 selected, QC-passing, existing FASTAs; 161 participants; no non-Longcycler row and no duplicated participant-timepoint key.
- **Selection bias report:** 516 Not_UTI and 16 UTI have selected passing Longcycler WGS; 49 Not_UTI and 2 UTI do not. The report correctly says genomics is conditional on available QC-passing WGS (`results/qc/qc_selection_bias_report.txt:4-11`).
- **Pairwise provenance design:** FASTA SHA-256, ordered pair inputs, cache signature, report SHA-256, tool command/version, and sidecar validation are implemented (`11_compare_strains_helpers.R:114-138,375-499`). The old mixed pair table has 963 unique unordered within-participant keys and no stale/missing SHA-bound cache row relative to its old 556 input manifest.
- **Current pair reconciliation:** filtering that table to exact current-manifest endpoints yields 893 rows, exactly the mathematically expected universe.
- **Longcycler adjacency rebuild:** 371 transitions from 139 participants; 362 retain old mixed adjacencies and 9 are newly adjacent after restriction; 140/371 meet the operational <=25 rule; 5/9 Not_UTI-to-UTI do so.
- **Current candidate selection:** all eight phenotype-switch candidates have a direct Longcycler pair, composite `Same` call, and <=25 SNPs; maximum 20. Candidate selection survives the graph audit, while locus annotation does not.
- **Existing association output is appropriately non-significant:** no FDR <0.05 result in either old mixed model table; warnings correctly flag low event count, singularity, and sparse/separation risk. Interpretation must remain exploratory and the tables must be regenerated.
- **AMR transparency:** script 29 explicitly records that true AMR integration was not performed rather than fabricating an AMR result.
- **Strict manifest verifier design:** `scripts/verify_longcycler_only_pipeline.R` checks manifest rows, assembler, key/path sets, pair count/endpoints, and transition count. It should become a non-writing pre-publication gate; the current frozen output state would fail.
- **Static validation:** 81 R files parsed; both shell runners pass `bash -n`; 105 test expectations pass with zero failures, errors, warnings, or skips; `git diff --check` passes.

## 6. Protocol reconciliation and estimand matrix

| Protocol element/research question | Implemented state | Classification | What may be claimed |
|---|---|---|---|
| UTI symptom set includes suprapubic pain | Detected but excluded from rule | **Defect** | Operational phenotype only; two-row sensitivity quantified, not adjudicated |
| Recent-onset/new symptoms | No field/rule | **Not verifiable / defect if called protocol-defined** | Cannot verify protocol phenotype |
| No alternative infectious focus | No field/rule | **Not verifiable / defect if called protocol-defined** | Cannot attribute nonspecific/systemic symptoms to UTI |
| Catheter-aware symptom rule | Separate systemic rule implemented | **Partially aligned** | Descriptive only; only four included catheter episodes |
| Separate collection-method strata | Method normalized and counts available; primary model pools them | **Undocumented adaptation** | Stratified descriptive counts; no stable catheter-specific effect |
| Specimen/culture handling | One primary >=10^3 support threshold plus coded fallbacks | **Undocumented adaptation requiring adjudication** | State exact implemented rule; all current operational UTIs also meet >=10^4 |
| Indeterminate phenotype | Code can label confidence indeterminate but binary outcome forces Not_UTI | **Design defect** | Frozen included data have no indeterminate row; future missingness is unsafe |
| Protocol ASB versus no-ASB/full cohort | `Not_UTI` combines bacteriuria and non-bacteriuria states | **Different estimand** | Do not call Not_UTI ASB |
| UTI incidence | Isolate-linked episode table, not verified full at-risk cohort | **Not estimable here** | No protocol incidence claim |
| Bacterial species dynamics | WGS analysis restricted to selected *E. coli* isolates | **Partial secondary analysis** | *E. coli*-conditional descriptions only |
| Relevant VF and AMR genes | VFDB presence available; dedicated AMR absent | **Partial / AMR not performed** | VF detection only, not expression/function; no AMR conclusion |
| Longitudinal VF/AMR change | VF adjacency analysis available; AMR absent; variants blocked | **Exploratory partial** | Current Longcycler VF stability after rerun; no AMR or mutation mechanism |
| Low-virulence *E. coli* ASB -> subsequent UTI | No validated low/high exposure; no full risk set/window/facility | **Not estimable here** | Hypothesis only |

### Clinical-context definitions

- [IDSA 2019 ASB guideline](https://www.idsociety.org/practice-guideline/asymptomatic-bacteriuria/) provides authoritative ASB context and cautions against treating bacteriuria without attributable UTI symptoms in older adults. It should define a labelled sensitivity/context state, not silently redefine the study phenotype.
- [Revised McGeer criteria](https://pmc.ncbi.nlm.nih.gov/articles/PMC3538836/) are surveillance definitions requiring clinical and microbiological evidence; they are suitable as a retrospective sensitivity definition.
- [Loeb minimum criteria](https://pubmed.ncbi.nlm.nih.gov/11232875/) concern the minimum criteria to initiate antibiotics in long-term care; they are not a substitute outcome definition for the protocol.

## 7. Current versus superseded artifact map

| Artifact/family | Status at freeze | Required handling |
|---|---|---|
| `results/clinical/status_map.csv` and clinical audit tables | Current operational phenotype | Retain; relabel as operational until P0-2 is adjudicated |
| `results/qc/canonical_assembly_selection.csv`, `analysis_assembly_manifest.csv` | Current Longcycler contract | Retain and content-hash; use as sole primary genomic input |
| `results/sensitivity/longcycler_only/longcycler_transitions.csv` | Numerically current for 532 policy; named/framed as sensitivity | Promote only after full rerun and graph fix; retain operational threshold caveat |
| `results/strain_compare/pairwise_metrics.csv` | Valid SHA-bound old mixed table; contains exact 893 current subset plus 70 obsolete rows | Archive mixed table; write fresh 893-row primary table |
| Core-SNP and Panaroo manifests/outputs | Stale mixed 556 despite GREEN text | Stop/withdraw as current; regenerate or mark report-only historical |
| `vf_pa_all.csv`, `vf_analysis_ready.csv`, model/module/score tables and related figures | Stale mixed 556 | Regenerate from 532 before primary use |
| `vf_transition_case_summary.csv` | Stale plus all-NA `vf_jaccard` defect | Withdraw and regenerate after code fix |
| `annotated_snps.csv`, `variant_annotation_detailed.csv`, mutation figures | Legacy-cache-derived and contig-unsafe | Withdraw named loci/mechanisms; regenerate from SHA-bound files |
| `results/strain_compare/archive/pre_sha256_cache_20260710T000634Z/` and old key-only caches | Historical invalid cache evidence | Archive only; never read by active workflow |
| `outputs/longcycler_only_methods_summary/` | Withdrawn values 116/371 and 7/9 | Archive with superseded banner; do not present |
| `results/thesis_audit/hamdi_thesis_discussion_brief_today.md` and its `.docx` companion | Superseded mixed values | Archive/regenerate |
| `README.md`, `docs/FINAL_SUMMARY.md`, `docs/DOCUMENTATION_INDEX.md` | Stale/withdrawn headline claims | Replace after validated rerun; retain history in archive only |
| `results/KEY_FINDINGS.md`, `results/ANALYSIS_README.md`, older report safety tables | Obsolete 579/12/557 generation | Archive; do not treat filename as authority |
| `results/summary/final_key_results_summary.md` | Mixed: clinical current, genomics stale | Do not cite as a whole; regenerate atomically from one run |
| May presentation families under `outputs/manual-20260527-current-review/` | Pre-SHA/current-policy history | Withdraw until regenerated |

## 8. Missing-provenance register

| Missing item | Why it matters | Status |
|---|---|---|
| Immutable analysis commit plus run ID/output directory | Prevents source and output generations from mixing | Missing |
| Full R lockfile and working package installer | Reconstructs statistical behavior | `renv.lock` and documented installer missing |
| Complete pinned Conda environment | Reconstructs bioinformatics commands | Current YAML incomplete/inconsistent |
| Per-output session/tool/database snapshot | Current shell versions do not prove generating versions | Missing for most outputs |
| Raw-read identifiers, coverage, chemistry, basecaller/model/version, polishing | Essential to assess long-read accuracy and assembly comparability | Not available/not verifiable |
| Assembly software/version/parameters and input-read hash per FASTA | Needed to reproduce each assembly | Not fully recorded in active result manifest |
| Implemented completeness/contamination/taxonomy QC | Configured but not computed by script 12a | Missing |
| Laboratory culture method, organism-selection process, specimen-specific threshold provenance | Needed to validate clinical phenotype and selection mechanism | Not fully available in reviewed files |
| Nursing-home/facility and provider metadata | Required for protocol nesting and facility confounding | Absent from analysis-ready table |
| GFF/Prokka input hash, version, and assembly-contig identity | Needed for locus annotation | Incomplete |
| Dedicated AMR tool/database/version output | Protocol objective includes AMR | Missing |
| Protocol-deviation/adaptation decision log | Distinguishes defects from approved secondary analyses | Missing |
| Prospective model/endpoint/multiplicity plan | Prevents outcome-driven selection | Missing |

Unavailable raw-read generation, basecalling, laboratory, nursing-home, and provider metadata are review limitations. They are not assumptions to fill in.

## 9. External evidence and decision mapping

| Source | Decision it informs | Agreement with current code |
|---|---|---|
| [YELLOW RoUTIne protocol](https://pmc.ncbi.nlm.nih.gov/articles/PMC11363575/) and [published PDF](https://pure.amsterdamumc.nl/files/142459523/Yellow-routine-prospective-cohort-study-protocol.pdf) | Prespecified phenotype, collection strata, repeated-measures/facility structure, ASB-to-subsequent-UTI question | **Partial/disagree:** current phenotype omits required elements and current model answers a different estimand |
| [Nursing-home <=25-SNP study](https://pmc.ncbi.nlm.nih.gov/articles/PMC9686610/) | Prior same-strain boundary | **Number adopted, method not equivalent:** current DNAdiff assembly totals require explicit operational labelling/validation |
| [Landman long-read benchmark](https://pmc.ncbi.nlm.nih.gov/articles/PMC11587635/) | Coverage and task-specific assembler choice | **Partial:** Longcycler is defensible for typing; universal Longcycler VF/AMR/plasmid use needs sensitivity/justification |
| [Parsnp/Harvest](https://pmc.ncbi.nlm.nih.gov/articles/PMC4262987/) and [official options](https://harvest.readthedocs.io/en/v1.0/content/parsnp/tutorial.html) | Core-genome alignment and recombination filtering | **Implementation/documentation mismatch:** stale/disconnected branch; `-x` comment wrong |
| [pyseer](https://pmc.ncbi.nlm.nih.gov/articles/PMC6289128/) | Microbial GWAS population structure | **Disagree:** current primary gene GLMM omits lineage/population structure |
| [Heinze & Schemper](https://pubmed.ncbi.nlm.nih.gov/12210625/) | Sparse/separated logistic regression | **Risk identified but not resolved:** current code flags separation yet accepts unstable/singular fits and permits ordinary-GLM fallback |
| [VFDB 2022](https://academic.oup.com/nar/article/50/D1/D912/6446532) | VF detection database | **Broadly aligned for screening:** cache/provenance and biological interpretation are insufficient |
| [Validated in-silico ExPEC/UPEC genotyping](https://journals.asm.org/doi/10.1128/jcm.01269-20) | Operational marker combinations and terminology | **Partial:** current one-component module presence must not be called a complete validated system |
| [AMRFinderPlus](https://pubmed.ncbi.nlm.nih.gov/34135355/) and [ResFinder 4.0](https://pmc.ncbi.nlm.nih.gov/articles/PMC7662176/) | Dedicated AMR determinant calling | **Not implemented:** current code correctly says AMR integration was not performed |
| [PlasmidFinder](https://pubmed.ncbi.nlm.nih.gov/24777092/) | Replicon detection/typing | **Aligned only for replicon evidence:** not proof of complete plasmids or gene location |
| [IDSA ASB guideline](https://academic.oup.com/cid/article/68/10/e83/5407612), [revised McGeer](https://pmc.ncbi.nlm.nih.gov/articles/PMC3538836/), [Loeb](https://pubmed.ncbi.nlm.nih.gov/11232875/) | Clinical context and labelled sensitivity definitions | **Must remain separate:** none should silently replace the protocol phenotype |

## 10. Validation performed

All checks were read-only with respect to pipeline code, data, and generated results. The existing verifier was inspected but not executed because it writes verification files; its checks were independently reproduced in memory.

| Check | Result |
|---|---|
| Parse reviewed R files | PASS: 81/81 |
| Shell syntax | PASS: both runners |
| Test suite | PASS: 3 files, 20 blocks, 105 expectations, 0 failures/errors/warnings/skips |
| `git diff --check` | PASS |
| Lint | 8,285 diagnostics; dominated by style plus 1,050 object-usage findings |
| Current dependency preflight | FAIL: Parsnp, Panaroo, Prokka, and PlasmidFinder unavailable in active environment |
| Clinical anchor | PASS: 583/166; 18 UTI and 565 Not_UTI |
| Current manifest | PASS internally: 532/161; Longcycler only; 16 UTI and 516 Not_UTI after status join |
| Prior plan anchors | REPRODUCED but reclassified: 556 VF rows; 963 mixed pairs; 394 mixed transitions; 138 at <=25 |
| Current all-pair denominator | PASS: expected 893; exact 893-row subset present in old pair table |
| Current adjacency denominator | PASS: 371 from 139; 362 retained plus 9 new adjacencies |
| Current operational same-strain summary | PASS: 140/371; Not_UTI-to-UTI 5/9 |
| Protocol suprapubic one-change counterfactual | 2 included episodes; conditional clinical 18->20 and current WGS 16->18 |
| Missing-WGS transition check | FAIL: 53 timeline episodes absent from current manifest; 63/418 adjacencies affected; all called replacement |
| Graph contradiction check | FAIL: 22 direct >25-SNP pairs in same components, maximum 139; 5 affected adjacent same-strain calls |
| Variant cache/contig check | FAIL: 8 current/8 legacy candidate files; one differs by 9 substitutions; 3 current pairs multi-contig |
| VF case Jaccard check | FAIL: 0/417 non-missing `vf_jaccard` despite 386 non-missing module Jaccards |
| Model output check | 0 FDR-significant; 6 singular; 37 sparse/separation; output stale under current manifest |
| AMR integration check | NOT PERFORMED, correctly disclosed by script 29 |
| Workspace preservation | PASS: excluding this README, the pre-existing Git-status fingerprint remains `03c87805a7e646dbbd1f486207b466102482afde493acde942a7dfe0608a95e2`; reviewed source and key-output hashes are unchanged |
| Confidentiality check | PASS: no participant identifier or row-level participant record appears in this README |

The full pipeline was not rerun because the approved review is findings-only, a rerun would overwrite generated results, and the active environment lacks required tools. Failures were recorded rather than repaired or suppressed.

## 11. Sequenced remediation roadmap

1. **Create an immutable analysis release.** Commit/label the intended source, create a new run directory, and prohibit in-place mixing with prior outputs.
2. **Adjudicate the clinical phenotype.** Resolve suprapubic pain, symptom recency, alternative focus, specimen-specific thresholds, source-row conflicts, and indeterminate states with protocol/clinical owners; version the result.
3. **Declare the estimands.** Separate the protocol ASB-to-subsequent-UTI estimand from cross-sectional isolate comparisons, longitudinal VF descriptions, and same-strain operational analyses.
4. **Declare task-specific assembly policy.** If 532 Longcycler is primary, justify it for typing and add prespecified Flye/assembly-method sensitivities for VF/AMR/plasmid questions; never use silent fallback.
5. **Repair provenance and preflight.** Pin Conda/R, record tool/database/read/assembly hashes, implement configured QC, fix true-force/output-completeness behavior, and make verifier checks non-writing until final publication.
6. **Regenerate upstream WGS branches.** Core, Panaroo, VF, MLST, and dedicated AMR/plasmid analyses must all consume the same declared manifest or an explicitly task-specific manifest.
7. **Repair pairwise/longitudinal logic.** Publish a fresh 893-row pair table, filter primary rows, represent missing WGS explicitly, use direct evidence for adjacent calls, and keep graph lineage context separate.
8. **Repair variants and runner order.** Persist current hashed SNP paths and contig identifiers, reorder dependencies, regenerate annotations, and manually validate any manuscript locus.
9. **Repair VF case/module logic.** Fix `vf_jaccard`, preserve zero-hit rows, fail on tool errors, bind caches to inputs/databases, and distinguish partial from complete systems.
10. **Redesign statistics before rerun.** Prespecify endpoints/testing families, event/facility/repeated-measures/population structure, sparse-outcome method, missingness, and sensitivity analyses.
11. **Regenerate every narrative/figure.** Build one claim registry and archive all obsolete “current/final” artifacts. No manual copy-forward of headline numbers.
12. **Independent release verification.** Run clean-room tests, manifest equality, stale-claim scanning, protocol counterfactuals, graph/variant fixtures, and a second scientific review before manuscript submission.

## 12. Publication-language guardrails

Until the roadmap is complete, use language like:

> In a Longcycler-selected subset of 532 *E. coli* isolate-linked episodes, 371 adjacent within-resident transitions were reconstructed. Under an operational assembly-to-assembly DNAdiff threshold of <=25 SNPs, 140 transitions met the threshold, including 5 of 9 Not_UTI-to-UTI transitions. These exploratory figures are conditional on culture, *E. coli* isolation, sequencing, QC, the implemented operational UTI definition, and available adjacent observations; the threshold is not method-equivalent to read-based core-SNP validation.

Do not use language like:

- “ASB protected residents from UTI.”
- “The same strain caused the UTI” without the operational-method qualifier.
- “A mutation in a named gene caused the phenotype switch.”
- “No virulence factor is associated with UTI.”
- “AMR was analyzed.”
- “Plasmids carried the genes” when only replicons were detected.
- “The analysis is publication-ready.”

## 13. Final audit conclusion

The source contains several strong foundations: explicit manifest validation, a repaired SHA-bound DNAdiff cache, reproducible clinical/Longcycler denominators, useful warning outputs, and passing helper tests. The current repository nevertheless fails the publication-grade requirement because source, outputs, summaries, and protocol estimands are not reconciled to one analysis state.

The highest-priority decision is not a statistical tweak. It is to freeze one phenotype and one declared assembly/analysis policy, regenerate all downstream artifacts from that state, and answer only estimands supported by the available cohort data. Once P0 items are closed and the verification tests pass, the supported exploratory Longcycler results can form a defensible secondary analysis; the protocol's main ASB-to-subsequent-UTI question still requires the full longitudinal cohort and facility metadata.
