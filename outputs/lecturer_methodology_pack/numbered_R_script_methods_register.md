# Audited register of numbered R scripts

This register explains the active root-level numbered R scripts and the main runner as they currently operate. It is a methods source of truth for lecturer and supervisor discussion, not an endorsement of every generated result.

The executed primary workflow is the **mixed canonical set of 556 selected assemblies**. The **532 selected Longcycler assemblies** form a separately labelled sensitivity analysis. SNP-derived persistence counts are intentionally not frozen in this register: they must be read from the fresh, provenance-validated rerun.

## Denominator and terminology anchors

- Clinical audit universe: 585 records from 167 participants.
- Primary clinical denominator: 583 episodes from 166 participants, comprising 18 UTI and 565 Not_UTI episodes.
- Active assembly-QC universe: 1,291 assembly records after primary/genomics curation. These include assembler alternatives and are not 1,291 independent biological episodes.
- Executed mixed canonical set: 556 selected QC-passing rows from 162 participants, comprising 532 Longcycler and 24 Flye fallback assemblies.
- Mixed VF/model set: 556 rows, 17 UTI and 539 Not_UTI, with 227 detected VFDB-derived binary gene columns.
- Longcycler sensitivity set: 532 rows from 161 participants, 16 UTI and 516 Not_UTI; 514 rows have a preferred MLST call, representing 80 ST labels.
- Mixed longitudinal set: 394 adjacent retained transitions from 144 participants.
- Longcycler sensitivity longitudinal set: 371 adjacent retained transitions from 139 participants, including nine Not_UTI-to-UTI transitions.
- A **dnadiff SNP difference** is an assembly-to-assembly MUMmer result. It is not a Parsnp core-genome distance and is not a wgMLST allele distance.
- The project rule of **<=25 dnadiff SNPs** is a predefined operational threshold for strong same-strain support. It is not presented as a universal biological boundary.
- Not_UTI means that the episode did not meet this study's combined culture-and-symptom UTI rule. It does not mean healthy, sterile urine or absence of bacteriuria.

## Infrastructure and clinical foundation

### `00_config.R`

- **Classification / runner status / unit:** Infrastructure; sourced by numbered scripts; no analysis unit.
- **Purpose / inputs:** Defines project paths, batch definitions, manual-curation helpers, output directories, the primary >=1e3 CFU/mL culture threshold and the operational <=25-SNP rule.
- **Filters, tools and tests:** Applies no biological filter or statistical test itself.
- **Repeated measures:** Not applicable.
- **Outputs and safe claim:** Supplies a shared configuration so numbered scripts use consistent paths and primary-analysis helpers.
- **Limitations:** The SNP threshold is described only as prior project practice, not externally calibrated here. Sourcing creates directories, and package/tool versions are not frozen.

### `00_input_snapshot.R`

- **Classification / runner status / unit:** Debug/audit utility; not called by the runner; unit is an input file preview.
- **Purpose / inputs:** Reads available VF, MLST, plasmid and clinical outputs, writes dimensions, names and small row/column previews, and lists existing dnadiff reports.
- **Filters, tools and tests:** Missing files are skipped; previews are deliberately truncated; no statistical test.
- **Repeated measures:** Not applicable.
- **Outputs and safe claim:** Writes under `results/within_person/debug/input_snapshot/`; supports input-shape debugging only.
- **Limitations:** It is not biological evidence, and its dimensions file is append-based.

### `00_make_assembly_metadata.r`

- **Classification / runner status / unit:** Required provenance builder; omitted by the runner; starts with expected isolates and expands to assembly rows.
- **Purpose / inputs:** Uses the overview workbook as the authoritative expected-isolate universe, supplements it with batch CSVs, discovers FASTAs, preserves Flye/Longcycler alternatives and writes exclusion audits.
- **Filters and denominator:** Unlinked/unexpected FASTAs and quarantined records are excluded from active episode analysis. The current metadata has 1,299 rows overall and 1,291 rows after primary/genomics curation.
- **Tools and tests:** Biostrings supplies lightweight contig, total-size and GC metrics when available; duplicate Assembly_ID checks; no inference.
- **Repeated measures:** Assembler alternatives remain explicit and must not be treated as independent episodes.
- **Outputs and safe claim:** Writes `assembly_metadata.csv`, two mirrored metadata outputs and QC/provenance tables. It can support the factual sample/FASTA provenance statement.
- **Limitations:** It does not describe sequencing, basecalling, read QC, assembler settings or polishing. Missing Biostrings causes lightweight metrics to be skipped.

### `00a_load_clean_clinical.R`

- **Classification / runner status / unit:** Main clinical preprocessing; Phase 0; raw clinical records harmonised to participant/timepoint fields.
- **Purpose / inputs:** Loads available batch 1-6 CSVs, reconciles names and types, separates symptom fields and retains the B1/B2 SnS fallback.
- **Filters and denominator:** Loads only configured batch files present on disk and warns for missing batches; it does not assign the final UTI outcome.
- **Tools and tests:** Tidyverse harmonisation; no inferential test.
- **Repeated measures:** Raw repeated records remain available for episode-level collapse in the next step.
- **Outputs and safe claim:** Writes `results/clinical/intermediate/clinical_merged.rds` and the detected symptom-column map; supports a reproducible data-harmonisation description.
- **Limitations:** The thesis still needs the provenance of batch CSV exports and an explanation for any missing batch file.

### `00b_classify_episodes.R`

- **Classification / runner status / unit:** Main clinical classification; Phase 0; one clinical episode per participant/timepoint.
- **Purpose / inputs:** Collapses harmonised records and creates the authoritative primary UTI/Not_UTI status while retaining legacy labels for audit.
- **Filters and denominator:** UTI requires culture support at >=1e3 CFU/mL and catheter-aware symptoms. Non-indwelling episodes require a local urinary symptom or flank pain plus a systemic symptom; indwelling episodes require a systemic symptom. All others are Not_UTI with a subgroup. Current primary denominator: 583 episodes/166 participants, 18 UTI and 565 Not_UTI.
- **Tools and tests:** Rule validation, duplicate-key checks and audit tables; no inferential test.
- **Repeated measures:** Multiple raw rows are collapsed before classification; participants can still contribute multiple episodes.
- **Outputs and safe claim:** `results/clinical/status_map.csv` plus classification, movement, symptom and CFU audits. This supports an exact operational outcome definition.
- **Limitations:** Not_UTI is heterogeneous and can include indeterminate evidence. It must not be described as healthy or bacteria-free.

### `00c_plot_clinical_summary.R`

- **Classification / runner status / unit:** Descriptive reporting; Phase 0; clinical episode, participant and assembly summaries.
- **Purpose / inputs:** Uses the current status map and metadata to visualise the denominator, rule components, reclassification and assembly context.
- **Filters, tools and tests:** Does not redefine status; descriptive counts and graphics only.
- **Repeated measures:** Participant counts are presented separately where appropriate; no repeated-measures model.
- **Outputs and safe claim:** Clinical plots plus waterfall, symptom, CFU, reclassification and inspection tables.
- **Limitations:** No association or causal claim is supported.

### `00d_derive_plot_timepoints.R`

- **Classification / runner status / unit:** Display-only preprocessing; Phase 0; clinical episode.
- **Purpose / inputs:** Derives poster-friendly timepoint labels and placement-confidence fields from dates and event labels.
- **Filters, tools and tests:** Preserves the primary clinical status; uncertain and unplaced event rows are audited; no statistical test.
- **Repeated measures:** Not applicable beyond retaining the episode structure.
- **Outputs and safe claim:** `status_map_with_poster_tp.csv` and placement audits; useful for explaining figure ordering.
- **Limitations:** Poster positions must never be treated as elapsed biological time or as an inferential time variable.

## Genomic screening, typing and QC

### `02_gene_presence_analysis.R`

- **Classification / runner status / unit:** Executed primary VF screen; Phase 1; one selected canonical assembly per participant/timepoint.
- **Purpose / inputs:** Screens selected FASTAs against VFDB and converts detected genes to an episode-level binary matrix.
- **Filters and denominator:** Applies primary/genomics curation and selected-canonical filtering. ABRicate thresholds are >=80% identity and >=80% coverage. Current mixed output: 556 rows/162 participants and 227 VF columns.
- **Tools and tests:** ABRicate VFDB; descriptive prevalence only.
- **Repeated measures:** One selected row per episode avoids counting assembler alternatives as independent observations.
- **Outputs and safe claim:** VF hit RDS, `vf_pa_all.csv`, prevalence summaries and plots. Safe wording is "VFDB genes detected under the stated thresholds in selected assemblies."
- **Limitations:** This is not a union across assemblers. Cache keys lack FASTA hashes, and ABRicate/VFDB version and date are not recorded. Non-detection is not proof of biological absence.

### `03_plotting.R`

- **Classification / runner status / unit:** Legacy exploratory visualisation; called only if `RUN_LEGACY_EXPLORATORY_PLOTS=1`; mixes episode, isolate, participant and pair units.
- **Purpose / inputs:** Generates older VF, phylogeny, epidemiology, transmission, nitrate and timeline plots from many optional outputs.
- **Filters, tools and tests:** Uses whichever optional inputs exist and includes legacy-output branches; mostly descriptive.
- **Repeated measures:** Units and clustering vary and are not governed by one model.
- **Outputs and safe claim:** Broad internal exploratory figures only.
- **Limitations:** Superseded by the structured VF workflow; heterogeneous denominators and legacy ASB wording make it unsuitable as a methods source.

### `04_gene_breakdown.R`

- **Classification / runner status / unit:** Gene curation plus exploratory focus-gene modelling; Phase 1b; gene hits and clinical episodes.
- **Purpose / inputs:** Creates/reuses the gene-category map, category and nitrate summaries, then attempts selected focus-gene models.
- **Filters and denominator:** Reuses curated mappings if present; otherwise uses ordered regex heuristics. Focus genes are joined to the clinical status table.
- **Tools and tests:** GLMM `Outcome ~ gene + timepoint + batch + (1|Participant_id)`, accepts singular fits, falls back to GLM, and applies BH over focus genes.
- **Repeated measures:** A resident random intercept is attempted only in the model subsection.
- **Outputs and safe claim:** Gene map, annotated/category/nitrate tables and exploratory focus-gene output. The map is safe as a transparent functional grouping aid.
- **Limitations:** Missing WGS/no-hit clinical rows can be filled as gene absent in the focus model. Regex categories are not validated curation, and sparse/fallback models are not primary inference.

### `05_gene_overview_plots.R`

- **Classification / runner status / unit:** Descriptive VF overview; Phase 1b; selected VF episode row.
- **Purpose / inputs:** Calculates prevalence and divides genes into >95% dataset-level core and <95% variable groups.
- **Filters, tools and tests:** Uses only canonical binary VF columns; no inferential test.
- **Repeated measures:** No adjustment because the output is descriptive.
- **Outputs and safe claim:** Core/variable lists and prevalence/heatmap figures; describes frequency among selected sequenced episodes.
- **Limitations:** Dataset-level core is not pangenome core, biological essentiality or evidence of protective colonisation.

### `06_MLST.R`

- **Classification / runner status / unit:** Executed MLST integration; Phase 1; canonical episode plus all-isolate provenance.
- **Purpose / inputs:** Uses provider SeqSphere MLST as primary and a labelled local fallback, retaining source fields.
- **Filters and denominator:** Provider calls require `PercGoodTargets >=95`. Mixed canonical output: 556 rows, 527 provider calls, 6 local fallbacks, 23 missing; 533 typed rows and 83 ST labels.
- **Tools and tests:** Provider/local source reconciliation and guardrail scripts; no phenotype test.
- **Repeated measures:** One canonical episode row is used downstream, but residents can repeat over time.
- **Outputs and safe claim:** Provider-preferred canonical/all tables and a source audit. ST can be used as lineage context with explicit provenance.
- **Limitations:** Provider assembler provenance often includes Flye. It is incorrect to call the Longcycler-selected subset "Longcycler-derived MLST." Scheme/version and provider manifest are still required, and ST does not prove strain identity.

### `07_explore_MLST.R`

- **Classification / runner status / unit:** Descriptive MLST exploration; Phase 1b; biological isolate after assembly-row collapse.
- **Purpose / inputs:** Uses provider-preferred all-isolate MLST when present, collapses to one row per Isolate_ID and plots frequency/top STs.
- **Filters and denominator:** Prefers complete/more-locus calls. Current all-isolate source has 1,291 assembly rows representing 593 isolates; the current frequency output contains 562 typed isolates and 90 STs.
- **Tools and tests:** Descriptive counts and plots only.
- **Repeated measures:** Multiple isolates from the same participant remain separate.
- **Outputs and safe claim:** ST frequency table, top-20 plots and an exploratory metadata join.
- **Limitations:** This is not the 556-row canonical episode denominator and is not a participant-level temporal trend analysis.

### `08_core_vs_plasmid.R`

- **Classification / runner status / unit:** Exploratory plasmid-lineage analysis; Phase 1b; selected canonical assembly/isolate episode.
- **Purpose / inputs:** Relates chromosomal ST to available pMLST schemes, or to PlasmidFinder replicons when pMLST schemes are absent.
- **Filters and denominator:** Screens selected canonical FASTAs; association screens use the five commonest STs and replicons present in >5% of joined rows.
- **Tools and tests:** `mlst`/pMLST or ABRicate PlasmidFinder; Fisher exact tests with BH.
- **Repeated measures:** Resident-level clustering is not used.
- **Outputs and safe claim:** ST frequencies, plasmid typing/replicon tables and exploratory ST-plasmid associations.
- **Limitations:** The header's "inferential-core" label is too strong. Tool availability changes the endpoint, versions are not frozen and episode-level Fisher tests violate independence when residents repeat.

### `09_inc_plasmid_network.R`

- **Classification / runner status / unit:** Executed exploratory plasmid screen; Phase 1b; selected assembly/isolate.
- **Purpose / inputs:** Runs PlasmidFinder on selected FASTAs and builds replicon co-occurrence and ST-replicon networks.
- **Filters and denominator:** Requires selected canonical existing FASTAs; network replicons must occur in at least three positive rows.
- **Tools and tests:** ABRicate PlasmidFinder and graph visualisation; no association model.
- **Repeated measures:** No participant clustering.
- **Outputs and safe claim:** Input manifest, long hits, P/A matrix and network PDFs; shows detected replicon patterns.
- **Limitations:** Cache keys are basename-only, database/settings provenance is incomplete, and zero-hit assemblies are absent from the matrix, which can conflate missing screening and zero detected replicons.

### `10_replicon_heatmap.R`

- **Classification / runner status / unit:** Descriptive plasmid visualisation; Phase 1b; replicon-matrix isolate row.
- **Purpose / inputs:** Creates a clustered replicon heatmap annotated with provider-preferred ST.
- **Filters, tools and tests:** Retains replicons present in at least five matrix rows, otherwise all; uses ComplexHeatmap or pheatmap; no inference.
- **Repeated measures:** No adjustment; descriptive only.
- **Outputs and safe claim:** Replicon heatmap PDF/PNG; visualises the matrix represented in the input.
- **Limitations:** Because zero-hit selected assemblies may be missing from the matrix, the heatmap denominator is not automatically all 556 selected episodes.

### `11_compare_strains_helpers.R`

- **Classification / runner status / unit:** Support/provenance code; sourced by pairwise scripts; pair endpoint/metric helper.
- **Purpose / inputs:** Loads assembly, status, MLST, VF and plasmid tables; resolves endpoints; calculates dnadiff, Mash and Jaccard metrics; applies a composite class.
- **Filters and denominator:** The repaired implementation must accept only selected-canonical QC-pass endpoints and reuse cache only when exact FASTA fingerprints match.
- **Tools and tests:** MUMmer dnadiff, Mash and set Jaccard; no standalone test.
- **Repeated measures:** Pair dependence is not modelled in the helper.
- **Outputs and safe claim:** No standalone output; supports reproducible metric generation when provenance validation passes.
- **Limitations:** The former implementation fell back to noncanonical assemblies and used SampleKey-only caches. The composite Same rule can pass missing accessory evidence and is not the strict SNP definition.

### `11_compare_strains.R`

- **Classification / runner status / unit:** Main comparative-genomics layer after cache repair; Phase 2; unordered within-participant assembly pair.
- **Purpose / inputs:** Calculates dnadiff identity/SNP, Mash, ST, VF and plasmid concordance for all eligible within-resident pairs.
- **Filters and denominator:** Primary rerun must use only current selected-canonical QC-pass FASTAs. The present 556-row selected set implies 963 unordered within-resident pairs. Strict support is `TotalSNPs <=25` and is reported separately from the composite class.
- **Tools and tests:** dnadiff, Mash and Jaccard; exploratory within/between Wilcoxon and status-group Kruskal-Wallis tests.
- **Repeated measures:** Pairs share residents and endpoints; the exploratory tests do not model this dependence.
- **Outputs and safe claim:** Fresh `pairwise_metrics.csv`, provenance, summaries and plots. Safe wording is "assembly-to-assembly dnadiff SNP differences with exact endpoint provenance."
- **Limitations:** dnadiff SNPs are not core-genome SNPs or wgMLST alleles. The <=25 rule is operational, not universal. The pre-repair 1,018-row output was invalid because stale reports and noncanonical endpoints were reused; strict results must be dynamic after rerun.

### `12a_wgs_qc.R`

- **Classification / runner status / unit:** Executed assembly QC and canonical selection; Phase 1; assembly row, then participant/timepoint selection.
- **Purpose / inputs:** Measures contiguity/size and chooses one passing assembly per episode, preferring Longcycler.
- **Filters and denominator:** Pass requires <=200 contigs, N50 >=20 kb, 4-6 Mb total size and no read error. Current active universe: 1,291 assembly records; selected mixed set: 556 rows/162 participants, 532 Longcycler and 24 Flye fallbacks.
- **Tools and tests:** seqinr FASTA metrics and an exploratory Fisher selection-by-status diagnostic.
- **Repeated measures:** One selected row per episode; Fisher does not cluster repeated episodes.
- **Outputs and safe claim:** QC summary, canonical selection, bias audits and QC plot. Supports the executed mixed denominator and the separate Longcycler sensitivity definition.
- **Limitations:** Completeness and contamination are not applied here. Passing is a contiguity/size screen, not proof of completeness or no contamination.

### `12b_core_snp.R`

- **Classification / runner status / unit:** Separate core-genome phylogeny branch; Phase 1; selected genome and cohort-wide pair distance.
- **Purpose / inputs:** Builds a Parsnp alignment, snp-dists matrix, neighbour-joining tree and a separate rule-based pair table.
- **Filters and denominator:** Targets the 556 selected QC-pass genomes. In this branch only, <=25 is Same, <=1000 Related and >1000 Different.
- **Tools and tests:** Parsnp `-c -r ! -x`, snp-dists and ape neighbour-joining; no phenotype model.
- **Repeated measures:** Not modelled.
- **Outputs and safe claim:** Core manifest, alignment, SNP distances, pair labels and tree; provides a separate phylogenetic context.
- **Limitations:** Script 11 does not consume these SNPs despite stale header wording. Reference selection is random; fingerprints use path/size/mtime rather than content; stale output can cause a warning and successful exit unless forced.

### `12c_panaroo.R`

- **Classification / runner status / unit:** Executed pangenome branch; Phase 1; selected genome/GFF.
- **Purpose / inputs:** Audits the GFF inventory, regenerates missing Prokka annotations and runs Panaroo on the selected set.
- **Filters and denominator:** Requires current metadata/QC and complete GFF coverage for selected QC-pass genomes; intended denominator is 556.
- **Tools and tests:** Prokka with compliant E. coli settings; Panaroo strict cleaning and invalid-gene removal; no association test.
- **Repeated measures:** Not applicable to pangenome construction.
- **Outputs and safe claim:** GFF/Panaroo manifests and pangenome outputs; describes accessory content conditional on this assembly/annotation workflow.
- **Limitations:** Prokka/Panaroo versions are not frozen and outputs are assembler/annotation dependent.

### `13_visualise_panaroo_selection.R`

- **Classification / runner status / unit:** QC visualisation/audit; Phase 1; assembly row and selected episode.
- **Purpose / inputs:** Contrasts QC pass, canonical selection, GFF availability and Panaroo inclusion across metadata strata.
- **Filters and denominator:** Distinguishes 1,291 active assembly records from 556 selected rows.
- **Tools and tests:** Descriptive plots and exploratory Fisher selection-by-status test.
- **Repeated measures:** Fisher does not cluster resident episodes.
- **Outputs and safe claim:** Selection CSVs, QC-bias output and plots; makes attrition visible.
- **Limitations:** A nonsignificant sparse Fisher result does not establish absence of selection bias.

## Modelling and longitudinal analysis

### `14_genotype_phenotype_model.R`

- **Classification / runner status / unit:** Exploratory adjusted association; Phase 2; VF-ready episode.
- **Purpose / inputs:** Screens genomic features and fits logistic participant-random-intercept models for UTI versus Not_UTI.
- **Filters and denominator:** Uses the 556-row mixed VF/model set: 17 UTI and 539 Not_UTI. Retains features at 5-95% prevalence; current output contains 119 Fisher screens and 50 selected models.
- **Tools and tests:** Fisher screening; GLMM `Outcome ~ feature + Timepoint + Batch + (1|Participant_id)`; singular fits accepted; GLM fallback; BH within the selected model set.
- **Repeated measures:** Resident intercept where GLMM succeeds; no clustering after GLM fallback.
- **Outputs and safe claim:** Denominator, univariable and model tables, warnings and plots. Supports exploratory adjusted estimates with explicit diagnostics.
- **Limitations:** Seventeen UTI rows are underpowered. Two-stage data-driven selection plus BH only on selected models is not confirmatory; ST is omitted. Current results have no FDR<0.05 findings, six singular fits and 37 sparse/separation flags.

### `15_longitudinal_patterns.R`

- **Classification / runner status / unit:** Unsupported primary persistence branch; invoked in Phase 3; clinical episode graph node.
- **Purpose / inputs:** Creates graph-component Strain_ID values from composite Same edges, timelines and switch candidates.
- **Filters and denominator:** Includes every status-map SampleKey as a graph vertex; uses composite `Classification == Same`; orders by date/fallback.
- **Tools and tests:** igraph connected components and descriptive transition counts.
- **Repeated measures:** No recurrent-event model.
- **Outputs and safe claim:** Timeline, transition, switch, persistence and swimmer outputs; at most exploratory after redesign.
- **Limitations:** No-WGS episodes receive singleton Strain_ID values, graph transitivity can merge contradictory pairs, composite Same is not strict SNP persistence, and the output uses the full 585-row audit status map rather than the primary genomic denominator.

### `16_within_host_evolution.R`

- **Classification / runner status / unit:** Unsupported downstream branch; Phase 3; script-15 switch pair.
- **Purpose / inputs:** Recomputes pair dnadiff metrics and lists VF/plasmid differences for switch candidates.
- **Filters, tools and tests:** Inherits script-15 candidates; dnadiff and set differences; descriptive only.
- **Repeated measures:** No clustering or inferential model.
- **Outputs and safe claim:** Evolution-event table and summary; no current primary claim.
- **Limitations:** Inherits invalid candidates and historical cache risks; VF selection can include metadata and missing plasmid rows can appear as empty profiles.

### `17_lineage_analysis.R`

- **Classification / runner status / unit:** Exploratory lineage analysis; Phase 2; typed VF-ready episode.
- **Purpose / inputs:** Calculates episode UTI proportion for common STs and tests each against all other STs.
- **Filters and denominator:** Uses typed rows and STs with at least five episodes.
- **Tools and tests:** Fisher exact tests with BH; simple plotted intervals.
- **Repeated measures:** Resident episodes are treated as independent.
- **Outputs and safe claim:** ST UTI-proportion profile and plot; exploratory over-representation only.
- **Limitations:** Use "proportion," not "risk." Sparse cells and repeated residents prevent confirmatory inference.

### `18_annotate_variants.R`

- **Classification / runner status / unit:** Unsupported variant parser; Phase 3; variant record from a switch pair.
- **Purpose / inputs:** Parses dnadiff `.snps` files to position/base/type fields.
- **Filters, tools and tests:** Uses script-15 candidates; no statistical test.
- **Repeated measures:** Not applicable.
- **Outputs and safe claim:** Raw annotated-SNP table/report only; no defensible mutation-location claim.
- **Limitations:** Reference/query contig identifiers are discarded, preventing safe mapping in multi-contig assemblies; candidate/cache provenance is inherited.

### `19_host_context.R`

- **Classification / runner status / unit:** Exploratory host-context extraction; Phase 3; switch endpoint.
- **Purpose / inputs:** Looks up collection method and raw symptom summaries for candidate endpoints.
- **Filters, tools and tests:** Exact-value heuristic over raw symptom columns; no model.
- **Repeated measures:** No adjustment.
- **Outputs and safe claim:** Host-context table for manual review only.
- **Limitations:** Inherits invalid candidates, may miss Dutch/differently encoded values, and lacks antibiotics and broader host factors.

### `20_variant_annotation_deep.R`

- **Classification / runner status / unit:** Unsupported mutation annotation; Phase 3; variant record.
- **Purpose / inputs:** Overlays script-18 positions on a selected/preferred Prokka GFF and assigns the first overlap.
- **Filters, tools and tests:** Restricts to target transitions/SNP rows; interval lookup only.
- **Repeated measures:** Not applicable.
- **Outputs and safe claim:** Detailed variant CSV; no current defensible gene-level mutation claim.
- **Limitations:** Upstream positions lack contig ID and GFF seqid is ignored, so a coordinate can match the wrong contig; To_Time grouping can also be mixed.

### `21_publication_figures.R`

- **Classification / runner status / unit:** Rendering of unsupported upstream results; Phase 3; timeline and variant rows.
- **Purpose / inputs:** Produces a swimmer plot and linear mutation map.
- **Filters, tools and tests:** Participants with multiple plotted timepoints; all available variants; plotting only.
- **Repeated measures:** No model.
- **Outputs and safe claim:** Two publication-named figures; no primary scientific claim.
- **Limitations:** Inherits scripts 15 and 20, uses hard-coded highlight IDs and places contig-local positions on an artificial 5 Mb line.

### `22_vf_build_analysis_dataset.R`

- **Classification / runner status / unit:** Executed central integration; Phase 2; selected participant/timepoint linked to a clinical episode.
- **Purpose / inputs:** Joins VF calls, primary status, gene categories and provider-preferred MLST; computes burden/category fields and the audited Uricult-to-UTI-N bridge.
- **Filters and denominator:** Applies primary/genomics curation. Bridge priority is date match then lowest UTI number, otherwise participant-only lowest number. Current mixed output: 556 rows/162 participants, 17 UTI and 539 Not_UTI.
- **Tools and tests:** Join, duplicate, freshness and unmatched-row diagnostics; no association test.
- **Repeated measures:** One selected row per episode; Uricult alternatives remain in a sensitivity audit.
- **Outputs and safe claim:** Canonical VF-ready/binary tables plus bridge/gap audits; defines the executed mixed model denominator.
- **Limitations:** Union-across-assembler comments are stale; MLST is not Longcycler-derived; participant-only bridging is an assumption; the 227 VFDB genes include unassigned/low-confidence features.

### `23_vf_cross_sectional.R`

- **Classification / runner status / unit:** Descriptive/exploratory VF comparison; Phase 4; VF-ready episode.
- **Purpose / inputs:** Summarises burden, prevalence and categories by UTI status in the full cohort and depth subsets.
- **Filters and denominator:** Runs `all`, >=2, >=3 and >=4 timepoint subsets; only primary-status rows enter comparisons.
- **Tools and tests:** Fisher gene/category and Wilcoxon category tests with BH where implemented.
- **Repeated measures:** Not modelled; tests are explicitly exploratory.
- **Outputs and safe claim:** Burden/prevalence/enrichment tables and plots; safe for descriptive medians and prevalence.
- **Limitations:** Stratified files contain overlapping repeated cohort subsets; filter `cohort == "all"` for overall counts. Sparse UTI, lineage and resident repetition limit inference.

### `24_vf_longitudinal_dynamics.R`

- **Classification / runner status / unit:** Main descriptive longitudinal layer after fresh pairwise rerun; Phase 4; adjacent retained VF-ready pair.
- **Purpose / inputs:** Rebuilds adjacent within-resident pairs, computes VF Jaccard/detected changes and layers SNP and ST context.
- **Filters and denominator:** Date-first ordering with labelled fallback. Mixed set: 394 transitions/144 participants. Separate Longcycler sensitivity: 371/139, including nine Not_UTI-to-UTI transitions.
- **Tools and tests:** Set Jaccard and descriptive summaries; no causal or recurrent-event model.
- **Repeated measures:** Residents can contribute multiple correlated transitions.
- **Outputs and safe claim:** Longitudinal transition and strain-context tables/plots; supports observed VF stability and separately labelled strict SNP evidence after repair.
- **Limitations:** Gained/lost means detected P/A change, not proven biological gain/loss. Adjacency is among retained samples. Stratified outputs duplicate pairs, and strict SNP results must be dynamic.

### `25_vf_lineage_vf_interaction.R`

- **Classification / runner status / unit:** Exploratory confounding diagnostic; Phase 4; typed VF-ready episode.
- **Purpose / inputs:** Tests whether burden varies across ST, within ST by status and whether ST composition differs by status.
- **Filters and denominator:** STs need >=5 episodes, relaxed to >=3 if none qualify.
- **Tools and tests:** Kruskal-Wallis, within-ST Wilcoxon and simulated Fisher with 10,000 replicates.
- **Repeated measures:** Residents are not clustered.
- **Outputs and safe claim:** Burden-by-ST tables, summary and plots; useful for showing lineage structure.
- **Limitations:** Nonsignificance cannot exclude confounding; sparse UTI and multiple test families limit power; script 14 does not automatically apply ST adjustment.

### `26_vf_define_gene_modules.R`

- **Classification / runner status / unit:** Supplementary rule-based curation; Phase 4; gene then episode module.
- **Purpose / inputs:** Assigns genes to ordered, first-match modules and builds module presence/count, confidence and review tables.
- **Filters, tools and tests:** Unassigned and low-confidence genes remain visible; no association test.
- **Repeated measures:** One row per episode; no inference.
- **Outputs and safe claim:** Module map, episode/summary/audit tables and plots; transparent analysis framework.
- **Limitations:** Hand-authored regex modules are not validated systems and must not be described as causal or predictive.

### `27_vf_score_framework.R`

- **Classification / runner status / unit:** Exploratory supplementary endpoints; Phase 4; VF-ready episode.
- **Purpose / inputs:** Creates burden, ExPEC-like marker and UPEC-system endpoints plus ordination and status/ST summaries.
- **Filters, tools and tests:** Current VF-ready/module rules; Fisher, Wilcoxon, Spearman, PCA and Jaccard PCoA; BH within output families.
- **Repeated measures:** Episode-level tests do not adjust for residents.
- **Outputs and safe claim:** Score/marker/correlation/ordination tables and plots; descriptive and hypothesis-generating only.
- **Limitations:** Not validated prediction scores; literature mapping requires citation; sparse UTI, multiplicity and repeated residents prevent confirmation.

### `28_vf_transition_case_studies.R`

- **Classification / runner status / unit:** Clinical-first descriptive transition audit; Phase 4; ordered clinical transition.
- **Purpose / inputs:** Retains all ordered clinical transitions, including missing genomics, and layers VF/module/score/ST/pairwise evidence.
- **Filters and denominator:** Date-first ordering with documented fallback; endpoint availability is explicit; Not_UTI-to-UTI does not require WGS.
- **Tools and tests:** Rule-based descriptive case classes; no causal model.
- **Repeated measures:** Residents may contribute several transitions; no recurrent-event method.
- **Outputs and safe claim:** Transition index/summary/change/notes tables and plots; provides a complete clinical-first inventory.
- **Limitations:** Denominator differs from VF-only script 24. Case classes do not identify causes and SNP fields require fresh pairwise provenance.

### `29_vf_amr_combined_profile.R`

- **Classification / runner status / unit:** VF-plus-plasmid exploratory analysis, not AMR; Phase 4; VF-ready episode.
- **Purpose / inputs:** Audits AMR-data availability and otherwise combines VF endpoints with replicon summaries.
- **Filters and denominator:** Current logic explicitly records that true AMR integration was not performed.
- **Tools and tests:** Descriptive summaries and Spearman VF-replicon correlation.
- **Repeated measures:** Not clustered.
- **Outputs and safe claim:** VF/plasmid tables, availability report and plots; supports plasmid context only.
- **Limitations:** Filenames containing AMR are misleading. Missing plasmid matches may become zero, and absence of screening output is not evidence of susceptibility.

### `30_vf_project_summary_tables.R`

- **Classification / runner status / unit:** Reporting/consolidation; Phase 4; inherits mixed units.
- **Purpose / inputs:** Collects generated clinical, WGS, MLST, VF, transition and sensitivity outputs into numbered tables and a markdown summary.
- **Filters, tools and tests:** File availability and selected freshness/denominator checks; no new core analysis.
- **Repeated measures:** Inherits each input's approach.
- **Outputs and safe claim:** Summary CSVs, markdown, optional workbook/RDS and QC log; useful as an index after validation.
- **Limitations:** "Manuscript-ready" is too strong. Unsafe or stale upstream claims propagate unless explicitly removed.

## Audits, robustness and final reporting

### `31_audit_fasta_usage.R`

- **Classification / runner status / unit:** Standalone FASTA provenance audit; not runner-invoked; FASTA/assembly file.
- **Purpose / inputs:** Checks disk FASTAs against metadata, QC, MLST, VF, GFF and Panaroo representation and flags conflicts/staleness.
- **Filters, tools and tests:** Recursive scan with classified exclusions; strict failure by default; no statistical test.
- **Repeated measures:** Not applicable.
- **Outputs and safe claim:** FASTA usage, summary, conflict, batch and GFF/Panaroo audits; supports sample-completeness checking.
- **Limitations:** File-audit rows are not biological denominators and the report must be rerun after file changes.

### `31_audit_uti_denominator_drop.R`

- **Classification / runner status / unit:** Standalone denominator audit; not runner-invoked; clinical UTI episode/attrition stage.
- **Purpose / inputs:** Traces UTI rows through metadata, FASTA, QC, VF and MLST layers without changing labels.
- **Filters, tools and tests:** Diagnostic only; no inference.
- **Repeated measures:** Not applicable.
- **Outputs and safe claim:** UTI cascade, missing-row, rule and denominator tables; supports factual attrition explanations when current.
- **Limitations:** Older generated copies contain obsolete counts and must not be cited without freshness verification/rerun.

### `32_uti_not_uti_diagnostic_stats.R`

- **Classification / runner status / unit:** Mixed-pipeline diagnostic/sensitivity; Phase 4; episode, resident or transition by endpoint.
- **Purpose / inputs:** Explains denominator/near-miss flow and estimates sparse-case uncertainty through bootstrap, leave-one, paired and transition diagnostics.
- **Filters and denominator:** Hard-validates 18/565 clinical and 17/539 mixed VF counts; uses selected scores/features.
- **Tools and tests:** Resident bootstrap B=1,000, exploratory Fisher, leave-one-UTI, paired sign/Wilcoxon and illustrative power grid.
- **Repeated measures:** Bootstrap resamples residents and paired summaries collapse within resident; Fisher/transition tests are less protected.
- **Outputs and safe claim:** Diagnostic tables/plots under results/audit and results/vf; quantifies mixed-analysis fragility.
- **Limitations:** Cannot be used unchanged for Longcycler sensitivity. Feature selection is partly data-driven, transition rows can repeat and power estimates are not prospective.

### `33_mechanism_first_addon.R`

- **Classification / runner status / unit:** Descriptive mechanism casebook; Phase 4; clinical transition case.
- **Purpose / inputs:** Combines clinical context, strain, VF/module, accessory, plasmid, variant and optional ResFinder evidence.
- **Filters and denominator:** Uses Not_UTI-to-UTI cases and selected canonical FASTAs for optional ResFinder screening at 80/80.
- **Tools and tests:** Rule-based descriptive buckets and optional ABRicate; no causal model.
- **Repeated measures:** Residents can contribute multiple cases; no recurrent-event adjustment.
- **Outputs and safe claim:** Casebook, summaries, validations, optional AMR tables and plots; structured descriptive evidence only after fresh inputs.
- **Limitations:** Does not establish mechanism or antibiotic effect. Cache is basename-based, variant input is unsafe, AMR genes are not antibiotic exposure and counts must be dynamic.

### `34_robustness_first_addon.R`

- **Classification / runner status / unit:** Mixed-pipeline robustness synthesis; Phase 4; episode, score or existing model output.
- **Purpose / inputs:** Compares the primary outcome rule with a near-miss-expanded sensitivity and consolidates QC, model, bootstrap, leave-one and power diagnostics.
- **Filters and denominator:** Hard-validates 583/18 clinical, 556/17 VF and the current near-miss count; alternative labels remain sensitivity-only.
- **Tools and tests:** Selected Wilcoxon/Fisher comparisons with BH plus upstream diagnostic summaries.
- **Repeated measures:** Handling varies; new episode tests are not fully clustered.
- **Outputs and safe claim:** Robustness tables, claim matrix, validation, summary and plots; shows sensitivity to the clinical rule.
- **Limitations:** Mixed analysis only and does not refit all feature GLMMs under the alternative rule; non-confirmatory.

### `35_final_figure_pack.R`

- **Classification / runner status / unit:** Rendering/validation only; end of Phase 4; inherited figure-specific units.
- **Purpose / inputs:** Renders final figures after mechanism, robustness and sensitivity validation passes.
- **Filters, tools and tests:** Rejects some legacy paths and requires PASS validations; no new inference.
- **Repeated measures:** Inherits upstream handling.
- **Outputs and safe claim:** Final-figure directory, manifest and validation; can render already validated descriptive results.
- **Limitations:** Contains hard-coded mixed denominators and interpretation text plus unsafe variant/AMR branches. It is not a methods source and needs regeneration after corrected inputs.

### `36_statistical_sensitivity_addon.R`

- **Classification / runner status / unit:** Exploratory targeted sensitivity; Phase 4; participant-collapsed status, paired resident, transition and episode model.
- **Purpose / inputs:** Adds resident-collapsed score, paired binary-feature, transition-module and score-GLMM sensitivity analyses.
- **Filters and denominator:** Uses four named score endpoints; its feature catalogue combines top previous results with explicit genes; GLMMs add batch, timepoint and optional collapsed ST group.
- **Tools and tests:** Bootstrap B=5,000, rank-sum/sign/signed-rank with BH, Fisher with BH and GLMM with resident intercept/GLM fallback.
- **Repeated measures:** Participant-collapsed/paired analyses address repetition; transition tests can repeat residents and GLM fallback loses clustering.
- **Outputs and safe claim:** Statistical-sensitivity tables, validation, summary and plots; targeted exploratory evidence.
- **Limitations:** "Prespecified" is inaccurate for data-selected features. Seventeen UTI rows with multiple factors can overfit, and ST grouping is hard-coded.

## Pipeline runner

### `RUN_COMPLETE_ANALYSIS.sh`

- **Classification / runner status / unit:** Manual pipeline orchestrator; inherited units.
- **Purpose / inputs:** Runs the clinical, assembly-to-results WGS, VF, pairwise, longitudinal, exploratory and reporting scripts in a fixed order.
- **Filters, tools and tests:** `set -e` and `set -u`; checks that some conda environment is active but not the intended environment identity; invokes each script's embedded filters/tests.
- **Repeated measures:** Inherited from individual scripts.
- **Outputs and safe claim:** Documents the actual execution order of the current assembly-to-results workflow.
- **Limitations:** It does not run the metadata builder, sequencing, basecalling, read QC/coverage, Longcycler/Flye assembly or polishing. Versions are not locked, unsupported scripts 15-21 are invoked, and the Parsnp branch can return success while leaving stale output.

## Optional numbered utilities outside the active root pipeline

These utilities are not part of the 46 active root numbered scripts and are not called by the main runner.

### `scripts/12d_wgs_badsize_rescue_screen.R`

- **Classification / runner status / unit:** Optional QC-rescue sensitivity; not runner-invoked; assembly candidate.
- **Purpose / inputs:** Screens assemblies excluded only for size with optional CheckM2, GUNC, QUAST, MLST, Mash and FCS-GX evidence, without changing canonical QC.
- **Filters and denominator:** Preserves primary canonical rows and adds only candidates meeting configured rescue evidence; default upper rescue size is 6.8 Mb.
- **Tools and tests:** Availability-dependent QC/taxonomy/distance tools; no phenotype model.
- **Repeated measures:** One sensitivity selection per episode.
- **Outputs and safe claim:** Rescue manifests, tool results, decisions and sensitivity selection; can support a separately labelled QC sensitivity.
- **Limitations:** Tool availability changes evidence, and rescued rows must never silently enter the primary denominator.

### `scripts/32_compare_primary_vs_rescue_vf.R`

- **Classification / runner status / unit:** Optional rescue comparison; not runner-invoked; primary/rescue and matched episode.
- **Purpose / inputs:** Compares canonical and rescue VF-ready denominators and matched VF profiles.
- **Filters, tools and tests:** Labels primary/rescue datasets and identifies rescue-only, primary-only and shared keys; descriptive only.
- **Repeated measures:** Matched episodes can be paired; no broader resident model.
- **Outputs and safe claim:** Rescue comparison tables/report; indicates whether the optional rescue changes calls or denominators.
- **Limitations:** Does not validate rescued assemblies or justify replacing the canonical analysis.

### `scripts/99_script_and_figure_index.R`

- **Classification / runner status / unit:** Documentation/meta utility; not runner-invoked; script file.
- **Purpose / inputs:** Parses first-100-line headers from root numbered and non-legacy utility scripts into an index.
- **Filters, tools and tests:** Filename-pattern and comment-tag heuristics; no statistical test.
- **Repeated measures:** Not applicable.
- **Outputs and safe claim:** `results/meta/script_and_figure_index.csv`; navigational aid only.
- **Limitations:** Header comments may be stale or incomplete and cannot replace this code-audited register.

## Excluded legacy scripts

Numbered files under `legacy/` and `legacy/unused_modules/` are explicitly excluded from the active methods register. They are archived or superseded and are not called by `RUN_COMPLETE_ANALYSIS.sh`. They should be mentioned only as historical provenance, not as part of the executed thesis method.

## Missing methods provenance to obtain

The scripts cannot supply the following details, so they must be recovered from laboratory records or the person who generated the upstream files:

- ONT platform, flow cell, library kit, basecaller version/model and read-QC/coverage criteria.
- Longcycler and Flye versions, commands, parameters and polishing steps.
- ABRicate and VFDB version/date.
- SeqSphere/provider MLST scheme version and source manifest.
- Prokka, Panaroo, Parsnp, MUMmer/dnadiff and Mash versions.
- Scientific citation or dataset-specific calibration for the operational <=25 dnadiff-SNP threshold.
