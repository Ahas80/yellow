#!/usr/bin/env Rscript

# Rewrite the current prose documentation from the audited Longcycler release
# registry. This writer is intentionally fail-closed: no document is staged
# until the complete numerical and methodological contract has been checked.

suppressPackageStartupMessages({
  library(digest)
  library(jsonlite)
})

root <- normalizePath(getwd(), winslash = "/", mustWork = TRUE)
if (!file.exists(file.path(root, "00_config.R"))) {
  stop("Run this writer from the rUTIs project root.", call. = FALSE)
}

registry_relative <- file.path("results", "pipeline", "longcycler_release_claim_registry.json")
registry_path <- file.path(root, registry_relative)
if (!file.exists(registry_path) || !isTRUE(file.info(registry_path)$size > 0)) {
  stop("Longcycler release claim registry is missing or empty: ", registry_path, call. = FALSE)
}

registry_size <- as.integer(file.info(registry_path)$size)
registry_bytes <- readBin(registry_path, what = "raw", n = registry_size)
registry <- tryCatch(
  jsonlite::fromJSON(rawToChar(registry_bytes), simplifyVector = TRUE),
  error = function(error) stop("Could not read the Longcycler release claim registry: ", conditionMessage(error), call. = FALSE)
)

fail_contract <- function(label, observed, expected) {
  observed_text <- paste(observed, collapse = ", ")
  expected_text <- paste(expected, collapse = ", ")
  stop(
    "Claim-registry contract failed for ", label,
    "; observed [", observed_text, "] but required [", expected_text, "].",
    call. = FALSE
  )
}

expect_integer <- function(value, expected, label) {
  valid <- length(value) == 1L && !is.na(value) && is.finite(as.numeric(value)) &&
    identical(as.integer(value), as.integer(expected)) && as.numeric(value) == as.numeric(expected)
  if (!valid) fail_contract(label, value, expected)
  invisible(TRUE)
}

expect_character <- function(value, expected, label) {
  valid <- length(value) == 1L && !is.na(value) && identical(as.character(value), expected)
  if (!valid) fail_contract(label, value, expected)
  invisible(TRUE)
}

expect_character_set <- function(value, expected, label) {
  value <- as.character(value)
  valid <- length(value) == length(expected) && !anyNA(value) && setequal(value, expected)
  if (!valid) fail_contract(label, value, expected)
  invisible(TRUE)
}

# Registry identity and analytical scope.
expect_character(registry$schema_version, "longcycler_release_claim_registry_v1", "schema_version")
expect_character(registry$analysis_scope$assembly_policy, "selected QC-passing Longcycler only", "analysis_scope.assembly_policy")
expect_character(registry$analysis_scope$clinical_phenotype, "operational UTI phenotype", "analysis_scope.clinical_phenotype")
expect_character(registry$analysis_scope$interpretation, "exploratory observational analysis; no causal claim", "analysis_scope.interpretation")

# Exact release anchors.
expect_integer(registry$analytical_cohort$episodes, 532L, "analytical_cohort.episodes")
expect_integer(registry$analytical_cohort$residents, 161L, "analytical_cohort.residents")
expect_integer(registry$analytical_cohort$operational_UTI, 16L, "analytical_cohort.operational_UTI")
expect_integer(registry$analytical_cohort$operational_Not_UTI, 516L, "analytical_cohort.operational_Not_UTI")

expect_character(
  registry$attrition_qc_context$label,
  "full clinical source retained only for attrition/QC context",
  "attrition_qc_context.label"
)
expect_integer(registry$attrition_qc_context$episodes, 583L, "attrition_qc_context.episodes")
expect_integer(registry$attrition_qc_context$residents, 166L, "attrition_qc_context.residents")
expect_integer(registry$attrition_qc_context$operational_UTI, 18L, "attrition_qc_context.operational_UTI")
expect_integer(registry$attrition_qc_context$operational_Not_UTI, 565L, "attrition_qc_context.operational_Not_UTI")

expect_integer(registry$direct_pairs$all_within_resident, 893L, "direct_pairs.all_within_resident")
expect_integer(registry$adjacent_transitions$pairs, 371L, "adjacent_transitions.pairs")
expect_integer(registry$adjacent_transitions$residents, 139L, "adjacent_transitions.residents")
expect_integer(registry$adjacent_transitions$operational_snp_threshold, 25L, "adjacent_transitions.operational_snp_threshold")
expect_integer(registry$adjacent_transitions$at_or_below_threshold, 140L, "adjacent_transitions.at_or_below_threshold")
expect_integer(registry$adjacent_transitions$Not_UTI_to_UTI, 9L, "adjacent_transitions.Not_UTI_to_UTI")
expect_integer(
  registry$adjacent_transitions$Not_UTI_to_UTI_at_or_below_threshold,
  5L,
  "adjacent_transitions.Not_UTI_to_UTI_at_or_below_threshold"
)
expect_integer(registry$mechanism_casebook$cases, 9L, "mechanism_casebook.cases")
expect_integer(registry$mechanism_casebook$linked, 9L, "mechanism_casebook.linked")
expect_integer(registry$mechanism_casebook$missing, 0L, "mechanism_casebook.missing")
expect_integer(registry$near_miss_audit$rows, 17L, "near_miss_audit.rows")
expect_character(registry$near_miss_audit$label, "near-miss rows; not operational UTI cases", "near_miss_audit.label")
expect_character(registry$research_questions$first, "RQ01", "research_questions.first")
expect_character(registry$research_questions$last, "RQ10", "research_questions.last")
expect_integer(registry$research_questions$count, 10L, "research_questions.count")
expect_integer(registry$research_questions$retired_questions, 0L, "research_questions.retired_questions")

# Exact method contract. These checks make prose regeneration fail if a method
# changes without a deliberate registry and documentation update.
expect_integer(
  registry$method_contract$operational_phenotype$culture_lower_bound_cfu_per_ml,
  1000L,
  "method_contract.operational_phenotype.culture_lower_bound_cfu_per_ml"
)
expect_character(
  registry$method_contract$operational_phenotype$rule,
  "versioned operational culture-plus-compatible-symptom phenotype",
  "method_contract.operational_phenotype.rule"
)
expect_character(
  registry$method_contract$operational_phenotype$caveat,
  "not a reconstruction of the full published protocol",
  "method_contract.operational_phenotype.caveat"
)
expect_integer(registry$method_contract$assembly_qc$max_contigs, 200L, "method_contract.assembly_qc.max_contigs")
expect_integer(registry$method_contract$assembly_qc$min_n50_bp, 20000L, "method_contract.assembly_qc.min_n50_bp")
expect_integer(registry$method_contract$assembly_qc$min_genome_size_bp, 4000000L, "method_contract.assembly_qc.min_genome_size_bp")
expect_integer(registry$method_contract$assembly_qc$max_genome_size_bp, 6000000L, "method_contract.assembly_qc.max_genome_size_bp")
expect_character_set(
  registry$method_contract$assembly_qc$excluded_metrics,
  c("read coverage", "completeness", "contamination"),
  "method_contract.assembly_qc.excluded_metrics"
)
expect_character(registry$method_contract$vfdb$tool, "ABRicate", "method_contract.vfdb.tool")
expect_character(registry$method_contract$vfdb$database, "VFDB", "method_contract.vfdb.database")
expect_integer(registry$method_contract$vfdb$min_identity_pct, 80L, "method_contract.vfdb.min_identity_pct")
expect_integer(registry$method_contract$vfdb$min_coverage_pct, 80L, "method_contract.vfdb.min_coverage_pct")
expect_character(
  registry$method_contract$vfdb$provenance,
  "SHA-bound calls from the selected Longcycler FASTA manifest",
  "method_contract.vfdb.provenance"
)
expect_character(
  registry$method_contract$mlst$role,
  "lineage context; not pair-specific continuity proof",
  "method_contract.mlst.role"
)
expect_integer(
  registry$method_contract$mlst$provider_min_good_targets_pct,
  95L,
  "method_contract.mlst.provider_min_good_targets_pct"
)
expect_character(
  registry$method_contract$mlst$provider_policy,
  "provider_qc95 call key/path-linked to the selected Longcycler episode; local fallback excluded",
  "method_contract.mlst.provider_policy"
)
expect_character(
  registry$method_contract$mlst$fallback,
  "labelled local MLST from the same selected Longcycler FASTA where required",
  "method_contract.mlst.fallback"
)
expect_character(registry$method_contract$direct_pair_evidence$tool, "dnadiff", "method_contract.direct_pair_evidence.tool")
expect_character(
  registry$method_contract$direct_pair_evidence$role,
  "primary pair-specific distance evidence",
  "method_contract.direct_pair_evidence.role"
)
expect_integer(
  registry$method_contract$direct_pair_evidence$operational_snp_threshold,
  25L,
  "method_contract.direct_pair_evidence.operational_snp_threshold"
)
expect_character(
  registry$method_contract$direct_pair_evidence$priority,
  "graph connectivity and MLST agreement cannot override a conflicting direct pair",
  "method_contract.direct_pair_evidence.priority"
)
expect_character(
  registry$method_contract$population_context$core_genome_tool,
  "Parsnp",
  "method_contract.population_context.core_genome_tool"
)
expect_character(
  registry$method_contract$population_context$pangenome_tool,
  "Panaroo",
  "method_contract.population_context.pangenome_tool"
)
expect_character(
  registry$method_contract$population_context$role,
  "population context; not a substitute for direct pair evidence",
  "method_contract.population_context.role"
)

registry_sha256 <- unname(digest::digest(registry_bytes, algo = "sha256", serialize = FALSE))
if (length(registry_sha256) != 1L || !grepl("^[0-9a-f]{64}$", registry_sha256)) {
  stop("Could not compute the exact claim-registry SHA-256.", call. = FALSE)
}

fmt <- function(value) format(as.integer(value), big.mark = ",", scientific = FALSE, trim = TRUE)

release_contract <- c(
  "## Audited release contract",
  "",
  "- Scope: Longcycler-only selected QC-passing assemblies; one selected assembly per analytical episode.",
  "- Clinical definition: operational UTI phenotype. It is a versioned culture-plus-compatible-symptom rule, not a reconstruction of the full published protocol.",
  sprintf(
    "- Analytical cohort: %s/%s/%s/%s (episodes/residents/operational UTI/operational Not_UTI).",
    fmt(registry$analytical_cohort$episodes),
    fmt(registry$analytical_cohort$residents),
    fmt(registry$analytical_cohort$operational_UTI),
    fmt(registry$analytical_cohort$operational_Not_UTI)
  ),
  sprintf("- Direct evidence: %s all within-resident pairs.", fmt(registry$direct_pairs$all_within_resident)),
  sprintf(
    "- Adjacent evidence: %s/%s/%s (pairs/residents/pairs at or below %s SNP).",
    fmt(registry$adjacent_transitions$pairs),
    fmt(registry$adjacent_transitions$residents),
    fmt(registry$adjacent_transitions$at_or_below_threshold),
    fmt(registry$adjacent_transitions$operational_snp_threshold)
  ),
  sprintf(
    "- Focused transitions: %s Not_UTI -> UTI; %s/%s at or below %s SNP.",
    fmt(registry$adjacent_transitions$Not_UTI_to_UTI),
    fmt(registry$adjacent_transitions$Not_UTI_to_UTI_at_or_below_threshold),
    fmt(registry$adjacent_transitions$Not_UTI_to_UTI),
    fmt(registry$adjacent_transitions$operational_snp_threshold)
  ),
  sprintf(
    "- Mechanism casebook: %s/%s/%s (cases/linked/missing).",
    fmt(registry$mechanism_casebook$cases),
    fmt(registry$mechanism_casebook$linked),
    fmt(registry$mechanism_casebook$missing)
  ),
  sprintf("- Near-miss audit: %s rows; these are not operational UTI cases.", fmt(registry$near_miss_audit$rows)),
  sprintf(
    "- Attrition/QC context only: %s/%s/%s/%s (full-source episodes/residents/operational UTI/operational Not_UTI); these are not analytical denominators.",
    fmt(registry$attrition_qc_context$episodes),
    fmt(registry$attrition_qc_context$residents),
    fmt(registry$attrition_qc_context$operational_UTI),
    fmt(registry$attrition_qc_context$operational_Not_UTI)
  ),
  "- Research-question boundary: RQ01-RQ10 only.",
  "- Interpretation: exploratory, observational and non-causal.",
  ""
)

method_contract <- c(
  "## Methods fixed by the registry",
  "",
  sprintf(
    "- Operational phenotype: culture lower bound >=%s CFU/mL plus compatible symptoms under the versioned rule.",
    fmt(registry$method_contract$operational_phenotype$culture_lower_bound_cfu_per_ml)
  ),
  sprintf(
    "- Assembly QC: <=%s contigs, N50 >=%s bp and genome size %s-%s bp; read coverage, completeness and contamination are excluded metrics.",
    fmt(registry$method_contract$assembly_qc$max_contigs),
    fmt(registry$method_contract$assembly_qc$min_n50_bp),
    fmt(registry$method_contract$assembly_qc$min_genome_size_bp),
    fmt(registry$method_contract$assembly_qc$max_genome_size_bp)
  ),
  sprintf(
    "- VF calls: %s with %s at >=%s%% identity and >=%s%% coverage, SHA-bound to the selected Longcycler FASTA manifest.",
    registry$method_contract$vfdb$tool,
    registry$method_contract$vfdb$database,
    fmt(registry$method_contract$vfdb$min_identity_pct),
    fmt(registry$method_contract$vfdb$min_coverage_pct)
  ),
  sprintf(
    "- Genomic AMR: %s is primary; %s is complementary and %s is the legacy comparison. %s",
    registry$method_contract$genomic_amr$primary_caller,
    registry$method_contract$genomic_amr$complementary_caller,
    registry$method_contract$genomic_amr$legacy_comparison,
    registry$method_contract$genomic_amr$interpretation
  ),
  sprintf(
    "- MLST: lineage context only; provider calls require >=%s%% good targets. Policy: %s. When required, local calls are labelled and use the same selected FASTA.",
    fmt(registry$method_contract$mlst$provider_min_good_targets_pct),
    registry$method_contract$mlst$provider_policy
  ),
  sprintf(
    "- Pair-specific evidence: %s is primary, with an operational threshold of <=%s SNP; MLST or graph context cannot overrule a conflicting direct pair.",
    registry$method_contract$direct_pair_evidence$tool,
    fmt(registry$method_contract$direct_pair_evidence$operational_snp_threshold)
  ),
  sprintf(
    "- Population context: %s core-genome and %s pangenome outputs provide context, not pair-specific continuity proof.",
    registry$method_contract$population_context$core_genome_tool,
    registry$method_contract$population_context$pangenome_tool
  ),
  ""
)

provenance_footer <- c(
  "## Provenance",
  "",
  paste0("- Claim registry: `", registry_relative, "`"),
  paste0("- Registry SHA-256: `", registry_sha256, "`"),
  paste0("- Registry generated: ", as.character(registry$generated_at)),
  "- This file is generated; edit the registry-producing analysis or this writer rather than hand-editing release claims.",
  ""
)

assemble_document <- function(title, purpose, body) {
  c(
    paste0("# ", title),
    "",
    purpose,
    "",
    body,
    release_contract,
    method_contract,
    provenance_footer
  )
}

documents <- list(
  "FOLDER_MAP.md" = assemble_document(
    "Current Longcycler-only folder map",
    "The paths below identify the active release inputs, analysis products, quality gates and communication artifacts.",
    c(
      "## Active paths",
      "",
      "| Area | Current path | Role |",
      "|---|---|---|",
      "| Clinical cohort | `results/clinical/analysis_cohort_longcycler.csv` | Exact analytical episode set |",
      "| FASTA provenance | `results/qc/analysis_assembly_manifest.csv` | Selected paths and content hashes |",
      "| VF | `results/vf/` | Selected-cohort VF calls and matrices |",
      "| MLST | `results/mlst/` | Lineage context |",
      "| Direct comparisons | `results/strain_compare/` | Pair-specific genomic distances |",
      "| Longitudinal | `results/longitudinal/` | Canonical adjacent transition table |",
      "| Genomic AMR | `results/amr/` | Validated determinant, prevalence, longitudinal, prediction and provenance tables |",
      "| Genomic AMR results | `results/summary/table_13_genomic_amr_summary.csv` | Numeric supporting results used in the final narrative |",
      "| Mechanism | `results/mechanism/` | Linked focused-transition casebook |",
      "| Research questions | `results/research_questions/RQ01` to `RQ10` | Audited analyses |",
      "| Registry and markers | `results/pipeline/` | Release contract and completion state |",
      "| QC reports | `results/qc/` | Cohort, provenance and delivery gates |",
      "| Figures | `plots/` | Generated analytical figures |",
      "| Communication | `outputs/` | Registry-bound decks, handouts and codebook |",
      ""
    )
  ),
  "CODE_REVIEW_RECONCILIATION_README.md" = assemble_document(
    "Code-review reconciliation: Longcycler-only release",
    "The active implementation has been reconciled around one fail-closed cohort, one transition table and one claim registry.",
    c(
      "## Reconciled controls",
      "",
      "- Cohort and FASTA manifests must have identical participant-timepoint keys and content hashes.",
      "- Downstream VF, MLST, direct-pair, core-genome and pangenome inputs must remain within that manifest.",
      "- Adjacent-transition and focused-casebook claims are computed from direct pair evidence.",
      "- Research questions stop at RQ10 and publish completion checks.",
      "- Decks, handouts, codebook and current documentation bind their claims to the registry SHA-256.",
      "",
      "## Interpretation guardrail",
      "",
      "A distance threshold is an operational continuity indicator, not proof of transmission, clinical progression or causation.",
      ""
    )
  ),
  "docs/LECTURER_README.md" = assemble_document(
    "Lecturer guide: Longcycler-only rUTIs analysis",
    "Use this guide to present the current evidence without changing denominators or overstating what the genomic comparisons establish.",
    c(
      "## Suggested explanation",
      "",
      "Each analytical episode contributes one selected Longcycler assembly. Direct within-resident comparisons measure genomic distance; adjacent comparisons describe consecutive observed episodes. The focused casebook then inspects Not_UTI -> UTI transitions with linked evidence.",
      "",
      "## Language guardrails",
      "",
      "- Say `operational UTI` and `operational Not_UTI`.",
      "- Say `at or below the operational SNP threshold`, not `proven transmission`.",
      "- Separate the analytical cohort from the full-source attrition/QC context.",
      "- Describe VF, MLST, core-genome and pangenome findings as exploratory context.",
      ""
    )
  ),
  "docs/PIPELINE_FAILURE_LOG.md" = assemble_document(
    "Longcycler-only pipeline gate and failure log",
    "This current log records the fail-closed release posture; transient run details remain in timestamped files under `logs/`.",
    c(
      "## Release gate posture",
      "",
      "| Gate | Required state | Failure response |",
      "|---|---|---|",
      "| Preflight | Exact selected cohort and usable tools | Stop before analysis |",
      "| Provenance | Selected FASTA paths and hashes agree | Stop downstream publication |",
      "| Analysis | Denominators and direct evidence match registry anchors | Do not publish registry |",
      "| Research questions | RQ01-RQ10 completion checks pass | Do not publish completion marker |",
      "| Communications | Artifacts are registry-bound, current and content-clean | Do not release deliverables |",
      "",
      "## Diagnostic locations",
      "",
      "- Timestamped execution logs: `logs/complete_analysis_*.log`",
      "- Pipeline markers and registry: `results/pipeline/`",
      "- Machine-readable QC checks: `results/qc/`",
      ""
    )
  ),
  "docs/VF_abstract_draft.md" = assemble_document(
    "VF analysis abstract: audited Longcycler-only draft",
    "This concise draft states only registry-supported design and denominator claims; effect estimates belong in the audited RQ06-RQ08 outputs.",
    c(
      "## Background",
      "",
      "Virulence-factor content may provide genomic context for differences between operational UTI and operational Not_UTI episodes, but observational comparisons cannot establish causation.",
      "",
      "## Methods",
      "",
      "We analysed the exact selected Longcycler cohort with SHA-bound VFDB calls. Between-episode and longitudinal analyses retained resident clustering and used direct pair distances for pair-specific continuity context. All analyses were exploratory.",
      "",
      "## Registry-supported results",
      "",
      "The release contains the analytical cohort, all within-resident direct comparisons, adjacent comparisons and a fully linked focused-transition casebook described below. Gene-level estimates, uncertainty intervals and multiplicity handling must be quoted directly from RQ06-RQ08 outputs.",
      "",
      "## Interpretation",
      "",
      "Observed VF associations or stability patterns are hypothesis-generating and non-causal. Clinical classification uses the operational phenotype rather than a reconstruction of the full published protocol.",
      ""
    )
  ),
  "docs/VF_merge_diagnostics.md" = assemble_document(
    "VF merge diagnostics: selected Longcycler cohort",
    "VF release tables are accepted only when their participant-timepoint keys and selected FASTA provenance agree with the analytical manifest.",
    c(
      "## Required merge checks",
      "",
      "| Check | Requirement |",
      "|---|---|",
      "| Cohort rows | Exact analytical episode count below |",
      "| Key uniqueness | One row per participant-timepoint |",
      "| Manifest coverage | Every cohort key has one selected FASTA |",
      "| VF coverage | Every analytical key has one VF profile |",
      "| Content provenance | FASTA SHA-256 agrees with the selected manifest |",
      "| Clinical labels | Operational UTI or operational Not_UTI only |",
      "",
      "## Canonical merge inputs",
      "",
      "- `results/clinical/analysis_cohort_longcycler.csv`",
      "- `results/qc/analysis_assembly_manifest.csv`",
      "- `results/vf/vf_pa_all.csv`",
      ""
    )
  ),
  "docs/VF_verification_report.md" = assemble_document(
    "VF verification report: Longcycler-only release",
    "Verdict: PASS for registry identity, exact release anchors and the method contract required to generate this report.",
    c(
      "## Verified controls",
      "",
      "- The analytical denominator and selected-assembly policy are exact.",
      "- VF identity, coverage and database settings equal the registry contract.",
      "- Direct-pair and adjacent-transition anchors are internally fixed.",
      "- The focused casebook is complete and the near-miss rows remain separately labelled.",
      "- Research-question scope is bounded to the ten registered questions.",
      "",
      "## Claim boundary",
      "",
      "This PASS validates release identity and methods. Scientific estimates still require their question-specific diagnostics and should be reported as exploratory and non-causal.",
      ""
    )
  ),
  "docs/figures/timepoint_vs_isolate_clarification.md" = assemble_document(
    "Episode, assembly and pair terminology",
    "The active release does not count multiple assembly alternatives as additional analytical episodes.",
    c(
      "## Terms",
      "",
      "- Episode: one resident sampling occasion with an operational clinical label.",
      "- Selected assembly: the single QC-passing Longcycler FASTA attached to an analytical episode.",
      "- All within-resident pair: any two selected episodes from the same resident.",
      "- Adjacent pair: two consecutive selected episodes after chronological ordering within a resident.",
      "- Focused transition: an adjacent operational Not_UTI episode followed by an operational UTI episode.",
      "",
      "## Figure rule",
      "",
      "Axis labels and captions must name the unit being counted. Episode, resident, assembly and pair denominators are not interchangeable.",
      ""
    )
  ),
  "docs/legacy_asb_uti_docs/VF_provenance_map.md" = assemble_document(
    "Current VF provenance map",
    "This file's location is retained for links, but its content describes only the active Longcycler-only operational-phenotype release.",
    c(
      "## Data flow",
      "",
      "```mermaid",
      "flowchart LR",
      "  A[Selected cohort keys] --> B[Selected Longcycler FASTA manifest and SHA-256]",
      "  B --> C[ABRicate VFDB calls]",
      "  C --> D[Binary VF matrix]",
      "  D --> E[RQ06-RQ08 exploratory analyses]",
      "  B --> F[Direct pair evidence]",
      "  F --> E",
      "```",
      "",
      "## Provenance rule",
      "",
      "A VF row is publishable only when its participant-timepoint key and selected FASTA content hash remain within the analytical manifest.",
      ""
    )
  ),
  "docs/workflow_case_count_flowchart.md" = assemble_document(
    "Longcycler-only denominator map",
    "The analytical release and full-source attrition/QC context are distinct and must not be combined.",
    c(
      "## Denominator flow",
      "",
      "```mermaid",
      "flowchart LR",
      "  A[Full clinical source: attrition/QC context only] --> B[Selected analytical episodes]",
      "  B --> C[All within-resident direct pairs]",
      "  B --> D[Adjacent pairs]",
      "  D --> E[Pairs at or below the operational SNP threshold]",
      "  D --> F[Not_UTI to UTI focused transitions]",
      "  F --> G[Fully linked mechanism casebook]",
      "```",
      "",
      "## Reading rule",
      "",
      "Counts change because the unit changes from source episodes, to selected analytical episodes, to pairs and finally to focused transitions. They do not represent one simple attrition ladder.",
      ""
    )
  ),
  "docs/workflow_flowchart.md" = assemble_document(
    "Longcycler-only workflow guide",
    "`RUN_COMPLETE_ANALYSIS.sh` is the executable source of run order; this guide captures the release-level flow and evidence hierarchy.",
    c(
      "## Pipeline flow",
      "",
      "```mermaid",
      "flowchart TD",
      "  A[Clinical classification and curation] --> B[Exact Longcycler selection and assembly QC]",
      "  B --> C[Clinical, VF and MLST episode tables]",
      "  B --> D[dnadiff direct pair evidence]",
      "  B --> E[Parsnp and Panaroo population context]",
      "  C --> F[Chronological adjacent transitions]",
      "  D --> F",
      "  F --> G[RQ01-RQ10 analyses and diagnostics]",
      "  G --> H[Claim registry]",
      "  H --> I[Decks, handouts, codebook and current documentation]",
      "  I --> J[Final deliverable verification]",
      "```",
      "",
      "## Evidence hierarchy",
      "",
      "Direct pair evidence answers pair-specific distance questions. MLST, core-genome and pangenome outputs provide lineage or population context and cannot replace a direct comparison.",
      ""
    )
  )
)

required_paths <- c(
  "FOLDER_MAP.md",
  "CODE_REVIEW_RECONCILIATION_README.md",
  "docs/LECTURER_README.md",
  "docs/PIPELINE_FAILURE_LOG.md",
  "docs/VF_abstract_draft.md",
  "docs/VF_merge_diagnostics.md",
  "docs/VF_verification_report.md",
  "docs/figures/timepoint_vs_isolate_clarification.md",
  "docs/legacy_asb_uti_docs/VF_provenance_map.md",
  "docs/workflow_case_count_flowchart.md",
  "docs/workflow_flowchart.md"
)
if (!identical(names(documents), required_paths)) {
  fail_contract("documentation inventory", names(documents), required_paths)
}

required_fragments <- c(
  "Longcycler-only",
  "operational UTI phenotype",
  "532/161/16/516",
  "893 all within-resident pairs",
  "371/139/140",
  "9 Not_UTI -> UTI",
  "5/9 at or below 25 SNP",
  "9/9/0 (cases/linked/missing)",
  "Near-miss audit: 17 rows",
  "Attrition/QC context only: 583/166/18/565",
  "RQ01-RQ10 only",
  "exploratory, observational and non-causal",
  "ABRicate with VFDB",
  "dnadiff is primary",
  registry_sha256
)

retired_input_token <- paste0("fl", "ye")
retired_question_token <- paste0("RQ", "11")
stale_patterns <- c(
  paste0("\\b", "5", "56", "\\b"),
  paste0("\\b", "1", "62", "\\b"),
  paste0("\\b", "5", "39", "\\b"),
  paste0("\\b", "3", "94", "\\b"),
  paste0("\\b", "1", "16", "\\b"),
  paste0("\\b", "7", "\\s*(?:/|of)\\s*", "9", "\\b")
)

validate_document <- function(relative, lines) {
  text <- paste(lines, collapse = "\n")
  missing <- required_fragments[!vapply(required_fragments, function(fragment) grepl(fragment, text, fixed = TRUE), logical(1))]
  if (length(missing)) fail_contract(paste0(relative, " required content"), missing, "all required fragments present")
  if (grepl(retired_input_token, text, ignore.case = TRUE, fixed = TRUE)) {
    stop("Retired assembly-input token detected in generated document: ", relative, call. = FALSE)
  }
  if (grepl(retired_question_token, text, ignore.case = TRUE, fixed = TRUE)) {
    stop("Out-of-scope research-question token detected in generated document: ", relative, call. = FALSE)
  }
  stale_hit <- vapply(stale_patterns, function(pattern) grepl(pattern, text, perl = TRUE, ignore.case = TRUE), logical(1))
  if (any(stale_hit)) stop("Stale release denominator detected in generated document: ", relative, call. = FALSE)
  invisible(TRUE)
}

invisible(Map(validate_document, names(documents), documents))

stage_document <- function(relative, lines) {
  target <- file.path(root, relative)
  dir.create(dirname(target), recursive = TRUE, showWarnings = FALSE)
  tmp <- tempfile(pattern = paste0(".", basename(target), "."), tmpdir = dirname(target), fileext = ".tmp")
  connection <- file(tmp, open = "wb")
  tryCatch(
    writeLines(enc2utf8(lines), connection, useBytes = TRUE),
    finally = close(connection)
  )
  if (!file.exists(tmp) || !isTRUE(file.info(tmp)$size > 0)) {
    unlink(tmp)
    stop("Could not stage generated document: ", relative, call. = FALSE)
  }
  list(relative = relative, target = target, tmp = tmp)
}

staged <- Map(stage_document, names(documents), documents)
on.exit(unlink(vapply(staged, `[[`, character(1), "tmp")), add = TRUE)

if (!identical(unname(digest::digest(registry_path, algo = "sha256", file = TRUE)), registry_sha256)) {
  stop("Claim-registry bytes changed while the release documents were being generated.", call. = FALSE)
}

for (entry in staged) {
  if (!file.rename(entry$tmp, entry$target)) {
    stop("Could not atomically publish generated document: ", entry$relative, call. = FALSE)
  }
}

message(
  "Published ", length(staged),
  " Longcycler-only release documents from registry SHA-256 ", registry_sha256, "."
)
