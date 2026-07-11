#!/usr/bin/env Rscript

# Independent audit of:
# Hamdi Hersi, "Longitudinal Dynamics of Asymptomatic Bacteriuria in
# Nursing Home Residents: Genomic Persistence and Diversity of Escherichia coli"
#
# This script deliberately does not attempt to reconstruct wgMLST, antibiotic,
# demographic, or host-factor results when their source data are unavailable.
# SNP results are reported as independent, non-equivalent corroboration only.

options(stringsAsFactors = FALSE, scipen = 999)

args <- commandArgs(trailingOnly = TRUE)

arg_value <- function(flag, default) {
  hit <- which(args == flag)
  if (length(hit) == 0) return(default)
  if (hit[[1]] == length(args)) stop("Missing value after ", flag)
  args[[hit[[1]] + 1]]
}

results_dir <- arg_value("--results-dir", "results")
out_dir <- arg_value("--out-dir", file.path(results_dir, "thesis_audit"))

paths <- list(
  clinical = file.path(results_dir, "clinical", "status_map.csv"),
  vf = file.path(results_dir, "vf", "vf_analysis_ready.csv"),
  transitions = file.path(results_dir, "vf", "vf_transitions_stratified.csv"),
  pairwise = file.path(results_dir, "strain_compare", "pairwise_metrics.csv"),
  assemblies = file.path(results_dir, "qc", "canonical_assembly_selection.csv")
)

missing_inputs <- names(paths)[!file.exists(unlist(paths))]
if (length(missing_inputs) > 0) {
  stop(
    "Required audit inputs are missing: ",
    paste(sprintf("%s (%s)", missing_inputs, unlist(paths[missing_inputs])), collapse = "; ")
  )
}

dir.create(out_dir, recursive = TRUE, showWarnings = FALSE)

read_anchor <- function(path) {
  read.csv(
    path,
    check.names = FALSE,
    na.strings = c("", "NA"),
    stringsAsFactors = FALSE
  )
}

clinical_all <- read_anchor(paths$clinical)
vf <- read_anchor(paths$vf)
transitions_all <- read_anchor(paths$transitions)
pairwise <- read_anchor(paths$pairwise)
assemblies_all <- read_anchor(paths$assemblies)

required_pairwise_provenance <- c(
  "Selected_canonical_A", "Selected_canonical_B", "QC_PASS_A", "QC_PASS_B",
  "within_participant", "Fasta_SHA256_A", "Fasta_SHA256_B",
  "dnadiff_report_path", "dnadiff_sidecar_path", "dnadiff_report_sha256",
  "dnadiff_cache_signature", "dnadiff_version"
)
missing_pairwise_fields <- setdiff(required_pairwise_provenance, names(pairwise))
if (length(missing_pairwise_fields) > 0) {
  stop(
    "Pairwise input lacks required cache-safe provenance fields: ",
    paste(missing_pairwise_fields, collapse = ", ")
  )
}
if (any(!pairwise$Selected_canonical_A %in% TRUE) ||
    any(!pairwise$Selected_canonical_B %in% TRUE) ||
    any(!pairwise$QC_PASS_A %in% TRUE) ||
    any(!pairwise$QC_PASS_B %in% TRUE) ||
    any(!pairwise$within_participant %in% TRUE)) {
  stop("Pairwise input contains a noncanonical, QC-failing or between-participant endpoint.")
}
nonempty_pairwise_fields <- c(
  "Fasta_SHA256_A", "Fasta_SHA256_B", "dnadiff_report_path",
  "dnadiff_sidecar_path", "dnadiff_report_sha256",
  "dnadiff_cache_signature", "dnadiff_version"
)
for (field in nonempty_pairwise_fields) {
  if (any(is.na(pairwise[[field]]) | !nzchar(as.character(pairwise[[field]])))) {
    stop("Pairwise provenance field `", field, "` contains a missing value.")
  }
}
if (any(!file.exists(pairwise$dnadiff_report_path)) ||
    any(!file.exists(pairwise$dnadiff_sidecar_path))) {
  stop("A pairwise dnadiff report or provenance sidecar referenced by the input is missing.")
}

if (!"cohort" %in% names(transitions_all)) {
  stop("Transition input lacks the required `cohort` field.")
}
transitions <- transitions_all[transitions_all$cohort == "all", , drop = FALSE]
if (nrow(transitions) == 0) stop("No transition rows remained after filtering `cohort == \"all\"`.")

clinical <- clinical_all[clinical_all$analysis_include_primary %in% TRUE, , drop = FALSE]
assemblies <- assemblies_all[assemblies_all$selected_canonical %in% TRUE, , drop = FALSE]

safe_unique_n <- function(x) length(unique(x[!is.na(x)]))
num <- function(x) suppressWarnings(as.numeric(x))
pct <- function(n, d, digits = 1) {
  if (is.na(d) || d == 0) return(NA_real_)
  round(100 * n / d, digits)
}
fmt_pct <- function(n, d, digits = 1) sprintf("%d/%d (%.1f%%)", n, d, pct(n, d, digits))
fmt_num <- function(x, digits = 1) format(round(x, digits), nsmall = digits, trim = TRUE)

# ---- Current independently computed anchors ---------------------------------

n_clinical <- nrow(clinical)
n_clinical_participants <- safe_unique_n(clinical$Participant_id)
n_clinical_all_participants <- safe_unique_n(clinical_all$Participant_id)
n_vf <- nrow(vf)
n_vf_participants <- safe_unique_n(vf$Participant_id)
n_batches <- safe_unique_n(vf$Batch)
n_st <- safe_unique_n(vf$ST)

tp_counts <- as.numeric(table(vf$Participant_id))
tp_quantiles <- quantile(tp_counts, probs = c(0, 0.25, 0.5, 0.75, 1), names = FALSE)

assembler_counts <- table(tolower(assemblies$assembler), useNA = "ifany")
n_flye <- if ("flye" %in% names(assembler_counts)) unname(assembler_counts[["flye"]]) else 0
n_longcycler <- if ("longcycler" %in% names(assembler_counts)) {
  unname(assembler_counts[["longcycler"]])
} else {
  0
}

st_snapshot <- function(tp, st) {
  x <- vf[vf$tp_lab == tp, , drop = FALSE]
  typed <- !is.na(x$ST)
  n_st_here <- sum(as.character(x$ST) == as.character(st), na.rm = TRUE)
  list(
    timepoint = tp,
    st = st,
    n = n_st_here,
    n_all = nrow(x),
    n_typed = sum(typed),
    pct_all = pct(n_st_here, nrow(x)),
    pct_typed = pct(n_st_here, sum(typed))
  )
}

st131_t0 <- st_snapshot("T0", 131)
st131_t6 <- st_snapshot("T6", 131)
st73_t0 <- st_snapshot("T0", 73)
st73_t6 <- st_snapshot("T6", 73)

n_transitions <- nrow(transitions)
n_transition_participants <- safe_unique_n(transitions$Participant_id)
strict_snp_thresholds <- unique(num(transitions$same_strain_snp_threshold))
strict_snp_thresholds <- strict_snp_thresholds[is.finite(strict_snp_thresholds)]
if (length(strict_snp_thresholds) != 1L) {
  stop(
    "Expected exactly one finite operational SNP threshold in the transition input; observed: ",
    paste(strict_snp_thresholds, collapse = ", ")
  )
}
strict_snp_threshold <- strict_snp_thresholds[[1]]
if (!identical(strict_snp_threshold, 25)) {
  stop("The audit is specified for an operational <=25-SNP threshold; observed ", strict_snp_threshold, ".")
}

is_strict_snp <- function(x) {
  x_num <- num(x)
  !is.na(x_num) & x_num <= strict_snp_threshold
}
is_above_snp <- function(x) {
  x_num <- num(x)
  !is.na(x_num) & x_num > strict_snp_threshold
}

n_strict_same <- sum(is_strict_snp(transitions$SNPs))
n_above_snp <- sum(is_above_snp(transitions$SNPs))
n_replacement <- sum(transitions$replacement_flag %in% TRUE, na.rm = TRUE)
n_not_confident_replacement <- sum(transitions$replacement_flag %in% FALSE, na.rm = TRUE)

assert_anchor <- function(actual, expected, label) {
  if (!identical(as.integer(actual), as.integer(expected))) {
    stop(
      "Current snapshot anchor failed for ", label, ": expected ",
      expected, ", observed ", actual,
      ". Review whether the authoritative June 2026 outputs changed."
    )
  }
}

assert_anchor(n_vf, 556L, "canonical WGS/VF rows")
assert_anchor(n_st, 83L, "unique typed STs")
assert_anchor(n_transitions, 394L, "cohort-all WGS-linked transitions")

context_strict_same <- sum(
  transitions$snp_strain_context == "Strong same strain",
  na.rm = TRUE
)
context_above_snp <- sum(
  transitions$snp_strain_context == "Above same-strain SNP threshold",
  na.rm = TRUE
)
if (n_strict_same != context_strict_same || n_above_snp != context_above_snp) {
  stop(
    "Transition SNP labels disagree with the current numeric <=", strict_snp_threshold,
    " rule: numeric strict/above = ", n_strict_same, "/", n_above_snp,
    "; labelled strict/above = ", context_strict_same, "/", context_above_snp, "."
  )
}

routine_to_event <- transitions[
  transitions$event_type_from == "Routine" &
    transitions$event_type_to == "UTI_event",
  ,
  drop = FALSE
]
routine_to_primary_uti <- routine_to_event[
  routine_to_event$status_from == "Not_UTI" &
    routine_to_event$status_to == "UTI",
  ,
  drop = FALSE
]
all_primary_notuti_uti <- transitions[
  transitions$status_from == "Not_UTI" &
    transitions$status_to == "UTI",
  ,
  drop = FALSE
]

transition_summary <- function(x) {
  list(
    n = nrow(x),
    participants = safe_unique_n(x$Participant_id),
    strict_same = sum(is_strict_snp(x$SNPs)),
    above_snp = sum(is_above_snp(x$SNPs)),
    replacement = sum(x$replacement_flag %in% TRUE, na.rm = TRUE),
    not_confident_replacement = sum(x$replacement_flag %in% FALSE, na.rm = TRUE)
  )
}

event_sum <- transition_summary(routine_to_event)
primary_event_sum <- transition_summary(routine_to_primary_uti)
primary_all_sum <- transition_summary(all_primary_notuti_uti)

# ---- Pre-routine -> UTI event -> post-routine triples -----------------------

pre_uti <- routine_to_event
post_uti <- transitions[
  transitions$event_type_from == "UTI_event" &
    transitions$event_type_to == "Routine",
  ,
  drop = FALSE
]

triples <- merge(
  pre_uti,
  post_uti,
  by.x = c("Participant_id", "tp_to"),
  by.y = c("Participant_id", "tp_from"),
  suffixes = c("_preuti", "_postuti")
)

find_pair <- function(pid, tp_a, tp_b) {
  z <- pairwise[
    pairwise$Participant_id_A == pid &
      pairwise$Participant_id_B == pid &
      (
        (pairwise$Timepoint_A == tp_a & pairwise$Timepoint_B == tp_b) |
          (pairwise$Timepoint_A == tp_b & pairwise$Timepoint_B == tp_a)
      ),
    ,
    drop = FALSE
  ]
  if (nrow(z) == 0) {
    return(data.frame(
      TotalSNPs = NA_real_,
      snp_strain_context = NA_character_,
      st_lineage_context = NA_character_,
      replacement_flag = NA
    ))
  }
  relevant <- unique(z[, c(
    "TotalSNPs", "snp_strain_context",
    "st_lineage_context", "replacement_flag"
  )])
  if (nrow(relevant) > 1) {
    stop(
      "Conflicting pairwise rows for participant ", pid,
      " and timepoints ", tp_a, "/", tp_b
    )
  }
  relevant
}

if (nrow(triples) > 0) {
  pre_post <- do.call(
    rbind,
    lapply(seq_len(nrow(triples)), function(i) {
      find_pair(
        triples$Participant_id[[i]],
        triples$tp_from[[i]],
        triples$tp_to_postuti[[i]]
      )
    })
  )
  triples$pre_post_snps <- pre_post$TotalSNPs
  triples$pre_post_snp_context <- pre_post$snp_strain_context
  triples$pre_post_st_context <- pre_post$st_lineage_context
  triples$pre_post_replacement_flag <- pre_post$replacement_flag
} else {
  triples$pre_post_snps <- numeric()
  triples$pre_post_snp_context <- character()
  triples$pre_post_st_context <- character()
  triples$pre_post_replacement_flag <- logical()
}

triple_summary <- function(x) {
  list(
    n = nrow(x),
    participants = safe_unique_n(x$Participant_id),
    uti_strict_same = sum(is_strict_snp(x$SNPs_preuti)),
    uti_above_snp = sum(is_above_snp(x$SNPs_preuti)),
    original_strict_return = sum(is_strict_snp(x$pre_post_snps)),
    original_not_confident_replaced = sum(
      x$pre_post_replacement_flag %in% FALSE,
      na.rm = TRUE
    ),
    original_replacement = sum(
      x$pre_post_replacement_flag %in% TRUE,
      na.rm = TRUE
    )
  )
}

triple_all_sum <- triple_summary(triples)
triples_primary <- triples[
  triples$status_from_preuti == "Not_UTI" &
    triples$status_to_preuti == "UTI",
  ,
  drop = FALSE
]
triple_primary_sum <- triple_summary(triples_primary)

# ---- Claim register ---------------------------------------------------------

allowed_statuses <- c(
  "Verified from available data",
  "Independently corroborated",
  "Discrepant",
  "Methodologically unsupported",
  "Not checkable with available data"
)

claims <- list()

add_claim <- function(
    id, section, claim, thesis_result, thesis_denominator, thesis_method,
    required_data, independent_result, independent_denominator,
    comparison_type, status, rationale, source_file, filter, formula) {
  if (!status %in% allowed_statuses) stop("Invalid audit status for ", id, ": ", status)
  claims[[length(claims) + 1]] <<- data.frame(
    claim_id = id,
    section = section,
    thesis_claim = claim,
    thesis_result = thesis_result,
    thesis_denominator = thesis_denominator,
    thesis_method = thesis_method,
    required_data = required_data,
    independent_result = independent_result,
    independent_denominator = independent_denominator,
    comparison_type = comparison_type,
    audit_status = status,
    rationale = rationale,
    source_file = source_file,
    filter = filter,
    formula = formula,
    stringsAsFactors = FALSE
  )
}

add_claim(
  "C01", "Study population", "The genomic analysis contained 556 E. coli isolates.",
  "556 isolates", "All analysed genomes", "Dataset count",
  "Canonical genome manifest", sprintf("%d canonical WGS/VF rows", n_vf), as.character(n_vf),
  "Exact count", "Verified from available data",
  "The independently generated canonical WGS/VF table has the same row count.",
  paths$vf, "No additional filter", "nrow(vf_analysis_ready)"
)
add_claim(
  "C02", "Study population", "The 556 isolates came from 167 residents.",
  "167 residents", "556 isolates", "Unique resident count",
  "Thesis resident-to-isolate linkage", sprintf(
    "%d WGS-ready participants; %d primary clinical participants",
    n_vf_participants, n_clinical_participants
  ), sprintf("%d WGS-ready rows; %d clinical rows", n_vf, n_clinical),
  "Denominator reconciliation", "Discrepant",
  paste0(
    "The full current clinical table contains ", n_clinical_all_participants,
    " IDs before primary exclusions, but the analysed canonical WGS table links ",
    n_vf, " genomes to ", n_vf_participants, " participants. The primary clinical table contains ",
    n_clinical_participants, " participants. The thesis must show how unlinked/unknown IDs were counted."
  ),
  paste(paths$vf, paths$clinical, sep = "; "),
  "WGS-ready rows; clinical analysis_include_primary == TRUE",
  "n_distinct(Participant_id) in each analysis universe"
)
add_claim(
  "C03", "Study population", "The isolates came from six sequencing batches.",
  "6 batches", "556 isolates", "Unique batch count",
  "Canonical WGS/VF table", sprintf("%d batches", n_batches), as.character(n_vf),
  "Exact count", "Verified from available data",
  "The canonical analysis contains batches 1 through 6.",
  paths$vf, "No additional filter", "n_distinct(Batch)"
)
add_claim(
  "C04", "Study population", "The median resident contributed 3 timepoints (IQR 2-4; range 1-9).",
  "Median 3; IQR 2-4; range 1-9", "167 residents", "Resident-level summary",
  "Resident-linked isolate table", sprintf(
    "Median %s; IQR %s-%s; range %s-%s",
    fmt_num(tp_quantiles[[3]], 0), fmt_num(tp_quantiles[[2]], 0),
    fmt_num(tp_quantiles[[4]], 0), fmt_num(tp_quantiles[[1]], 0),
    fmt_num(tp_quantiles[[5]], 0)
  ), sprintf("%d WGS-ready participants", n_vf_participants),
  "Independent non-equivalent denominator", "Independently corroborated",
  "The distribution is identical in the current WGS-ready universe, although that universe contains 162 rather than 167 linked participants.",
  paths$vf, "No additional filter", "Per-participant row counts; quartiles and range"
)

baseline_claims <- list(
  c("C05", "Median age", "86.5 years (IQR 82-91)"),
  c("C06", "Female residents", "152/167 (91.0%)"),
  c("C07", "Psychogeriatric and somatic ward distribution", "82/167 (49.1%) and 85/167 (50.9%)"),
  c("C08", "Dementia prevalence", "87/167 (52.1%)"),
  c("C09", "Urinary incontinence prevalence", "109/167 (65.3%)"),
  c("C10", "Charlson Comorbidity Index", "Median 2 (IQR 1-4)"),
  c("C11", "ADL score", "Median 8 (IQR 7-10)"),
  c("C12", "Any antibiotic exposure during follow-up", "100/166 (60.2%)"),
  c("C13", "Residents with at least one recorded UTI", "89/167 (53.3%)")
)
for (x in baseline_claims) {
  add_claim(
    x[[1]], "Study population", paste0(x[[2]], "."),
    x[[3]], "Thesis cohort", "Descriptive summary",
    "Castor demographic/clinical extract used by the thesis", "Not recalculated", "Unavailable",
    "Unavailable source data", "Not checkable with available data",
    "The current analysis anchors do not contain the necessary resident-level demographic, comorbidity, ADL, antibiotic, or complete Castor UTI-history variables.",
    "Required thesis/Castor source absent", "None", "Not applicable"
  )
}
event_labeled_uti_participants <- safe_unique_n(
  vf$Participant_id[vf$Event_type == "UTI_event"]
)
primary_uti_participants <- safe_unique_n(
  vf$Participant_id[vf$UTI_Status == "UTI"]
)
add_claim(
  "C14", "Study population", "29/167 residents had at least one sequenced UTI isolate.",
  "29/167 (17.4%)", "Thesis cohort", "Resident-level sequenced UTI count",
  "Thesis UTI case definition and isolate list", sprintf(
    "%d participants have a sequenced UTI-event-labelled isolate; %d have a sequenced isolate meeting the current primary UTI definition",
    event_labeled_uti_participants, primary_uti_participants
  ), sprintf("%d WGS-ready participants", n_vf_participants),
  "Independent definition sensitivity", "Discrepant",
  "The event-label count is close to the thesis value, but the current symptom/culture-based UTI definition yields a much smaller participant count. The thesis must state whether 'sequenced UTI' means event label or clinically classified UTI.",
  paths$vf, "Participant has Event_type == UTI_event; separately UTI_Status == UTI",
  "n_distinct(Participant_id) under each UTI definition"
)

add_claim(
  "C15", "Strain carriage", "There were 368 consecutive within-resident isolate pairs.",
  "368 pairs", "167 residents", "Chronological consecutive-pair construction",
  "Thesis wgMLST pair table and ordering code", sprintf(
    "%d current consecutive WGS-linked pairs from %d participants",
    n_transitions, n_transition_participants
  ), as.character(n_transitions),
  "Independent denominator check", "Discrepant",
  "The current independently ordered WGS-linked transition table contains 394 pairs. Different linkage, exclusions, or ordering may explain the gap, but the thesis needs an attrition table.",
  paths$transitions, "cohort == \"all\"", "nrow(filtered transitions)"
)
add_claim(
  "C16", "Strain carriage", "231/368 pairs (62.8%) were persistent at <=25 wgMLST alleles.",
  "231/368 (62.8%)", "368 consecutive pairs", "wgMLST <=25 alleles",
  "Thesis wgMLST pair-distance table", sprintf(
    "Strict SNP same: %s; not confidently replaced: %s",
    fmt_pct(n_strict_same, n_transitions),
    fmt_pct(n_not_confident_replacement, n_transitions)
  ), as.character(n_transitions),
  "Directional SNP corroboration only", "Not checkable with available data",
  "wgMLST allele distance and assembly-to-assembly dnadiff SNP evidence are not interchangeable. The current data show that the estimated persistence fraction is highly definition-sensitive.",
  paths$transitions, "cohort == \"all\"",
  "Strict: SNPs <=25; broad: replacement_flag == FALSE"
)
add_claim(
  "C17", "Strain carriage", "26 pairs formed a 26-100 allele grey zone and all had the same ST and fewer than 60 SNPs.",
  "26/368 (7.1%)", "368 consecutive pairs", "wgMLST, MLST and SNP cross-check",
  "Thesis wgMLST pair table plus exact SNP output", "Not recalculated", "Unavailable",
  "Unavailable wgMLST source", "Not checkable with available data",
  "The independent pipeline does not contain the thesis wgMLST distances or a mapping of those 26 pairs. Same ST and <60 SNPs alone do not establish same-strain identity under a 25-SNP rule.",
  "Required thesis wgMLST pair mapping absent", "None", "Not applicable"
)
add_claim(
  "C18", "Strain carriage", "111/368 pairs (30.2%) were replacements above 100 wgMLST alleles.",
  "111/368 (30.2%)", "368 consecutive pairs", "wgMLST >100 alleles",
  "Thesis wgMLST pair-distance table", sprintf(
    "Current independent replacement flag: %s; above 25-SNP threshold: %s",
    fmt_pct(n_replacement, n_transitions), fmt_pct(n_above_snp, n_transitions)
  ), as.character(n_transitions),
  "Directional SNP corroboration only", "Not checkable with available data",
  "The current replacement flag and SNP threshold do not reproduce a >100-wgMLST-allele definition.",
  paths$transitions, "cohort == \"all\"",
  "replacement_flag == TRUE; separately SNPs >25"
)
add_claim(
  "C19", "Strain carriage", "The five trajectory groups contained 67, 23, 19, 16 and 11 residents.",
  "136 residents classified; stable persistence 67 (49.3%)", "Residents with >=2 pairs",
  "Resident trajectory algorithm", "Thesis pair table and trajectory code",
  "Not recalculated because thesis handling of grey-zone pairs is unresolved", "Unavailable",
  "Definition unavailable", "Not checkable with available data",
  "The trajectory counts cannot be validated until it is clear whether 26-100 allele pairs were treated as replacements, uncertain pairs, or excluded.",
  "Required thesis trajectory code absent", "None", "Not applicable"
)
add_claim(
  "C20", "Lineages", "There were 83 unique sequence types.",
  "83 STs", "556 isolates", "Achtman MLST",
  "Current canonical ST table", sprintf("%d unique typed STs", n_st), as.character(n_vf),
  "Exact count", "Verified from available data",
  "The independent canonical table contains the same number of distinct typed STs.",
  paths$vf, "!is.na(ST)", "n_distinct(ST)"
)
add_claim(
  "C21", "Lineages", "ST131 declined from 18.2% at T0 to 5.5% at T6.",
  "18.2% to 5.5%", "T0 and T6 isolates", "Timepoint prevalence",
  "Current canonical ST table", sprintf(
    "ST131: T0 %s of all / %.1f%% typed; T6 %s of all / %.1f%% typed",
    fmt_pct(st131_t0$n, st131_t0$n_all), st131_t0$pct_typed,
    fmt_pct(st131_t6$n, st131_t6$n_all), st131_t6$pct_typed
  ), sprintf("T0 n=%d; T6 n=%d", st131_t0$n_all, st131_t6$n_all),
  "Independent trend check", "Independently corroborated",
  "The direction and approximate magnitude are supported, but exact percentages differ and depend on whether missing STs are included.",
  paths$vf, "tp_lab in {T0,T6}", "ST count / all rows; ST count / typed rows"
)
add_claim(
  "C22", "Lineages", "ST73 increased from 6.8% at T0 to 16.4% at T6.",
  "6.8% to 16.4%", "T0 and T6 isolates", "Timepoint prevalence",
  "Current canonical ST table", sprintf(
    "ST73: T0 %s of all / %.1f%% typed; T6 %s of all / %.1f%% typed",
    fmt_pct(st73_t0$n, st73_t0$n_all), st73_t0$pct_typed,
    fmt_pct(st73_t6$n, st73_t6$n_all), st73_t6$pct_typed
  ), sprintf("T0 n=%d; T6 n=%d", st73_t0$n_all, st73_t6$n_all),
  "Independent trend check", "Independently corroborated",
  "The independent analysis supports an increase, but not the exact reported percentages.",
  paths$vf, "tp_lab in {T0,T6}", "ST count / all rows; ST count / typed rows"
)

drift_claims <- list(
  c("C23", "39 residents with stable persistence and T0 were analysed", "39 residents"),
  c("C24", "Median distance rose from 3.5 alleles at T1 to 15.5 at T6", "3.5 to 15.5 alleles"),
  c("C25", "Persistent strains accumulated approximately one allele per month", "Approximately 1 allele/month"),
  c("C26", "No stable-persistence resident exceeded 25 alleles", "0 residents exceeded threshold")
)
for (x in drift_claims) {
  add_claim(
    x[[1]], "Within-host drift", paste0(x[[2]], "."),
    x[[3]], "Stable-persistence residents", "T0-referenced wgMLST trajectory",
    "Thesis wgMLST profiles, dates and stable-persistence inclusion list",
    "Not recalculated", "Unavailable",
    "Unavailable wgMLST source", "Not checkable with available data",
    "Current assembly-to-assembly dnadiff SNP results cannot substitute for the thesis T0-referenced wgMLST profiles.",
    "Required thesis wgMLST longitudinal table absent", "None", "Not applicable"
  )
}
add_claim(
  "C27", "Within-host drift", "The one-allele-per-month statement is a supported evolutionary rate.",
  "Approximately 1 allele/month", "39 selected residents", "Medians across nominal timepoints",
  "Individual dates and repeated-measures rate model", "No regression estimate or uncertainty is reported", "39 thesis-selected residents",
  "Method review", "Methodologically unsupported",
  "A change between endpoint medians is not a rate estimate. A time-based repeated-measures model with uncertainty is needed, and selection on stable persistence can induce circularity.",
  "Thesis Methods 2.3 and Results 3.3", "Textual methodology audit",
  "Requires allele_distance ~ elapsed_time with resident-level dependence"
)

recurrence_claims <- list(
  c("C28", "27 residents showed re-emergence of the original strain", "27 residents"),
  c("C29", "79% of recurrence cases were genomic relapses", "79%"),
  c("C30", "The median displacement-to-re-emergence gap was five timepoints", "5 timepoints")
)
for (x in recurrence_claims) {
  add_claim(
    x[[1]], "Recurrence", paste0(x[[2]], "."),
    x[[3]], "Residents with a claimed recurrence", "wgMLST recurrence classification",
    "Thesis recurrence case list and wgMLST distances", "Not recalculated", "Unavailable",
    "Unavailable wgMLST source", "Not checkable with available data",
    "The required original/displacing/returning isolate mapping is not available in the current outputs.",
    "Required thesis recurrence table absent", "None", "Not applicable"
  )
}
add_claim(
  "C31", "Recurrence", "The recurrence description is internally coherent.",
  "Remaining recurrence cases were reinfections because the original strain did not return",
  "27 residents selected for re-emergence", "Narrative definition",
  "No additional data required", "Selection says the original strain re-emerged; later wording says it did not return", "Thesis text",
  "Internal consistency check", "Methodologically unsupported",
  "The reported selection criterion and explanation of the remaining cases contradict each other; the recurrence categories require correction and explicit case definitions.",
  "Thesis Results 3.3", "Textual methodology audit", "Internal logical consistency"
)

add_claim(
  "C32", "ASB to UTI", "There were 34 ASB-to-UTI pairs across 26 residents.",
  "34 pairs; 26 residents", "Routine ASB followed by UTI", "Chronological wgMLST pairs",
  "Thesis transition list and clinical definitions", sprintf(
    "Event-label analysis: %d routine-to-UTI-event pairs in %d participants; current primary definition: %d pairs",
    event_sum$n, event_sum$participants, primary_event_sum$n
  ), sprintf("%d current WGS-linked transitions", n_transitions),
  "Independent denominator check", "Discrepant",
  "Neither the thesis-event-label analysis nor the current clinical definition reproduces 34 WGS-linked pairs. A row-level attrition reconciliation is required.",
  paths$transitions, "cohort == \"all\"; Routine -> UTI_event; primary UTI sensitivity",
  "Count chronological transition rows under each clinical definition"
)
add_claim(
  "C33", "ASB to UTI", "16/34 UTIs used the same strain and 18/34 (52.9%) used a new strain.",
  "16 same; 18 new", "34 ASB-to-UTI pairs", "wgMLST <=25 vs >25 alleles",
  "Thesis wgMLST transition table", sprintf(
    paste0(
      "Event-label analysis: strict SNP same %s, above SNP threshold %s, confident replacement %s; ",
      "current primary UTI analysis: strict SNP same %s, above threshold %s, confident replacement %s"
    ),
    fmt_pct(event_sum$strict_same, event_sum$n),
    fmt_pct(event_sum$above_snp, event_sum$n),
    fmt_pct(event_sum$replacement, event_sum$n),
    fmt_pct(primary_event_sum$strict_same, primary_event_sum$n),
    fmt_pct(primary_event_sum$above_snp, primary_event_sum$n),
    fmt_pct(primary_event_sum$replacement, primary_event_sum$n)
  ), sprintf("Event labels n=%d; primary UTI n=%d", event_sum$n, primary_event_sum$n),
  "Independent sensitivity analysis", "Discrepant",
  "The conclusion is sensitive to both UTI definition and strain definition. Under the current primary UTI definition, most linked routine-to-UTI transitions have <=25 SNPs; therefore the headline 'majority new strain' claim is not independently corroborated.",
  paths$transitions, "cohort == \"all\"; two clinical definitions",
  "Separate strict SNP context from replacement_flag"
)
add_claim(
  "C34", "ASB to UTI", "The separate transition-rate analysis used 336 ASB-ASB, 26 ASB-UTI and 6 UTI-UTI pairs.",
  "336 + 26 + 6 = 368 pairs", "All 368 consecutive pairs", "Transition-type tabulation",
  "Thesis pair-level transition table", sprintf(
    "Current primary transitions: Not_UTI->Not_UTI %d; Not_UTI->UTI %d; UTI->Not_UTI %d; UTI->UTI %d",
    sum(transitions$transition_type == "Not_UTI\u2192Not_UTI"),
    sum(transitions$transition_type == "Not_UTI\u2192UTI"),
    sum(transitions$transition_type == "UTI\u2192Not_UTI"),
    sum(transitions$transition_type == "UTI\u2192UTI")
  ), as.character(n_transitions),
  "Internal and independent denominator check", "Discrepant",
  "The thesis separately reports 34 ASB-to-UTI pairs but uses 26 in this complete 368-pair partition. It also omits an explicit UTI-to-ASB category despite analysing post-UTI samples.",
  paths$transitions, "cohort == \"all\"", "Table of status_from by status_to"
)
add_claim(
  "C35", "Post-UTI return", "The original ASB strain returned in 25/37 cases (67.6%).",
  "25/37 (67.6%)", "37 post-UTI pairs", "Pre-ASB strain compared with post-UTI isolate",
  "Thesis three-isolate case table and wgMLST distances", sprintf(
    paste0(
      "Current complete pre-routine/UTI-event/post-routine triples: n=%d; ",
      "strict SNP return %s; not confidently replaced %s. ",
      "Primary-UTI triples: n=%d; strict return %s; not confidently replaced %s"
    ),
    triple_all_sum$n,
    fmt_pct(triple_all_sum$original_strict_return, triple_all_sum$n),
    fmt_pct(triple_all_sum$original_not_confident_replaced, triple_all_sum$n),
    triple_primary_sum$n,
    fmt_pct(triple_primary_sum$original_strict_return, triple_primary_sum$n),
    fmt_pct(triple_primary_sum$original_not_confident_replaced, triple_primary_sum$n)
  ), sprintf("All event-label triples n=%d; primary UTI triples n=%d", triple_all_sum$n, triple_primary_sum$n),
  "Independent three-isolate sensitivity analysis", "Discrepant",
  "The current complete triple analysis does not reproduce 37 evaluable cases or a 67.6% return estimate. The estimate remains method-sensitive because SNP and wgMLST distances differ.",
  paste(paths$transitions, paths$pairwise, sep = "; "),
  "cohort == \"all\"; matched Routine -> UTI_event -> Routine triples",
  "Compare pre-routine isolate directly with post-UTI routine isolate"
)
add_claim(
  "C36", "Post-UTI return", "Return was 94% after same-strain UTI and 48% after new-strain UTI.",
  "94% and 48%", "37 post-UTI cases", "Stratified wgMLST return analysis",
  "Thesis three-isolate case table", sprintf(
    "Current event-label triples: pre-to-UTI strict same %s; pre-to-post strict return %s",
    fmt_pct(triple_all_sum$uti_strict_same, triple_all_sum$n),
    fmt_pct(triple_all_sum$original_strict_return, triple_all_sum$n)
  ), sprintf("%d complete current triples", triple_all_sum$n),
  "Independent three-isolate sensitivity analysis", "Not checkable with available data",
  "The current sample is smaller and uses SNP evidence. Exact stratum-specific wgMLST return estimates require the thesis case table.",
  paste(paths$transitions, paths$pairwise, sep = "; "),
  "Matched triples; no substitution of SNP for wgMLST",
  "Stratify direct pre-to-post comparison by pre-to-UTI strain relation"
)

antibiotic_claims <- list(
  c("C37", "Antibiotic linkage coverage", "337/368 pairs (90.8%)"),
  c("C38", "Replacement rates with versus without recent antibiotics", "56.0% versus 31.3%"),
  c("C39", "Antibiotics and replacement odds", "OR 2.78; 95% CI 1.49-5.20; p=0.001"),
  c("C40", "Replacement rates by six individual antibiotics", "43.7%, 40.7%, 39.1%, 37.3%, 32.0%, 26.2%"),
  c("C41", "Antibiotics and subsequent sequenced UTI", "24/100 versus 5/66; OR 3.83; 95% CI 1.33-13.60; p=0.006"),
  c("C42", "Host factors and replacement", "No statistically significant associations")
)
for (x in antibiotic_claims) {
  add_claim(
    x[[1]], "Antibiotics and host factors", paste0(x[[2]], "."),
    x[[3]], "Thesis antibiotic/host analysis universe", "Castor-linked analysis",
    "Castor exposure, host-factor and outcome extract plus thesis code",
    "Not recalculated", "Unavailable",
    "Unavailable source data", "Not checkable with available data",
    "The current pipeline explicitly records that antibiotic exposure data are missing from its inputs; proxies were not used.",
    "Required thesis/Castor source absent", "None", "Not applicable"
  )
}
add_claim(
  "C43", "Cox model", "A Cox model included 368 pairs and 137 replacement events and estimated HR 2.13.",
  "HR 2.13; 95% CI 1.40-3.26; p=0.0005; n=368; 137 events",
  "368 consecutive pairs", "Cox proportional hazards model",
  "Model-ready data, formula, time variables, missing-data handling and diagnostics",
  "Not recalculated; thesis states antibiotic data existed for only 337 pairs", "337 exposed-status pairs reported",
  "Internal methods audit", "Methodologically unsupported",
  "The reported n=368 is unreconciled with 31 pairs lacking antibiotic exposure. The thesis does not specify start/stop times, censoring, recurrent-event handling, resident clustering, proportional-hazards checks, or how 137 events relate to 111 >100-allele replacements.",
  "Thesis Methods 2.4-2.5 and Results 3.5", "Textual methodology audit",
  "Requires explicit survival data structure and model formula"
)
add_claim(
  "C44", "Statistical dependence", "Pair-level Fisher, Wilcoxon and Cox analyses adequately handle repeated pairs from the same resident.",
  "No clustering method reported", "Repeated pair-level observations", "Pair-level inferential tests",
  "Analysis code and resident-cluster specification", "Current independent analysis confirms repeated transitions per participant", sprintf("%d transitions from %d participants", n_transitions, n_transition_participants),
  "Method review", "Methodologically unsupported",
  "Pairs from the same resident are not independent. The thesis does not report clustered standard errors, mixed models, GEE, frailty, or another repeated-measures method.",
  paths$transitions, "cohort == \"all\"", "Compare number of rows with number of unique participants"
)
add_claim(
  "C45", "Multiple testing", "The p<0.05 threshold is adequate across all antibiotic and host-factor comparisons.",
  "Unadjusted p<0.05", "Multiple exposures and host factors", "Several hypothesis tests",
  "Complete test family and prespecified analysis plan", "No multiple-testing correction is described", "Thesis text",
  "Method review", "Methodologically unsupported",
  "Testing several antibiotics and host characteristics without a defined primary contrast or multiplicity strategy inflates false-positive risk.",
  "Thesis Methods 2.5 and Results 3.5", "Textual methodology audit",
  "Requires prespecified test families and adjusted or clearly exploratory inference"
)
add_claim(
  "C46", "Assembly methodology", "Flye was selected because it had higher Nanopore assembly accuracy than Longcycler.",
  "Flye preferred for wgMLST", "556 isolates", "Cited assembler benchmark",
  "Exact basecalling, coverage, polishing, assembly QC and scheme provenance", sprintf(
    "Current canonical pipeline selected %d Longcycler and %d Flye assemblies",
    n_longcycler, n_flye
  ), as.character(nrow(assemblies)),
  "Citation and independent-pipeline audit", "Methodologically unsupported",
  "Landman et al. report that Longcycler and Miniasm were best for wgMLST allele calling, whereas Flye performed better for AMR/replicon detection. For E. coli, Flye's 95th-percentile wgMLST discrepancy was 26.1 alleles, exceeding the thesis's 25-allele threshold.",
  paste(paths$assemblies, "Landman et al. 2024", sep = "; "),
  "selected_canonical == TRUE", "Count selected assembler; compare with cited benchmark"
)
add_claim(
  "C47", "wgMLST scheme", "The cited E. coli wgMLST scheme contains 4,512 loci.",
  "4,512 loci", "E. coli wgMLST scheme", "SeqSphere v9.0.1",
  "Scheme name, version/export date and target list", "The cited Landman paper reports 4,503 loci", "Published scheme description",
  "Citation audit", "Discrepant",
  "This may reflect a later scheme revision, but the thesis must identify the exact scheme version rather than relying on a mismatched citation.",
  "Landman et al. 2024; thesis Methods 2.2", "Textual citation audit", "Compare reported scheme size and version"
)
add_claim(
  "C48", "Causal interpretation", "Antibiotic exposure was the primary driver of strain replacement and subsequent UTI.",
  "Primary driver; independent effect", "Observational cohort", "Associational models",
  "Time-resolved exposure/outcome data and a causal analysis plan", "Not established by the available independent data", "Unavailable antibiotic data",
  "Causal-inference review", "Methodologically unsupported",
  "Confounding by indication, reverse causality, exposure timing, repeated measures and treatment of UTI itself are not adequately resolved. Association cannot establish a primary causal driver.",
  "Thesis Results/Discussion/Conclusion", "Textual causal-claim audit",
  "Requires temporally ordered exposure, confounder control and sensitivity analyses"
)
add_claim(
  "C49", "Causal interpretation", "The findings demonstrate protective bladder colonisation.",
  "Supports/demonstrates protective colonisation", "Observational genomic cohort", "Interpretation of persistence and UTI transitions",
  "Appropriate causal comparator or intervention", "Some persistence and lineage trends are independently corroborated, but protection is not tested directly", "Current genomic outputs",
  "Causal-inference review", "Methodologically unsupported",
  "Persistence and strain return are compatible with the hypothesis but do not demonstrate that colonisation prevents UTI.",
  "Thesis Discussion/Conclusion", "Textual causal-claim audit",
  "Distinguish compatibility with a hypothesis from causal demonstration"
)
add_claim(
  "C50", "Overall conclusion", "Stable carriage is common and lineage composition changes over time.",
  "Persistence common; ST131 down and ST73 up", "Longitudinal genomic cohort", "Descriptive genomic analysis",
  "Current transition and ST tables", sprintf(
    "Not confidently replaced %s; ST131 decreased and ST73 increased in the current data",
    fmt_pct(n_not_confident_replacement, n_transitions)
  ), sprintf("%d transitions; %d genomes", n_transitions, n_vf),
  "Independent descriptive corroboration", "Independently corroborated",
  "The broad descriptive direction is supported. Exact persistence percentages remain method-dependent and should not be called replicated.",
  paste(paths$transitions, paths$vf, sep = "; "),
  "cohort == \"all\"; canonical ST rows", "Descriptive counts only"
)

claim_matrix <- do.call(rbind, claims)

if (any(is.na(claim_matrix$audit_status)) || any(!claim_matrix$audit_status %in% allowed_statuses)) {
  stop("Every claim must have exactly one valid audit status.")
}
if (anyDuplicated(claim_matrix$claim_id)) stop("Duplicate claim IDs detected.")
trace_fields <- c("source_file", "filter", "formula")
for (field in trace_fields) {
  if (any(is.na(claim_matrix[[field]]) | !nzchar(claim_matrix[[field]]))) {
    stop("Every claim must have a non-empty `", field, "` value.")
  }
}

claim_path <- file.path(out_dir, "claim_matrix.csv")
write.csv(claim_matrix, claim_path, row.names = FALSE, na = "")

# ---- Markdown report --------------------------------------------------------

status_counts <- table(factor(claim_matrix$audit_status, levels = allowed_statuses))
status_line <- function(status) sprintf("%s: **%d**", status, status_counts[[status]])

md <- c(
  "# Independent Scientific Audit of Hamdi Hersi's Thesis",
  "",
  sprintf("**Generated:** %s", format(Sys.Date(), "%Y-%m-%d")),
  "",
  "## Executive verdict",
  "",
  paste0(
    "The available evidence supports several basic cohort and lineage facts, but it does **not** permit a full validation of the thesis's wgMLST, antibiotic, demographic, or host-factor analyses. ",
    "The strongest defensible conclusion is therefore **partly corroborated but not reproducible from the materials currently available**. ",
    "Several denominator and method inconsistencies require resolution before the headline causal conclusions can be treated as sound."
  ),
  "",
  "This is a scientific-validity and reproducibility audit. It does not assess intent and does not make an allegation of fabrication or misconduct.",
  "",
  "Claim classifications:",
  "",
  paste0("- ", vapply(allowed_statuses, status_line, character(1))),
  "",
  "## Scope and guardrails",
  "",
  "- Current June 2026 generated outputs were treated as the independent analysis anchors.",
  "- SNP and wgMLST allele distances were kept separate; SNP evidence is directional corroboration, not an exact replication.",
  "- No demographic, antibiotic, host-factor or wgMLST result was reconstructed from a proxy.",
  "- Transition analyses use only `cohort == \"all\"`; stratified duplicate cohorts were excluded.",
  "- The thesis appears to be a tracked draft. Editorial annotations were not treated as evidence of analytical invalidity.",
  "",
  "## Independently computed anchors",
  "",
  "| Measure | Current result | Source/filter |",
  "|---|---:|---|",
  sprintf("| Primary clinical episodes | %d from %d participants | `status_map.csv`; `analysis_include_primary == TRUE` |", n_clinical, n_clinical_participants),
  sprintf("| Canonical WGS/VF episodes | %d from %d participants | `vf_analysis_ready.csv` |", n_vf, n_vf_participants),
  sprintf("| Sequencing batches | %d | Canonical WGS/VF rows |", n_batches),
  sprintf("| Canonical selected assemblers | %d Longcycler; %d Flye | `selected_canonical == TRUE` |", n_longcycler, n_flye),
  sprintf("| Unique typed STs | %d | Non-missing `ST` |", n_st),
  sprintf("| Consecutive WGS-linked transitions | %d from %d participants | `cohort == \"all\"` |", n_transitions, n_transition_participants),
  sprintf("| Strict dnadiff SNP support (operational) | %s | numeric `SNPs <= %s`; label checked for agreement |", fmt_pct(n_strict_same, n_transitions), fmt_num(strict_snp_threshold, 0)),
  sprintf("| Confident replacement flag | %s | `replacement_flag == TRUE` |", fmt_pct(n_replacement, n_transitions)),
  sprintf("| Not confidently replaced | %s | `replacement_flag == FALSE` |", fmt_pct(n_not_confident_replacement, n_transitions)),
  "",
  "The strict SNP and broader replacement results answer different questions. They must not be collapsed into one persistence estimate.",
  "",
  "## Findings that are supported",
  "",
  sprintf("- The isolate count (**%d**), batch count (**%d**) and number of unique STs (**%d**) match the thesis.", n_vf, n_batches, n_st),
  sprintf(
    "- ST131 decreases from T0 (%s of all isolates; %.1f%% of typed isolates) to T6 (%s; %.1f%% typed).",
    fmt_pct(st131_t0$n, st131_t0$n_all), st131_t0$pct_typed,
    fmt_pct(st131_t6$n, st131_t6$n_all), st131_t6$pct_typed
  ),
  sprintf(
    "- ST73 increases from T0 (%s of all isolates; %.1f%% of typed isolates) to T6 (%s; %.1f%% typed).",
    fmt_pct(st73_t0$n, st73_t0$n_all), st73_t0$pct_typed,
    fmt_pct(st73_t6$n, st73_t6$n_all), st73_t6$pct_typed
  ),
  sprintf(
    "- The WGS-ready timepoint distribution is median %s, IQR %s-%s, range %s-%s, matching the reported descriptive distribution.",
    fmt_num(tp_quantiles[[3]], 0), fmt_num(tp_quantiles[[2]], 0),
    fmt_num(tp_quantiles[[4]], 0), fmt_num(tp_quantiles[[1]], 0),
    fmt_num(tp_quantiles[[5]], 0)
  ),
  "- These agreements corroborate dataset identity and broad lineage trends; they do not validate the thesis's wgMLST or antibiotic models.",
  "",
  "## Important denominator and result discrepancies",
  "",
  sprintf(
    "- **Participants:** the thesis states 167 residents contributed 556 isolates. The current primary clinical universe has %d participants, while the 556 canonical WGS rows link to %d participants.",
    n_clinical_participants, n_vf_participants
  ),
  sprintf(
    "- **Consecutive pairs:** the thesis reports 368; the current WGS-linked table contains %d after the required `cohort == \"all\"` filter.",
    n_transitions
  ),
  "- **Grey-zone handling:** 231 persistent + 26 grey-zone + 111 >100-allele replacements equals 368, but later analyses report 137 replacement events. The latter equals 111 + 26, suggesting that grey-zone pairs may have been counted as replacement in the models while being presented separately in the descriptive results.",
  "- **Antibiotic model denominator:** exposure is reported for 337 pairs, but the Cox model is reported with n=368. Missing-exposure handling is not described.",
  "- **ASB-to-UTI denominator:** 34 pairs across 26 residents are reported in one analysis, while a separate complete transition partition uses only 26 ASB-to-UTI pairs.",
  "- **Recurrence wording:** cases are introduced as residents whose original strain re-emerged, but the remaining cases are then described as cases where the original strain did not return.",
  "",
  "## Independent UTI-transition sensitivity analyses",
  "",
  "### Thesis-style event labels",
  "",
  sprintf(
    "There are **%d** WGS-linked Routine→UTI-event transitions from **%d** participants. Strict SNP evidence classifies %s as same strain and %s as above the 25-SNP threshold; only %s carry a confident replacement flag.",
    event_sum$n, event_sum$participants,
    fmt_pct(event_sum$strict_same, event_sum$n),
    fmt_pct(event_sum$above_snp, event_sum$n),
    fmt_pct(event_sum$replacement, event_sum$n)
  ),
  "",
  "### Current primary clinical definition",
  "",
  sprintf(
    "Only **%d** Routine→UTI-event transitions also satisfy the current Not_UTI→UTI definition. Of these, %s have strict same-strain SNP evidence, %s exceed 25 SNPs and %s carry a confident replacement flag.",
    primary_event_sum$n,
    fmt_pct(primary_event_sum$strict_same, primary_event_sum$n),
    fmt_pct(primary_event_sum$above_snp, primary_event_sum$n),
    fmt_pct(primary_event_sum$replacement, primary_event_sum$n)
  ),
  "",
  "Accordingly, the thesis conclusion that most UTIs were caused by a new strain is **not robustly corroborated** by the current primary clinical definition. The difference could arise from clinical definitions, wgMLST/SNP differences, assembly selection or pair inclusion; the thesis row-level transition table is required to distinguish these explanations.",
  "",
  "### Direct pre-UTI to post-UTI return",
  "",
  sprintf(
    "The current data contain **%d** complete Routine→UTI-event→Routine triples. Direct comparison of the pre-UTI and post-UTI routine isolates finds strict SNP-defined return in %s and no confident replacement in %s.",
    triple_all_sum$n,
    fmt_pct(triple_all_sum$original_strict_return, triple_all_sum$n),
    fmt_pct(triple_all_sum$original_not_confident_replaced, triple_all_sum$n)
  ),
  sprintf(
    "Restricting to the **%d** triples whose UTI event satisfies the current primary definition gives strict return in %s and no confident replacement in %s.",
    triple_primary_sum$n,
    fmt_pct(triple_primary_sum$original_strict_return, triple_primary_sum$n),
    fmt_pct(triple_primary_sum$original_not_confident_replaced, triple_primary_sum$n)
  ),
  "",
  "These are independent sensitivity analyses, not exact replications of the thesis's 37 wgMLST-defined post-UTI pairs.",
  "",
  "## Methodological assessment",
  "",
  "### 1. Assembly choice and wgMLST threshold — major concern",
  "",
  paste0(
    "The Methods state that Flye was chosen because it was more accurate than Longcycler. ",
    "The cited Landman et al. study instead reports that Longcycler and Miniasm were best for wgMLST allele calling, while Flye performed better for AMR genes and replicons. ",
    "For E. coli, Flye's 95th-percentile discrepancy relative to Illumina was **26.1 wgMLST alleles**, which exceeds the thesis's **25-allele** strain threshold. ",
    "This is especially relevant to the 26 grey-zone pairs. The current independent pipeline selected ",
    n_longcycler, " Longcycler and ", n_flye, " Flye assemblies."
  ),
  "",
  "The cited paper describes a 4,503-locus E. coli scheme, whereas the thesis reports 4,512 loci. A later scheme revision is possible, but the exact SeqSphere scheme name, version, target export and date must be supplied.",
  "",
  "Primary reference: [Landman et al., 2024](https://pure.eur.nl/ws/portalfiles/portal/187878763/Landman_et_al_2024_Final.pdf).",
  "",
  "### 2. Repeated observations — major concern",
  "",
  paste0(
    "The current independent dataset has ", n_transitions, " transitions from ",
    n_transition_participants, " participants, confirming substantial within-resident dependence. ",
    "The thesis describes pair-level Fisher, Wilcoxon and Cox analyses but does not report mixed models, GEE, clustered standard errors, frailty or another repeated-measures correction."
  ),
  "",
  "### 3. Cox model specification — major concern",
  "",
  "The thesis does not define the survival time scale, start/stop structure, censoring, recurrent-event framework, participant clustering, proportional-hazards diagnostics or missing-exposure handling. Reporting 368 model rows despite antibiotic data for 337 pairs requires explicit reconciliation.",
  "",
  "### 4. Antibiotic causal claims — major concern",
  "",
  "The observational analyses do not adequately resolve confounding by indication, reverse causality or whether an antibiotic was prescribed for the UTI outcome itself. The results may describe an association, but phrases such as “primary driver,” “increased UTI risk,” and “strengthened causal inference” are not justified by the described methods.",
  "",
  "### 5. Genetic drift rate — moderate concern",
  "",
  "Calling the change in timepoint medians “one allele per month” is not a formal rate estimate. The analysis conditions on stable persistence and needs an elapsed-time repeated-measures model with uncertainty. The statement that no selected stable resident crossed the identity threshold is also partly dependent on how stable persistence was defined.",
  "",
  "### 6. Multiple testing and covariate selection — moderate concern",
  "",
  "Several antibiotics and host characteristics were tested without a stated multiplicity strategy. Excluding host covariates because their univariable tests were non-significant does not establish that the antibiotic effect is independent of those factors.",
  "",
  "## Evidence required before stronger validation",
  "",
  "1. The isolate-level and pair-level wgMLST table, including comparable loci, missing loci and exact scheme version.",
  "2. The complete analysis code and package/session information.",
  "3. The assembly manifest, basecaller, polishing status, coverage and QC metrics for every isolate.",
  "4. The row-level ASB/UTI transition and recurrence case tables.",
  "5. The Castor antibiotic extract with prescription indication, start/end dates, drug and linkage rules.",
  "6. Exact statistical model formulas, model-ready datasets, resident clustering method and diagnostics.",
  "7. The supplementary tables referenced by the analysis; the supplied draft's Supplementary Materials section is empty.",
  "",
  "## Bottom line",
  "",
  "The thesis is **not fully verifiable from the available materials**. Basic dataset size, batch count, ST diversity and the directions of the ST131/ST73 trends are supported. The exact persistence, recurrence, ASB-to-UTI and antibiotic estimates remain unverified, and several methodological inconsistencies materially weaken the causal conclusions. The appropriate next step is targeted author clarification and release of the listed reproducibility materials—not an inference about intent.",
  "",
  "## Reproducibility",
  "",
  sprintf("- Claim matrix: `%s`", claim_path),
  sprintf("- Clinical anchor: `%s`", paths$clinical),
  sprintf("- WGS/VF anchor: `%s`", paths$vf),
  sprintf("- Transition anchor: `%s` filtered to `cohort == \"all\"`", paths$transitions),
  sprintf("- Pairwise anchor: `%s`", paths$pairwise),
  sprintf("- Canonical assembly anchor: `%s` filtered to `selected_canonical == TRUE`", paths$assemblies),
  "",
  "Run:",
  "",
  "```bash",
  "Rscript scripts/audit_hamdi_thesis.R",
  "```"
)

report_path <- file.path(out_dir, "hamdi_thesis_scientific_audit.md")
writeLines(md, report_path, useBytes = TRUE)

cat("Audit complete\n")
cat("Claim matrix:", claim_path, "\n")
cat("Report:", report_path, "\n")
