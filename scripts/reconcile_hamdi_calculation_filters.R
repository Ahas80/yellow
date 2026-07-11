#!/usr/bin/env Rscript

# Reconcile how Hamdi Hersi's thesis appears to calculate key results with the
# current rUTI/YELLOW RoUTIne analysis filters. This is a descriptive audit:
# it does not recompute wgMLST, antibiotic, demographic or host-factor results
# from proxies.

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
dir.create(out_dir, recursive = TRUE, showWarnings = FALSE)

paths <- list(
  clinical = file.path(results_dir, "clinical", "status_map.csv"),
  vf = file.path(results_dir, "vf", "vf_analysis_ready.csv"),
  transitions = file.path(results_dir, "vf", "vf_transitions_stratified.csv"),
  pairwise = file.path(results_dir, "strain_compare", "pairwise_metrics.csv"),
  assemblies = file.path(results_dir, "qc", "canonical_assembly_selection.csv"),
  cohort_flow = file.path(results_dir, "summary", "table_01_cohort_episode_flow.csv")
)

missing <- names(paths)[!file.exists(unlist(paths))]
if (length(missing) > 0) {
  stop("Missing required inputs: ", paste(sprintf("%s=%s", missing, unlist(paths[missing])), collapse = "; "))
}

read_anchor <- function(path) {
  read.csv(path, check.names = FALSE, na.strings = c("", "NA"), stringsAsFactors = FALSE)
}

clinical_all <- read_anchor(paths$clinical)
clinical <- clinical_all[clinical_all$analysis_include_primary %in% TRUE, , drop = FALSE]
vf <- read_anchor(paths$vf)
transitions_all <- read_anchor(paths$transitions)
transitions <- transitions_all[transitions_all$cohort == "all", , drop = FALSE]
pairwise <- read_anchor(paths$pairwise)
assemblies_all <- read_anchor(paths$assemblies)
assemblies <- assemblies_all[assemblies_all$selected_canonical %in% TRUE, , drop = FALSE]
cohort_flow <- read_anchor(paths$cohort_flow)

safe_unique_n <- function(x) length(unique(x[!is.na(x)]))
pct <- function(n, d, digits = 1) if (d == 0 || is.na(d)) NA_real_ else round(100 * n / d, digits)
fmt_pct <- function(n, d) sprintf("%d/%d (%.1f%%)", n, d, pct(n, d))
fmt_range <- function(x) {
  x <- as.numeric(x)
  sprintf("%s-%s", min(x, na.rm = TRUE), max(x, na.rm = TRUE))
}
collapse_counts <- function(tab) {
  if (length(tab) == 0) return("none")
  paste(paste0(names(tab), "=", as.integer(tab)), collapse = "; ")
}
fmt_vf_event_status <- function(tab) {
  dn <- dimnames(tab)
  out <- character()
  for (status in dn[[1]]) {
    for (event in dn[[2]]) {
      out <- c(out, sprintf("%s/%s=%d", status, event, as.integer(tab[status, event])))
    }
  }
  paste(out, collapse = "; ")
}
fmt_transition_table <- function(tab) {
  dn <- dimnames(tab)
  out <- character()
  for (from in dn[[1]]) {
    for (to in dn[[2]]) {
      out <- c(out, sprintf("%s->%s=%d", from, to, as.integer(tab[from, to])))
    }
  }
  paste(out, collapse = "; ")
}

tp_counts <- as.integer(table(vf$Participant_id))
tp_iqr <- as.numeric(quantile(tp_counts, probs = c(0.25, 0.5, 0.75), names = FALSE))
tp_distribution <- sprintf(
  "median %s; IQR %s-%s; range %s",
  tp_iqr[[2]], tp_iqr[[1]], tp_iqr[[3]], fmt_range(tp_counts)
)

assembler_counts <- table(tolower(assemblies$assembler), useNA = "ifany")
asm_text <- collapse_counts(assembler_counts)

st_snapshot <- function(tp, st) {
  x <- vf[vf$tp_lab == tp, , drop = FALSE]
  typed <- !is.na(x$ST)
  n_st <- sum(as.character(x$ST) == as.character(st), na.rm = TRUE)
  sprintf(
    "%s ST%s: %d/%d all (%.1f%%); %d/%d typed (%.1f%%)",
    tp, st, n_st, nrow(x), pct(n_st, nrow(x)),
    n_st, sum(typed), pct(n_st, sum(typed))
  )
}

status_tab <- table(clinical$UTI_Status, useNA = "ifany")
vf_status_tab <- table(vf$UTI_Status, useNA = "ifany")
vf_event_status <- with(vf, table(UTI_Status, Event_type, useNA = "ifany"))
transition_status_tab <- with(transitions, table(status_from, status_to, useNA = "ifany"))

n_strict_same <- sum(transitions$snp_strain_context == "Strong same strain", na.rm = TRUE)
n_above_snp <- sum(transitions$snp_strain_context == "Above same-strain SNP threshold", na.rm = TRUE)
n_replacement <- sum(transitions$replacement_flag %in% TRUE, na.rm = TRUE)
n_same_lineage_not_same_snp <- sum(
  transitions$pair_interpretation == "Same lineage, not same strain by SNP",
  na.rm = TRUE
)

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
primary_notuti_to_uti <- transitions[
  transitions$status_from == "Not_UTI" &
    transitions$status_to == "UTI",
  ,
  drop = FALSE
]

transition_summary <- function(x) {
  sprintf(
    "%d pairs from %d participants; strict SNP same %s; above 25 SNPs %s; replacement likely %s",
    nrow(x), safe_unique_n(x$Participant_id),
    fmt_pct(sum(x$snp_strain_context == "Strong same strain", na.rm = TRUE), nrow(x)),
    fmt_pct(sum(x$snp_strain_context == "Above same-strain SNP threshold", na.rm = TRUE), nrow(x)),
    fmt_pct(sum(x$replacement_flag %in% TRUE, na.rm = TRUE), nrow(x))
  )
}

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
    return(data.frame(TotalSNPs = NA_real_, snp_strain_context = NA_character_, replacement_flag = NA))
  }
  z[1, c("TotalSNPs", "snp_strain_context", "replacement_flag"), drop = FALSE]
}

post_uti <- transitions[
  transitions$event_type_from == "UTI_event" &
    transitions$event_type_to == "Routine",
  ,
  drop = FALSE
]
triples <- merge(
  routine_to_event,
  post_uti,
  by.x = c("Participant_id", "tp_to"),
  by.y = c("Participant_id", "tp_from"),
  suffixes = c("_preuti", "_postuti")
)
if (nrow(triples) > 0) {
  pre_post <- do.call(rbind, lapply(seq_len(nrow(triples)), function(i) {
    find_pair(triples$Participant_id[[i]], triples$tp_from[[i]], triples$tp_to_postuti[[i]])
  }))
  triples$pre_post_snp_context <- pre_post$snp_strain_context
  triples$pre_post_replacement_flag <- pre_post$replacement_flag
} else {
  triples$pre_post_snp_context <- character()
  triples$pre_post_replacement_flag <- logical()
}
primary_triples <- triples[
  triples$status_from_preuti == "Not_UTI" &
    triples$status_to_preuti == "UTI",
  ,
  drop = FALSE
]
triple_summary <- function(x) {
  sprintf(
    "%d triples; pre/post strict SNP return %s; no confident replacement %s",
    nrow(x),
    fmt_pct(sum(x$pre_post_snp_context == "Strong same strain", na.rm = TRUE), nrow(x)),
    fmt_pct(sum(x$pre_post_replacement_flag %in% FALSE, na.rm = TRUE), nrow(x))
  )
}

flow_lookup <- function(stage) {
  row <- cohort_flow[cohort_flow$stage_name == stage, , drop = FALSE]
  if (nrow(row) == 0) return("not available")
  sprintf("%s episodes; %s participants", row$n_episodes[[1]], row$n_participants[[1]])
}

rows <- list()
add_row <- function(
  topic,
  thesis_location,
  thesis_calculation,
  thesis_filter,
  thesis_denominator,
  current_calculation,
  current_filter,
  current_denominator,
  compatibility,
  rationale,
  checkable
) {
  rows[[length(rows) + 1]] <<- data.frame(
    topic = topic,
    thesis_location = thesis_location,
    thesis_calculation = thesis_calculation,
    thesis_filter = thesis_filter,
    thesis_denominator = thesis_denominator,
    current_calculation = current_calculation,
    current_filter = current_filter,
    current_denominator = current_denominator,
    compatibility = compatibility,
    rationale = rationale,
    checkable = checkable,
    stringsAsFactors = FALSE
  )
}

add_row(
  "Cohort and isolates",
  "Methods 2.1; Results 3.1",
  "Counts residents with at least one MALDI-TOF E. coli isolate at >=1e3 CFU/mL submitted for WGS.",
  "Genomic resident/isolate universe; routine T0-T6 treated as ASB in prose.",
  "556 isolates from 167 residents; six batches.",
  sprintf(
    "%d VF/WGS rows from %d linked participants; clinical primary universe is %d rows from %d participants; full status map has %d rows from %d participants.",
    nrow(vf), safe_unique_n(vf$Participant_id),
    nrow(clinical), safe_unique_n(clinical$Participant_id),
    nrow(clinical_all), safe_unique_n(clinical_all$Participant_id)
  ),
  "vf_analysis_ready.csv already filtered to primary/genomics-expected rows; status_map.csv primary uses analysis_include_primary == TRUE.",
  sprintf("VF: %d rows/%d participants; primary clinical: %d rows/%d participants.", nrow(vf), safe_unique_n(vf$Participant_id), nrow(clinical), safe_unique_n(clinical$Participant_id)),
  "Partly compatible, denominator labels differ",
  "The 556 isolate count matches, but the participant denominator changes depending on whether unlinked/unknown and primary clinical exclusions are included.",
  "Partly checkable"
)

add_row(
  "Participant linkage",
  "Methods 2.2",
  "States genomic and Castor IDs were harmonised; 166/167 linked, one UNKNOWN.",
  "All genomic residents before final clinical/VF-ready linkage.",
  "167 genomic residents; 166 linked; one UNKNOWN.",
  sprintf("Current primary clinical table has %d participants; VF-ready table links the 556 rows to %d participants.", safe_unique_n(clinical$Participant_id), safe_unique_n(vf$Participant_id)),
  "Current outputs do not retain her UNKNOWN genomic-resident row in the VF-ready participant count.",
  sprintf("Clinical primary participants=%d; VF participants=%d.", safe_unique_n(clinical$Participant_id), safe_unique_n(vf$Participant_id)),
  "Needs attrition reconciliation",
  "Her own linkage statement is compatible with 166 clinically linked residents, but not enough to explain why the current WGS-ready count is 162.",
  "Partly checkable"
)

add_row(
  "Primary UTI/ASB definition",
  "Methods 2.1 and 3.4",
  "Routine scheduled samples are described as ASB; UTI episodes are UTI-1 through UTI-5 event samples.",
  "Event-label definition: Routine = ASB, UTI-* = UTI episode.",
  "Not stated as a row-level clinical rule in the thesis.",
  sprintf("Primary clinical status: %s. VF event/status cross-tab: %s.", collapse_counts(status_tab), fmt_vf_event_status(vf_event_status)),
  "Current primary UTI requires culture support plus compatible symptoms; Event_type is separate from UTI_Status.",
  sprintf("VF status: %s.", collapse_counts(vf_status_tab)),
  "Not apples-to-apples unless event labels are analysed separately",
  "Many UTI_event rows are Not_UTI under the current primary definition, so her ASB->UTI event analysis and your Not_UTI->UTI clinical analysis answer different questions.",
  "Checkable with current outputs"
)

add_row(
  "Assembly and wgMLST source",
  "Methods 2.2; Limitations",
  "Uses Flye assemblies, SeqSphere v9.0.1, E. coli wgMLST scheme of 4,512 loci, <=25 allele same-strain threshold, and excludes pairs with <100 comparable loci.",
  "RIVM wgMLST table from Flye assemblies.",
  "556 isolate allele profiles; pair denominators after comparable-loci filtering.",
  sprintf("Current canonical assembly selection is %s; current transition audit uses SNP/pairwise context, not wgMLST allele calls.", asm_text),
  "Current pipeline selected canonical assemblies, mostly Longcycler, and classifies strain context using <=25 SNPs plus ST context.",
  sprintf("%d selected canonical assemblies.", nrow(assemblies)),
  "Methodologically non-equivalent",
  "SNP and wgMLST are directional corroboration only. Her exact allele-distance results need the wgMLST export, comparable-loci counts and scheme version.",
  "Not exactly checkable"
)

add_row(
  "Consecutive pairs",
  "Methods 2.3; Results 3.2",
  "Pairs each isolate with the next available isolate from the subsequent/next available timepoint.",
  "Within-resident chronological pairs after wgMLST comparable-loci exclusions.",
  "368 consecutive pairs across 167 residents.",
  sprintf("%d cohort=='all' WGS/VF-linked consecutive transitions from %d participants.", nrow(transitions), safe_unique_n(transitions$Participant_id)),
  "vf_transitions_stratified.csv filtered to cohort == 'all'; ordered by Collection_Date when available, fallback timepoint order otherwise.",
  sprintf("%d transitions from %d participants.", nrow(transitions), safe_unique_n(transitions$Participant_id)),
  "Plausible filter difference but unreconciled",
  "Both analyses use next-observed within-resident pairs, but her wgMLST comparable-loci exclusion and resident universe are unavailable in current outputs.",
  "Partly checkable"
)

add_row(
  "Persistence/replacement categories",
  "Methods 2.3; Results 3.2; Discussion",
  "Methods define persistent as <=25 alleles and replacement as >25; Results split 26-100 as grey zone and >100 as replacement.",
  "wgMLST allele distance categories: <=25, 26-100, >100.",
  "231 persistent, 26 grey-zone, 111 replacement; all sum to 368.",
  sprintf("Strict SNP same %s; above 25 SNPs %s; replacement likely %s; same-lineage-not-same-strain-by-SNP %s.", fmt_pct(n_strict_same, nrow(transitions)), fmt_pct(n_above_snp, nrow(transitions)), fmt_pct(n_replacement, nrow(transitions)), fmt_pct(n_same_lineage_not_same_snp, nrow(transitions))),
  "Current strain context is <=25 SNPs for strict same strain; replacement_flag only when replacement is likely, not every >25-SNP same-ST pair.",
  sprintf("%d transitions.", nrow(transitions)),
  "Internally inconsistent in thesis; non-equivalent to current SNP audit",
  "The thesis toggles between >25 and >100 replacement definitions. Your analysis intentionally separates strict same strain from broader not-confidently-replaced categories.",
  "Partly checkable"
)

add_row(
  "Grey-zone handling",
  "Results 3.2; Discussion; Cox result",
  "Reports 26 pairs at 26-100 alleles, same ST and <60 SNPs; later says grey-zone pairs were classed as replacement; Cox has 137 events.",
  "Grey-zone pairs described separately but apparently counted as replacements in some analyses.",
  "26 grey-zone; 111 >100 replacement; 137 model events.",
  sprintf("Current above-threshold same-ST/same-lineage by SNP category is %d pairs; replacement likely category is %d pairs.", n_same_lineage_not_same_snp, n_replacement),
  "Current does not force same-ST >25 SNP pairs into replacement_flag.",
  sprintf("%d transitions.", nrow(transitions)),
  "Does not make sense as written",
  "111 + 26 = 137, so the model event count strongly suggests grey-zone pairs were treated as events after being presented as ambiguous.",
  "Needs thesis row-level table"
)

add_row(
  "Trajectory groups",
  "Methods 2.3; Results 3.2",
  "Residents with >=2 consecutive pairs are classified into stable persistence, multiple replacement, late replacement, single replacement, early replacement then stable.",
  "Requires the thesis replacement coding per pair, including grey-zone decisions.",
  "136 residents; group sizes 67, 23, 19, 16, 11.",
  sprintf("Current longitudinal transition universe has %s; transition participants=%d.", flow_lookup("Repeated-measures VF longitudinal subset"), safe_unique_n(transitions$Participant_id)),
  "Current outputs do not define the same five trajectory groups because wgMLST/grey-zone classifications are unavailable.",
  sprintf("%d transition participants.", safe_unique_n(transitions$Participant_id)),
  "Not directly checkable",
  "The trajectory denominators depend on her pair count and replacement coding; current data can audit transition counts but not her five-group assignment.",
  "Not checkable with current outputs"
)

add_row(
  "ST temporal trends",
  "Results 3.2",
  "Counts ST prevalence at T0 and T6 among isolates.",
  "All 556 isolates; denominator for percentages is not explicitly stated as all isolates or typed isolates.",
  "ST131 18.2% to 5.5%; ST73 6.8% to 16.4%; 83 STs.",
  sprintf("%d unique typed STs. %s; %s; %s; %s.", safe_unique_n(vf$ST), st_snapshot("T0", 131), st_snapshot("T6", 131), st_snapshot("T0", 73), st_snapshot("T6", 73)),
  "Current counts use vf_analysis_ready.csv; both all-isolate and typed-isolate denominators are reported.",
  sprintf("%d VF rows.", nrow(vf)),
  "Broadly compatible",
  "Direction and ST diversity match. Exact percentages differ modestly, likely due to typed-vs-all denominators or slightly different ST source/export.",
  "Checkable"
)

add_row(
  "Within-host genetic drift",
  "Methods 2.3; Results 3.3",
  "For residents with stable persistence and T0, computes wgMLST distance from T0 to later timepoints and describes median change as ~one allele/month.",
  "Subset selected on stable persistence plus T0 availability.",
  "39 residents; median T1 3.5 alleles and T6 15.5; none exceed 25.",
  "Current audit does not recompute T0-referenced wgMLST drift; SNP pairwise outputs are not a replacement.",
  "Current pairwise SNP metrics are not wgMLST allele trajectories and not restricted to her stable-persistence subset.",
  "Not applicable.",
  "Not directly checkable; rate claim unsupported",
  "The subset definition is circular for claims about staying under the same-strain threshold, and endpoint medians are not a formal time-based rate model.",
  "Not checkable with current outputs"
)

add_row(
  "Recurrence/re-emergence",
  "Methods 2.3; Results 3.3",
  "Defines recurrence as original strain re-emerging after displacement; returning strain <=25 alleles is relapse, otherwise reinfection.",
  "Requires original/displacing/returning isolate mapping.",
  "27 residents; 79% relapse; median 5 timepoints.",
  "Current outputs do not contain the original/displacing/returning wgMLST case table.",
  "Current pairwise/triple logic can test some pre/post returns but not her recurrence case definitions.",
  "Not applicable.",
  "Not checkable and wording inconsistent",
  "The paper says residents showed re-emergence, then says remaining cases are those where the original strain did not return.",
  "Not checkable with current outputs"
)

add_row(
  "ASB/routine to UTI event",
  "Methods 2.3; Results 3.4",
  "Identifies consecutive pairs where routine ASB timepoint is followed by UTI episode; compares wgMLST distance.",
  "Event-label transition: Routine -> UTI episode.",
  "34 ASB->UTI pairs across 26 residents; 16 same strain, 18 new strain.",
  paste(
    "Event-label sensitivity:", transition_summary(routine_to_event),
    "Primary clinical sensitivity:", transition_summary(routine_to_primary_uti),
    "All primary Not_UTI->UTI:", transition_summary(primary_notuti_to_uti)
  ),
  "Current audit separately evaluates Routine->UTI_event and current Not_UTI->UTI definitions.",
  sprintf("Routine->UTI_event n=%d; Routine->primary-UTI n=%d; all Not_UTI->UTI n=%d.", nrow(routine_to_event), nrow(routine_to_primary_uti), nrow(primary_notuti_to_uti)),
  "Definition-sensitive discrepancy",
  "Her conclusion can only be compared to event-label analyses. Under the current primary UTI rule, most evaluable Routine->UTI transitions are strict SNP same strain.",
  "Partly checkable"
)

add_row(
  "Post-UTI return of original strain",
  "Results 3.4",
  "Examines post-UTI pairs and asks whether original ASB strain returned after UTI/treatment.",
  "Requires pre-ASB, UTI, and post-UTI isolate mapping plus wgMLST distances.",
  "37 post-UTI pairs; original strain returned in 25/37; 94% after same-strain UTI and 48% after new-strain UTI.",
  sprintf("Current complete triples: event-label %s; primary-UTI %s.", triple_summary(triples), triple_summary(primary_triples)),
  "Current triple audit uses direct pre-routine to post-routine SNP context and separates event-label from primary UTI definitions.",
  sprintf("Event-label triples n=%d; primary-UTI triples n=%d.", nrow(triples), nrow(primary_triples)),
  "Not reproduced; method-sensitive",
  "The current complete-triple denominator is smaller and SNP-based. Exact replication requires her row-level pre/UTI/post wgMLST table.",
  "Partly checkable"
)

add_row(
  "Transition-type replacement rates",
  "Results 3.4",
  "Compares replacement rates across ASB->ASB, ASB->UTI and UTI->UTI transition categories.",
  "Uses thesis ASB/UTI labels across the 368 consecutive pairs.",
  "ASB->ASB n=336, ASB->UTI n=26, UTI->UTI n=6.",
  sprintf("Current primary status transitions: %s.", fmt_transition_table(transition_status_tab)),
  "Current uses UTI_Status transition_type after filtering cohort == 'all'.",
  sprintf("%d transitions.", nrow(transitions)),
  "Discrepant denominator/partition",
  "Her 336+26+6 partition omits UTI->ASB/Not_UTI despite later post-UTI analyses. Current primary table has Not_UTI->Not_UTI, Not_UTI->UTI, UTI->Not_UTI and UTI->UTI.",
  "Checkable as a mismatch"
)

add_row(
  "Antibiotic exposure linkage",
  "Methods 2.4; Results 3.5",
  "For T1-T6, records any antibiotic in previous 3 months and links exposure at Tn to pair T(n-1)->Tn.",
  "Pair-level Castor exposure joined to consecutive pairs.",
  "337/368 pairs have antibiotic exposure.",
  "Current pipeline does not contain Castor antibiotic exposure, indication, timing or dose fields.",
  "No proxy antibiotic variables were used in the audit.",
  "Not applicable.",
  "Not checkable; model denominator problem",
  "The exposure definition is plausible, but reporting n=337 exposure coverage and n=368 Cox model rows needs missing-exposure handling.",
  "Not checkable with current outputs"
)

add_row(
  "Antibiotics and replacement OR",
  "Results 3.5",
  "Compares replacement among pair-level exposed versus unexposed windows using Fisher's exact test.",
  "Pair-level exposed/unexposed among pairs with antibiotic linkage.",
  "56.0% replacement exposed vs 31.3% unexposed; OR 2.78, p=0.001.",
  "Not recalculated.",
  "Antibiotic data unavailable; wgMLST replacement coding unavailable.",
  "Not applicable.",
  "Not checkable; causal wording too strong",
  "Even if association is correct, confounding by indication and reverse causality remain unless exposure timing and indication are modelled.",
  "Not checkable with current outputs"
)

add_row(
  "Antibiotic type analysis",
  "Methods 2.4; Results 3.5",
  "Classifies residents by whether they received each antibiotic at any time; reports mean replacement rate per drug.",
  "Resident-level ever-exposure by drug, while no-antibiotic reference is reported as pairs.",
  "Amox-clav n=45 residents; nitrofurantoin n=31; fosfomycin n=18; ciprofloxacin n=17; no-antibiotic n=140 pairs.",
  "Not recalculated.",
  "Drug-specific exposure data unavailable; row unit appears mixed between residents and pairs.",
  "Not applicable.",
  "Methodologically unclear",
  "The text mixes resident-level ever exposure with pair-level replacement rates/reference counts. It needs the exact aggregation formula.",
  "Not checkable with current outputs"
)

add_row(
  "Antibiotics and subsequent sequenced UTI",
  "Results 3.5",
  "Classifies residents as antibiotic exposed or unexposed and compares development of at least one sequenced UTI episode.",
  "Resident-level ever-exposure; outcome is sequenced UTI episode.",
  "100 exposed residents, 66 unexposed; 24.0% vs 7.6%; OR 3.83, p=0.006.",
  sprintf("Current VF data: %d participants have UTI_event-labelled isolates; %d have primary-UTI VF isolates.", safe_unique_n(vf$Participant_id[vf$Event_type == "UTI_event"]), safe_unique_n(vf$Participant_id[vf$UTI_Status == "UTI"])),
  "Current audit has no antibiotic data and distinguishes event-labelled UTI from primary clinical UTI.",
  sprintf("Event-labelled UTI participants=%d; primary-UTI participants=%d.", safe_unique_n(vf$Participant_id[vf$Event_type == "UTI_event"]), safe_unique_n(vf$Participant_id[vf$UTI_Status == "UTI"])),
  "Outcome definition-sensitive; exposure not checkable",
  "The thesis outcome likely uses UTI_event/sequenced UTI labels. Your stricter primary-UTI outcome is much smaller.",
  "Partly checkable for outcome denominator only"
)

add_row(
  "Host factors",
  "Methods 2.4; Results 3.5",
  "Tests age, ward type, dementia, incontinence, Charlson and ADL against resident-level replacement rate.",
  "Resident-level host covariates and replacement-rate summary.",
  "No significant associations reported.",
  "Not recalculated.",
  "Demographic/host factor data unavailable in current generated anchors.",
  "Not applicable.",
  "Not checkable",
  "Non-significant univariable tests do not prove antibiotic effects are independent of host factors.",
  "Not checkable with current outputs"
)

add_row(
  "Cox model",
  "Methods 2.5; Results 3.5",
  "Fits Cox proportional-hazards model with consecutive pairs as the unit of analysis.",
  "Pair-level model, antibiotic exposure in preceding 3 months; sequence type apparently adjusted.",
  "n=368 pairs; 137 replacement events; HR 2.13.",
  sprintf("Current transitions have %d pairs from %d participants; replacement likely events=%d. No antibiotic exposure/time-to-event fields are available.", nrow(transitions), safe_unique_n(transitions$Participant_id), n_replacement),
  "Current audit can count repeated pair structure but cannot fit her Cox model.",
  sprintf("%d transitions/%d participants.", nrow(transitions), safe_unique_n(transitions$Participant_id)),
  "Methodologically unsupported as reported",
  "A Cox model needs time origin, survival time, censoring, recurrent-event handling, participant clustering, PH diagnostics and missing-exposure handling.",
  "Not checkable with current outputs"
)

add_row(
  "Repeated observations",
  "Methods 2.5",
  "Uses pair-level Fisher/Wilcoxon/Kruskal/Spearman/Cox tests.",
  "Repeated pairs nested within residents.",
  "368 pairs across 167 residents.",
  sprintf("Current transition table has %d pairs from %d participants, confirming repeated within-resident observations.", nrow(transitions), safe_unique_n(transitions$Participant_id)),
  "Current audit treats pair dependence as a design issue; no inferential re-analysis is run without the model-ready data.",
  sprintf("Median transitions per participant with transitions: %.1f.", median(as.integer(table(transitions$Participant_id)))),
  "Does not make sense for inference unless clustered/recurrent methods were used",
  "Pair rows are not independent. The thesis does not report clustering, frailty, mixed models, GEE or robust SEs.",
  "Checkable as a design concern"
)

add_row(
  "Interpretive causal claims",
  "Abstract; Discussion; Conclusion",
  "States antibiotics are a primary driver and findings support/provide evidence for protective colonisation.",
  "Inference from observational associations plus genomic transitions.",
  "No causal design described.",
  "Current independent data corroborate persistence and lineage trends but not antibiotic causality or direct protection.",
  "Current audit separates descriptive reproducibility from causal interpretation.",
  "Not applicable.",
  "Overstated",
  "The descriptive observations are plausible, but causal phrases require stronger design or careful qualification.",
  "Scientific interpretation"
)

recon <- do.call(rbind, rows)

csv_path <- file.path(out_dir, "calculation_filter_reconciliation.csv")
write.csv(recon, csv_path, row.names = FALSE, na = "")

md_path <- file.path(out_dir, "hamdi_thesis_calculation_filter_reconciliation.md")
lines <- c(
  "# Hamdi thesis calculation and filter reconciliation",
  "",
  sprintf("Generated: %s", format(Sys.Date())),
  "",
  "## Executive read",
  "",
  "The thesis calculations are easiest to understand as using a broader event-label genomic universe: scheduled routine samples are treated as ASB, UTI-* samples are treated as UTI events, and strain identity is based on wgMLST allele distances from a SeqSphere export. The current pipeline uses a stricter primary clinical `UTI_Status`, explicit `analysis_include_primary` and `genomics_expected_include` filters, canonical assembly selection, and SNP/ST-based strain-context labels. Therefore, several apparent conflicts are probably denominator/definition conflicts rather than direct computational contradictions.",
  "",
  "However, some items do not make sense as written. The largest are the grey-zone/replacement switch, the antibiotic model denominator, the ASB-to-UTI denominators, and the Cox model/repeated-observation specification.",
  "",
  "## Current filter anchors",
  "",
  sprintf("- Clinical status map: %d rows from %d participants; primary included: %d rows from %d participants.", nrow(clinical_all), safe_unique_n(clinical_all$Participant_id), nrow(clinical), safe_unique_n(clinical$Participant_id)),
  sprintf("- Primary clinical statuses: %s.", collapse_counts(status_tab)),
  sprintf("- VF/WGS-ready rows: %d rows from %d participants; statuses: %s.", nrow(vf), safe_unique_n(vf$Participant_id), collapse_counts(vf_status_tab)),
  sprintf("- Repeated-measures VF subset: %s.", flow_lookup("Repeated-measures VF longitudinal subset")),
  sprintf("- Consecutive WGS/VF transitions: %d rows from %d participants after `cohort == 'all'`.", nrow(transitions), safe_unique_n(transitions$Participant_id)),
  sprintf("- Current strain context: strict SNP same strain %s; replacement likely %s; same lineage but not same strain by SNP %s.", fmt_pct(n_strict_same, nrow(transitions)), fmt_pct(n_replacement, nrow(transitions)), fmt_pct(n_same_lineage_not_same_snp, nrow(transitions))),
  "",
  "## Key reconciliation points",
  "",
  "1. `Routine ASB` in the thesis is not the same object as your current `Not_UTI`. The thesis often uses event labels; your pipeline uses culture-plus-symptom clinical classification.",
  "2. `UTI_event` in your VF table is not always primary `UTI`. That explains why thesis-style Routine->UTI-event transitions and current Not_UTI->UTI transitions give different denominators and different same-strain proportions.",
  "3. Her 368 consecutive-pair denominator may reflect wgMLST comparable-loci filtering, but this cannot be verified without the row-level wgMLST pair table.",
  "4. Her persistence/replacement text is internally inconsistent: Methods imply >25 alleles is replacement, Results use >100 for replacement with 26-100 grey zone, and later analyses appear to count 111+26=137 replacement events.",
  "5. Antibiotic claims cannot be checked from current outputs and should not be recreated from proxies.",
  "6. The Cox model is under-specified because repeated resident-level observations, missing antibiotic exposure, time origin, censoring and proportional hazards diagnostics are not documented.",
  "",
  "## Calculation-by-calculation table",
  ""
)

for (i in seq_len(nrow(recon))) {
  lines <- c(
    lines,
    sprintf("### %02d. %s", i, recon$topic[[i]]),
    "",
    sprintf("- Thesis calculation: %s", recon$thesis_calculation[[i]]),
    sprintf("- Thesis filter/denominator: %s Denominator: %s", recon$thesis_filter[[i]], recon$thesis_denominator[[i]]),
    sprintf("- Current calculation: %s", recon$current_calculation[[i]]),
    sprintf("- Current filter/denominator: %s Denominator: %s", recon$current_filter[[i]], recon$current_denominator[[i]]),
    sprintf("- Compatibility: %s", recon$compatibility[[i]]),
    sprintf("- Rationale: %s", recon$rationale[[i]]),
    ""
  )
}

lines <- c(
  lines,
  "## What to ask Hamdi/supervisor for",
  "",
  "- Row-level wgMLST isolate table with comparable-loci counts, missing loci, scheme name/version/date, and allele distances.",
  "- Consecutive-pair table showing inclusion/exclusion rules and whether the 26 grey-zone pairs were excluded, retained as ambiguous, or counted as replacement.",
  "- Row-level Routine/ASB -> UTI and pre-UTI -> post-UTI case tables.",
  "- Castor antibiotic extract with indication, start/end dates, drug, exposure windows and missingness.",
  "- Cox model formula, survival time scale, censoring rules, clustering/recurrent-event method and diagnostics.",
  "- Supplementary tables referenced by the draft."
)

writeLines(lines, md_path)

cat("Wrote ", csv_path, "\n", sep = "")
cat("Wrote ", md_path, "\n", sep = "")
