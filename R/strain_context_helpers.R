# Shared same-strain/replacement helpers for longitudinal VF interpretation.

if (!exists("%||%", mode = "function")) {
  `%||%` <- function(x, y) {
    if (is.null(x) || length(x) == 0 || all(is.na(x))) y else x
  }
}

strain_snp_threshold <- function() {
  if (exists("SAME_STRAIN_SNP_THRESHOLD", inherits = TRUE)) {
    as.integer(get("SAME_STRAIN_SNP_THRESHOLD", inherits = TRUE))
  } else {
    25L
  }
}

normalise_strain_st_label <- function(x) {
  x <- as.character(x)
  x <- trimws(x)
  unknown <- c("", "-", "ST-", "NA", "N/A", "UNKNOWN", "UNK", "NT",
               "NON-TYPABLE", "NONTYPABLE", "NOT TYPED")
  x[toupper(x) %in% unknown] <- NA_character_
  x
}

normalise_strain_tp_label <- function(x) {
  if (exists("normalise_timepoint_preserve_events", mode = "function")) {
    normalise_timepoint_preserve_events(x)
  } else {
    trimws(as.character(x))
  }
}

prepare_pairwise_for_strain_context <- function(pairwise) {
  if (is.null(pairwise) || nrow(pairwise) == 0) return(pairwise)
  pairwise$Participant_id_A <- as.character(pairwise$Participant_id_A)
  pairwise$Participant_id_B <- as.character(pairwise$Participant_id_B)
  pairwise$Timepoint_A <- normalise_strain_tp_label(pairwise$Timepoint_A)
  pairwise$Timepoint_B <- normalise_strain_tp_label(pairwise$Timepoint_B)
  pairwise
}

match_pairwise_context <- function(pairwise, pid, tp_from, tp_to) {
  if (is.null(pairwise) || nrow(pairwise) == 0) return(NULL)
  pid <- as.character(pid)
  tp_from <- normalise_strain_tp_label(tp_from)
  tp_to <- normalise_strain_tp_label(tp_to)
  matched <- pairwise[
    (pairwise$Participant_id_A == pid & pairwise$Participant_id_B == pid &
       pairwise$Timepoint_A == tp_from & pairwise$Timepoint_B == tp_to) |
      (pairwise$Participant_id_A == pid & pairwise$Participant_id_B == pid &
         pairwise$Timepoint_A == tp_to & pairwise$Timepoint_B == tp_from),
    ,
    drop = FALSE
  ]
  if (nrow(matched) == 0) return(NULL)
  matched[1, , drop = FALSE]
}

match_evolution_context <- function(evol, pid, tp_from, tp_to) {
  if (is.null(evol) || nrow(evol) == 0) return(NULL)
  pid <- as.character(pid)
  tp_from <- normalise_strain_tp_label(tp_from)
  tp_to <- normalise_strain_tp_label(tp_to)
  matched <- evol[
    evol$Participant_id == pid &
      ((evol$From_Time == tp_from & evol$To_Time == tp_to) |
         (evol$From_Time == tp_to & evol$To_Time == tp_from)),
    ,
    drop = FALSE
  ]
  if (nrow(matched) == 0) return(NULL)
  matched[1, , drop = FALSE]
}

classify_strain_context <- function(has_vf_pair = TRUE,
                                    ST_from = NA_character_,
                                    ST_to = NA_character_,
                                    SNPs = NA_real_,
                                    AvgIdentity = NA_real_,
                                    Pairwise_Classification = NA_character_,
                                    Pairwise_RuleUsed = NA_character_,
                                    threshold = strain_snp_threshold()) {
  ST_from <- normalise_strain_st_label(ST_from)[1]
  ST_to <- normalise_strain_st_label(ST_to)[1]
  st_known_pair <- !is.na(ST_from) && !is.na(ST_to)
  same_ST <- if (st_known_pair) ST_from == ST_to else NA
  snps_known <- !is.na(SNPs)
  classification <- as.character(Pairwise_Classification %||% NA_character_)[1]
  threshold <- as.integer(threshold)

  snp_context <- if (!snps_known) {
    "Missing SNP evidence"
  } else if (SNPs <= threshold) {
    "Strong same strain"
  } else {
    "Above same-strain SNP threshold"
  }

  st_context <- if (!st_known_pair) {
    "Missing ST evidence"
  } else if (same_ST) {
    "Same ST"
  } else {
    "Different ST"
  }

  has_pairwise_classification <- !is.na(classification) && nzchar(classification)
  pair_interpretation <- dplyr::case_when(
    !isTRUE(has_vf_pair) ~ "Missing strain metrics",
    snp_context == "Strong same strain" && st_context == "Different ST" ~
      "Conflict: SNP same-strain but ST differs",
    snp_context == "Strong same strain" ~
      "Strong same strain",
    snp_context != "Strong same strain" &&
      (classification == "Different" || st_context == "Different ST") ~
      "Replacement likely",
    snp_context == "Above same-strain SNP threshold" && st_context == "Same ST" ~
      "Same lineage, not same strain by SNP",
    snp_context == "Missing SNP evidence" && st_context == "Same ST" ~
      "ST-consistent, SNP missing",
    snp_context == "Missing SNP evidence" && st_context == "Missing ST evidence" &&
      !has_pairwise_classification ~
      "Missing strain metrics",
    snp_context == "Above same-strain SNP threshold" ~
      "Above same-strain SNP threshold",
    snp_context == "Missing SNP evidence" ~
      "Missing SNP evidence",
    TRUE ~ "Missing strain metrics"
  )

  note <- sprintf(
    "SNP context: %s (same-strain threshold <=%d); ST context: %s; pairwise classification: %s.",
    snp_context,
    threshold,
    st_context,
    classification %||% "missing"
  )
  if (!isTRUE(has_vf_pair)) {
    note <- "One or both endpoints lack WGS/VF data"
  }

  pair_levels <- c(
    "Strong same strain",
    "Conflict: SNP same-strain but ST differs",
    "Same lineage, not same strain by SNP",
    "ST-consistent, SNP missing",
    "Above same-strain SNP threshold",
    "Missing SNP evidence",
    "Replacement likely",
    "Missing strain metrics"
  )

  tibble::tibble(
    ST_from = ST_from,
    ST_to = ST_to,
    same_ST = same_ST,
    SNPs = as.numeric(SNPs),
    AvgIdentity = as.numeric(AvgIdentity),
    Pairwise_Classification = classification,
    Pairwise_RuleUsed = as.character(Pairwise_RuleUsed %||% NA_character_)[1],
    snp_strain_context = factor(
      snp_context,
      levels = c("Strong same strain", "Above same-strain SNP threshold",
                 "Missing SNP evidence")
    ),
    st_lineage_context = factor(
      st_context,
      levels = c("Same ST", "Different ST", "Missing ST evidence")
    ),
    pair_interpretation = factor(pair_interpretation, levels = pair_levels),
    same_strain_evidence = as.character(pair_interpretation),
    strain_context_level = factor(
      pair_interpretation,
      levels = pair_levels
    ),
    replacement_flag = pair_interpretation == "Replacement likely",
    strain_context_note = note,
    same_strain_snp_threshold = threshold
  )
}

lookup_strain_context <- function(pairwise,
                                  pid,
                                  tp_from,
                                  tp_to,
                                  ST_from = NA_character_,
                                  ST_to = NA_character_,
                                  has_vf_pair = TRUE,
                                  evol = NULL,
                                  threshold = strain_snp_threshold()) {
  snps <- NA_real_
  avg_id <- NA_real_
  classification <- NA_character_
  rule_used <- NA_character_
  vf_jaccard_pairwise <- NA_real_
  inc_jaccard <- NA_real_

  ev <- match_evolution_context(evol, pid, tp_from, tp_to)
  if (!is.null(ev)) {
    if ("SNPs" %in% names(ev)) snps <- suppressWarnings(as.numeric(ev$SNPs[1]))
    if ("AvgIdentity" %in% names(ev)) avg_id <- suppressWarnings(as.numeric(ev$AvgIdentity[1]))
  }

  pw <- match_pairwise_context(pairwise, pid, tp_from, tp_to)
  if (!is.null(pw)) {
    if (is.na(snps) && "TotalSNPs" %in% names(pw)) snps <- suppressWarnings(as.numeric(pw$TotalSNPs[1]))
    if (is.na(avg_id) && "AvgIdentity" %in% names(pw)) avg_id <- suppressWarnings(as.numeric(pw$AvgIdentity[1]))
    if ("Classification" %in% names(pw)) classification <- as.character(pw$Classification[1])
    if ("RuleUsed" %in% names(pw)) rule_used <- as.character(pw$RuleUsed[1])
    if ("VF_Jaccard" %in% names(pw)) vf_jaccard_pairwise <- suppressWarnings(as.numeric(pw$VF_Jaccard[1]))
    if ("Inc_Jaccard" %in% names(pw)) inc_jaccard <- suppressWarnings(as.numeric(pw$Inc_Jaccard[1]))
  }

  classify_strain_context(
    has_vf_pair = has_vf_pair,
    ST_from = ST_from,
    ST_to = ST_to,
    SNPs = snps,
    AvgIdentity = avg_id,
    Pairwise_Classification = classification,
    Pairwise_RuleUsed = rule_used,
    threshold = threshold
  ) |>
    dplyr::mutate(
      VF_Jaccard_pairwise = vf_jaccard_pairwise,
      Inc_Jaccard = inc_jaccard,
      Classification = Pairwise_Classification,
      RuleUsed = Pairwise_RuleUsed
    )
}
