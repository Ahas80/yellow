# Helpers for Script 29 predicted-plasmid localization.
#
# All classifications here are assembly based. "Predicted linkage" means that
# two calls were placed in the same MOB-recon plasmid bin; it does not establish
# circularity, transfer, transmission, or causality.

plasmid_atomic_write_csv <- function(x, path) {
  dir.create(dirname(path), recursive = TRUE, showWarnings = FALSE)
  tmp <- tempfile(
    pattern = paste0(".", basename(path), "."),
    tmpdir = dirname(path), fileext = ".tmp"
  )
  readr::write_csv(x, tmp, na = "")
  if (file.exists(path) && !file.remove(path)) {
    stop("Could not replace plasmid mechanism output: ", path)
  }
  if (!file.rename(tmp, path)) {
    stop("Could not publish plasmid mechanism output: ", path)
  }
  invisible(path)
}

plasmid_split_set <- function(x) {
  if (length(x) != 1L || is.na(x) || !nzchar(trimws(x))) {
    return(character())
  }
  z <- trimws(unlist(strsplit(x, ";", fixed = TRUE), use.names = FALSE))
  sort(unique(z[nzchar(z) & z != "-"]))
}

plasmid_collapse_set <- function(x) {
  x <- sort(unique(as.character(x[!is.na(x) & nzchar(trimws(x))])))
  if (length(x)) paste(x, collapse = ";") else ""
}

plasmid_set_metrics <- function(from, to) {
  from <- plasmid_split_set(from)
  to <- plasmid_split_set(to)
  union <- union(from, to)
  intersection <- intersect(from, to)
  list(
    jaccard = if (!length(union)) 1 else length(intersection) / length(union),
    both_empty = !length(from) && !length(to),
    gained = setdiff(to, from),
    lost = setdiff(from, to)
  )
}

plasmid_sequence_table <- function(path, id_name) {
  sequences <- Biostrings::readDNAStringSet(path)
  ids <- sub("[[:space:]].*$", "", names(sequences))
  if (any(!nzchar(ids)) || anyDuplicated(ids)) {
    stop("FASTA contains empty or duplicated identifiers: ", path)
  }
  hashes <- vapply(
    as.character(sequences),
    digest::digest, character(1),
    algo = "sha256", serialize = FALSE
  )
  out <- tibble::tibble(
    contig_id = ids,
    sequence_sha256 = hashes,
    sequence_length = as.integer(Biostrings::width(sequences)),
    fasta_order = seq_along(ids)
  )
  key <- paste(out$sequence_sha256, out$sequence_length, sep = "|")
  out$duplicate_occurrence <- ave(
    seq_along(key), key, FUN = function(i) seq_along(i)
  )
  names(out)[names(out) == "contig_id"] <- id_name
  out
}

plasmid_build_prokka_contig_map <- function(manifest) {
  maps <- lapply(seq_len(nrow(manifest)), function(i) {
    original <- plasmid_sequence_table(
      manifest$fasta_path[[i]], "original_contig_id"
    )
    prokka <- plasmid_sequence_table(
      manifest$fna_path[[i]], "prokka_contig_id"
    )
    by <- c(
      "sequence_sha256", "sequence_length", "duplicate_occurrence"
    )
    mapped <- prokka %>%
      dplyr::left_join(
        original %>%
          dplyr::select(
            dplyr::all_of(by), original_contig_id,
            original_fasta_order = fasta_order
          ),
        by = by, relationship = "one-to-one"
      ) %>%
      dplyr::mutate(
        Isolate_ID = as.character(manifest$Isolate_ID[[i]]),
        Assembly_ID = as.character(manifest$Assembly_ID[[i]]),
        Participant_id = as.character(manifest$Participant_id[[i]]),
        tp_lab = as.character(manifest$tp_lab[[i]]),
        fasta_path = as.character(manifest$fasta_path[[i]]),
        prokka_fna_path = as.character(manifest$fna_path[[i]])
      ) %>%
      dplyr::relocate(
        Isolate_ID, Assembly_ID, Participant_id, tp_lab,
        fasta_path, prokka_fna_path, prokka_contig_id,
        original_contig_id
      )
    if (
      nrow(mapped) != nrow(original) ||
        anyNA(mapped$original_contig_id) ||
        anyDuplicated(mapped$prokka_contig_id) ||
        anyDuplicated(mapped$original_contig_id) ||
        !setequal(mapped$original_contig_id, original$original_contig_id)
    ) {
      stop(
        "Prokka-to-original contig sequence map is not full one-to-one for ",
        manifest$Isolate_ID[[i]], call. = FALSE
      )
    }
    mapped
  })
  dplyr::bind_rows(maps)
}

plasmid_add_mob_location <- function(calls, mob) {
  calls %>%
    dplyr::left_join(
      mob %>%
        dplyr::select(
          Isolate_ID,
          original_contig_id = contig_id,
          mob_molecule_type = molecule_type,
          mob_primary_cluster = primary_cluster_id,
          mob_secondary_cluster = secondary_cluster_id,
          mob_assignment_confidence = assignment_confidence,
          mob_uncertainty_reason = uncertainty_reason,
          mob_contig_size = input_size
        ),
      by = c("Isolate_ID", "original_contig_id"),
      relationship = "many-to-one"
    ) %>%
    dplyr::mutate(
      localization = dplyr::case_when(
        mob_molecule_type == "predicted_plasmid" ~ "predicted_plasmid",
        mob_molecule_type == "predicted_chromosome" ~ "chromosome",
        TRUE ~ "ambiguous_or_unassigned"
      ),
      localization_status = dplyr::case_when(
        is.na(reported_contig_id) | !nzchar(reported_contig_id) ~
          "reported_contig_missing",
        is.na(original_contig_id) | !nzchar(original_contig_id) ~
          "original_contig_mapping_unresolved",
        is.na(mob_molecule_type) ~ "original_contig_absent_from_mob_report",
        mob_molecule_type == "unassigned" ~ "explicitly_unassigned_by_mob",
        TRUE ~ "mapped"
      ),
      predicted_plasmid_bin_id = dplyr::if_else(
        localization == "predicted_plasmid" &
          !is.na(mob_primary_cluster) &
          nzchar(mob_primary_cluster) &
          mob_primary_cluster != "-",
        paste(Isolate_ID, mob_primary_cluster, sep = "::"),
        NA_character_
      ),
      predicted_linkage_interpretation = dplyr::if_else(
        localization == "predicted_plasmid" &
          !is.na(predicted_plasmid_bin_id),
        "placed in a MOB predicted plasmid bin; same-bin calls are predicted linkage only",
        "no predicted same-plasmid-bin linkage assigned"
      ),
      interpretation_scope =
        "assembly-based localization; not confirmed circularity, transfer, transmission, phenotype, or causality"
    )
}

plasmid_profile_changes <- function(transitions, episode_profiles) {
  episode_profiles <- episode_profiles %>%
    dplyr::mutate(
      Participant_id = as.character(Participant_id),
      tp_lab = as.character(tp_lab)
    )
  profile_columns <- c(
    "replicon_markers", "mob_primary_clusters",
    "predicted_plasmid_bin_ids", "plasmid_binned_vf_genes",
    "plasmid_binned_informative_amr_genes",
    "plasmid_binned_informative_vf_amr_features",
    "predicted_linkage_groups"
  )
  from <- episode_profiles %>%
    dplyr::select(
      Participant_id, tp_from = tp_lab,
      predicted_plasmid_count_from = predicted_plasmid_count,
      mob_high_confidence_profile_from = mob_high_confidence_profile,
      dplyr::all_of(profile_columns)
    ) %>%
    dplyr::rename_with(
      ~ paste0(.x, "_from"), dplyr::all_of(profile_columns)
    )
  to <- episode_profiles %>%
    dplyr::select(
      Participant_id, tp_to = tp_lab,
      predicted_plasmid_count_to = predicted_plasmid_count,
      mob_high_confidence_profile_to = mob_high_confidence_profile,
      dplyr::all_of(profile_columns)
    ) %>%
    dplyr::rename_with(
      ~ paste0(.x, "_to"), dplyr::all_of(profile_columns)
    )
  out <- transitions %>%
    dplyr::mutate(
      Participant_id = as.character(Participant_id),
      tp_from = as.character(tp_from),
      tp_to = as.character(tp_to)
    ) %>%
    dplyr::left_join(
      from, by = c("Participant_id", "tp_from"),
      relationship = "many-to-one"
    ) %>%
    dplyr::left_join(
      to, by = c("Participant_id", "tp_to"),
      relationship = "many-to-one"
    )
  if (anyNA(out$predicted_plasmid_count_from) ||
      anyNA(out$predicted_plasmid_count_to)) {
    stop("A transition endpoint lacks a mechanism profile.", call. = FALSE)
  }
  changes <- lapply(seq_len(nrow(out)), function(i) {
    rep <- plasmid_set_metrics(
      out$replicon_markers_from[[i]], out$replicon_markers_to[[i]]
    )
    mob <- plasmid_set_metrics(
      out$mob_primary_clusters_from[[i]],
      out$mob_primary_clusters_to[[i]]
    )
    bins <- plasmid_set_metrics(
      out$predicted_plasmid_bin_ids_from[[i]],
      out$predicted_plasmid_bin_ids_to[[i]]
    )
    vf <- plasmid_set_metrics(
      out$plasmid_binned_vf_genes_from[[i]],
      out$plasmid_binned_vf_genes_to[[i]]
    )
    amr <- plasmid_set_metrics(
      out$plasmid_binned_informative_amr_genes_from[[i]],
      out$plasmid_binned_informative_amr_genes_to[[i]]
    )
    combined <- plasmid_set_metrics(
      out$plasmid_binned_informative_vf_amr_features_from[[i]],
      out$plasmid_binned_informative_vf_amr_features_to[[i]]
    )
    tibble::tibble(
      replicon_jaccard = rep$jaccard,
      replicon_both_empty = rep$both_empty,
      replicons_gained = plasmid_collapse_set(rep$gained),
      replicons_lost = plasmid_collapse_set(rep$lost),
      any_replicon_profile_change =
        length(rep$gained) + length(rep$lost) > 0L,
      mob_cluster_jaccard = mob$jaccard,
      mob_cluster_both_empty = mob$both_empty,
      mob_clusters_gained = plasmid_collapse_set(mob$gained),
      mob_clusters_lost = plasmid_collapse_set(mob$lost),
      any_mob_cluster_change =
        length(mob$gained) + length(mob$lost) > 0L,
      predicted_plasmid_bins_gained =
        plasmid_collapse_set(bins$gained),
      predicted_plasmid_bins_lost =
        plasmid_collapse_set(bins$lost),
      predicted_plasmid_count_difference =
        out$predicted_plasmid_count_to[[i]] -
        out$predicted_plasmid_count_from[[i]],
      plasmid_binned_vf_genes_gained =
        plasmid_collapse_set(vf$gained),
      plasmid_binned_vf_genes_lost =
        plasmid_collapse_set(vf$lost),
      plasmid_binned_informative_amr_genes_gained =
        plasmid_collapse_set(amr$gained),
      plasmid_binned_informative_amr_genes_lost =
        plasmid_collapse_set(amr$lost),
      plasmid_binned_informative_features_gained =
        plasmid_collapse_set(combined$gained),
      plasmid_binned_informative_features_lost =
        plasmid_collapse_set(combined$lost),
      mob_high_confidence_profiles_both =
        isTRUE(out$mob_high_confidence_profile_from[[i]]) &&
        isTRUE(out$mob_high_confidence_profile_to[[i]]),
      mechanism_interpretation =
        "descriptive predicted-bin change; no physical transfer or causal UTI inference"
    )
  })
  dplyr::bind_cols(out, dplyr::bind_rows(changes))
}

run_plasmid_gene_localization <- function(
    root, amr_analysis,
    output_root = file.path(root, "results", "plasmids", "mob_suite")) {
  required_packages <- c(
    "Biostrings", "digest", "dplyr", "readr", "tibble", "tidyr"
  )
  absent <- required_packages[
    !vapply(required_packages, requireNamespace, logical(1), quietly = TRUE)
  ]
  if (length(absent)) {
    stop(
      "Missing package(s) for predicted-plasmid localization: ",
      paste(absent, collapse = ", ")
    )
  }
  required_files <- c(
    mob_marker = file.path(output_root, "RUN_COMPLETE.txt"),
    mob_contigs = file.path(output_root, "contig_assignments.csv"),
    mob_profiles = file.path(output_root, "episode_plasmid_profiles.csv"),
    plasmidfinder_hits = file.path(
      root, "results", "plasmids", "plasmidfinder_hits_long.csv"
    ),
    plasmidfinder_pa = file.path(
      root, "results", "plasmids",
      "plasmidfinder_presence_absence.csv"
    ),
    plasmidfinder_catalog = file.path(
      root, "results", "plasmids",
      "plasmidfinder_replicon_catalog.csv"
    ),
    virulencefinder_marker = file.path(
      root, "results", "virulencefinder_cge_3_2_1", "RUN_COMPLETE.txt"
    ),
    virulencefinder_hits = file.path(
      root, "results", "virulencefinder_cge_3_2_1", "hits_long.csv"
    )
  )
  missing <- required_files[!file.exists(required_files)]
  if (length(missing)) {
    stop(
      "Predicted-plasmid localization inputs are incomplete: ",
      paste(missing, collapse = ", "), call. = FALSE
    )
  }
  manifest <- amr_analysis$manifest %>%
    dplyr::mutate(
      Isolate_ID = as.character(Isolate_ID),
      Assembly_ID = as.character(Assembly_ID),
      Participant_id = as.character(Participant_id),
      tp_lab = as.character(tp_lab),
      episode_key = paste(Participant_id, tp_lab, sep = "||")
    )
  if (nrow(manifest) != 532L || anyDuplicated(manifest$episode_key)) {
    stop("Localization manifest is not exactly 532 episodes.", call. = FALSE)
  }

  mob <- readr::read_csv(
    required_files[["mob_contigs"]], show_col_types = FALSE
  ) %>%
    dplyr::mutate(Isolate_ID = as.character(Isolate_ID))
  mob_profiles <- readr::read_csv(
    required_files[["mob_profiles"]], show_col_types = FALSE
  ) %>%
    dplyr::mutate(
      Isolate_ID = as.character(Isolate_ID),
      Participant_id = as.character(Participant_id),
      tp_lab = as.character(tp_lab)
    )
  if (
    dplyr::n_distinct(mob$Isolate_ID) != 532L ||
      anyDuplicated(mob[c("Isolate_ID", "contig_id")]) ||
      nrow(mob_profiles) != 532L ||
      anyDuplicated(mob_profiles$Isolate_ID)
  ) {
    stop("MOB contig/profile denominator contract failed.", call. = FALSE)
  }

  prokka_map <- plasmid_build_prokka_contig_map(manifest)
  plasmid_atomic_write_csv(
    prokka_map, file.path(output_root, "prokka_original_contig_map.csv")
  )

  vf <- readr::read_csv(
    required_files[["virulencefinder_hits"]],
    show_col_types = FALSE
  ) %>%
    dplyr::filter(profile == "web_default_id90_cov60") %>%
    dplyr::semi_join(
      manifest %>% dplyr::select(episode_key),
      by = "episode_key"
    ) %>%
    dplyr::mutate(
      Participant_id = as.character(Participant_id),
      tp_lab = as.character(tp_lab),
      Isolate_ID = manifest$Isolate_ID[
        match(episode_key, manifest$episode_key)
      ],
      call_id = paste0("VF_CGE_", dplyr::row_number()),
      evidence_source = "VF",
      caller = "CGE_VirulenceFinder_3.2.1_web_default_id90_cov60",
      determinant_type = "virulence_gene",
      gene = as.character(gene_family),
      reported_contig_id = as.character(query_id),
      original_contig_id = reported_contig_id,
      start = as.numeric(query_start),
      end = as.numeric(query_end),
      identity = as.numeric(identity_pct),
      coverage = as.numeric(coverage_pct),
      background_flag = FALSE,
      analysis_informative = TRUE,
      coordinate_mapping_method =
        "CGE query_id is the original selected Longcycler contig ID"
    ) %>%
    dplyr::select(
      call_id, Participant_id, tp_lab, episode_key, Assembly_ID,
      Isolate_ID, evidence_source, caller, determinant_type, gene,
      reported_contig_id, original_contig_id, start, end,
      identity, coverage, background_flag, analysis_informative,
      coordinate_mapping_method
    )

  amr <- amr_analysis$long %>%
    dplyr::mutate(
      Participant_id = as.character(Participant_id),
      tp_lab = as.character(tp_lab),
      episode_key = paste(Participant_id, tp_lab, sep = "||")
    ) %>%
    dplyr::left_join(
      manifest %>%
        dplyr::select(
          episode_key, Isolate_ID,
          manifest_Assembly_ID = Assembly_ID
        ),
      by = "episode_key", relationship = "many-to-one"
    ) %>%
    dplyr::mutate(
      call_id = paste0("AMR_", dplyr::row_number()),
      evidence_source = "AMR",
      gene = as.character(normalized_symbol),
      reported_contig_id = as.character(contig),
      original_contig_id = dplyr::if_else(
        caller == "amrfinderplus", NA_character_, reported_contig_id
      ),
      analysis_informative =
        caller == "amrfinderplus" &
        determinant_type == "acquired_gene" &
        !background_flag,
      coordinate_mapping_method = dplyr::if_else(
        caller == "amrfinderplus",
        "Prokka contig mapped by sequence SHA-256, length, and duplicate occurrence",
        "caller used the original selected Longcycler FASTA"
      ),
      Assembly_ID = dplyr::coalesce(
        as.character(Assembly_ID), as.character(manifest_Assembly_ID)
      )
    ) %>%
    dplyr::left_join(
      prokka_map %>%
        dplyr::select(
          Isolate_ID, prokka_contig_id,
          prokka_original_contig_id = original_contig_id
        ),
      by = c(
        "Isolate_ID",
        "reported_contig_id" = "prokka_contig_id"
      ),
      relationship = "many-to-one"
    ) %>%
    dplyr::mutate(
      original_contig_id = dplyr::if_else(
        caller == "amrfinderplus",
        prokka_original_contig_id, original_contig_id
      )
    ) %>%
    dplyr::select(
      call_id, Participant_id, tp_lab, episode_key, Assembly_ID,
      Isolate_ID, evidence_source, caller, determinant_type, gene,
      reported_contig_id, original_contig_id, start, end,
      identity, coverage, background_flag, analysis_informative,
      coordinate_mapping_method
    )

  pf_raw <- readr::read_csv(
    required_files[["plasmidfinder_hits"]],
    show_col_types = FALSE
  )
  required_pf_columns <- c(
    "Isolate_ID", "Assembly_ID", "GENE", "SEQUENCE",
    "start", "end", "identity", "coverage"
  )
  missing_pf_columns <- setdiff(required_pf_columns, names(pf_raw))
  if (length(missing_pf_columns)) {
    stop(
      "Canonical PlasmidFinder hits lack coordinate column(s): ",
      paste(missing_pf_columns, collapse = ", "),
      ". Rerun 09_inc_plasmid_network.R to derive coordinate-bearing hits ",
      "from the hash-bound raw ABRicate caches.",
      call. = FALSE
    )
  }
  pf <- pf_raw %>%
    dplyr::mutate(
      Isolate_ID = as.character(Isolate_ID),
      episode_key = manifest$episode_key[
        match(Isolate_ID, manifest$Isolate_ID)
      ],
      Participant_id = manifest$Participant_id[
        match(Isolate_ID, manifest$Isolate_ID)
      ],
      tp_lab = manifest$tp_lab[match(Isolate_ID, manifest$Isolate_ID)],
      call_id = paste0("PF_", dplyr::row_number()),
      evidence_source = "PlasmidFinder",
      caller = "ABRicate_PlasmidFinder_primary_80_80",
      determinant_type = "replicon_marker",
      gene = as.character(GENE),
      reported_contig_id = as.character(SEQUENCE),
      original_contig_id = reported_contig_id,
      start = as.numeric(start),
      end = as.numeric(end),
      identity = as.numeric(identity),
      coverage = as.numeric(coverage),
      background_flag = FALSE,
      analysis_informative = FALSE,
      coordinate_mapping_method =
        "ABRicate SEQUENCE is the original selected Longcycler contig ID"
    ) %>%
    dplyr::select(
      call_id, Participant_id, tp_lab, episode_key, Assembly_ID,
      Isolate_ID, evidence_source, caller, determinant_type, gene,
      reported_contig_id, original_contig_id, start, end,
      identity, coverage, background_flag, analysis_informative,
      coordinate_mapping_method
    )

  locations <- dplyr::bind_rows(vf, amr, pf) %>%
    plasmid_add_mob_location(mob) %>%
    dplyr::arrange(
      Participant_id, tp_lab, evidence_source, gene, call_id
    )
  if (
    anyNA(locations$localization) ||
      anyNA(locations$localization_status)
  ) {
    stop(
      "Every VF/AMR/replicon call must be mapped or explicitly ambiguous.",
      call. = FALSE
    )
  }
  plasmid_atomic_write_csv(
    locations,
    file.path(output_root, "plasmid_gene_locations_long.csv")
  )

  linkage <- locations %>%
    dplyr::filter(
      localization == "predicted_plasmid",
      !is.na(predicted_plasmid_bin_id)
    ) %>%
    dplyr::group_by(
      Participant_id, tp_lab, episode_key, Isolate_ID,
      predicted_plasmid_bin_id, mob_primary_cluster
    ) %>%
    dplyr::summarise(
      vf_genes = plasmid_collapse_set(
        gene[evidence_source == "VF"]
      ),
      informative_amr_genes = plasmid_collapse_set(
        gene[evidence_source == "AMR" & analysis_informative]
      ),
      replicon_markers = plasmid_collapse_set(
        gene[evidence_source == "PlasmidFinder"]
      ),
      n_calls = dplyr::n(),
      linkage_statement =
        "calls share one MOB predicted plasmid bin; predicted linkage only",
      .groups = "drop"
    )
  plasmid_atomic_write_csv(
    linkage,
    file.path(output_root, "predicted_plasmid_linkage_groups.csv")
  )

  pf_pa <- readr::read_csv(
    required_files[["plasmidfinder_pa"]], show_col_types = FALSE
  )
  pf_catalog <- readr::read_csv(
    required_files[["plasmidfinder_catalog"]], show_col_types = FALSE
  )
  feature_columns <- sort(unique(as.character(pf_catalog$GENE)))
  if (
    nrow(pf_pa) != 532L ||
      !all(feature_columns %in% names(pf_pa))
  ) {
    stop("Corrected PlasmidFinder matrix contract failed.", call. = FALSE)
  }
  replicon_profiles <- pf_pa %>%
    tidyr::pivot_longer(
      dplyr::all_of(feature_columns),
      names_to = "replicon", values_to = "present"
    ) %>%
    dplyr::filter(as.numeric(present) > 0) %>%
    dplyr::group_by(Isolate_ID) %>%
    dplyr::summarise(
      replicon_markers = plasmid_collapse_set(replicon),
      .groups = "drop"
    )

  location_profile <- locations %>%
    dplyr::group_by(
      Participant_id, tp_lab, episode_key, Isolate_ID
    ) %>%
    dplyr::summarise(
      plasmid_binned_vf_genes = plasmid_collapse_set(
        gene[
          evidence_source == "VF" &
            localization == "predicted_plasmid"
        ]
      ),
      plasmid_binned_informative_amr_genes =
        plasmid_collapse_set(
          gene[
            evidence_source == "AMR" &
              analysis_informative &
              localization == "predicted_plasmid"
          ]
        ),
      plasmid_binned_informative_vf_amr_features =
        plasmid_collapse_set(c(
          paste0(
            "VF:", gene[
              evidence_source == "VF" &
                localization == "predicted_plasmid"
            ]
          ),
          paste0(
            "AMR:", gene[
              evidence_source == "AMR" &
                analysis_informative &
                localization == "predicted_plasmid"
            ]
          )
        )),
      chromosome_binned_informative_vf_amr_features =
        plasmid_collapse_set(c(
          paste0(
            "VF:", gene[
              evidence_source == "VF" &
                localization == "chromosome"
            ]
          ),
          paste0(
            "AMR:", gene[
              evidence_source == "AMR" &
                analysis_informative &
                localization == "chromosome"
            ]
          )
        )),
      ambiguous_or_unassigned_informative_call_count =
        sum(
          analysis_informative &
            localization == "ambiguous_or_unassigned"
        ),
      .groups = "drop"
    )
  bin_profiles <- linkage %>%
    dplyr::group_by(Isolate_ID) %>%
    dplyr::summarise(
      predicted_plasmid_bin_ids =
        plasmid_collapse_set(predicted_plasmid_bin_id),
      predicted_linkage_groups = plasmid_collapse_set(
        paste0(
          predicted_plasmid_bin_id, "[VF=", vf_genes,
          "|AMR=", informative_amr_genes,
          "|replicon=", replicon_markers, "]"
        )
      ),
      .groups = "drop"
    )
  episode_profiles <- manifest %>%
    dplyr::select(
      Participant_id, tp_lab, episode_key, Assembly_ID, Isolate_ID,
      fasta_path, fasta_sha256
    ) %>%
    dplyr::left_join(
      mob_profiles %>%
        dplyr::select(
          Isolate_ID, predicted_plasmid_count,
          predicted_plasmid_bp, mob_primary_clusters,
          mob_high_confidence_profile
        ),
      by = "Isolate_ID", relationship = "one-to-one"
    ) %>%
    dplyr::left_join(
      replicon_profiles, by = "Isolate_ID",
      relationship = "one-to-one"
    ) %>%
    dplyr::left_join(
      location_profile,
      by = c("Participant_id", "tp_lab", "episode_key", "Isolate_ID"),
      relationship = "one-to-one"
    ) %>%
    dplyr::left_join(
      bin_profiles, by = "Isolate_ID", relationship = "one-to-one"
    ) %>%
    dplyr::mutate(
      dplyr::across(
        c(
          replicon_markers, mob_primary_clusters,
          plasmid_binned_vf_genes,
          plasmid_binned_informative_amr_genes,
          plasmid_binned_informative_vf_amr_features,
          chromosome_binned_informative_vf_amr_features,
          predicted_plasmid_bin_ids, predicted_linkage_groups
        ),
        ~ tidyr::replace_na(as.character(.x), "")
      ),
      ambiguous_or_unassigned_informative_call_count =
        tidyr::replace_na(
          ambiguous_or_unassigned_informative_call_count, 0L
        ),
      plasmid_binned_vf_burden = lengths(lapply(
        plasmid_binned_vf_genes, plasmid_split_set
      )),
      plasmid_binned_informative_amr_burden = lengths(lapply(
        plasmid_binned_informative_amr_genes, plasmid_split_set
      )),
      plasmid_binned_informative_vf_amr_burden = lengths(lapply(
        plasmid_binned_informative_vf_amr_features,
        plasmid_split_set
      )),
      interpretation_scope =
        "CGE VF 90/60 plus primary informative AMRFinder genes placed in assembly-based MOB bins"
    )
  if (nrow(episode_profiles) != 532L ||
      anyDuplicated(episode_profiles$episode_key) ||
      anyNA(episode_profiles$predicted_plasmid_count)) {
    stop("Episode mechanism profile contract failed.", call. = FALSE)
  }
  plasmid_atomic_write_csv(
    episode_profiles,
    file.path(output_root, "episode_mechanism_profiles.csv")
  )

  adjacent <- plasmid_profile_changes(
    amr_analysis$transitions, episode_profiles
  )
  focused <- adjacent %>%
    dplyr::filter(status_from == "Not_UTI", status_to == "UTI")
  if (nrow(adjacent) != 371L || nrow(focused) != 9L) {
    stop("Plasmid pair denominator contract is not exactly 371/9.")
  }
  plasmid_atomic_write_csv(
    adjacent,
    file.path(output_root, "adjacent_pair_plasmid_metrics_371.csv")
  )
  plasmid_atomic_write_csv(
    focused,
    file.path(output_root, "not_uti_to_uti_plasmid_metrics_9.csv")
  )

  pilot <- locations %>%
    dplyr::filter(Isolate_ID == "2510C119001-1")
  pilot_mob <- mob %>%
    dplyr::filter(Isolate_ID == "2510C119001-1")
  pilot_mob_incf <- pilot_mob %>%
    dplyr::filter(
      molecule_type == "predicted_plasmid",
      primary_cluster_id == "AA179"
    )
  pilot_mob_incf_types <- sort(unique(trimws(unlist(
    strsplit(
      pilot_mob_incf$`rep_type(s)`[
        !is.na(pilot_mob_incf$`rep_type(s)`)
      ],
      ",", fixed = TRUE
    ),
    use.names = FALSE
  ))))
  pilot_pf <- pilot %>%
    dplyr::filter(
      evidence_source == "PlasmidFinder",
      gene %in% c(
        "IncFIB(AP001918)_1", "IncFIC(FII)_1", "IncFIA_1"
      )
    )
  pilot_vf <- pilot %>% dplyr::filter(evidence_source == "VF")
  pilot_amr <- pilot %>%
    dplyr::filter(
      evidence_source == "AMR", caller == "amrfinderplus"
    )
  pilot_checks <- tibble::tribble(
    ~check, ~expected, ~observed, ~pass,
    "pilot_contigs_accounted", "9", as.character(nrow(pilot_mob)),
    nrow(pilot_mob) == 9L,
    "pilot_chromosome_contigs", "3",
    as.character(sum(pilot_mob$molecule_type == "predicted_chromosome")),
    sum(pilot_mob$molecule_type == "predicted_chromosome") == 3L,
    "pilot_predicted_plasmid_contigs", "6",
    as.character(sum(pilot_mob$molecule_type == "predicted_plasmid")),
    sum(pilot_mob$molecule_type == "predicted_plasmid") == 6L,
    "pilot_PlasmidFinder_IncF_markers_same_bin", "one AA179 bin",
    plasmid_collapse_set(pilot_pf$mob_primary_cluster),
    nrow(pilot_pf) == 3L &&
      dplyr::n_distinct(pilot_pf$predicted_plasmid_bin_id) == 1L &&
      unique(pilot_pf$mob_primary_cluster) == "AA179",
    "pilot_MOB_IncF_replicons_same_bin", "IncFIB;IncFIC;IncFII on AA179",
    paste0(
      plasmid_collapse_set(pilot_mob_incf_types),
      " on ", plasmid_collapse_set(pilot_mob_incf$primary_cluster_id)
    ),
    nrow(pilot_mob_incf) == 1L &&
      all(c("IncFIB", "IncFIC", "IncFII") %in% pilot_mob_incf_types),
    "pilot_vf_predicted_plasmid", "14",
    as.character(sum(pilot_vf$localization == "predicted_plasmid")),
    sum(pilot_vf$localization == "predicted_plasmid") == 14L,
    "pilot_vf_chromosome", "35",
    as.character(sum(pilot_vf$localization == "chromosome")),
    sum(pilot_vf$localization == "chromosome") == 35L,
    "pilot_blaTEM1_AA179", "predicted plasmid cluster AA179",
    plasmid_collapse_set(
      pilot_amr$mob_primary_cluster[
        grepl("^blaTEM-1", pilot_amr$gene)
      ]
    ),
    any(
      grepl("^blaTEM-1", pilot_amr$gene) &
        pilot_amr$localization == "predicted_plasmid" &
        pilot_amr$mob_primary_cluster == "AA179"
    ),
    "pilot_cyaA_S352T_chromosome", "chromosome",
    plasmid_collapse_set(
      pilot_amr$localization[
        grepl("^cyaA[:_]?S352T$", pilot_amr$gene)
      ]
    ),
    any(
      grepl("^cyaA[:_]?S352T$", pilot_amr$gene) &
        pilot_amr$localization == "chromosome"
    )
  )
  if (!all(pilot_checks$pass)) {
    failed <- pilot_checks$check[!pilot_checks$pass]
    stop(
      "Plasmid localization pilot fixture failed: ",
      paste(failed, collapse = ", "), call. = FALSE
    )
  }

  validation <- dplyr::bind_rows(
    pilot_checks %>%
      dplyr::transmute(
        check, expected, observed, pass,
        detail = "integration fixture 2510C119001-1"
      ),
    tibble::tibble(
      check = c(
        "episode_profiles", "adjacent_pairs", "focused_transitions",
        "all_calls_explicitly_localized", "prokka_original_map_one_to_one"
      ),
      expected = c("532", "371", "9", "all", "all selected contigs"),
      observed = as.character(c(
        nrow(episode_profiles), nrow(adjacent), nrow(focused),
        sum(!is.na(locations$localization)),
        nrow(prokka_map)
      )),
      pass = c(
        nrow(episode_profiles) == 532L,
        nrow(adjacent) == 371L,
        nrow(focused) == 9L,
        all(!is.na(locations$localization)),
        all(!is.na(prokka_map$original_contig_id))
      ),
      detail = c(
        "exact selected cohort", "exact adjacent pairs",
        "descriptive Not_UTI-to-UTI cases",
        "mapped or explicit ambiguous/unassigned",
        "sequence SHA-256 + length + deterministic duplicate occurrence"
      )
    )
  )
  plasmid_atomic_write_csv(
    validation,
    file.path(output_root, "plasmid_gene_location_validation.csv")
  )
  if (!all(validation$pass)) {
    stop("Predicted-plasmid localization validation failed.")
  }
  list(
    locations = locations, episodes = episode_profiles,
    adjacent = adjacent, focused = focused,
    linkage = linkage, validation = validation
  )
}
