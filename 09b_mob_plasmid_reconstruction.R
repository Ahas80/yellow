#!/usr/bin/env Rscript
# ==============================================================================
# 09b_mob_plasmid_reconstruction.R
# ==============================================================================
#
# Assembly-first MOB-suite reconstruction for all 532 selected Longcycler
# assemblies. Outputs are predicted plasmid bins and typing context; they do
# not establish circularity, transfer, transmission, or causal mechanism.
# ==============================================================================

source("00_config.R")
suppressPackageStartupMessages({
  library(Biostrings)
  library(dplyr)
  library(readr)
  library(tidyr)
  library(purrr)
  library(furrr)
  library(stringr)
  library(tibble)
})

SCHEMA_VERSION <- "mob_recon_longcycler_v1"
EXPECTED_EPISODES <- 532L
MOB_SUITE_VERSION <- "3.1.9"
BLAST_VERSION <- "2.15.0"
MASH_VERSION <- "2.3"
ZENODO_RECORD <- "10304948"
EXPECTED_DB_HASHES <- c(
  "ncbi_plasmid_full_seqs.fas" = "939ca451e97f9a5e6ea3e38286d5d11c620aafdefe693705c2611f4ed98355d5",
  "clusters.txt" = "2811a3f1af5a632d985f6e64b6c3be46ec249acc4bf8f76b0bc3bdf569937e11",
  "rep.dna.fas" = "fdc10d866d8fdc3c75db0f945c8b52797392abb7b7bf389c5de36a2429500eb0",
  "mob.proteins.faa" = "e5db0e07f9e94f7252f8adceaba4355851331b1ea4e60f51ce90af97490c47a9",
  "mpf.proteins.faa" = "fc6a9f78465271826120659e46c6876aa1b0f17baa0b806a525104b2e41f12fa",
  "orit.fas" = "821b2d39ff960f9c05dd467a3f9694d108b2e45c368923d21dfbf9a7383b921d",
  "repetitive.dna.fas" = "bf039b8cf672ffed353ea35d8088e9567ea9d16191003bb10a5ec5f2157d3bf7"
)

MOB_ROOT <- file.path(DIR_PLASMIDS, "mob_suite")
CACHE_ROOT <- file.path(MOB_ROOT, "cache")
dir.create(MOB_ROOT, recursive = TRUE, showWarnings = FALSE)
dir.create(CACHE_ROOT, recursive = TRUE, showWarnings = FALSE)

runtime_prefix <- Sys.getenv(
  "MOB_RUNTIME_PREFIX",
  file.path(DIR_ROOT, "data", "mob_suite_runtime", "env")
)
mob_recon_bin <- file.path(runtime_prefix, "bin", "mob_recon")
python_bin <- file.path(runtime_prefix, "bin", "python")
blast_bin <- file.path(runtime_prefix, "bin", "blastn")
mash_bin <- file.path(runtime_prefix, "bin", "mash")
required_bins <- c(mob_recon_bin, python_bin, blast_bin, mash_bin)
if (!all(file.exists(required_bins))) {
  stop(
    "Pinned MOB-suite runtime is missing. Run scripts/setup_mob_suite_runtime.sh first.",
    call. = FALSE
  )
}

command_text <- function(command, args) {
  res <- processx::run(command, args, error_on_status = FALSE, stderr_to_stdout = TRUE)
  if (res$status != 0L) stop("Version command failed: ", command, call. = FALSE)
  trimws(res$stdout)
}
mob_version <- command_text(mob_recon_bin, "--version")
blast_version <- command_text(blast_bin, "-version")
mash_version <- command_text(mash_bin, "--version")
if (!str_detect(mob_version, fixed(MOB_SUITE_VERSION)) ||
    !str_detect(blast_version, fixed(BLAST_VERSION)) ||
    !identical(mash_version, MASH_VERSION)) {
  stop("Pinned MOB-suite/BLAST/Mash version contract failed.", call. = FALSE)
}

db_dir <- command_text(
  python_bin,
  c(
    "-c",
    "from pathlib import Path; import mob_suite; print(Path(mob_suite.__file__).resolve().parent / 'databases')"
  )
)
if (!dir.exists(db_dir)) stop("MOB-suite database directory is missing.", call. = FALSE)

sha256_file <- function(path) {
  digest::digest(path, algo = "sha256", file = TRUE, serialize = FALSE)
}
archive_superseded_mob_outputs <- function() {
  marker <- file.path(MOB_ROOT, "RUN_COMPLETE.txt")
  marker_text <- if (file.exists(marker)) {
    paste(readLines(marker, warn = FALSE), collapse = "\n")
  } else {
    ""
  }
  if (grepl(
    paste0("schema=", SCHEMA_VERSION), marker_text, fixed = TRUE
  )) {
    return(invisible(NULL))
  }
  names_to_archive <- c(
    "RUN_COMPLETE.txt", "database_manifest.csv", "run_manifest.csv",
    "sample_status.csv", "contig_assignments.csv", "plasmids_long.csv",
    "episode_plasmid_profiles.csv", "prokka_original_contig_map.csv",
    "plasmid_gene_locations_long.csv",
    "predicted_plasmid_linkage_groups.csv",
    "episode_mechanism_profiles.csv",
    "adjacent_pair_plasmid_metrics_371.csv",
    "not_uti_to_uti_plasmid_metrics_9.csv",
    "plasmid_gene_location_validation.csv"
  )
  candidates <- file.path(MOB_ROOT, names_to_archive)
  candidates <- candidates[file.exists(candidates)]
  if (!length(candidates)) return(invisible(NULL))
  stamp <- format(Sys.time(), "%Y%m%dT%H%M%SZ", tz = "UTC")
  archive_dir <- file.path(
    MOB_ROOT, "archive",
    paste0("superseded_pre_", SCHEMA_VERSION, "_", stamp)
  )
  dir.create(archive_dir, recursive = TRUE, showWarnings = FALSE)
  destinations <- file.path(archive_dir, basename(candidates))
  moved <- mapply(file.rename, candidates, destinations)
  if (!all(moved)) {
    stop(
      "Could not archive every superseded MOB output: ",
      paste(basename(candidates[!moved]), collapse = ", "),
      call. = FALSE
    )
  }
  writeLines(
    c(
      paste0("Superseded generated MOB outputs archived at ", stamp),
      paste0("Replacement schema: ", SCHEMA_VERSION),
      "Selected FASTAs, MOB databases, resumable caches and unrelated user files were not moved.",
      paste0("Archived: ", paste(basename(candidates), collapse = ";"))
    ),
    file.path(archive_dir, "ARCHIVE_MANIFEST.txt")
  )
  invisible(archive_dir)
}
observed_db_hashes <- vapply(
  file.path(db_dir, names(EXPECTED_DB_HASHES)),
  sha256_file,
  character(1)
)
names(observed_db_hashes) <- names(EXPECTED_DB_HASHES)
if (!identical(observed_db_hashes, EXPECTED_DB_HASHES)) {
  bad <- names(EXPECTED_DB_HASHES)[observed_db_hashes != EXPECTED_DB_HASHES]
  stop("Pinned MOB-suite database hash mismatch: ", paste(bad, collapse = ", "), call. = FALSE)
}

database_manifest <- tibble(
  schema_version = SCHEMA_VERSION,
  mob_suite_version = mob_version,
  blast_version = str_split(blast_version, "\n", simplify = TRUE)[1L],
  mash_version = mash_version,
  zenodo_record = ZENODO_RECORD,
  database_file = names(EXPECTED_DB_HASHES),
  database_sha256 = unname(EXPECTED_DB_HASHES),
  database_path = normalizePath(
    file.path(db_dir, names(EXPECTED_DB_HASHES)),
    winslash = "/", mustWork = TRUE
  )
)

manifest <- load_analysis_assemblies(
  FILE_ANALYSIS_ASSEMBLY_MANIFEST, require_files = TRUE
) %>%
  mutate(
    Isolate_ID = as.character(Isolate_ID),
    Assembly_ID = as.character(Assembly_ID),
    Participant_id = as.character(Participant_id),
    tp_lab = as.character(tp_lab),
    fasta_path = normalizePath(full_path, winslash = "/", mustWork = TRUE),
    fasta_sha256 = vapply(fasta_path, sha256_file, character(1))
  ) %>%
  distinct(Isolate_ID, .keep_all = TRUE)
if (nrow(manifest) != EXPECTED_EPISODES ||
    anyDuplicated(manifest$Isolate_ID) ||
    anyDuplicated(manifest$fasta_path)) {
  stop("MOB reconstruction requires the exact 532 selected Longcycler assemblies.", call. = FALSE)
}
archive_superseded_mob_outputs()
if (file.exists(file.path(MOB_ROOT, "RUN_COMPLETE.txt"))) {
  unlink(file.path(MOB_ROOT, "RUN_COMPLETE.txt"))
}
write_csv(database_manifest, file.path(MOB_ROOT, "database_manifest.csv"))

run_overhang <- identical(Sys.getenv("MOB_RUN_OVERHANG", "1"), "1")
threads <- suppressWarnings(as.integer(Sys.getenv("MOB_THREADS_PER_SAMPLE", "1")))
workers <- suppressWarnings(as.integer(Sys.getenv("MOB_WORKERS", "4")))
if (!is.finite(threads) || threads < 1L) threads <- 1L
if (!is.finite(workers) || workers < 1L) workers <- 1L
workers <- min(workers, nrow(manifest))

db_signature <- paste(names(EXPECTED_DB_HASHES), EXPECTED_DB_HASHES, collapse = "\n")
command_signature <- paste(
  SCHEMA_VERSION,
  mob_version,
  blast_version,
  mash_version,
  db_signature,
  paste0("threads=", threads),
  paste0("run_overhang=", run_overhang),
  "min_length=1000",
  "min_rep_ident=80",
  "min_rep_cov=80",
  "min_con_ident=80",
  "min_con_cov=60",
  "max_plasmid_size=450000",
  sep = "\n"
)

read_input_contigs <- function(fasta) {
  seqs <- Biostrings::readDNAStringSet(fasta)
  ids <- sub("[[:space:]].*$", "", names(seqs))
  if (any(!nzchar(ids)) || anyDuplicated(ids)) {
    stop("Input FASTA has empty or duplicate contig identifiers: ", fasta, call. = FALSE)
  }
  gc <- rowSums(Biostrings::letterFrequency(
    seqs, letters = c("G", "C"), as.prob = TRUE
  ))
  tibble(
    contig_id = ids,
    input_size = width(seqs),
    input_gc = as.numeric(gc)
  )
}

validate_cached_result <- function(result_dir, fasta) {
  required <- file.path(result_dir, c("contig_report.txt", "RUN_COMPLETE.txt"))
  if (!all(file.exists(required))) return(FALSE)
  contigs <- read_input_contigs(fasta)
  report <- tryCatch(
    read_tsv(file.path(result_dir, "contig_report.txt"), show_col_types = FALSE),
    error = function(e) NULL
  )
  if (is.null(report) || !"contig_id" %in% names(report) || anyDuplicated(report$contig_id)) {
    return(FALSE)
  }
  all(report$contig_id %in% contigs$contig_id)
}

run_one <- function(isolate_id, fasta, fasta_sha256) {
  cache_key <- digest::digest(
    paste(fasta_sha256, command_signature, sep = "\n"),
    algo = "sha256", serialize = FALSE
  )
  safe_id <- str_replace_all(isolate_id, "[^A-Za-z0-9_.-]", "_")
  result_dir <- file.path(CACHE_ROOT, paste0(safe_id, "__", substr(cache_key, 1L, 24L)))
  reused <- validate_cached_result(result_dir, fasta)
  if (!reused) {
    tmp_root <- tempfile(pattern = paste0(".", safe_id, "."), tmpdir = CACHE_ROOT)
    dir.create(tmp_root, recursive = TRUE, showWarnings = FALSE)
    on.exit(unlink(tmp_root, recursive = TRUE, force = TRUE), add = TRUE)
    out_dir <- file.path(tmp_root, "recon")
    args <- c(
      "-i", fasta,
      "-o", out_dir,
      "-s", isolate_id,
      "-n", as.character(threads),
      "-f"
    )
    if (run_overhang) args <- c(args, "-c")
    env <- c(PATH = paste(file.path(runtime_prefix, "bin"), Sys.getenv("PATH"), sep = .Platform$path.sep))
    res <- processx::run(
      mob_recon_bin, args,
      env = env,
      error_on_status = FALSE,
      stderr_to_stdout = FALSE
    )
    writeLines(res$stdout, file.path(tmp_root, "mob_recon.stdout.log"))
    writeLines(res$stderr, file.path(tmp_root, "mob_recon.stderr.log"))
    if (res$status != 0L ||
        !file.exists(file.path(out_dir, "contig_report.txt"))) {
      stop(
        "MOB-recon failed for ", isolate_id, " (exit ", res$status, "): ",
        trimws(tail(str_split(res$stderr, "\n")[[1L]], 5L)),
        call. = FALSE
      )
    }
    if (!file.exists(file.path(out_dir, "mobtyper_results.txt"))) {
      writeLines(
        c(
          "MOB-recon completed successfully without a plasmid biomarker report.",
          "This is a valid zero-predicted-plasmid profile; contig_report.txt remains authoritative for contig accounting."
        ),
        file.path(out_dir, "NO_PREDICTED_PLASMIDS.txt")
      )
    }
    writeLines(
      c(
        "MOB-recon sample: PASS",
        paste0("schema=", SCHEMA_VERSION),
        paste0("cache_key=", cache_key),
        paste0("Isolate_ID=", isolate_id),
        paste0("fasta_path=", normalizePath(fasta, winslash = "/", mustWork = TRUE)),
        paste0("fasta_sha256=", fasta_sha256),
        paste0("mob_suite_version=", mob_version),
        paste0("database_signature_sha256=", digest::digest(db_signature, "sha256", serialize = FALSE))
      ),
      file.path(out_dir, "RUN_COMPLETE.txt")
    )
    file.copy(file.path(tmp_root, "mob_recon.stdout.log"), out_dir, overwrite = TRUE)
    file.copy(file.path(tmp_root, "mob_recon.stderr.log"), out_dir, overwrite = TRUE)
    if (dir.exists(result_dir)) {
      stop("Refusing to overwrite an invalid versioned MOB cache: ", result_dir, call. = FALSE)
    }
    if (!file.rename(out_dir, result_dir)) {
      stop("Could not atomically publish MOB result for ", isolate_id, call. = FALSE)
    }
    if (!validate_cached_result(result_dir, fasta)) {
      stop("Published MOB result failed validation for ", isolate_id, call. = FALSE)
    }
  }
  tibble(
    Isolate_ID = isolate_id,
    fasta_path = normalizePath(fasta, winslash = "/", mustWork = TRUE),
    fasta_sha256 = fasta_sha256,
    cache_key = cache_key,
    result_dir = normalizePath(result_dir, winslash = "/", mustWork = TRUE),
    cache_reused = reused,
    mobtyper_report_present = file.exists(
      file.path(result_dir, "mobtyper_results.txt")
    ),
    zero_predicted_plasmid_profile = !file.exists(
      file.path(result_dir, "mobtyper_results.txt")
    ),
    status = "complete"
  )
}

future::plan(future::multisession, workers = workers)
sample_status <- future_pmap_dfr(
  list(manifest$Isolate_ID, manifest$fasta_path, manifest$fasta_sha256),
  run_one,
  .progress = interactive(),
  .options = furrr::furrr_options(seed = TRUE)
)
future::plan(future::sequential)

if (nrow(sample_status) != EXPECTED_EPISODES ||
    any(sample_status$status != "complete") ||
    !setequal(sample_status$Isolate_ID, manifest$Isolate_ID)) {
  stop("MOB-suite did not complete the exact 532-episode cohort.", call. = FALSE)
}

parse_one_result <- function(i) {
  status_row <- sample_status[i, ]
  metadata <- manifest[match(status_row$Isolate_ID, manifest$Isolate_ID), ]
  input_contigs <- read_input_contigs(status_row$fasta_path)
  report <- read_tsv(
    file.path(status_row$result_dir, "contig_report.txt"),
    show_col_types = FALSE,
    col_types = cols(.default = "c")
  )
  if (anyDuplicated(report$contig_id) || any(!report$contig_id %in% input_contigs$contig_id)) {
    stop("MOB contig report violates input-contig identity for ", status_row$Isolate_ID, call. = FALSE)
  }
  report <- input_contigs %>%
    left_join(report, by = "contig_id") %>%
    mutate(
      Isolate_ID = status_row$Isolate_ID,
      Assembly_ID = metadata$Assembly_ID,
      Participant_id = metadata$Participant_id,
      tp_lab = metadata$tp_lab,
      fasta_path = status_row$fasta_path,
      fasta_sha256 = status_row$fasta_sha256,
      cache_key = status_row$cache_key,
      molecule_type = case_when(
        molecule_type == "plasmid" ~ "predicted_plasmid",
        molecule_type == "chromosome" ~ "predicted_chromosome",
        TRUE ~ "unassigned"
      ),
      filtering_reason = case_when(
        !is.na(filtering_reason) & nzchar(filtering_reason) ~ filtering_reason,
        molecule_type == "unassigned" & input_size < 1000L ~ "below_mob_min_length",
        molecule_type == "unassigned" ~ "not_classified_by_mob_recon",
        TRUE ~ NA_character_
      )
    ) %>%
    relocate(
      Isolate_ID, Assembly_ID, Participant_id, tp_lab,
      fasta_path, fasta_sha256, cache_key, contig_id,
      input_size, input_gc, molecule_type
    )
  if (nrow(report) != nrow(input_contigs) ||
      anyDuplicated(report$contig_id) ||
      any(is.na(report$molecule_type))) {
    stop("Not every input contig received exactly one explicit MOB assignment.", call. = FALSE)
  }

  typer_path <- file.path(status_row$result_dir, "mobtyper_results.txt")
  typer <- if (file.exists(typer_path) && file.size(typer_path) > 0L) {
    read_tsv(
      typer_path,
      show_col_types = FALSE,
      col_types = cols(.default = "c")
    )
  } else {
    tibble()
  }
  if (nrow(typer)) {
    typer <- typer %>%
      mutate(
        Isolate_ID = status_row$Isolate_ID,
        Assembly_ID = metadata$Assembly_ID,
        Participant_id = metadata$Participant_id,
        tp_lab = metadata$tp_lab,
        fasta_path = status_row$fasta_path,
        fasta_sha256 = status_row$fasta_sha256,
        cache_key = status_row$cache_key
      ) %>%
      relocate(
        Isolate_ID, Assembly_ID, Participant_id, tp_lab,
        fasta_path, fasta_sha256, cache_key
      )
  }
  list(contigs = report, plasmids = typer)
}

parsed <- map(seq_len(nrow(sample_status)), parse_one_result)
contig_assignments <- map_dfr(parsed, "contigs")
plasmids_long <- map_dfr(parsed, "plasmids")
if (n_distinct(contig_assignments$Isolate_ID) != EXPECTED_EPISODES ||
    anyDuplicated(contig_assignments[c("Isolate_ID", "contig_id")])) {
  stop("Combined MOB contig assignment contract failed.", call. = FALSE)
}

valid_cluster <- function(x) {
  !is.na(x) & nzchar(trimws(x)) & trimws(x) != "-"
}

if (nrow(plasmids_long)) {
  plasmids_long <- plasmids_long %>%
    mutate(
      mash_neighbor_distance_numeric = suppressWarnings(
        as.numeric(mash_neighbor_distance)
      ),
      high_confidence_assignment =
        valid_cluster(primary_cluster_id) &
        is.finite(mash_neighbor_distance_numeric) &
        mash_neighbor_distance_numeric <= 0.025,
      uncertainty_reason = case_when(
        high_confidence_assignment ~ "",
        !valid_cluster(primary_cluster_id) ~ "missing_primary_cluster",
        !is.finite(mash_neighbor_distance_numeric) ~
          "missing_nearest_reference_distance",
        TRUE ~ "nearest_reference_distance_above_0.025"
      ),
      interpretation =
        "assembly-based predicted plasmid bin; not confirmed circular or transferred"
    )
}

contig_assignments <- contig_assignments %>%
  mutate(
    mash_neighbor_distance_numeric = suppressWarnings(
      as.numeric(mash_neighbor_distance)
    ),
    assignment_confidence = case_when(
      molecule_type == "unassigned" ~ "unassigned",
      molecule_type == "predicted_chromosome" ~ "mob_predicted_chromosome",
      valid_cluster(primary_cluster_id) &
        is.finite(mash_neighbor_distance_numeric) &
        mash_neighbor_distance_numeric <= 0.025 ~ "high_confidence_plasmid_bin",
      TRUE ~ "lower_confidence_plasmid_bin"
    ),
    uncertainty_reason = case_when(
      molecule_type == "unassigned" ~ coalesce(
        filtering_reason, "not_classified_by_mob_recon"
      ),
      molecule_type == "predicted_chromosome" ~
        "assembly_based_chromosome_prediction",
      assignment_confidence == "high_confidence_plasmid_bin" ~ "",
      !valid_cluster(primary_cluster_id) ~ "missing_primary_cluster",
      !is.finite(mash_neighbor_distance_numeric) ~
        "missing_nearest_reference_distance",
      TRUE ~ "nearest_reference_distance_above_0.025"
    )
  )

collapse_values <- function(x) {
  values <- unique(trimws(unlist(str_split(x[!is.na(x) & nzchar(x)], ","))))
  values <- sort(values[nzchar(values) & values != "-"])
  if (length(values)) paste(values, collapse = ";") else ""
}

plasmid_summary <- if (nrow(plasmids_long)) {
  plasmids_long %>%
    group_by(Isolate_ID) %>%
    summarise(
      predicted_plasmid_count = n(),
      predicted_plasmid_bp = sum(suppressWarnings(as.numeric(size)), na.rm = TRUE),
      high_confidence_predicted_plasmid_count =
        sum(high_confidence_assignment, na.rm = TRUE),
      mob_primary_clusters = collapse_values(primary_cluster_id),
      mob_secondary_clusters = collapse_values(secondary_cluster_id),
      mob_replicon_types = collapse_values(`rep_type(s)`),
      mob_relaxase_types = collapse_values(`relaxase_type(s)`),
      mob_mpf_types = collapse_values(mpf_type),
      mob_orit_types = collapse_values(`orit_type(s)`),
      predicted_mobility_types = collapse_values(predicted_mobility),
      .groups = "drop"
    )
} else {
  tibble(Isolate_ID = character())
}

contig_summary <- contig_assignments %>%
  group_by(Isolate_ID) %>%
  summarise(
    input_contig_count = n(),
    predicted_plasmid_contigs = sum(molecule_type == "predicted_plasmid"),
    predicted_chromosome_contigs = sum(molecule_type == "predicted_chromosome"),
    unassigned_contigs = sum(molecule_type == "unassigned"),
    .groups = "drop"
  )

episode_profiles <- manifest %>%
  select(Isolate_ID, Assembly_ID, Participant_id, tp_lab, fasta_path, fasta_sha256) %>%
  left_join(contig_summary, by = "Isolate_ID") %>%
  left_join(plasmid_summary, by = "Isolate_ID") %>%
  mutate(
    predicted_plasmid_count = replace_na(predicted_plasmid_count, 0L),
    predicted_plasmid_bp = replace_na(predicted_plasmid_bp, 0),
    high_confidence_predicted_plasmid_count = replace_na(
      high_confidence_predicted_plasmid_count, 0L
    ),
    across(
      c(
        mob_primary_clusters, mob_secondary_clusters, mob_replicon_types,
        mob_relaxase_types, mob_mpf_types, mob_orit_types, predicted_mobility_types
      ),
      ~ replace_na(.x, "")
    ),
    mob_high_confidence_profile =
      unassigned_contigs == 0L &
      predicted_plasmid_count == high_confidence_predicted_plasmid_count,
    mob_call_status = "complete",
    assembly_only_prediction = TRUE
  )

run_manifest <- sample_status %>%
  left_join(
    manifest %>% select(Isolate_ID, Assembly_ID, Participant_id, tp_lab),
    by = "Isolate_ID"
  ) %>%
  mutate(
    schema_version = SCHEMA_VERSION,
    mob_suite_version = mob_version,
    blast_version = str_split(blast_version, "\n", simplify = TRUE)[1L],
    mash_version = mash_version,
    zenodo_record = ZENODO_RECORD,
    database_signature_sha256 = digest::digest(db_signature, "sha256", serialize = FALSE),
    workers = workers,
    threads_per_sample = threads,
    run_overhang = run_overhang,
    interpretation = "assembly-based predicted plasmid bins; not transfer or transmission evidence"
  )

write_csv(run_manifest, file.path(MOB_ROOT, "run_manifest.csv"))
write_csv(sample_status, file.path(MOB_ROOT, "sample_status.csv"))
write_csv(contig_assignments, file.path(MOB_ROOT, "contig_assignments.csv"))
write_csv(plasmids_long, file.path(MOB_ROOT, "plasmids_long.csv"))
write_csv(episode_profiles, file.path(MOB_ROOT, "episode_plasmid_profiles.csv"))

writeLines(
  c(
    "MOB-suite cohort reconstruction: PASS",
    paste0("schema=", SCHEMA_VERSION),
    paste0("episodes=", nrow(episode_profiles)),
    paste0("input_contigs=", nrow(contig_assignments)),
    paste0("predicted_plasmids=", nrow(plasmids_long)),
    paste0("mob_suite_version=", mob_version),
    paste0("database_zenodo_record=", ZENODO_RECORD),
    paste0("database_signature_sha256=", digest::digest(db_signature, "sha256", serialize = FALSE)),
    "interpretation=assembly-based predicted plasmid bins; no circularity, transfer, transmission, or causal claim"
  ),
  file.path(MOB_ROOT, "RUN_COMPLETE.txt")
)
msg("✓ MOB-suite reconstruction complete for all 532 selected assemblies.")
