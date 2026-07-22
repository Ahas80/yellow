#!/usr/bin/env Rscript

source("00_config.R")
source("R/wgs_helpers.R")
source("R/amr_helpers.R")

suppressPackageStartupMessages({
  library(dplyr)
  library(readr)
  library(tibble)
})

checks <- list()
add_check <- function(check, pass, observed, required) {
  checks[[length(checks) + 1L]] <<- tibble(
    check = check,
    pass = isTRUE(pass),
    observed = as.character(observed),
    required = as.character(required)
  )
}

required_packages <- c(
  "Biostrings", "ComplexHeatmap", "ComplexUpset", "RColorBrewer", "ape",
  "broom", "broom.mixed", "digest", "dplyr", "forcats", "fs", "furrr",
  "future", "future.apply", "ggalluvial", "ggplot2", "ggraph", "ggrepel",
  "glue", "gridExtra", "igraph", "jsonlite", "lme4", "lubridate",
  "optparse", "patchwork", "pheatmap", "processx", "purrr", "ragg",
  "randomForest", "readr", "scales", "seqinr", "stringr", "tibble",
  "tidyr", "tidyverse", "writexl"
)
package_ok <- vapply(required_packages, requireNamespace, logical(1), quietly = TRUE)
add_check(
  "R packages",
  all(package_ok),
  if (all(package_ok)) length(package_ok) else paste(names(package_ok)[!package_ok], collapse = ";"),
  paste0(length(required_packages), " installed")
)

assembly_manifest <- tryCatch(
  load_analysis_assemblies(
    FILE_ANALYSIS_ASSEMBLY_MANIFEST, require_files = TRUE
  ),
  error = function(e) {
    attr(e, "manifest_error") <- TRUE
    e
  }
)
if (inherits(assembly_manifest, "error")) {
  add_check(
    "selected Longcycler FASTA manifest", FALSE,
    conditionMessage(assembly_manifest),
    "532 unique existing Longcycler FASTAs with matching SHA-256"
  )
} else {
  assembly_manifest <- assembly_manifest %>%
    mutate(
      Isolate_ID = as.character(Isolate_ID),
      full_path = normalizePath(full_path, winslash = "/", mustWork = TRUE),
      fasta_sha256 = tolower(as.character(fasta_sha256))
    )
  observed_fasta_hashes <- vapply(
    assembly_manifest$full_path,
    digest::digest,
    character(1),
    algo = "sha256",
    file = TRUE,
    serialize = FALSE
  )
  manifest_ok <-
    nrow(assembly_manifest) == 532L &&
    !anyDuplicated(assembly_manifest$Isolate_ID) &&
    !anyDuplicated(assembly_manifest$full_path) &&
    all(tolower(assembly_manifest$assembler) == "longcycler") &&
    all(nzchar(assembly_manifest$fasta_sha256)) &&
    identical(
      unname(tolower(observed_fasta_hashes)),
      unname(assembly_manifest$fasta_sha256)
    )
  add_check(
    "selected Longcycler FASTA manifest", manifest_ok,
    paste0(
      nrow(assembly_manifest), " rows; ",
      sum(tolower(observed_fasta_hashes) ==
            assembly_manifest$fasta_sha256),
      " matching content hashes"
    ),
    "532 unique existing Longcycler FASTAs with matching SHA-256"
  )
}

amr_prefix <- Sys.getenv(
  "AMR_RUNTIME_PREFIX", file.path(DIR_ROOT, "data", "amr_runtime", "env")
)
mob_prefix <- Sys.getenv(
  "MOB_RUNTIME_PREFIX",
  file.path(DIR_ROOT, "data", "mob_suite_runtime", "env")
)
required_tools <- c(
  parsnp = unname(Sys.which("parsnp")),
  harvesttools = unname(Sys.which("harvesttools")),
  `snp-dists` = unname(Sys.which("snp-dists")),
  panaroo = unname(find_prokka_bin("panaroo")),
  prokka = unname(find_prokka_bin("prokka")),
  abricate = unname(Sys.which("abricate")),
  amrfinder = Sys.getenv(
    "AMRFINDER_BIN", file.path(amr_prefix, "bin", "amrfinder")
  ),
  mob_recon = file.path(mob_prefix, "bin", "mob_recon"),
  mob_blastn = file.path(mob_prefix, "bin", "blastn"),
  mob_mash = file.path(mob_prefix, "bin", "mash"),
  resfinder_python = Sys.getenv(
    "RESFINDER_PYTHON", file.path(amr_prefix, "bin", "python")
  ),
  kma = Sys.getenv("KMA_BIN", file.path(amr_prefix, "bin", "kma")),
  mlst = unname(Sys.which("mlst")),
  dnadiff = unname(Sys.which("dnadiff")),
  nucmer = unname(Sys.which("nucmer")),
  `show-snps` = unname(Sys.which("show-snps")),
  mash = unname(Sys.which("mash")),
  java = unname(Sys.which("java"))
)
tool_ok <- nzchar(required_tools) & file.exists(required_tools)
add_check(
  "executables",
  all(tool_ok),
  if (all(tool_ok)) paste(names(required_tools), collapse = ";") else paste(names(required_tools)[!tool_ok], collapse = ";"),
  "all required tools discoverable"
)

mob_version <- if (file.exists(required_tools[["mob_recon"]])) {
  paste(
    system2(required_tools[["mob_recon"]], "--version",
            stdout = TRUE, stderr = TRUE),
    collapse = " "
  )
} else ""
mob_blast_version <- if (file.exists(required_tools[["mob_blastn"]])) {
  paste(
    system2(required_tools[["mob_blastn"]], "-version",
            stdout = TRUE, stderr = TRUE),
    collapse = " "
  )
} else ""
mob_mash_version <- if (file.exists(required_tools[["mob_mash"]])) {
  paste(
    system2(required_tools[["mob_mash"]], "--version",
            stdout = TRUE, stderr = TRUE),
    collapse = " "
  )
} else ""
add_check("MOB-suite version", grepl("3\\.1\\.9", mob_version),
          mob_version, "3.1.9")
add_check("MOB runtime BLAST version", grepl("2\\.15\\.0", mob_blast_version),
          mob_blast_version, "2.15.0")
add_check("MOB runtime Mash version", identical(trimws(mob_mash_version), "2.3"),
          mob_mash_version, "2.3")

mob_db_dir <- if (file.exists(file.path(mob_prefix, "bin", "python"))) {
  out <- system2(
    file.path(mob_prefix, "bin", "python"),
    c(
      "-c",
      shQuote(
        "from pathlib import Path; import mob_suite; print(Path(mob_suite.__file__).resolve().parent / 'databases')"
      )
    ),
    stdout = TRUE, stderr = TRUE
  )
  if (length(out) == 1L && dir.exists(out)) out else ""
} else ""
mob_expected_hashes <- c(
  "ncbi_plasmid_full_seqs.fas" = "939ca451e97f9a5e6ea3e38286d5d11c620aafdefe693705c2611f4ed98355d5",
  "clusters.txt" = "2811a3f1af5a632d985f6e64b6c3be46ec249acc4bf8f76b0bc3bdf569937e11",
  "rep.dna.fas" = "fdc10d866d8fdc3c75db0f945c8b52797392abb7b7bf389c5de36a2429500eb0",
  "mob.proteins.faa" = "e5db0e07f9e94f7252f8adceaba4355851331b1ea4e60f51ce90af97490c47a9",
  "mpf.proteins.faa" = "fc6a9f78465271826120659e46c6876aa1b0f17baa0b806a525104b2e41f12fa",
  "orit.fas" = "821b2d39ff960f9c05dd467a3f9694d108b2e45c368923d21dfbf9a7383b921d",
  "repetitive.dna.fas" = "bf039b8cf672ffed353ea35d8088e9567ea9d16191003bb10a5ec5f2157d3bf7"
)
mob_db_paths <- file.path(mob_db_dir, names(mob_expected_hashes))
mob_db_exists <- nzchar(mob_db_dir) && all(file.exists(mob_db_paths))
mob_observed_hashes <- if (mob_db_exists) {
  vapply(
    mob_db_paths, digest::digest, character(1),
    algo = "sha256", file = TRUE, serialize = FALSE
  )
} else {
  rep(NA_character_, length(mob_expected_hashes))
}
add_check(
  "MOB-suite pinned database snapshot",
  mob_db_exists && identical(unname(mob_observed_hashes), unname(mob_expected_hashes)),
  if (mob_db_exists) {
    paste(names(mob_expected_hashes)[mob_observed_hashes == mob_expected_hashes],
          collapse = ";")
  } else "missing",
  "seven core files from Zenodo record 10304948 with exact SHA-256"
)

abricate_list <- tryCatch(system2(required_tools[["abricate"]], "--list", stdout = TRUE, stderr = TRUE),
                         error = function(e) character())
required_databases <- c("vfdb", "plasmidfinder", "resfinder")
database_names <- tolower(sub("[[:space:]].*$", "", trimws(abricate_list)))
database_ok <- required_databases %in% database_names
add_check(
  "ABRicate databases",
  all(database_ok),
  paste(required_databases[database_ok], collapse = ";"),
  paste(required_databases, collapse = ";")
)

amr_database_root <- Sys.getenv(
  "AMR_DATABASE_ROOT", file.path(DIR_ROOT, "data", "amr_runtime", "databases")
)
amrfinder_db_version <- Sys.getenv(
  "AMRFINDER_DB_VERSION", AMR_AMRFINDER_DB_VERSION
)
amrfinder_db <- Sys.getenv(
  "AMRFINDER_DB",
  file.path(
    Sys.getenv(
      "AMRFINDER_DB_ROOT", file.path(amr_database_root, "amrfinderplus")
    ),
    amrfinder_db_version
  )
)
amr_database_dirs <- file.path(
  amr_database_root, c("resfinder_db", "pointfinder_db")
)
add_check(
  "ResFinder databases",
  all(dir.exists(amr_database_dirs)),
  paste(basename(amr_database_dirs[dir.exists(amr_database_dirs)]), collapse = ";"),
  "resfinder_db and pointfinder_db"
)
add_check(
  "AMRFinderPlus project-local database",
  dir.exists(amrfinder_db) &&
    file.exists(file.path(amrfinder_db, "database_format_version.txt")) &&
    startsWith(
      normalizePath(amrfinder_db, winslash = "/", mustWork = FALSE),
      normalizePath(
        file.path(DIR_ROOT, "data", "amr_runtime", "databases"),
        winslash = "/", mustWork = TRUE
      )
    ),
  amrfinder_db,
  paste0("project-local database version ", AMR_AMRFINDER_DB_VERSION)
)
amrfinder_version <- if (file.exists(required_tools[["amrfinder"]])) {
  amr_command_version(required_tools[["amrfinder"]], "--version")
} else ""
resfinder_version <- if (file.exists(required_tools[["resfinder_python"]])) {
  amr_command_version(
    required_tools[["resfinder_python"]],
    c("-c", "import importlib.metadata as m; print(m.version('resfinder'))")
  )
} else ""
add_check(
  "AMRFinderPlus version", grepl("4\\.2\\.7", amrfinder_version),
  amrfinder_version, "4.2.7"
)
add_check(
  "ResFinder version", grepl("4\\.7\\.2", resfinder_version),
  resfinder_version, "4.7.2"
)
amrfinder_db_report <- if (file.exists(required_tools[["amrfinder"]]) &&
                            dir.exists(amrfinder_db)) {
  amr_command_version(
    required_tools[["amrfinder"]], c("-d", amrfinder_db, "-V")
  )
} else ""
add_check(
  "AMRFinderPlus database version",
  grepl(AMR_AMRFINDER_DB_VERSION, amrfinder_db_report, fixed = TRUE),
  amrfinder_db_report, AMR_AMRFINDER_DB_VERSION
)
resfinder_commit <- if (dir.exists(amr_database_dirs[[1L]])) {
  amr_git_head(amr_database_dirs[[1L]])
} else ""
pointfinder_commit <- if (dir.exists(amr_database_dirs[[2L]])) {
  amr_git_head(amr_database_dirs[[2L]])
} else ""
add_check(
  "ResFinder database commit",
  identical(resfinder_commit, AMR_RESFINDER_DB_COMMIT),
  resfinder_commit, AMR_RESFINDER_DB_COMMIT
)
add_check(
  "PointFinder database commit",
  identical(pointfinder_commit, AMR_POINTFINDER_DB_COMMIT),
  pointfinder_commit, AMR_POINTFINDER_DB_COMMIT
)

df_out <- system2("df", c("-Pk", "."), stdout = TRUE)
df_fields <- strsplit(trimws(tail(df_out, 1L)), "[[:space:]]+")[[1L]]
free_gib <- suppressWarnings(as.numeric(df_fields[[4L]]) / 1024^2)
min_free_gib <- suppressWarnings(as.numeric(Sys.getenv("MIN_FREE_GIB", "60")))
if (!is.finite(min_free_gib)) min_free_gib <- 60
add_check("free disk", is.finite(free_gib) && free_gib >= min_free_gib,
          sprintf("%.1f GiB", free_gib), sprintf(">= %.1f GiB", min_free_gib))

add_check("analysis assembler", identical(tolower(ANALYSIS_ASSEMBLER), "longcycler"),
          ANALYSIS_ASSEMBLER, "longcycler")
add_check("assembler fallback", identical(ALLOW_ASSEMBLER_FALLBACK, FALSE),
          ALLOW_ASSEMBLER_FALLBACK, "FALSE")
add_check("GFF matching", Sys.getenv("GFF_ALLOW_SUBSTRING_MATCH", "0") != "1",
          Sys.getenv("GFF_ALLOW_SUBSTRING_MATCH", "0"), "not 1")
add_check("legacy exploratory plots", Sys.getenv("RUN_LEGACY_EXPLORATORY_PLOTS", "0") != "1",
          Sys.getenv("RUN_LEGACY_EXPLORATORY_PLOTS", "0"), "not 1")

rq_scripts <- file.path("scripts", "research_questions", c(
  "run_all.R", "run_rq01_05.R", "run_rq06_08.R", "run_rq09_10.R"
))
add_check("RQ01-RQ10 scripts", all(file.exists(rq_scripts)),
          paste(basename(rq_scripts[file.exists(rq_scripts)]), collapse = ";"),
          paste(basename(rq_scripts), collapse = ";"))

numbered_plasmid_scripts <- c(
  "08_core_vs_plasmid.R", "09_inc_plasmid_network.R",
  "09b_mob_plasmid_reconstruction.R", "10_replicon_heatmap.R"
)
add_check(
  "numbered plasmid scripts",
  all(file.exists(numbered_plasmid_scripts)),
  paste(numbered_plasmid_scripts[file.exists(numbered_plasmid_scripts)],
        collapse = ";"),
  paste(numbered_plasmid_scripts, collapse = ";")
)

if (identical(Sys.getenv("REQUIRE_COMPLETED_AMR", "0"), "1")) {
  amr_marker <- file.path(DIR_RESULTS, "amr", "RUN_COMPLETE.txt")
  amr_marker_text <- if (file.exists(amr_marker)) {
    paste(readLines(amr_marker, warn = FALSE), collapse = "\n")
  } else ""
  required_amr_marker_lines <- c(
    "status=complete", "episodes=532", "residents=161",
    "adjacent_pairs=371", "focused_transitions=9"
  )
  add_check(
    "completed genomic-AMR gate",
    file.exists(amr_marker) &&
      all(vapply(
        required_amr_marker_lines,
        grepl,
        logical(1),
        x = amr_marker_text,
        fixed = TRUE
      )),
    if (file.exists(amr_marker)) amr_marker_text else "missing",
    "completed 532-episode Script-29 marker"
  )
}

requested_boot <- Sys.getenv("RQ_BOOTSTRAP_REPS", "")
requested_perm <- Sys.getenv("RQ_PERMUTATIONS", "")
requested_amr_boot <- Sys.getenv("AMR_BOOTSTRAP_REPS", "")
add_check("RQ bootstrap contract", !nzchar(requested_boot) || requested_boot == "10000",
          ifelse(nzchar(requested_boot), requested_boot, "unset"), "unset or 10000")
add_check("RQ permutation contract", !nzchar(requested_perm) || requested_perm == "10000",
          ifelse(nzchar(requested_perm), requested_perm, "unset"), "unset or 10000")
add_check("AMR bootstrap contract", !nzchar(requested_amr_boot) || requested_amr_boot == "10000",
          ifelse(nzchar(requested_amr_boot), requested_amr_boot, "unset"), "unset or 10000")

result <- bind_rows(checks)
dir.create(file.path("results", "pipeline"), recursive = TRUE, showWarnings = FALSE)
write_csv(result, file.path("results", "pipeline", "preflight_checks.csv"))

failed <- result %>% filter(!pass)
if (nrow(failed)) {
  print(failed)
  stop("Complete-analysis preflight failed: ", paste(failed$check, collapse = ", "), call. = FALSE)
}
message("Complete-analysis preflight passed: ", nrow(result), " checks.")
