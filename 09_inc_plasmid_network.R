#!/usr/bin/env Rscript
# ==============================================================================
# 09_inc_plasmid_network.R
# ==============================================================================
#
# Canonical PlasmidFinder replicon screen for the exact 532-episode Longcycler
# cohort. ABRicate is run once at an 80% identity / 60% coverage floor; the
# primary 80/80 and stringent 90/90 profiles are derived from that immutable
# result. Presence/absence features are exact PlasmidFinder GENE labels.
#
# These are replicon-marker calls. They are not reconstructed plasmids and do
# not establish physical linkage, circularity, transfer, or transmission.
# ==============================================================================

source("00_config.R")
suppressPackageStartupMessages({
  library(dplyr)
  library(readr)
  library(tidyr)
  library(stringr)
  library(furrr)
  library(fs)
  library(igraph)
  library(ggraph)
  library(ggplot2)
  library(purrr)
  library(tibble)
})

SCHEMA_VERSION <- "plasmidfinder_gene_profile_v2"
SCAN_MIN_ID <- 80
SCAN_MIN_COV <- 60
PRIMARY_MIN_ID <- 80
PRIMARY_MIN_COV <- 80
STRINGENT_MIN_ID <- 90
STRINGENT_MIN_COV <- 90
EXPECTED_EPISODES <- 532L
EXPECTED_PRIMARY_HITS <- 1257L
EXPECTED_PRIMARY_GENES <- 42L
EXPECTED_POSITIVE_EPISODES <- 422L
EXPECTED_NO_HIT_EPISODES <- 110L

abricate_bin <- unname(Sys.which("abricate"))
if (!nzchar(abricate_bin) || !file.exists(abricate_bin)) {
  stop("ABRicate is required for canonical PlasmidFinder screening.", call. = FALSE)
}

ensure_dir(DIR_PLASMIDS)
ensure_dir(DIR_PLOTS_PLASMIDS)
DIR_CACHE <- file.path(DIR_PLASMIDS, "abricate_cache_gene_v2")
ensure_dir(DIR_CACHE)

sha256_file <- function(path) {
  digest::digest(path, algo = "sha256", file = TRUE, serialize = FALSE)
}

atomic_write_lines <- function(lines, path) {
  tmp <- tempfile(pattern = paste0(".", basename(path), "."), tmpdir = dirname(path))
  on.exit(unlink(tmp), add = TRUE)
  writeLines(lines, tmp, useBytes = TRUE)
  if (!file.rename(tmp, path)) stop("Could not atomically publish ", path, call. = FALSE)
  invisible(path)
}

atomic_write_csv <- function(x, path) {
  tmp <- tempfile(pattern = paste0(".", basename(path), "."), tmpdir = dirname(path), fileext = ".csv")
  on.exit(unlink(tmp), add = TRUE)
  write_csv(x, tmp, na = "")
  if (!file.rename(tmp, path)) stop("Could not atomically publish ", path, call. = FALSE)
  invisible(path)
}

archive_superseded_plasmid_outputs <- function() {
  prior_marker <- file.path(
    DIR_PLASMIDS, "PLASMIDFINDER_RUN_COMPLETE.txt"
  )
  prior_text <- if (file.exists(prior_marker)) {
    paste(readLines(prior_marker, warn = FALSE), collapse = "\n")
  } else {
    ""
  }
  if (grepl(
    paste0("schema=", SCHEMA_VERSION), prior_text, fixed = TRUE
  )) {
    return(invisible(NULL))
  }

  result_candidates <- list.files(
    DIR_PLASMIDS,
    pattern = "^(PLASMIDFINDER_RUN_COMPLETE|plasmidfinder_|replicon_|ST_replicon).*(csv|tsv|txt|rds)$",
    full.names = TRUE,
    ignore.case = FALSE
  )
  plot_candidates <- list.files(
    DIR_PLOTS_PLASMIDS,
    pattern = "^(replicon_cooccurrence|ST_vs_replicon_network|replicon_heatmap)\\.(png|pdf)$",
    full.names = TRUE,
    ignore.case = FALSE
  )
  candidates <- c(result_candidates, plot_candidates)
  candidates <- candidates[file.exists(candidates) & !dir.exists(candidates)]
  if (!length(candidates)) return(invisible(NULL))

  stamp <- format(Sys.time(), "%Y%m%dT%H%M%SZ", tz = "UTC")
  archive_root <- file.path(
    DIR_PLASMIDS, "archive", paste0("superseded_pre_", SCHEMA_VERSION, "_", stamp)
  )
  result_archive <- file.path(archive_root, "results")
  plot_archive <- file.path(archive_root, "plots")
  dir.create(result_archive, recursive = TRUE, showWarnings = FALSE)
  dir.create(plot_archive, recursive = TRUE, showWarnings = FALSE)
  destinations <- vapply(candidates, function(path) {
    target_dir <- if (startsWith(
      normalizePath(path, winslash = "/", mustWork = TRUE),
      normalizePath(DIR_PLOTS_PLASMIDS, winslash = "/", mustWork = TRUE)
    )) {
      plot_archive
    } else {
      result_archive
    }
    file.path(target_dir, basename(path))
  }, character(1))
  moved <- mapply(file.rename, candidates, destinations)
  if (!all(moved)) {
    stop(
      "Could not archive every superseded generated plasmid output: ",
      paste(basename(candidates[!moved]), collapse = ", "),
      call. = FALSE
    )
  }
  writeLines(
    c(
      paste0("Superseded generated plasmid outputs archived at ", stamp),
      paste0("Replacement schema: ", SCHEMA_VERSION),
      "Source FASTAs, provider inputs, databases and unrelated user files were not moved.",
      paste0("Archived: ", paste(basename(candidates), collapse = ";"))
    ),
    file.path(archive_root, "ARCHIVE_MANIFEST.txt")
  )
  invisible(archive_root)
}

command_version <- function(command, args = "--version") {
  res <- processx::run(command, args, error_on_status = FALSE, stderr_to_stdout = TRUE)
  trimws(paste(res$stdout, collapse = " "))
}

abricate_version <- command_version(abricate_bin)
abricate_help <- processx::run(
  abricate_bin, "--help", error_on_status = FALSE, stderr_to_stdout = TRUE
)$stdout
datadir <- Sys.getenv("ABRICATE_DATADIR", "")
if (!nzchar(datadir)) {
  datadir_line <- str_subset(
    str_split(abricate_help, "\n", simplify = FALSE)[[1L]],
    "--datadir"
  )[1L]
  datadir <- str_match(datadir_line, "\\[([^]]+)\\]")[, 2L]
}
if (is.na(datadir) || !nzchar(datadir) || !dir.exists(datadir)) {
  candidates <- Sys.glob(
    file.path(
      dirname(dirname(abricate_bin)), "Cellar", "abricate",
      "*", "libexec", "db"
    )
  )
  candidates <- candidates[dir.exists(candidates)]
  if (length(candidates)) datadir <- candidates[[1L]]
}
if (is.na(datadir) || !nzchar(datadir) || !dir.exists(datadir)) {
  stop("Could not resolve the active ABRicate database directory.", call. = FALSE)
}
datadir <- normalizePath(datadir, winslash = "/", mustWork = TRUE)
db_sequence_file <- file.path(datadir, "plasmidfinder", "sequences")
if (!file.exists(db_sequence_file)) {
  stop("PlasmidFinder database sequences are missing: ", db_sequence_file, call. = FALSE)
}

abricate_list_raw <- processx::run(
  abricate_bin, "--list", error_on_status = TRUE, stderr_to_stdout = TRUE
)$stdout
abricate_list <- read_tsv(I(abricate_list_raw), show_col_types = FALSE, progress = FALSE)
db_row <- abricate_list %>%
  filter(tolower(DATABASE) == "plasmidfinder")
if (nrow(db_row) != 1L) {
  stop("ABRicate --list did not return exactly one plasmidfinder database row.", call. = FALSE)
}

db_sha256 <- sha256_file(db_sequence_file)
database_manifest <- tibble(
  schema_version = SCHEMA_VERSION,
  caller = normalizePath(abricate_bin, winslash = "/", mustWork = TRUE),
  caller_version = abricate_version,
  database = "plasmidfinder",
  database_sequences = as.integer(db_row$SEQUENCES),
  database_type = as.character(db_row$DBTYPE),
  database_date = as.character(db_row$DATE),
  database_path = normalizePath(db_sequence_file, winslash = "/", mustWork = TRUE),
  database_sha256 = db_sha256,
  scan_min_identity = SCAN_MIN_ID,
  scan_min_coverage = SCAN_MIN_COV,
  primary_min_identity = PRIMARY_MIN_ID,
  primary_min_coverage = PRIMARY_MIN_COV,
  stringent_min_identity = STRINGENT_MIN_ID,
  stringent_min_coverage = STRINGENT_MIN_COV
)
fasta_manifest <- load_analysis_assemblies(
  FILE_ANALYSIS_ASSEMBLY_MANIFEST, require_files = TRUE
) %>%
  mutate(
    Isolate_ID = as.character(Isolate_ID),
    Assembly_ID = as.character(Assembly_ID),
    Participant_id = as.character(Participant_id),
    tp_lab = as.character(tp_lab),
    full_path = normalizePath(full_path, winslash = "/", mustWork = TRUE),
    fasta_sha256 = vapply(full_path, sha256_file, character(1))
  ) %>%
  distinct(Isolate_ID, .keep_all = TRUE)

if (nrow(fasta_manifest) != EXPECTED_EPISODES ||
    anyDuplicated(fasta_manifest$Isolate_ID) ||
    anyDuplicated(fasta_manifest$full_path)) {
  stop("PlasmidFinder input must be the exact 532-row selected Longcycler cohort.", call. = FALSE)
}
archive_superseded_plasmid_outputs()
if (file.exists(file.path(
  DIR_PLASMIDS, "PLASMIDFINDER_RUN_COMPLETE.txt"
))) {
  unlink(file.path(DIR_PLASMIDS, "PLASMIDFINDER_RUN_COMPLETE.txt"))
}
atomic_write_csv(
  database_manifest,
  file.path(DIR_PLASMIDS, "plasmidfinder_database_manifest.csv")
)
atomic_write_csv(fasta_manifest, file.path(DIR_PLASMIDS, "plasmidfinder_input_manifest.csv"))

call_signature <- paste(
  SCHEMA_VERSION,
  normalizePath(abricate_bin, winslash = "/", mustWork = TRUE),
  abricate_version,
  db_sha256,
  paste0("minid=", SCAN_MIN_ID),
  paste0("mincov=", SCAN_MIN_COV),
  "threads=1",
  sep = "\n"
)

run_abricate <- function(fasta, fasta_sha256, isolate_id) {
  fasta_norm <- normalizePath(fasta, winslash = "/", mustWork = TRUE)
  cache_key <- digest::digest(
    paste(fasta_sha256, call_signature, sep = "\n"),
    algo = "sha256", serialize = FALSE
  )
  stem <- paste0(
    str_replace_all(isolate_id, "[^A-Za-z0-9_.-]", "_"),
    "__", substr(cache_key, 1L, 24L)
  )
  out <- file.path(DIR_CACHE, paste0(stem, ".tsv"))
  sidecar <- file.path(DIR_CACHE, paste0(stem, ".provenance.csv"))

  expected <- tibble(
    schema_version = SCHEMA_VERSION,
    cache_key = cache_key,
    Isolate_ID = isolate_id,
    fasta_path = fasta_norm,
    fasta_sha256 = fasta_sha256,
    caller_path = normalizePath(abricate_bin, winslash = "/", mustWork = TRUE),
    caller_version = abricate_version,
    database = "plasmidfinder",
    database_sha256 = db_sha256,
    min_identity = as.character(SCAN_MIN_ID),
    min_coverage = as.character(SCAN_MIN_COV),
    threads = "1",
    status = "complete"
  )

  cache_ok <- FALSE
  if (file.exists(out) && file.exists(sidecar)) {
    observed <- tryCatch(
      read_csv(sidecar, show_col_types = FALSE, col_types = cols(.default = "c")),
      error = function(e) NULL
    )
    cache_ok <- !is.null(observed) && nrow(observed) == 1L &&
      identical(names(observed), names(expected)) &&
      identical(as.list(observed[1, ]), as.list(expected[1, ]))
  }

  if (!cache_ok) {
    args <- c(
      "--quiet", "--threads", "1", "--db", "plasmidfinder",
      "--minid", as.character(SCAN_MIN_ID),
      "--mincov", as.character(SCAN_MIN_COV),
      fasta_norm
    )
    res <- processx::run(
      abricate_bin, args,
      error_on_status = FALSE,
      stderr_to_stdout = FALSE
    )
    if (res$status != 0L) {
      stop(
        "PlasmidFinder failed for ", isolate_id, " (exit ", res$status, "): ",
        trimws(res$stderr), call. = FALSE
      )
    }
    atomic_write_lines(res$stdout, out)
    atomic_write_csv(expected, sidecar)
  }

  tibble(
    Isolate_ID = isolate_id,
    cache_key = cache_key,
    tab_file = normalizePath(out, winslash = "/", mustWork = TRUE),
    provenance_file = normalizePath(sidecar, winslash = "/", mustWork = TRUE),
    cache_reused = cache_ok
  )
}

workers <- suppressWarnings(as.integer(Sys.getenv("PLASMIDFINDER_WORKERS", as.character(CORES_USE))))
if (!is.finite(workers) || workers < 1L) workers <- 1L
workers <- min(workers, nrow(fasta_manifest))
future::plan(future::multisession, workers = workers)
call_files <- future_pmap_dfr(
  list(fasta_manifest$full_path, fasta_manifest$fasta_sha256, fasta_manifest$Isolate_ID),
  run_abricate,
  .progress = interactive(),
  .options = furrr::furrr_options(seed = TRUE)
)
future::plan(future::sequential)

read_one_tab <- function(tab_file, metadata) {
  lines <- readLines(tab_file, warn = FALSE)
  if (!length(lines) || !any(nzchar(trimws(lines)))) return(tibble())
  dat <- read_tsv(
    I(paste(lines, collapse = "\n")),
    col_types = cols(.default = "c"),
    show_col_types = FALSE,
    progress = FALSE
  )
  required <- c(
    "ACCESSION", "%IDENTITY", "%COVERAGE", "SEQUENCE", "GENE",
    "START", "END", "STRAND"
  )
  if (!all(required %in% names(dat))) {
    stop("ABRicate result lacks required columns: ", tab_file, call. = FALSE)
  }
  dat %>%
    transmute(
      Isolate_ID = metadata$Isolate_ID,
      Assembly_ID = metadata$Assembly_ID,
      Participant_id = metadata$Participant_id,
      tp_lab = metadata$tp_lab,
      fasta_path = metadata$full_path,
      fasta_sha256 = metadata$fasta_sha256,
      accession = ACCESSION,
      identity = parse_number(`%IDENTITY`),
      coverage = parse_number(`%COVERAGE`),
      SEQUENCE = SEQUENCE,
      GENE = GENE,
      start = parse_number(START),
      end = parse_number(END),
      strand = as.character(STRAND)
    )
}

tab_manifest <- fasta_manifest %>%
  left_join(call_files, by = "Isolate_ID")
if (any(is.na(tab_manifest$tab_file))) {
  stop("At least one selected assembly lacks a completed PlasmidFinder call.", call. = FALSE)
}

floor_hits <- map_dfr(seq_len(nrow(tab_manifest)), function(i) {
  read_one_tab(tab_manifest$tab_file[[i]], tab_manifest[i, ])
}) %>%
  filter(
    is.finite(identity), is.finite(coverage),
    identity >= SCAN_MIN_ID, coverage >= SCAN_MIN_COV,
    !is.na(GENE), nzchar(GENE)
  )

primary_hits <- floor_hits %>%
  filter(identity >= PRIMARY_MIN_ID, coverage >= PRIMARY_MIN_COV)
stringent_hits <- floor_hits %>%
  filter(identity >= STRINGENT_MIN_ID, coverage >= STRINGENT_MIN_COV)

atomic_write_csv(floor_hits, file.path(DIR_PLASMIDS, "plasmidfinder_hits_long_id80_cov60.csv"))
atomic_write_csv(primary_hits, file.path(DIR_PLASMIDS, "plasmidfinder_hits_long.csv"))
atomic_write_csv(
  primary_hits %>% mutate(raw_hit_id = paste0("PF_", row_number()), .before = 1L),
  file.path(DIR_PLASMIDS, "plasmidfinder_hits_primary_raw.csv")
)
atomic_write_csv(stringent_hits, file.path(DIR_PLASMIDS, "plasmidfinder_hits_long_id90_cov90.csv"))

tier_counts <- bind_rows(
  floor_hits %>% count(Isolate_ID, name = "n_hits") %>% mutate(tier = "id80_cov60"),
  primary_hits %>% count(Isolate_ID, name = "n_hits") %>% mutate(tier = "id80_cov80"),
  stringent_hits %>% count(Isolate_ID, name = "n_hits") %>% mutate(tier = "id90_cov90")
) %>%
  complete(
    Isolate_ID = fasta_manifest$Isolate_ID,
    tier = c("id80_cov60", "id80_cov80", "id90_cov90"),
    fill = list(n_hits = 0L)
  ) %>%
  arrange(match(Isolate_ID, fasta_manifest$Isolate_ID), tier)
atomic_write_csv(tier_counts, file.path(DIR_PLASMIDS, "plasmidfinder_threshold_concordance.csv"))

primary_call_counts <- primary_hits %>%
  count(Isolate_ID, name = "n_primary_hits")
run_manifest <- tab_manifest %>%
  left_join(primary_call_counts, by = "Isolate_ID") %>%
  mutate(
    n_primary_hits = replace_na(n_primary_hits, 0L),
    call_status = "complete",
    no_hit = n_primary_hits == 0L,
    schema_version = SCHEMA_VERSION,
    caller_version = abricate_version,
    database_sha256 = db_sha256,
    scan_min_identity = SCAN_MIN_ID,
    scan_min_coverage = SCAN_MIN_COV,
    primary_min_identity = PRIMARY_MIN_ID,
    primary_min_coverage = PRIMARY_MIN_COV
  )
atomic_write_csv(run_manifest, file.path(DIR_PLASMIDS, "plasmidfinder_run_manifest.csv"))

matrix_hits <- primary_hits %>%
  distinct(Isolate_ID, GENE) %>%
  mutate(present = 1L) %>%
  pivot_wider(names_from = GENE, values_from = present, values_fill = 0L)
pa_matrix <- fasta_manifest %>%
  distinct(Isolate_ID) %>%
  left_join(matrix_hits, by = "Isolate_ID") %>%
  mutate(across(-Isolate_ID, ~ replace_na(as.integer(.x), 0L)))

feature_cols <- setdiff(names(pa_matrix), "Isolate_ID")
if (!length(feature_cols) ||
    any(feature_cols %in% unique(primary_hits$accession)) ||
    any(!vapply(pa_matrix[feature_cols], function(x) all(x %in% c(0L, 1L)), logical(1)))) {
  stop("Canonical plasmid matrix is not an exact binary GENE-label matrix.", call. = FALSE)
}
atomic_write_csv(pa_matrix, file.path(DIR_PLASMIDS, "plasmidfinder_presence_absence.csv"))

replicon_catalog <- primary_hits %>%
  distinct(GENE, accession, Isolate_ID) %>%
  count(GENE, accession, name = "n_episodes") %>%
  arrange(GENE, accession)
atomic_write_csv(replicon_catalog, file.path(DIR_PLASMIDS, "plasmidfinder_replicon_catalog.csv"))

ap001918_genes <- primary_hits %>%
  filter(accession == "AP001918") %>%
  distinct(GENE) %>%
  pull(GENE)
if (length(ap001918_genes) < 3L ||
    !all(c("IncFIB(AP001918)_1", "IncFIC(FII)_1", "IncFIA_1") %in% ap001918_genes)) {
  stop("AP001918 did not resolve to the expected separate IncF replicon labels.", call. = FALSE)
}

positive_episodes <- n_distinct(primary_hits$Isolate_ID)
no_hit_episodes <- nrow(fasta_manifest) - positive_episodes
if (nrow(primary_hits) != EXPECTED_PRIMARY_HITS ||
    n_distinct(primary_hits$GENE) != EXPECTED_PRIMARY_GENES ||
    positive_episodes != EXPECTED_POSITIVE_EPISODES ||
    no_hit_episodes != EXPECTED_NO_HIT_EPISODES) {
  stop(
    "Pinned primary PlasmidFinder baseline changed: observed ",
    nrow(primary_hits), " hits, ", n_distinct(primary_hits$GENE), " genes, ",
    positive_episodes, " positive and ", no_hit_episodes,
    " no-hit episodes; expected 1257/42/422/110.",
    call. = FALSE
  )
}

num_mat <- pa_matrix %>% select(-Isolate_ID) %>% as.matrix()
keep <- which(colSums(num_mat) >= 3L)
if (length(keep) > 1L) {
  mat <- num_mat[, keep, drop = FALSE]
  coocc <- t(mat) %*% mat
  diag(coocc) <- 0L
  coocc[lower.tri(coocc)] <- 0L
  edges <- as.data.frame(as.table(coocc), stringsAsFactors = FALSE) %>%
    filter(Freq > 0L) %>%
    transmute(from = as.character(Var1), to = as.character(Var2), weight = as.integer(Freq))
  atomic_write_csv(edges, file.path(DIR_PLASMIDS, "replicon_cooccurrence_edges.csv"))
  graph <- graph_from_data_frame(edges, directed = FALSE)
  set.seed(20260714)
  grDevices::pdf(file.path(DIR_PLOTS_PLASMIDS, "replicon_cooccurrence.pdf"), width = 7, height = 7)
  print(
    ggraph(graph, layout = "nicely") +
      geom_edge_link(aes(edge_alpha = weight), edge_width = 0.9) +
      scale_edge_alpha_continuous(range = c(0.25, 0.9), name = "Co-occurring episodes") +
      geom_node_point(size = 4) +
      geom_node_text(aes(label = name), repel = TRUE, seed = 20260714, max.overlaps = Inf, size = 3) +
      theme_void() +
      labs(
        title = "PlasmidFinder replicon-marker co-occurrence",
        subtitle = "Exact GENE labels detected in at least three selected episodes",
        caption = "Edges denote episode-level co-occurrence. Marker co-occurrence does not prove that replicons occupy the same physical plasmid or that transfer occurred."
      )
  )
  grDevices::dev.off()
}

if (!file.exists(FILE_MLST_CANONICAL)) stop("Missing canonical MLST input.", call. = FALSE)
mlst <- read_csv(FILE_MLST_CANONICAL, show_col_types = FALSE) %>%
  select(any_of(c("Isolate_ID", "ST", "ST_source")))
bipartite_edges <- primary_hits %>%
  distinct(Isolate_ID, GENE) %>%
  left_join(mlst, by = "Isolate_ID") %>%
  filter(!is.na(ST), nzchar(as.character(ST))) %>%
  count(ST, GENE, name = "weight")
atomic_write_csv(bipartite_edges, file.path(DIR_PLASMIDS, "ST_replicon_descriptive_edges.csv"))

if (nrow(bipartite_edges) > 0L) {
  edge_graph <- bipartite_edges %>%
    transmute(from = paste0("ST:", ST), to = GENE, weight)
  graph <- graph_from_data_frame(edge_graph, directed = FALSE)
  V(graph)$node_type <- ifelse(startsWith(V(graph)$name, "ST:"), "Sequence type", "Replicon marker")
  set.seed(20260714)
  grDevices::pdf(file.path(DIR_PLOTS_PLASMIDS, "ST_vs_replicon_network.pdf"), width = 8, height = 6)
  print(
    ggraph(graph, layout = "fr") +
      geom_edge_link(aes(edge_alpha = weight), edge_width = 0.8) +
      scale_edge_alpha_continuous(range = c(0.2, 0.85), name = "Episodes") +
      geom_node_point(aes(color = node_type, shape = node_type), size = 4) +
      scale_color_manual(
        values = c("Replicon marker" = "#0072B2", "Sequence type" = "#E69F00"),
        name = "Node type"
      ) +
      scale_shape_manual(
        values = c("Replicon marker" = 16, "Sequence type" = 17),
        name = "Node type"
      ) +
      geom_node_text(aes(label = name), repel = TRUE, seed = 20260714, max.overlaps = Inf, size = 3) +
      theme_void() +
      labs(
        title = "Sequence-type and replicon-marker co-occurrence",
        caption = "Descriptive episode counts only. Repeated episodes are not independent, and co-occurrence does not demonstrate physical linkage or transmission."
      )
  )
  grDevices::dev.off()
}

writeLines(
  c(
    "PlasmidFinder canonical screen: PASS",
    paste0("schema=", SCHEMA_VERSION),
    paste0("episodes=", nrow(fasta_manifest)),
    paste0("primary_hits=", nrow(primary_hits)),
    paste0("primary_gene_labels=", n_distinct(primary_hits$GENE)),
    paste0("positive_episodes=", positive_episodes),
    paste0("valid_no_hit_episodes=", no_hit_episodes),
    paste0("database_sha256=", db_sha256),
    "interpretation=replicon markers only; not reconstructed plasmids or transfer evidence"
  ),
  file.path(DIR_PLASMIDS, "PLASMIDFINDER_RUN_COMPLETE.txt")
)
msg("✓ Canonical gene-level PlasmidFinder analysis complete.")
