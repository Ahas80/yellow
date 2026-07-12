# R/wgs_helpers.R
# Shared helper functions for WGS pipeline scripts (12a-12e)

suppressPackageStartupMessages({
    library(dplyr)
    library(readr)
    library(stringr)
    library(tibble)
    library(tidyr)
    library(data.table)
    library(ggplot2)
    library(jsonlite)
    library(fs)
    library(processx)
})

# --- Environment Setup ---
# Add common conda paths to PATH
conda_paths <- c(
    file.path(Sys.getenv("HOME"), "miniconda3/envs/asm-snp-x86/bin"),
    file.path(Sys.getenv("HOME"), "miniconda3/envs/wgs/bin"),
    file.path(Sys.getenv("HOME"), "miniconda3/envs/panaroo/bin")
)
curr_path <- Sys.getenv("PATH")
Sys.setenv(PATH = paste(c(conda_paths, curr_path), collapse = .Platform$path.sep))

# --- Configuration & Paths ---
# Directories (fixed by project structure)
# These should be consistent across all scripts
ASM_DIR <- "ont-yellow-routine-fastas"
READS_DIR <- "data/reads"
BASE_OUT <- "results/wgs"
PLOT_DIR <- "results/plots"

# Ensure directories exist
ensure_dirs <- function() {
    dirs <- c(
        file.path(BASE_OUT, "logs"),
        file.path(BASE_OUT, "core"),
        file.path(BASE_OUT, "kmer"),
        file.path(BASE_OUT, "pangenome"),
        file.path(BASE_OUT, "sv"),
        file.path(BASE_OUT, "vaf"),
        file.path(BASE_OUT, "plasmids"),
        file.path(BASE_OUT, "model"),
        file.path(BASE_OUT, "masks"),
        file.path(BASE_OUT, "reports"),
        file.path(BASE_OUT, "sensitivity"),
        file.path(BASE_OUT, "qc"),
        PLOT_DIR
    )
    for (d in dirs) dir.create(d, recursive = TRUE, showWarnings = FALSE)
}

# --- Logging & Messaging ---
msg <- function(fmt, ...) {
    cat(sprintf("[%s] %s\n", format(Sys.time(), "%Y-%m-%d %H:%M:%S"), sprintf(fmt, ...)))
    flush.console()
}

# --- Command Execution ---
run_cmd <- function(cmd, args = character(), log_file = NULL, env = character(), echo = FALSE, fail_ok = FALSE) {
    if (is.null(log_file)) {
        log_file <- file.path(
            BASE_OUT, "logs",
            paste0(gsub("[^A-Za-z0-9_]+", "_", basename(cmd)), "_", as.integer(Sys.time()), ".log")
        )
    }
    dir.create(dirname(log_file), recursive = TRUE, showWarnings = FALSE)

    if (echo) msg("CMD: %s %s", cmd, paste(vapply(args, as.character, ""), collapse = " "))

    status <- suppressWarnings(
        system2(command = cmd, args = args, stdout = log_file, stderr = log_file, env = env)
    )
    status <- as.integer(if (is.null(status)) 0L else status)

    if (status != 0L && !fail_ok) {
        msg("⚠️  Command failed (status=%d). See %s", status, log_file)
    }
    invisible(list(status = status, log = log_file))
}

check_bin <- function(bin_candidates) {
    cand <- unique(c(
        file.path(Sys.getenv("CONDA_PREFIX"), "bin", bin_candidates),
        file.path(Sys.getenv("HOME"), "miniconda3/envs/asm-snp-x86/bin", bin_candidates),
        file.path(Sys.getenv("HOME"), "opt/miniconda3/envs/asm-snp-x86/bin", bin_candidates),
        file.path(Sys.getenv("HOME"), "miniconda3/envs/wgs/bin", bin_candidates),
        file.path(Sys.getenv("HOME"), "miniconda3/envs/panaroo/bin", bin_candidates),
        bin_candidates
    ))
    hits <- Sys.which(cand)
    p <- hits[hits != ""]
    if (length(p)) {
        return(unname(p[1]))
    }
    ""
}

tool_version <- function(bin, args = c("--version")) {
    if (bin == "") {
        return(NA_character_)
    }
    call <- function(a) {
        out <- suppressWarnings(tryCatch(system2(bin, a, stdout = TRUE, stderr = TRUE),
            error = function(e) character()
        ))
        if (length(out)) paste(head(out, 3), collapse = " | ") else NA_character_
    }
    v <- call(args)
    if (is.na(v)) v <- call(c("--help"))
    v
}

# --- File I/O Helpers ---
safe_read_csv <- function(p) {
    if (file.exists(p)) suppressMessages(readr::read_csv(p, show_col_types = FALSE)) else NULL
}

safe_read_tsv <- function(p) {
    if (file.exists(p)) suppressMessages(readr::read_tsv(p, show_col_types = FALSE)) else NULL
}

safe_copy <- function(from, to, overwrite = TRUE) {
    nf <- tryCatch(normalizePath(from, winslash = "/", mustWork = FALSE), error = function(e) from)
    nt <- tryCatch(normalizePath(to, winslash = "/", mustWork = FALSE), error = function(e) to)
    if (identical(nf, nt)) {
        return(TRUE)
    }
    dir.create(dirname(to), recursive = TRUE, showWarnings = FALSE)
    ok <- base::file.copy(from, to, overwrite = overwrite)
    if (!ok) warning("safe_copy: failed to copy '", from, "' -> '", to, "'")
    ok
}

write_plot <- function(p, path, w = 8, h = 6, dpi = 300) {
    if (inherits(p, "ggplot")) try(ggsave(path, p, width = w, height = h, dpi = dpi), silent = TRUE)
}

# --- Data Processing Helpers ---
tp_norm <- function(x) {
    x <- as.character(x)
    is_uricult <- str_detect(x, regex("uricult", ignore_case = TRUE))
    tp_num <- suppressWarnings(as.integer(str_extract(x, "\\d+")))
    tp_lab <- case_when(
        is_uricult ~ "Uricult",
        !is.na(tp_num) ~ paste0("T", tp_num),
        TRUE ~ "Unscheduled"
    )
    tibble(tp_lab = tp_lab, tp_num = tp_num)
}

sanitize_id <- function(x) {
    x <- gsub("[^A-Za-z0-9._-]", "_", x)
    ifelse(grepl("^-", x), paste0("X", x), x)
}

parse_pid_tp <- function(paths) {
    if (length(paths) == 0) {
        return(tibble(Participant_id = character(), Timepoint = character()))
    }
    stems <- tools::file_path_sans_ext(basename(paths))
    PID <- stringr::str_extract(stems, "^[A-Za-z0-9]+")
    TP1 <- stringr::str_extract(stems, "(?i)(T\\d+|Uricult|U\\d+)")
    num <- suppressWarnings(as.integer(stringr::str_extract(stems, "(?<=_)\\d+|(?<=-)\\d+")))
    TP <- ifelse(is.na(TP1) & !is.na(num), paste0("T", num), TP1)
    tibble(Participant_id = PID, Timepoint = TP)
}

fasta_stats <- function(fa) {
    if (is.na(fa) || !nzchar(fa) || !file.exists(fa)) {
        return(tibble(n_contigs = NA_integer_, n50 = NA_real_, total_bp = NA_real_))
    }
    con <- if (grepl("\\.gz$", fa, ignore.case = TRUE)) gzfile(fa, open = "rt") else file(fa, "rt")
    on.exit(try(close(con), silent = TRUE), add = TRUE)
    L <- readLines(con, warn = FALSE)

    sel <- which(startsWith(L, ">"))
    lens <- integer()
    if (length(sel)) {
        starts <- c(sel, length(L) + 1L)
        for (i in seq_along(sel)) {
            seq_lines <- L[(sel[i] + 1):(starts[i + 1] - 1)]
            seq_lines <- seq_lines[seq_lines != ""]
            lens <- c(lens, sum(nchar(seq_lines), na.rm = TRUE))
        }
    } else {
        lens <- sum(nchar(L[!startsWith(L, ">")]), na.rm = TRUE)
    }
    lens <- sort(as.integer(lens), decreasing = TRUE)
    total <- sum(lens)
    n50 <- if (total > 0) lens[min(which(cumsum(lens) >= total / 2))] else NA_real_
    tibble(n_contigs = length(lens), n50 = as.numeric(n50), total_bp = as.numeric(total))
}

# --- Progress Tracking ---
read_progress <- function(pf) {
    if (!file.exists(pf)) {
        return(tibble::tibble(stage = character(), pid = character(), status = character(), t = as.POSIXct(character())))
    }
    sz <- suppressWarnings(file.size(pf))
    if (is.na(sz) || sz == 0) {
        return(tibble::tibble(stage = character(), pid = character(), status = character(), t = as.POSIXct(character())))
    }
    out <- data.table::fread(pf, nThread = 1, showProgress = FALSE)
    tibble::as_tibble(out)
}

write_progress <- function(df, pf) {
    dir.create(dirname(pf), recursive = TRUE, showWarnings = FALSE)
    tmp <- paste0(pf, ".", Sys.getpid(), ".tmp")
    readr::write_csv(df, tmp)
    ok <- suppressWarnings(try(file.rename(tmp, pf), silent = TRUE))
    if (!isTRUE(ok)) {
        try(file.copy(tmp, pf, overwrite = TRUE), silent = TRUE)
        try(unlink(tmp), silent = TRUE)
    }
}

progress_tick <- function(stage, pid = NA_character_, status = "tick") {
    if (!is.character(pid)) pid <- ifelse(is.na(pid), NA_character_, as.character(pid))

    pf <- file.path(BASE_OUT, "reports", "progress.csv")
    df <- read_progress(pf)
    if (nrow(df)) {
        df$stage <- as.character(df$stage)
        df$pid <- as.character(df$pid)
        df$status <- as.character(df$status)
    } else {
        df <- tibble::tibble(stage = character(), pid = character(), status = character(), t = as.POSIXct(character()))
    }
    df <- dplyr::bind_rows(
        df,
        tibble(
            stage  = as.character(stage),
            pid    = as.character(pid),
            status = as.character(status),
            t      = Sys.time()
        )
    )
    write_progress(df, pf)

    # tiny progress plot
    plotf <- file.path(PLOT_DIR, "_progress.png")
    summ <- df %>% count(stage, status, name = "n")
    p <- ggplot(summ, aes(stage, n, fill = status)) +
        geom_col() +
        coord_flip() +
        labs(title = "Pipeline progress (rolling)", x = NULL, y = "Events") +
        theme_minimal(base_size = 10)
    write_plot(p, plotf, 6, 4, 150)
}

# --- Sample Discovery ---
discover_samples <- function(asm_dir = ASM_DIR, reads_dir = READS_DIR, pids = "ALL") {
    msg("Loading Longcycler-only analysis assemblies and discovering reads…")

    if (!exists("load_analysis_assemblies", mode = "function")) {
        stop("Longcycler analysis-manifest helper is unavailable; source 00_config.R first.")
    }
    analysis_manifest <- load_analysis_assemblies(FILE_ANALYSIS_ASSEMBLY_MANIFEST, require_files = TRUE)
    asm_files <- analysis_manifest$full_path
    reads_R1 <- list.files(reads_dir, pattern = "_R1\\.(fastq|fq)(\\.gz)?$", full.names = TRUE)
    reads_R2 <- gsub("_R1\\.(fastq|fq)(\\.gz)?$", "_R2.\\1\\2", reads_R1, perl = TRUE)
    paired_ok <- file.exists(reads_R2)
    reads_R1 <- reads_R1[paired_ok]
    reads_R2 <- reads_R2[paired_ok]

    if (length(asm_files) > 0) {
        samples_from_asm <- tibble(
            assembly = asm_files,
            SampleID = tools::file_path_sans_ext(basename(asm_files))
        )
        samples <- samples_from_asm %>%
            bind_cols(parse_pid_tp(samples_from_asm$assembly)) %>%
            mutate(Timepoint = dplyr::coalesce(Timepoint, "Unscheduled"))
    } else if (length(reads_R1) > 0) {
        samples_from_reads <- tibble(
            R1 = reads_R1,
            R2 = reads_R2,
            SampleID = sub("_R1\\.(fastq|fq)(\\.gz)?$", "", basename(reads_R1))
        )
        samples <- samples_from_reads %>%
            bind_cols(parse_pid_tp(samples_from_reads$SampleID)) %>%
            mutate(
                Timepoint = dplyr::coalesce(Timepoint, "Unscheduled"),
                assembly = NA_character_
            )
    } else {
        samples <- tibble(
            assembly = character(), SampleID = character(),
            Participant_id = character(), Timepoint = character()
        )
    }

    # Attach reads by SampleID if present
    if (nrow(samples) > 0 && length(reads_R1) > 0) {
        reads_tbl <- tibble(R1 = reads_R1, R2 = reads_R2) %>%
            mutate(SampleID = sub("_R1\\.(fastq|fq)(\\.gz)?$", "", basename(R1)))
        samples <- samples %>% left_join(reads_tbl, by = "SampleID")
    } else {
        if (!"R1" %in% names(samples)) samples$R1 <- NA_character_
        if (!"R2" %in% names(samples)) samples$R2 <- NA_character_
    }

    # Attach metadata from the validated Longcycler-only manifest.
    if (nrow(analysis_manifest) > 0 && "assembly" %in% names(samples)) {
        assembly_raw <- analysis_manifest
        file_col <- intersect(c("file_name", "file", "assembly", "fasta", "fasta_file"), names(assembly_raw))
        if (length(file_col)) {
            sid_col <- intersect(c("Isolate_ID", "SampleID", "Sample_ID", "IsolateID"), names(assembly_raw))
            sid_col <- if (length(sid_col)) sid_col[1] else NULL

            assembly_df <- assembly_raw %>%
                dplyr::mutate(
                    assembly = if ("full_path" %in% names(assembly_raw)) {
                        normalizePath(.data$full_path, winslash = "/", mustWork = FALSE)
                    } else {
                        normalizePath(file.path(asm_dir, .data[[file_col[1]]]), winslash = "/", mustWork = FALSE)
                    },
                    SampleID_meta = if (!is.null(sid_col)) as.character(.data[[sid_col]]) else tools::file_path_sans_ext(.data[[file_col[1]]]),
                    Participant_id_meta = if ("Participant_id" %in% names(.)) as.character(.data[["Participant_id"]]) else NA_character_,
                    Timepoint_meta = if ("Timepoint" %in% names(.)) as.character(.data[["Timepoint"]]) else NA_character_
                ) %>%
                dplyr::select(assembly, SampleID_meta, Participant_id_meta, Timepoint_meta)

            samples <- samples %>%
                dplyr::left_join(assembly_df, by = "assembly") %>%
                dplyr::mutate(
                    SampleID       = dplyr::coalesce(SampleID_meta, SampleID),
                    Participant_id = dplyr::coalesce(Participant_id_meta, Participant_id),
                    Timepoint      = dplyr::coalesce(Timepoint_meta, Timepoint)
                ) %>%
                dplyr::select(-SampleID_meta, -Participant_id_meta, -Timepoint_meta)
        }
    }

    samples <- samples %>% distinct(SampleID, .keep_all = TRUE)

    # Normalize timepoints
    if (nrow(samples) > 0) {
        tp <- tp_norm(samples$Timepoint)
        samples <- samples %>% mutate(tp_lab = tp$tp_lab, tp_num = tp$tp_num)
    }

    # PID restriction
    if (!identical(toupper(pids), "ALL") && nrow(samples) > 0) {
        keep_pids <- strsplit(pids, ",", fixed = TRUE)[[1]] |> trimws()
        samples <- samples %>% filter(Participant_id %in% keep_pids)
    }

    samples <- samples %>% mutate(
        SampleID       = as.character(SampleID),
        Participant_id = as.character(Participant_id),
        Timepoint      = as.character(Timepoint)
    )

    # Sanitize IDs
    samples$SampleID <- sanitize_id(samples$SampleID)

    return(samples)
}

# --- Restored Functions for Modular Pipeline (12a-12e) ---

# 1. QC Configuration
get_qc_config <- function() {
    list(
        MAX_CONTIGS = 200,
        MIN_N50 = 20000,
        MIN_GENOME_SIZE = 4.0e6,
        MAX_GENOME_SIZE = 6.0e6,
        MIN_COMPLETENESS = 95.0,
        MAX_CONTAMINATION = 5.0
    )
}

# 2. Logging Adapters
log_info <- function(...) {
    msg <- paste0(...)
    message(format(Sys.time(), "[%Y-%m-%d %H:%M:%S] [INFO] "), msg)
}

log_warn <- function(...) {
    msg <- paste0(...)
    message(format(Sys.time(), "[%Y-%m-%d %H:%M:%S] [WARN] "), msg)
}

log_error <- function(...) {
    msg <- paste0(...)
    message(format(Sys.time(), "[%Y-%m-%d %H:%M:%S] [ERROR] "), msg)
}

# 3. Data Loading
load_qc_summary <- function(qc_file = NULL) {
    if (is.null(qc_file)) {
        # Try to find it in standard location
        # Assuming DIR_RESULTS is available from 00_config.R
        if (exists("DIR_RESULTS")) {
            qc_file <- file.path(DIR_RESULTS, "wgs", "qc_summary.csv")
        } else {
            qc_file <- "results/wgs/qc_summary.csv"
        }
    }

    if (!file.exists(qc_file)) {
        stop("QC summary file not found: ", qc_file, "\nPlease run 12a_wgs_qc.R first.")
    }

    readr::read_csv(qc_file, show_col_types = FALSE)
}

get_valid_genomes <- function(pid, qc_df) {
    required <- c("Participant_id", "QC_PASS", "selected_canonical", "full_path")
    missing <- setdiff(required, names(qc_df))
    if (length(missing)) stop("get_valid_genomes requires Longcycler manifest columns: ", paste(missing, collapse = ", "))
    candidate <- qc_df %>%
        dplyr::filter(
            as.character(Participant_id) == as.character(pid),
            QC_PASS %in% TRUE,
            selected_canonical %in% TRUE
        )
    assert_analysis_assembly_manifest(
        candidate,
        context = paste0("get_valid_genomes participant ", pid),
        require_selected = TRUE,
        require_qc = TRUE,
        require_files = TRUE,
        require_unique_episode = TRUE
    )
    candidate$full_path
}

# 4. Tool Checks
check_wgs_tool <- function(tool_name) {
    path <- Sys.which(tool_name)
    if (path == "") {
        warning(sprintf("Tool '%s' not found in PATH. WGS steps requiring it will fail.", tool_name))
        return(FALSE)
    }
    return(TRUE)
}
