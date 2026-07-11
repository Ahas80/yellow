#!/usr/bin/env Rscript
# 12d_sv_analysis.R
# Step D: Structural Variant Analysis with Optimization

source("R/wgs_helpers.R")
library(optparse)
library(future.apply)

option_list <- list(
    make_option("--threads", type = "integer", default = 8, help = "Threads"),
    make_option("--pids", type = "character", default = "ALL", help = "Comma-separated PIDs"),
    make_option("--strict", action = "store_true", default = FALSE, help = "Only run SV on 'Same Strain' or close pairs"),
    make_option("--mash_thresh", type = "double", default = 0.001, help = "Mash distance threshold for strict mode")
)
opt <- parse_args(OptionParser(option_list = option_list))

ensure_dirs()
BIN_nucmer <- check_bin("nucmer")
BIN_showdiff <- check_bin("show-diff")
BIN_showcoords <- check_bin("show-coords")

if (BIN_nucmer == "") stop("Nucmer not found.")

samples <- discover_samples(pids = opt$pids)
all_pids <- unique(na.omit(samples$Participant_id))

sv_dir <- file.path(BASE_OUT, "sv")
dir.create(sv_dir, recursive = TRUE, showWarnings = FALSE)

# Parallel setup
workers <- max(1L, floor(opt$threads / 2))
future::plan(future::multisession, workers = workers)

msg("Running SV analysis (Strict Mode: %s)...", opt$strict)

# Load Mash distances if strict mode
mash_dist <- NULL
if (opt$strict) {
    mash_file <- file.path(BASE_OUT, "kmer", "mash_pairs_long.csv")
    if (file.exists(mash_file)) {
        mash_dist <- read_csv(mash_file, show_col_types = FALSE)
    } else {
        warning("Mash distances not found. Strict mode may fall back to all pairs.")
    }
}

process_pid_sv <- function(pid) {
    dfp <- samples %>% filter(Participant_id == pid)
    if (nrow(dfp) < 2) {
        return(NULL)
    }

    comb <- t(combn(dfp$SampleID, 2))
    res_list <- list()

    for (i in seq_len(nrow(comb))) {
        a <- comb[i, 1]
        b <- comb[i, 2]

        # Optimization: Check strict criteria
        if (opt$strict && !is.null(mash_dist)) {
            d <- mash_dist %>%
                filter((A == a & B == b) | (A == b & B == a)) %>%
                pull(Mash_distance)
            if (length(d) == 0 || d > opt$mash_thresh) next
        }

        fa <- dfp$assembly[dfp$SampleID == a][1]
        fb <- dfp$assembly[dfp$SampleID == b][1]
        if (is.na(fa) || is.na(fb) || !file.exists(fa) || !file.exists(fb)) next

        out_prefix <- file.path(sv_dir, paste0(pid, "_", a, "_vs_", b))
        dir.create(dirname(out_prefix), recursive = TRUE, showWarnings = FALSE)
        delta <- paste0(out_prefix, ".delta")

        if (!file.exists(delta)) {
            run_cmd(BIN_nucmer, c("--mum", "-p", out_prefix, fa, fb), fail_ok = TRUE)
        }

        if (file.exists(delta)) {
            # show-diff
            diff_file <- paste0(out_prefix, ".diff")
            run_cmd(BIN_showdiff, c("-q", "-r", delta), log_file = diff_file, fail_ok = TRUE)

            ins <- del <- inv <- trans <- 0L
            ins_bp <- del_bp <- inv_bp <- trans_bp <- 0

            if (file.exists(diff_file) && file.size(diff_file) > 0) {
                df <- suppressWarnings(read_tsv(diff_file, col_names = FALSE, show_col_types = FALSE))
                if (nrow(df)) {
                    # Parsing logic adapted from 12_wgs_exact_compare.R
                    type_col <- which(sapply(df, function(c) all(nchar(unique(na.omit(c))) <= 2)))
                    if (length(type_col)) {
                        tcol <- type_col[1]
                        num_cols <- which(sapply(df, function(x) all(grepl("^\\d+$", as.character(x)))))
                        lcol <- if (length(num_cols)) num_cols[length(num_cols)] else NA_integer_
                        types <- df[[tcol]]
                        lens <- if (!is.na(lcol)) as.numeric(df[[lcol]]) else rep(1, nrow(df))
                        ins <- sum(types == "I", na.rm = TRUE)
                        del <- sum(types == "D", na.rm = TRUE)
                        inv <- sum(types == "R", na.rm = TRUE)
                        trans <- sum(types == "T", na.rm = TRUE)
                        ins_bp <- sum(lens[types == "I"], na.rm = TRUE)
                        del_bp <- sum(lens[types == "D"], na.rm = TRUE)
                        inv_bp <- sum(lens[types == "R"], na.rm = TRUE)
                        trans_bp <- sum(lens[types == "T"], na.rm = TRUE)
                    }
                }
            } else if (BIN_showcoords != "") {
                # Fallback to show-coords for inversions
                coords <- paste0(out_prefix, ".coords")
                run_cmd(BIN_showcoords, c("-TrH", delta), log_file = coords, fail_ok = TRUE)
                if (file.exists(coords) && file.size(coords) > 0) {
                    co <- suppressMessages(read_tsv(coords, col_names = FALSE, show_col_types = FALSE))
                    if (nrow(co) >= 1) {
                        # S1 E1 S2 E2
                        s1 <- as.numeric(co[[1]])
                        e1 <- as.numeric(co[[2]])
                        s2 <- as.numeric(co[[3]])
                        e2 <- as.numeric(co[[4]])
                        ori <- sign(e1 - s1) * sign(e2 - s2)
                        inv <- sum(ori < 0, na.rm = TRUE)
                    }
                }
            }

            res_list[[length(res_list) + 1]] <- tibble(
                Participant_id = pid, A = a, B = b,
                insertions = ins, deletions = del, inversions = inv, translocations = trans,
                ins_bp = ins_bp, del_bp = del_bp, inv_bp = inv_bp, trans_bp = trans_bp,
                SV_bp_total = sum(ins_bp, del_bp, inv_bp, trans_bp, na.rm = TRUE)
            )
        }
    }
    bind_rows(res_list)
}

sv_results <- future_lapply(all_pids, process_pid_sv, future.seed = TRUE)
sv_all <- bind_rows(sv_results)

if (nrow(sv_all) > 0) {
    write_csv(sv_all, file.path(sv_dir, "sv_pairs.csv"))
    msg("SV analysis complete. Processed %d pairs.", nrow(sv_all))
} else {
    msg("No SV pairs processed.")
}
