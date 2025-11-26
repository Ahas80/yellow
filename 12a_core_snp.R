#!/usr/bin/env Rscript
# 12a_core_snp.R
# Step A: Core SNP Analysis using minimap2, samtools, bcftools, and Gubbins

source("R/wgs_helpers.R")
library(optparse)
library(future.apply)

option_list <- list(
    make_option("--threads", type = "integer", default = 8, help = "Threads (default 8)"),
    make_option("--pids", type = "character", default = "ALL", help = "Comma-separated PIDs"),
    make_option("--force", action = "store_true", default = FALSE, help = "Force rebuild"),
    make_option("--min_mapq", type = "integer", default = 0, help = "Min mapping quality"),
    make_option("--min_qual", type = "integer", default = 20, help = "Min base quality"),
    make_option("--mpileup_max_depth", type = "integer", default = 8000, help = "Max depth for mpileup"),
    make_option("--min_depth", type = "integer", default = 10, help = "Min depth for callable sites")
)
opt <- parse_args(OptionParser(option_list = option_list))

ensure_dirs()
samples <- discover_samples(pids = opt$pids)
all_pids <- unique(na.omit(samples$Participant_id))

# Tools
BIN_minimap2 <- check_bin("minimap2")
BIN_samtools <- check_bin("samtools")
BIN_bcftools <- check_bin("bcftools")
BIN_gubbins <- check_bin(c("gubbins", "run_gubbins.py"))
BIN_snpdists <- check_bin("snp-dists")
BIN_mosdepth <- check_bin("mosdepth")

if (any(c(BIN_minimap2, BIN_samtools, BIN_bcftools) == "")) {
    stop("Missing required tools: minimap2, samtools, or bcftools.")
}

# Parallel setup
workers <- max(1L, min(length(all_pids), floor(opt$threads / 2L)))
future::plan(future::multisession, workers = workers)
per_tool_threads <- max(1L, floor(opt$threads / workers))
msg("Parallel plan: %d workers, %d threads per tool.", workers, per_tool_threads)

# --- Helper Functions ---

pick_ref_for_pid <- function(df_pid) {
    dfp <- df_pid %>%
        filter(!is.na(assembly) & file.exists(assembly)) %>%
        mutate(stats = purrr::map(assembly, fasta_stats)) %>%
        tidyr::unnest(stats) %>%
        arrange(is.na(tp_num), tp_num, desc(n50), desc(total_bp))

    if (nrow(dfp) == 0) {
        return(NULL)
    }
    tibble(ref_fa = dfp$assembly[1], ref_sample = dfp$SampleID[1])
}

make_bam_for_sample <- function(ref_fa, srow) {
    s <- srow[1, ]
    bam_dir <- file.path(BASE_OUT, "bam")
    out_bam <- file.path(bam_dir, s$Participant_id, paste0(s$SampleID, ".markdup.bam"))
    out_bai <- paste0(out_bam, ".bai")
    dir.create(dirname(out_bam), recursive = TRUE, showWarnings = FALSE)

    if (file.exists(out_bam) && file.exists(out_bai)) {
        return(out_bam)
    }

    sam_sort <- file.path(bam_dir, s$Participant_id, paste0(s$SampleID, ".sorted.bam"))
    tmp_prefix <- file.path(bam_dir, s$Participant_id, paste0(s$SampleID, ".tmp"))
    dir.create(dirname(tmp_prefix), recursive = TRUE, showWarnings = FALSE)

    rg <- sprintf("@RG\\tID:%s\\tSM:%s", s$SampleID, s$SampleID)
    rg_sam <- sprintf("ID:%s\tSM:%s", s$SampleID, s$SampleID)

    # Prefer reads if available
    use_reads <- !is.na(s$R1) && !is.na(s$R2) && file.exists(s$R1) && file.exists(s$R2)

    if (use_reads) {
        mm2_args <- c("-t", as.character(per_tool_threads), "-ax", "sr", "-R", rg, ref_fa, s$R1, s$R2)
    } else {
        if (is.na(s$assembly) || !file.exists(s$assembly)) {
            return(NA_character_)
        }
        mm2_args <- c("-t", as.character(per_tool_threads), "-ax", "asm10", ref_fa, s$assembly)
    }

    sort_threads <- as.integer(max(1L, min(per_tool_threads, 4L)))
    pipe_cmd <- sprintf(
        "%s %s | %s sort -@ %d -T %s -O BAM -o %s -",
        shQuote(BIN_minimap2), paste(vapply(mm2_args, shQuote, character(1)), collapse = " "),
        shQuote(BIN_samtools), sort_threads, shQuote(tmp_prefix), shQuote(sam_sort)
    )

    run_cmd("/bin/bash", c("-c", sprintf("set -o pipefail; %s", pipe_cmd)),
        log_file = file.path(BASE_OUT, "logs", paste0("map_sort_", s$SampleID, ".log")),
        fail_ok = FALSE
    )

    if (!file.exists(sam_sort)) {
        return(NA_character_)
    }

    run_cmd(BIN_samtools, c("markdup", "-@", as.character(per_tool_threads), "-s", sam_sort, out_bam), fail_ok = TRUE)

    if (!file.exists(out_bam) || file.size(out_bam) == 0) {
        safe_copy(sam_sort, out_bam, overwrite = TRUE)
    }

    tmp_bam <- paste0(out_bam, ".rg.bam")
    run_cmd(BIN_samtools, c("addreplacerg", "-r", rg_sam, "-o", tmp_bam, out_bam), fail_ok = FALSE)
    file.rename(tmp_bam, out_bam)
    run_cmd(BIN_samtools, c("index", "-@", as.character(min(4L, per_tool_threads)), out_bam), fail_ok = FALSE)
    unlink(sam_sort)

    if (file.exists(out_bam) && file.exists(out_bai)) out_bam else NA_character_
}

callable_bp_mosdepth <- function(bam, min_depth = opt$min_depth) {
    if (BIN_mosdepth == "") {
        return(NA_real_)
    }
    prefix <- file.path(dirname(bam), paste0(basename(bam), ".mosdepth"))
    if (!file.exists(paste0(prefix, ".global.dist.txt"))) {
        run_cmd(BIN_mosdepth, c("-t", as.character(min(4L, per_tool_threads)), "--no-per-base", "--fast-mode", "-n", prefix, bam),
            log_file = file.path(BASE_OUT, "logs", paste0("mosdepth_", basename(bam), ".log")), fail_ok = TRUE
        )
    }
    p <- paste0(prefix, ".global.dist.txt")
    if (!file.exists(p)) {
        return(NA_real_)
    }
    df <- tryCatch(readr::read_tsv(p, col_names = FALSE, show_col_types = FALSE), error = function(e) NULL)
    if (is.null(df)) {
        return(NA_real_)
    }
    if (ncol(df) == 2) names(df) <- c("depth", "count") else names(df)[1:3] <- c("chrom", "depth", "count")
    df$depth <- suppressWarnings(as.numeric(df$depth))
    df$count <- suppressWarnings(as.numeric(df$count))
    sum(df$count[df$depth >= min_depth], na.rm = TRUE)
}

# --- Main Analysis Loop ---
A_results <- future_lapply(all_pids, function(pid) {
    df_pid <- samples %>% filter(Participant_id == pid)
    if (!nrow(df_pid)) {
        return(NULL)
    }

    ref <- pick_ref_for_pid(df_pid)
    if (is.null(ref)) {
        msg("PID %s: no reference found.", pid)
        return(NULL)
    }
    ref_fa <- ref$ref_fa
    msg("PID %s: reference=%s", pid, basename(ref_fa))

    if (!file.exists(paste0(ref_fa, ".fai"))) run_cmd(BIN_samtools, c("faidx", ref_fa))

    bams <- lapply(split(df_pid, df_pid$SampleID), function(s) make_bam_for_sample(ref_fa, s))
    bam_vec <- unlist(bams, use.names = TRUE)
    bam_vec <- bam_vec[!is.na(bam_vec) & file.exists(bam_vec)]

    if (!length(bam_vec)) {
        return(NULL)
    }

    vcf_out <- file.path(BASE_OUT, "core", paste0("joint_", pid, ".vcf.gz"))
    bam_list_file <- file.path(BASE_OUT, "core", paste0("bam_list_", pid, ".txt"))
    writeLines(unname(bam_vec), bam_list_file)

    if (opt$force || !file.exists(vcf_out)) {
        mp <- sprintf(
            "%s mpileup -f %s -a AD,DP -q %d -Q %d -d %d -Ou -b %s",
            shQuote(BIN_bcftools), shQuote(ref_fa), opt$min_mapq, opt$min_qual, opt$mpileup_max_depth, shQuote(bam_list_file)
        )
        cl <- sprintf("%s call -m --ploidy 1 --threads %d -Oz -o %s", shQuote(BIN_bcftools), per_tool_threads, shQuote(vcf_out))
        run_cmd("/bin/bash", c("-c", sprintf("set -o pipefail; %s | %s", mp, cl)), fail_ok = FALSE)
        run_cmd(BIN_bcftools, c("index", "-t", vcf_out))
    }

    # Consensus
    cons_dir <- file.path(BASE_OUT, "core", paste0("consensus_", pid))
    dir.create(cons_dir, showWarnings = FALSE, recursive = TRUE)
    cons_fastas <- c()

    smps_in_vcf <- system2(BIN_bcftools, c("query", "-l", vcf_out), stdout = TRUE)
    vcf_name_map <- setNames(smps_in_vcf, sub("\\.markdup\\.bam$", "", basename(smps_in_vcf)))

    for (sid in df_pid$SampleID) {
        cons_fa <- file.path(cons_dir, paste0(sid, ".fa"))
        sname <- vcf_name_map[[sid]]
        if (is.null(sname)) next

        if (opt$force || !file.exists(cons_fa)) {
            tmp_vcf <- file.path(cons_dir, paste0("tmp_", sid, ".vcf.gz"))
            run_cmd(BIN_bcftools, c("view", "-s", sname, "-@", as.character(per_tool_threads), "-Oz", "-o", tmp_vcf, vcf_out), fail_ok = TRUE)
            run_cmd(BIN_bcftools, c("index", "-t", tmp_vcf), fail_ok = TRUE)
            cmd <- sprintf(
                "%s consensus -f %s -s %s %s > %s",
                shQuote(BIN_bcftools), shQuote(ref_fa), shQuote(sname), shQuote(tmp_vcf), shQuote(cons_fa)
            )
            run_cmd("/bin/sh", c("-c", cmd), fail_ok = TRUE)
            unlink(c(tmp_vcf, paste0(tmp_vcf, ".tbi")))
        }
        if (file.exists(cons_fa) && file.size(cons_fa) > 5) cons_fastas <- c(cons_fastas, cons_fa)
    }

    # Multi-FASTA
    multi_fa <- file.path(BASE_OUT, "core", paste0("multi_", pid, ".fa"))
    if (length(cons_fastas) > 0) {
        tmp_multi <- paste0(multi_fa, ".tmp")
        if (file.exists(tmp_multi)) unlink(tmp_multi)
        for (fa in cons_fastas) {
            L <- readLines(fa, warn = FALSE)
            header <- paste0(">", tools::file_path_sans_ext(basename(fa)))
            seq_str <- paste(L[!startsWith(L, ">")], collapse = "")
            cat(header, "\n", seq_str, "\n", file = tmp_multi, sep = "", append = TRUE)
        }
        file.rename(tmp_multi, multi_fa)
    }

    # Gubbins
    core_align <- file.path(BASE_OUT, "core", paste0("core_alignment_", pid, ".fa"))
    if (BIN_gubbins != "" && file.exists(multi_fa)) {
        outdir <- file.path(BASE_OUT, "core", paste0("gubbins_", pid))
        dir.create(outdir, showWarnings = FALSE, recursive = TRUE)

        gub_cmd <- sprintf("%s --threads %d --prefix %s %s", shQuote(BIN_gubbins), per_tool_threads, shQuote(paste0("run_", pid)), shQuote(multi_fa))
        run_cmd("/bin/sh", c("-c", sprintf("cd %s && %s", shQuote(outdir), gub_cmd)), fail_ok = TRUE)

        cand_align <- list.files(outdir, pattern = "filtered_polymorphic_sites\\.fa(sta)?$", full.names = TRUE)
        if (length(cand_align)) safe_copy(cand_align[1], core_align, overwrite = TRUE)
    }

    if (!file.exists(core_align) && file.exists(multi_fa)) safe_copy(multi_fa, core_align, overwrite = TRUE)

    # SNP Dists
    dist_tsv <- file.path(BASE_OUT, "core", paste0("core_snp_matrix_", pid, ".tsv"))
    if (BIN_snpdists != "" && file.exists(core_align)) {
        run_cmd("/bin/sh", c("-c", sprintf("%s -c %s > %s", shQuote(BIN_snpdists), shQuote(core_align), shQuote(dist_tsv))), fail_ok = TRUE)
    }

    return(list(pid = pid, vcf = vcf_out, dist = dist_tsv))
}, future.seed = TRUE)

msg("Core SNP analysis complete.")
