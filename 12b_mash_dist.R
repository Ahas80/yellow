#!/usr/bin/env Rscript
# 12b_mash_dist.R
# Step B: Mash Distance Analysis

source("R/wgs_helpers.R")
library(optparse)

option_list <- list(
    make_option("--threads", type = "integer", default = 8, help = "Threads"),
    make_option("--force", action = "store_true", default = FALSE, help = "Force rebuild")
)
opt <- parse_args(OptionParser(option_list = option_list))

ensure_dirs()
BIN_mash <- check_bin("mash")
if (BIN_mash == "") stop("Mash not found.")

samples <- discover_samples() # All samples
asm_tbl <- samples %>% filter(!is.na(assembly) & file.exists(assembly))

kmer_dir <- file.path(BASE_OUT, "kmer")
sk_dir <- file.path(kmer_dir, "sketches")
dir.create(sk_dir, recursive = TRUE, showWarnings = FALSE)

msg("Sketching %d assemblies...", nrow(asm_tbl))

# Sketch
for (i in seq_len(nrow(asm_tbl))) {
    sname <- asm_tbl$SampleID[i]
    out_msh <- file.path(sk_dir, paste0(sname, ".msh"))
    if (file.exists(out_msh) && !opt$force) next
    run_cmd(BIN_mash, c("sketch", "-p", "1", "-o", out_msh, asm_tbl$assembly[i]))
}

# Paste
paste_out <- file.path(kmer_dir, "mash_paste.msh")
msh_files <- list.files(sk_dir, pattern = "\\.msh$", full.names = TRUE)
if (length(msh_files) > 0) {
    if (opt$force || !file.exists(paste_out)) {
        # Write manifest to avoid arg length limits
        manifest <- file.path(kmer_dir, "msh_manifest.txt")
        writeLines(msh_files, manifest)
        # Use xargs or similar if list is huge, but mash paste takes files
        # We'll use the manifest approach if mash supports it, or just pass args
        # Original script used processx with args.
        run_cmd(BIN_mash, c("paste", paste_out, msh_files))
    }
}

# Dist
dist_out <- file.path(kmer_dir, "mash_all_vs_all.tab")
if (file.exists(paste_out) && (opt$force || !file.exists(dist_out))) {
    run_cmd(BIN_mash, c("dist", "-p", as.character(opt$threads), paste_out, paste_out),
        log_file = dist_out
    ) # Capture stdout to file
}

msg("Mash analysis complete.")
