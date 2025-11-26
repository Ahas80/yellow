#!/usr/bin/env Rscript
# 12_wgs_runner.R
# Master script to run the modular WGS pipeline (12a-12e) sequentially.
# Usage: Rscript 12_wgs_runner.R [--threads 8] [--force]

library(optparse)

option_list <- list(
    make_option("--threads", type = "integer", default = 8, help = "Threads to use"),
    make_option("--force", action = "store_true", default = FALSE, help = "Force rebuild of all steps"),
    make_option("--stage", type = "character", default = "all", help = "Run specific stage (all, a, b, c, d, e)")
)
opt <- parse_args(OptionParser(option_list = option_list))

run_step <- function(script, description) {
    cat(sprintf("\n=== Running %s (%s) ===\n", script, description))
    cmd <- sprintf("Rscript %s --threads %d %s", script, opt$threads, if (opt$force) "--force" else "")
    res <- system(cmd)
    if (res != 0) {
        stop(sprintf("Step %s failed with exit code %d", script, res))
    }
    cat(sprintf("=== %s completed successfully ===\n", script))
}

stages <- list(
    a = list(script = "12a_core_snp.R", desc = "Core SNP Analysis"),
    b = list(script = "12b_mash_dist.R", desc = "Mash Distance"),
    c = list(script = "12c_panaroo.R", desc = "Panaroo Pangenome"),
    d = list(script = "12d_sv_analysis.R", desc = "Structural Variants"),
    e = list(script = "12e_report.R", desc = "Reporting")
)

to_run <- if (opt$stage == "all") names(stages) else tolower(opt$stage)

for (s in to_run) {
    if (s %in% names(stages)) {
        run_step(stages[[s]]$script, stages[[s]]$desc)
    } else {
        warning(sprintf("Unknown stage: %s", s))
    }
}

cat("\nAll requested stages completed.\n")
