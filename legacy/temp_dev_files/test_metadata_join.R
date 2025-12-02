source("R/wgs_helpers.R")

# Mock options
opt <- list(pids = "ALL")

# Run discovery
samples <- discover_samples()

print(paste("Total samples found:", nrow(samples)))
if (nrow(samples) > 0) {
    print("First few samples:")
    print(head(samples))
}

# Check for PR0010
pr_samples <- samples %>% filter(grepl("PR0010", SampleID))

if (nrow(pr_samples) > 0) {
    print("Checking PR0010 samples...")
    print(head(pr_samples %>% select(SampleID, Participant_id, Timepoint)))

    # Check if Participant_id is "PR0010" (BAD) or something else (GOOD)
    bad_count <- sum(pr_samples$Participant_id == "PR0010")
    good_count <- sum(pr_samples$Participant_id != "PR0010" & !is.na(pr_samples$Participant_id))

    cat(sprintf("\nTotal PR0010 batch samples: %d\n", nrow(pr_samples)))
    cat(sprintf("Samples with Participant_id = 'PR0010' (BAD): %d\n", bad_count))
    cat(sprintf("Samples with Participant_id != 'PR0010' (GOOD): %d\n", good_count))

    if (bad_count > 0) {
        stop("TEST FAILED: Some samples still have batch ID as Participant_id!")
    } else {
        print("TEST PASSED: All PR0010 batch samples have correct Participant IDs.")
    }
} else {
    print("No PR0010 samples found to test.")
}
