library(tidyverse)
library(lubridate)

df <- read_csv("assembly_metadata.csv", show_col_types = FALSE)

# Convert date
df <- df %>%
    mutate(date_parsed = dmy(Collection_Date))

# Check distribution
cat("Latest dates in metadata:\n")
print(head(sort(df$date_parsed, decreasing = TRUE), 20))

# Check Batch for late dates (e.g. after Sept 2024)
late_samples <- df %>%
    filter(date_parsed > ymd("2024-09-01"))

if (nrow(late_samples) > 0) {
    cat("\nBatch labels for samples after Sept 1, 2024:\n")
    print(table(late_samples$Batch, useNA = "ifany"))

    cat("\nFirst few late samples:\n")
    print(late_samples %>% select(Isolate_ID, Collection_Date, Batch) %>% head())
} else {
    cat("\nNo samples found after Sept 1, 2024.\n")
}
