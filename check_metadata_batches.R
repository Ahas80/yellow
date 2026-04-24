library(tidyverse)
# Check assembly metadata if it exists
if (file.exists("assembly_metadata.csv")) {
    df <- read_csv("assembly_metadata.csv", show_col_types = FALSE)
    print("Columns in assembly_metadata.csv:")
    print(colnames(df))

    # Look for columns that might contain batch info
    possible_cols <- grep("batch|run|date", colnames(df), ignore.case = TRUE, value = TRUE)
    if (length(possible_cols) > 0) {
        print("Possible batch columns:")
        print(possible_cols)
        for (col in possible_cols) {
            print(paste("Unique values in", col, ":"))
            print(table(df[[col]], useNA = "ifany"))
        }
    } else {
        print("No obvious batch columns found in assembly_metadata.csv")
        # Peek at Sample IDs, sometimes batch is encoded there
        print("Sample ID sample (first 10):")
        print(head(df$sample_id, 10))
    }
} else {
    print("assembly_metadata.csv not found in root")
}
