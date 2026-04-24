library(readxl)

files <- c(
  "data/inputs/list Ecoli YELLOW sequencing - batch 1 - uitgebreid.xlsx",
  "data/inputs/list Ecoli YELLOW sequencing - batch 2 - uitgebreid - correctie.xlsx",
  "data/inputs/Lijst Ecoli YELLOW sequencing - batch 3 - UPDATE (definitief).xlsx"
)

for (f in files) {
  cat("\n========================================\n")
  cat("FILE:", basename(f), "\n")
  if (file.exists(f)) {
    df <- read_excel(f, n_max = 50)
    
    # Print all column names to see what we're dealing with
    cat("ALL COLUMNS:\n")
    print(names(df))
    
    cat("\nMATCHING COLUMNS (time/tijd/point/tp):\n")
    cols <- names(df)
    tp_cols <- cols[grepl("time|tp|tijd|point", tolower(cols))]
    print(tp_cols)
    
    if (length(tp_cols) > 0) {
      col <- tp_cols[1]
      cat("\nFIRST FEW VALUES IN '", col, "':\n")
      print(head(as.character(df[[col]]), 15))
    }
  } else {
    cat("File not found!\n")
  }
}
