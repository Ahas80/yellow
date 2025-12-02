library(readr)
library(dplyr)

f <- "results/longitudinal/annotated_snps.csv"
if (!file.exists(f)) stop("File not found")

df <- read_csv(f, show_col_types = FALSE)
print(colnames(df))
print(head(df))

# Try the select
df2 <- df %>% select(Participant_id, From_Time)
print(head(df2))
