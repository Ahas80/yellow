df <- read.csv("results/clinical/status_map.csv")
print(colnames(df))
print(sort(unique(df$Batch)))
