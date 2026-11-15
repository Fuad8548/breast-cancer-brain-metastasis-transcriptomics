# Milestone 1: Data Acquisition
library(GEOquery)
library(tidyverse)

if (!exists("gse_data")) {
    gse <- getGEO("GSE25066", GSEMatrix = TRUE)
    gse_data <- gse[[1]]
}
metadata <- pData(gse_data)

cat("--- Scanning columns for treatment outcomes ---\n")

# Create a clean, un-flooded summary table of all characteristics columns
summary_list <- list()

for (col in colnames(metadata)) {
    if (startsWith(col, "characteristics_ch1")) {
        # Get a sample value from this column to see what it contains
        sample_val <- head(metadata[[col]], 1)
        cat(paste0("Column [", col, "] contains data like: ", sample_val, "\n"))
    }
}

# Let's save a clean version of just the clinical traits to an Excel-readable CSV
clinical_subset <- metadata %>%
    select(geo_accession, title, starts_with("characteristics_ch1"))

write.csv(clinical_subset, "data/GSE25066_Clinical_Data_Look.csv", row.names = FALSE)
cat("\n[SUCCESS] Saved clean clinical matrix to 'data/GSE25066_Clinical_Data_Look.csv'!\n")
