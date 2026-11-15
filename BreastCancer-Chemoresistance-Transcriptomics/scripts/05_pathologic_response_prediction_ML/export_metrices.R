# Final R Export Step for Python Integration
library(GEOquery)

cat("--- Step 1: Re-instantiating GSE25066 Data Environment ---\n")
gse <- getGEO("GSE25066", GSEMatrix = TRUE)
gse_data <- gse[[1]]

# 1. Slice expression matrix to contain ONLY our significant 2,938 genes
ml_expression_matrix <- exprs(gse_data)[rownames(significant_genes), colnames(expression_aligned)]

# 2. Transpose so patients are rows and genes are columns
ml_expression_matrix <- t(ml_expression_matrix)

# 3. Export data matrices to CSV
write.csv(ml_expression_matrix, "05_GSE25066_ML_Features.csv", row.names = TRUE)

# 4. Export binary labels (1 = Resistant/RD, 0 = Sensitive/pCR)
ml_labels <- data.frame(
  Sample_ID = rownames(meta_aligned),
  Label = ifelse(meta_aligned$Response == "Resistant (RD)", 1, 0)
)
write.csv(ml_labels, "05_GSE25066_ML_Labels.csv", row.names = FALSE)

cat("[SUCCESS] Data structures successfully exported for the Python tech stack!\n")


# Add to the end of export_metrics.R
library(AnnotationDbi)
library(hgu133a.db)

probe_map <- data.frame(
  Probe_ID = colnames(ml_expression_matrix),
  Gene_Symbol = mapIds(
    hgu133a.db,
    keys = colnames(ml_expression_matrix),
    column = "SYMBOL",
    keytype = "PROBEID",
    multiVals = "first"
  ),
  stringsAsFactors = FALSE
)
probe_map$Gene_Symbol[is.na(probe_map$Gene_Symbol)] <- probe_map$Probe_ID[is.na(probe_map$Gene_Symbol)]
write.csv(probe_map, "GSE25066_Probe_to_Gene_Map.csv", row.names = FALSE)
cat("[SUCCESS] Probe-to-Gene annotation map exported!\n")
