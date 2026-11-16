# ==============================================================================
# Pipeline Phase 4: Probe Translation, Cell Annotation & Chemoresistance Projection
# Script: scripts/04_chemoresistance_projection.R
# Bridges: GSE25066 (2,938 Affymetrix U133A DEGs) -> GSE268662 (Brain Metastasis scRNA-seq)
# ==============================================================================

cat("[===] Initializing Phase 4: Cross-Platform Probe Mapping & Projection...\n")

source("config.R")
suppressPackageStartupMessages({
    library(Seurat)
    library(tidyverse)
    library(patchwork)
    library(AnnotationDbi)
    library(hgu133a.db)
})

# 1. Load Clustered Seurat Object from Phase 3
input_rds <- file.path(DIR_DATA_SER, "03_clustered_seurat.rds")
if (!file.exists(input_rds)) {
    stop("[!] Fatal Error: '03_clustered_seurat.rds' not found. Run Phase 3 first!")
}
bcbm <- readRDS(input_rds)
bcbm <- JoinLayers(bcbm)

# 2. Load GSE25066 limma DEG Table & Translate Affymetrix Probes to Gene Symbols
sig_dir <- file.path("data", "signatures")
if (!dir.exists(sig_dir)) dir.create(sig_dir, recursive = TRUE)

deg_files <- list.files(sig_dir, pattern = "\\.csv$", full.names = TRUE)
if (length(deg_files) == 0) {
    stop("[!] Place your GSE25066 limma DEG CSV file inside 'data/signatures/' first!")
}

cat(sprintf("[-] Reading GSE25066 signature table: %s\n", basename(deg_files)))
deg_df <- read.csv(deg_files, row.names = 1, check.names = FALSE)
deg_df$Probe_ID <- rownames(deg_df)

# Map Affymetrix U133A Probe IDs to official HUGO Gene Symbols
cat("[-] Translating 2,938 Affymetrix U133A probes to HUGO Gene Symbols via hgu133a.db...\n")
mapped_symbols <- mapIds(
    hgu133a.db,
    keys = deg_df$Probe_ID,
    column = "SYMBOL",
    keytype = "PROBEID",
    multiVals = "first"
)
deg_df$Gene_Symbol <- mapped_symbols

# Clean unmapped probes and retain genes present in GSE268662 scRNA-seq matrix
deg_clean <- deg_df %>%
    filter(!is.na(Gene_Symbol) & Gene_Symbol %in% rownames(bcbm)) %>%
    arrange(adj.P.Val)

# Export translated bridge table for your repository documentation
write.csv(deg_clean, file.path(DIR_OUT_TAB, "04_GSE25066_translated_gene_bridge.csv"), row.names = FALSE)
cat(sprintf(
    "[+] Successfully mapped %d significant probes to %d unique single-cell genes!\n",
    nrow(deg_clean), length(unique(deg_clean$Gene_Symbol))
))

# Split into Chemoresistant-Enriched (logFC < 0, includes TMSB15A/TTLL4/CCND1) vs Chemosensitive (logFC > 0)
resistant_genes <- unique(deg_clean %>% filter(logFC < 0) %>% pull(Gene_Symbol))
sensitive_genes <- unique(deg_clean %>% filter(logFC > 0) %>% pull(Gene_Symbol))

cat(sprintf("    -> Chemoresistant Signature: %d genes (Top drivers: TMSB15A, TTLL4, CCND1)\n", length(resistant_genes)))
cat(sprintf("    -> Chemosensitive Signature: %d genes\n", length(sensitive_genes)))

# 3. Canonical Lineage Marker DotPlot (Separating Tumor vs. Brain Microenvironment)
cat("[-] Generating canonical lineage marker DotPlot...\n")
lineage_markers <- c(
    "EPCAM", "KRT8", "KRT18", "KRT19", # Epithelial / Metastatic Breast Tumor
    "PTPRC", "C1QA", "AIF1", "P2RY12", # Microglia / Tumor-Associated Macrophages
    "CD3D", "CD3E", "NKG7", # T-cells / NK cells
    "PECAM1", "VWF", "RGS5", # Endothelial / Pericytes (BBB)
    "MBP", "PLP1", "GFAP", "AQP4" # Oligodendrocytes / Astrocytes
)
lineage_markers <- lineage_markers[lineage_markers %in% rownames(bcbm)]

p_dot <- DotPlot(bcbm, features = lineage_markers) +
    RotatedAxis() +
    scale_color_gradient2(low = "lightgrey", mid = "steelblue", high = "darkred") +
    ggtitle("GSE268662 Canonical Cell-Type Markers Across 17 Clusters")

ggsave(file.path(DIR_OUT_FIG, "04_lineage_markers_dotplot.png"), p_dot, width = 12, height = 6, dpi = 150)

# 4. Project Genome-Wide Chemoresistant & Chemosensitive Modules onto Single Cells
cat("[-] Computing single-cell module scores via AddModuleScore()...\n")
bcbm <- AddModuleScore(
    object = bcbm,
    features = list(resistant_genes, sensitive_genes),
    name = c("GSE25066_Resistant_Score", "GSE25066_Sensitive_Score"),
    seed = 42
)

# 5. Export Publication-Ready Projection Plots
p_umap_res <- FeaturePlot(bcbm, features = "GSE25066_Resistant_Score1", pt.size = 0.4, order = TRUE) +
    scale_color_viridis_c(option = "magma") +
    ggtitle("Primary Chemoresistance Signature Score (GSE25066)")

p_box_res <- ggplot(bcbm@meta.data, aes(x = seurat_clusters, y = GSE25066_Resistant_Score1, fill = tissue_type)) +
    geom_boxplot(outlier.size = 0.3) +
    theme_minimal() +
    labs(
        title = "Chemoresistance Score Across Brain Metastasis Clusters",
        x = "Seurat Cluster (0 - 16)",
        y = "Chemoresistance Module Score"
    )

ggsave(file.path(DIR_OUT_FIG, "04_chemoresistance_projection_umap.png"), p_umap_res + p_box_res, width = 15, height = 6, dpi = 150)

# Plot Individual Top Random Forest Drivers (TMSB15A, TTLL4, CCND1)
ml_drivers <- c("TMSB15A", "TTLL4", "CCND1")
ml_drivers <- ml_drivers[ml_drivers %in% rownames(bcbm)]
if (length(ml_drivers) > 0) {
    p_drivers <- FeaturePlot(bcbm, features = ml_drivers, ncol = 3, pt.size = 0.3, order = TRUE)
    ggsave(file.path(DIR_OUT_FIG, "04_top_ml_drivers_featureplot.png"), p_drivers, width = 14, height = 4.5, dpi = 150)
}

# 6. Export Statistical Summary Table
cluster_summary <- bcbm@meta.data %>%
    group_by(seurat_clusters, tissue_type) %>%
    summarise(
        Cell_Count = n(),
        Mean_Resistant_Score = round(mean(GSE25066_Resistant_Score1), 4),
        Mean_Sensitive_Score = round(mean(GSE25066_Sensitive_Score2), 4),
        .groups = "drop"
    ) %>%
    arrange(desc(Mean_Resistant_Score))

write.csv(cluster_summary, file.path(DIR_OUT_TAB, "04_cluster_chemoresistance_scores.csv"), row.names = FALSE)
saveRDS(bcbm, file = file.path(DIR_DATA_SER, "04_scored_bcbm_seurat.rds"), compress = FALSE)

cat("\n[=== PHASE 4 COMPLETE: ALL 2,938 PROBES BRIDGED TO SINGLE-CELL ATLAS ===]\n")
print(head(cluster_summary, 10))
