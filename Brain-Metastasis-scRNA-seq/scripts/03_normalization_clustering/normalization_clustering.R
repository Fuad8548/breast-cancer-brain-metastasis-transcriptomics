#!/usr/bin/env Rscript

# ==============================================================================
# Pipeline Phase 3: Low-Memory Normalization, Harmony Integration & UMAP
# Script: 03_normalization_clustering.R
# Target Cohort: GSE268662 (13,186 cells across 9 samples - 8GB RAM Optimized)
# ==============================================================================

cat("[===] Initializing Phase 3: Lightweight Normalization & Clustering...\n")

source("config.R")
suppressPackageStartupMessages({
    library(Seurat)
    library(tidyverse)
    library(patchwork)
})

has_harmony <- requireNamespace("harmony", quietly = TRUE)

# 1. Load Filtered Merged Object
input_rds <- file.path(DIR_DATA_SER, "01_filtered_merged_seurat.rds")
bcbm_seurat <- readRDS(input_rds)

# Annotate Tissue Type from Sample ID (_T = Tumor Lesion, _N = Adjacent/Normal)
bcbm_seurat$tissue_type <- ifelse(grepl("T$", bcbm_seurat$sample_id), "Brain_Metastasis_Tumor", "Adjacent_Brain_Tissue")
cat("[-] Tissue breakdown:\n")
print(table(bcbm_seurat$tissue_type))

# 2. Sparse Log-Normalization & Top 2,000 Variable Features
cat("[-] Running sparse LogNormalize and FindVariableFeatures...\n")
bcbm_seurat <- NormalizeData(bcbm_seurat, normalization.method = "LogNormalize", scale.factor = 10000, verbose = FALSE)
bcbm_seurat <- FindVariableFeatures(bcbm_seurat, selection.method = "vst", nfeatures = 2000, verbose = FALSE)

# 3. Scale ONLY Variable Features (8GB RAM Safe)
cat("[-] Scaling top 2,000 variable features...\n")
bcbm_seurat <- ScaleData(bcbm_seurat, features = VariableFeatures(bcbm_seurat), vars.to.regress = "percent.mt", verbose = FALSE)
gc()

# 4. Run PCA
cat("[-] Running PCA...\n")
bcbm_seurat <- RunPCA(bcbm_seurat, features = VariableFeatures(bcbm_seurat), npcs = TARGET_PCS, verbose = FALSE)

# 5. Batch Correction (Harmony) & UMAP Dimensionality Reduction
if (has_harmony) {
    cat("[-] Running Harmony batch correction across 9 samples...\n")
    library(harmony)
    bcbm_seurat <- RunHarmony(bcbm_seurat, group.by.vars = "sample_id", verbose = FALSE)
    reduction_use <- "harmony"
} else {
    cat("[!] 'harmony' package not installed; using standard PCA UMAP.\n")
    reduction_use <- "pca"
}

cat("[-] Computing SNN graph, clusters, and UMAP...\n")
bcbm_seurat <- FindNeighbors(bcbm_seurat, reduction = reduction_use, dims = 1:20, verbose = FALSE)
bcbm_seurat <- FindClusters(bcbm_seurat, resolution = 0.5, verbose = FALSE)
bcbm_seurat <- RunUMAP(bcbm_seurat, reduction = reduction_use, dims = 1:20, verbose = FALSE)

# 6. Export Publication-Ready UMAP Figures
p_cluster <- DimPlot(bcbm_seurat, reduction = "umap", label = TRUE, pt.size = 0.3) +
    ggtitle("GSE268662 Cell Clusters")

p_tissue <- DimPlot(bcbm_seurat, reduction = "umap", group.by = "tissue_type", pt.size = 0.3) +
    ggtitle("Tumor vs. Adjacent Brain Tissue")

p_sample <- DimPlot(bcbm_seurat, reduction = "umap", group.by = "sample_id", pt.size = 0.3) +
    ggtitle("Sample Integration (9 Samples)")

png(file.path(DIR_OUT_FIG, "03_umap_clusters_and_tissue.png"), width = 1600, height = 550, res = 110)
print(p_cluster + p_tissue + p_sample)
dev.off()

# 7. Save Clustered Object (Uncompressed for Speed)
saveRDS(bcbm_seurat, file = file.path(DIR_DATA_SER, "03_clustered_seurat.rds"), compress = FALSE)

cat("\n[=== PHASE 3 COMPLETE ===]\n")
cat(sprintf(
    "[+] Generated %d clusters across %d cells. Figure saved to outputs/figures/.\n",
    length(unique(Idents(bcbm_seurat))), ncol(bcbm_seurat)
))
