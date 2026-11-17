# ==============================================================================
# Pipeline Phase 4 (Complete Unified Script):
# Probe Translation, Lineage Annotation & Chemoresistance Projection
# Script: scripts/04_chemoresistance_projection.R
# Bridges: GSE25066 (2,938 Affymetrix U133A DEGs & Top RF Drivers) -> GSE268662
# ==============================================================================

cat(
  "[===] Initializing Phase 4: Unified Annotation & Chemoresistance Projection...\n"
)

source("config.R")
suppressPackageStartupMessages({
  library(Seurat)
  library(tidyverse)
  library(patchwork)
  library(ggrepel)
  library(AnnotationDbi)
  library(hgu133a.db)
})

# ------------------------------------------------------------------------------
# STEP 1: Load Clustered Seurat Object (Phase 3) & Enforce 11 Macro-Clusters (0-10)
# ------------------------------------------------------------------------------
input_rds <- file.path(DIR_DATA_SER, "03_clustered_seurat.rds")
if (!file.exists(input_rds)) {
  stop(
    "[!] Fatal Error: '03_clustered_seurat.rds' not found. Run Phase 3 first!"
  )
}
bcbm <- readRDS(input_rds)
bcbm <- JoinLayers(bcbm)

if ("RNA_snn_res.0.2" %in% colnames(bcbm@meta.data)) {
  Idents(bcbm) <- "RNA_snn_res.0.2"
  bcbm$seurat_clusters <- bcbm$RNA_snn_res.0.2
} else {
  reduction_use <- ifelse(
    "harmony" %in% names(bcbm@reductions),
    "harmony",
    "pca"
  )
  bcbm <- FindNeighbors(
    bcbm,
    reduction = reduction_use,
    dims = 1:20,
    verbose = FALSE
  )
  bcbm <- FindClusters(bcbm, resolution = 0.2, verbose = FALSE)
}
n_clusters <- length(unique(Idents(bcbm)))
cat(sprintf(
  "[+] Active Macro-Clusters: %d (IDs: 0 to %d)\n",
  n_clusters,
  n_clusters - 1
))

# ------------------------------------------------------------------------------
# STEP 2: Canonical Marker DotPlot & Biological Cell-Type Annotation
# ------------------------------------------------------------------------------
lineage_markers <- c(
  "EPCAM",
  "KRT8",
  "KRT18",
  "KRT19", # Epithelial / Metastatic Breast Tumor
  "PTPRC",
  "C1QA",
  "AIF1",
  "P2RY12", # Microglia / Tumor-Associated Macrophages
  "CD3D",
  "CD3E",
  "NKG7", # T-cells / NK cells
  "PECAM1",
  "VWF",
  "RGS5", # Endothelial / Pericytes (BBB)
  "MBP",
  "PLP1",
  "GFAP",
  "AQP4" # Oligodendrocytes / Astrocytes
)
lineage_markers <- lineage_markers[lineage_markers %in% rownames(bcbm)]

p_dot <- DotPlot(bcbm, features = lineage_markers) +
  RotatedAxis() +
  scale_color_gradient2(
    low = "lightgrey",
    mid = "steelblue",
    high = "darkred"
  ) +
  ggtitle(sprintf(
    "GSE268662 Canonical Cell-Type Markers Across %d Macro-Clusters",
    n_clusters
  ))

ggsave(
  file.path(DIR_OUT_FIG, "04_lineage_markers_dotplot.png"),
  p_dot,
  width = 12,
  height = 5.5,
  dpi = 300,
  bg = "white"
)

# Assign clean biological names & marker-proof labels to Clusters 0-10
cluster_annotations <- c(
  "0" = "T_NK_Lymphocytes",
  "1" = "BBB_Endothelial_1",
  "2" = "Oligodendrocytes",
  "3" = "Malignant_Epithelial_1",
  "4" = "Microglia_TAMs",
  "5" = "Pericytes_Vascular",
  "6" = "Resident_Myeloid_Glia",
  "7" = "Malignant_Plastic_EMT",
  "8" = "Malignant_Plastic_EMT",
  "9" = "Malignant_Epithelial_2",
  "10" = "BBB_Endothelial_2"
)

marker_proof_labels <- c(
  "Malignant_Epithelial_1" = "Malignant Epithelial 1\n(EPCAM+, KRT18+, KRT19+)",
  "Malignant_Epithelial_2" = "Malignant Epithelial 2\n(KRT18+, KRT8+)",
  "Malignant_Plastic_EMT" = "Malignant Plastic / EMT\n(Low-KRT Tumor State)",
  "Microglia_TAMs" = "Microglia / TAMs\n(C1QA+, AIF1+, P2RY12+)",
  "Resident_Myeloid_Glia" = "Resident Myeloid / Glia\n(PTPRC+, AIF1+)",
  "T_NK_Lymphocytes" = "T & NK Lymphocytes\n(CD3D+, CD3E+, NKG7+)",
  "BBB_Endothelial_1" = "BBB Endothelial 1\n(PECAM1+, VWF+)",
  "BBB_Endothelial_2" = "BBB Capillary (Normal)\n(PECAM1+, VWF+)",
  "Pericytes_Vascular" = "Vascular Pericytes\n(RGS5+)",
  "Oligodendrocytes" = "Oligodendrocytes (CNS)\n(MBP+, PLP1+)"
)

bcbm$cell_type <- unname(cluster_annotations[as.character(Idents(bcbm))])
bcbm$cell_type_with_markers <- unname(marker_proof_labels[bcbm$cell_type])
Idents(bcbm) <- "cell_type"

# Export Self-Explanatory Annotated UMAP (Cell Type + Canonical Marker Proof)
umap_df <- cbind(Embeddings(bcbm, "umap"), bcbm@meta.data)
centroids <- umap_df |>
  group_by(cell_type_with_markers) |>
  summarise(umap_1 = median(umap_1), umap_2 = median(umap_2), .groups = "drop")

p_annot_proof <- ggplot(
  umap_df,
  aes(x = umap_1, y = umap_2, color = cell_type)
) +
  geom_point(size = 0.4, alpha = 0.65) +
  geom_label_repel(
    data = centroids,
    aes(x = umap_1, y = umap_2, label = cell_type_with_markers),
    fill = "white",
    color = "black",
    fontface = "bold",
    size = 3.1,
    box.padding = 0.6,
    point.padding = 0.3,
    force = 5,
    max.overlaps = Inf,
    min.segment.length = 0,
    segment.color = "grey25",
    inherit.aes = FALSE
  ) +
  theme_bw(base_size = 11) +
  labs(
    title = "GSE268662 Annotated Cellular Ecosystem & Canonical Marker Proof (13,186 Cells)",
    x = "UMAP 1",
    y = "UMAP 2"
  ) +
  theme(legend.position = "none", plot.title = element_text(face = "bold"))

ggsave(
  file.path(DIR_OUT_FIG, "04B_annotated_umap_with_markers.png"),
  p_annot_proof,
  width = 12,
  height = 7.5,
  dpi = 300,
  bg = "white"
)

# ------------------------------------------------------------------------------
# STEP 3: Translate GSE25066 Probes -> Gene Symbols & Auto-Verify logFC Direction
# ------------------------------------------------------------------------------
sig_dir <- file.path("data", "signatures")
deg_files <- list.files(sig_dir, pattern = "\\.csv$", full.names = TRUE)
if (length(deg_files) == 0) {
  stop(
    "[!] Place your GSE25066 limma DEG CSV file inside 'data/signatures/' first!"
  )
}

deg_df <- read.csv(deg_files, row.names = 1, check.names = FALSE)
deg_df$Probe_ID <- rownames(deg_df)
deg_df$Gene_Symbol <- mapIds(
  hgu133a.db,
  keys = deg_df$Probe_ID,
  column = "SYMBOL",
  keytype = "PROBEID",
  multiVals = "first"
)

deg_clean <- deg_df |>
  filter(!is.na(Gene_Symbol) & Gene_Symbol %in% rownames(bcbm)) |>
  arrange(adj.P.Val)

write.csv(
  deg_clean,
  file.path(DIR_OUT_TAB, "04_GSE25066_translated_gene_bridge.csv"),
  row.names = FALSE
)

# Auto-detect logFC sign of Top RF Chemoresistance Drivers (TMSB15A, TTLL4, EPHB3)
rf_check <- deg_clean |> filter(Gene_Symbol %in% c("TMSB15A", "TTLL4", "EPHB3"))
rd_is_positive <- median(rf_check$logFC, na.rm = TRUE) > 0

if (rd_is_positive) {
  resistant_genes <- unique(deg_clean |> filter(logFC > 0) |> pull(Gene_Symbol))
  sensitive_genes <- unique(deg_clean |> filter(logFC < 0) |> pull(Gene_Symbol))
} else {
  resistant_genes <- unique(deg_clean |> filter(logFC < 0) |> pull(Gene_Symbol))
  sensitive_genes <- unique(deg_clean |> filter(logFC > 0) |> pull(Gene_Symbol))
}

# ------------------------------------------------------------------------------
# STEP 4: Compute Single-Cell Module Scores & Standardize (Z-Score)
# ------------------------------------------------------------------------------
bcbm <- AddModuleScore(
  object = bcbm,
  features = list(resistant_genes, sensitive_genes),
  name = c("Verified_Resistant_Score", "Verified_Sensitive_Score"),
  seed = 42
)
bcbm$Resistance_Z <- as.numeric(scale(bcbm$Verified_Resistant_Score1))
bcbm$Sensitivity_Z <- as.numeric(scale(bcbm$Verified_Sensitive_Score2))

# Figure: Standardized Chemoresistance Score on UMAP + Named Cell-Type Boxplot
p_umap_z <- FeaturePlot(
  bcbm,
  features = "Resistance_Z",
  pt.size = 0.4,
  order = TRUE
) +
  scale_color_viridis_c(option = "inferno", name = "Resistance\n(Z-Score)") +
  ggtitle("A. Primary Chemoresistance Signature on UMAP (Bright = High)") +
  theme_bw(base_size = 11) +
  theme(plot.title = element_text(face = "bold"))

plot_df <- bcbm@meta.data |>
  group_by(cell_type, tissue_type) |>
  filter(n() >= 20) |>
  ungroup()

p_box_z <- ggplot(
  plot_df,
  aes(
    x = reorder(cell_type, Resistance_Z, FUN = median),
    y = Resistance_Z,
    fill = tissue_type
  )
) +
  geom_boxplot(outlier.size = 0.25) +
  geom_hline(
    yintercept = 0,
    linetype = "dashed",
    color = "black",
    linewidth = 0.6
  ) +
  scale_fill_manual(
    name = "Tissue Origin",
    values = c(
      "Adjacent_Brain_Tissue" = "#E41A1C",
      "Brain_Metastasis_Tumor" = "#1F78B4"
    ),
    labels = c("Adjacent Normal Brain", "Brain Metastasis Tumor")
  ) +
  coord_flip() +
  theme_bw(base_size = 11) +
  labs(
    title = "B. Chemoresistance Score Ranked by Annotated Lineage",
    subtitle = "0 = Cohort Average | Above 0 = Enriched | Below 0 = Depleted",
    x = NULL,
    y = "Standardized Chemoresistance Score (Z-Score)"
  ) +
  theme(legend.position = "bottom", plot.title = element_text(face = "bold"))

ggsave(
  file.path(DIR_OUT_FIG, "04_chemoresistance_by_celltype.png"),
  p_umap_z + p_box_z,
  width = 15.5,
  height = 6.5,
  dpi = 300,
  bg = "white"
)

# Figure: High-Contrast 2x2 Grid of Top Random Forest Gini Drivers
ml_drivers <- c("TMSB15A", "TTLL4", "EPHB3", "PGK1")
driver_titles <- c(
  "TMSB15A" = "#1 RF Driver: TMSB15A (Actin Motility Switch)",
  "TTLL4" = "#2 RF Driver: TTLL4 (Tubulin Taxane-Evasion)",
  "EPHB3" = "#3 RF Driver: EPHB3 (Survival Signaling)",
  "PGK1" = "#10 RF Sentinel: PGK1 (Hypoxic Glycolysis)"
)
driver_plots <- lapply(ml_drivers[ml_drivers %in% rownames(bcbm)], function(g) {
  FeaturePlot(bcbm, features = g, pt.size = 0.35, order = TRUE) +
    scale_color_gradientn(
      colors = c("grey88", "#FDB863", "#E66101", "#B2182B"),
      name = "Expr"
    ) +
    ggtitle(driver_titles[g]) +
    theme_bw(base_size = 10) +
    theme(plot.title = element_text(face = "bold", size = 10.5))
})
ggsave(
  file.path(DIR_OUT_FIG, "04_top_ml_drivers_highcontrast.png"),
  wrap_plots(driver_plots, ncol = 2),
  width = 12,
  height = 9.5,
  dpi = 300,
  bg = "white"
)

# ------------------------------------------------------------------------------
# STEP 5: Export Named Cell-Type Summary Table & Annotated Seurat Object
# ------------------------------------------------------------------------------
celltype_table <- bcbm@meta.data |>
  group_by(Cell_Type = cell_type, Tissue_Origin = tissue_type) |>
  summarise(
    Cell_Count = n(),
    Raw_Resistant_Score = round(mean(Verified_Resistant_Score1), 4),
    Resistance_Z_Score = round(mean(Resistance_Z), 3),
    Sensitivity_Z_Score = round(mean(Sensitivity_Z), 3),
    .groups = "drop"
  ) |>
  filter(Cell_Count >= 20) |>
  arrange(desc(Resistance_Z_Score))
write.csv(
  celltype_table,
  file.path(DIR_OUT_TAB, "04_celltype_chemoresistance_scores.csv"),
  row.names = FALSE
)
saveRDS(
  bcbm,
  file = file.path(DIR_DATA_SER, "04_annotated_bcbm_seurat.rds"),
  compress = FALSE
)
cat("\n[=== PHASE 4 COMPLETE ===]\n")
print(celltype_table, n = 25)
