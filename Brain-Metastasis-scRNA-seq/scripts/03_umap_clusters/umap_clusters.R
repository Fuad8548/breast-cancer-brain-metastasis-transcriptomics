# ==============================================================================
# Publication-Grade UMAP Upgrade: Direct Patient Labels on Map + Composition Bar
# ==============================================================================
source("config.R")
suppressPackageStartupMessages({
  library(Seurat)
  library(tidyverse)
  library(patchwork)
  library(ggrepel)
})

bcbm_seurat <- readRDS(file.path(DIR_DATA_SER, "03_clustered_seurat.rds"))
reduction_use <- ifelse(
  "harmony" %in% names(bcbm_seurat@reductions),
  "harmony",
  "pca"
)

# 1. Re-cluster at Macro-Resolution (0.2 -> ~8-10 clean clusters)
bcbm_seurat <- FindNeighbors(
  bcbm_seurat,
  reduction = reduction_use,
  dims = 1:20,
  verbose = FALSE
)
bcbm_seurat <- FindClusters(bcbm_seurat, resolution = 0.2, verbose = FALSE)

# 2. Extract UMAP coordinates & clean Sample IDs (BN1, BT, etc.)
umap_df <- cbind(Embeddings(bcbm_seurat, "umap"), bcbm_seurat@meta.data)
umap_df$Macro_Cluster <- factor(Idents(bcbm_seurat))
umap_df$Sample_Clean <- gsub("^GSM[0-9]+_", "", umap_df$sample_id)

# Anchor each patient label to the median of its LARGEST cluster (never empty space!)
sample_centroids <- umap_df %>%
  group_by(Sample_Clean, tissue_type, Macro_Cluster) %>%
  summarise(
    n_cells = n(),
    umap_1 = median(umap_1),
    umap_2 = median(umap_2),
    .groups = "drop"
  ) |>
  group_by(Sample_Clean) |>
  slice_max(order_by = n_cells, n = 1, with_ties = FALSE) |>
  ungroup()

# 3. Panel A: Broad Macro-Clusters
p_macro <- DimPlot(
  bcbm_seurat,
  reduction = "umap",
  label = TRUE,
  label.box = TRUE,
  repel = TRUE,
  pt.size = 0.3
) +
  ggtitle("A. Broad Cellular Clusters (Res = 0.2)") +
  theme_bw(base_size = 11) +
  theme(legend.position = "none", plot.title = element_text(face = "bold"))

# 4. Panel B: Every badge now lands directly on a real cell cluster!
p_direct_labels <- ggplot(
  umap_df,
  aes(x = umap_1, y = umap_2, color = tissue_type)
) +
  geom_point(size = 0.35, alpha = 0.55) +
  scale_color_manual(
    name = NULL,
    values = c(
      "Adjacent_Brain_Tissue" = "#E41A1C",
      "Brain_Metastasis_Tumor" = "#1F78B4"
    ),
    labels = c(
      "Adjacent Normal Brain (BN1, BN2, DN, EN)",
      "Brain Metastasis Tumor (BT, DT, ET, FT, GT)"
    )
  ) +
  geom_label_repel(
    data = sample_centroids,
    aes(x = umap_1, y = umap_2, label = Sample_Clean, fill = tissue_type),
    color = "white",
    fontface = "bold",
    size = 3.5,
    box.padding = 0.45,
    point.padding = 0.25,
    min.segment.length = 0,
    segment.color = "grey25",
    show.legend = FALSE
  ) +
  scale_fill_manual(
    values = c(
      "Adjacent_Brain_Tissue" = "#B2182B",
      "Brain_Metastasis_Tumor" = "#08519C"
    )
  ) +
  guides(color = guide_legend(override.aes = list(size = 3.5, alpha = 1))) +
  theme_bw(base_size = 11) +
  labs(
    title = "B. Dominant Sample Localization (Red = Normal | Blue = Tumor)",
    x = "UMAP 1",
    y = "UMAP 2"
  ) +
  theme(legend.position = "bottom", plot.title = element_text(face = "bold"))

# 5. Panel C: Exact Tissue Proportion per Cluster (Proves Tumor vs. Host Islands)
prop_df <- umap_df %>%
  count(Macro_Cluster, tissue_type) %>%
  group_by(Macro_Cluster) %>%
  mutate(Pct = n / sum(n) * 100)

p_prop <- ggplot(prop_df, aes(x = Macro_Cluster, y = Pct, fill = tissue_type)) +
  geom_col(width = 0.72) +
  scale_fill_manual(
    values = c(
      "Adjacent_Brain_Tissue" = "#E41A1C",
      "Brain_Metastasis_Tumor" = "#1F78B4"
    ),
    labels = c("Normal Brain", "Brain Tumor")
  ) +
  coord_flip() +
  theme_bw(base_size = 11) +
  labs(
    title = "C. Tissue Composition per Cluster",
    x = "Macro Cluster",
    y = "Percentage of Cells (%)",
    fill = NULL
  ) +
  theme(legend.position = "bottom", plot.title = element_text(face = "bold"))

# 6. Combine & Save 3-Panel Publication Figure
final_umap_fig <- (p_macro | p_direct_labels | p_prop) +
  plot_layout(widths = c(1, 1.35, 0.85))
ggsave(
  file.path(DIR_OUT_FIG, "03_umap_direct_labels_and_composition.png"),
  plot = final_umap_fig,
  width = 17,
  height = 6,
  dpi = 300,
  bg = "white"
)
cat("[SUCCESS] Saved '03_umap_direct_labels_and_composition.png'!\n")
