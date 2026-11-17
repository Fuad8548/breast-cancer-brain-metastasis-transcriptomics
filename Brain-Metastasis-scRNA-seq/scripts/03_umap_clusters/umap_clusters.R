# ==============================================================================
# Unified 2-Panel UMAP: Clusters (0-10) + Sample Badges on One Map + Composition
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

# 1. Macro-Clustering (Resolution = 0.2 -> 11 clean clusters: 0-10)
bcbm_seurat <- FindNeighbors(
  bcbm_seurat,
  reduction = reduction_use,
  dims = 1:20,
  verbose = FALSE
)
bcbm_seurat <- FindClusters(bcbm_seurat, resolution = 0.2, verbose = FALSE)

umap_df <- cbind(Embeddings(bcbm_seurat, "umap"), bcbm_seurat@meta.data)
umap_df$Macro_Cluster <- factor(Idents(bcbm_seurat))
umap_df$Sample_Clean <- gsub("^GSM[0-9]+_", "", umap_df$sample_id)

# 2. Compute cluster and sample centroids for labeling
all_badges <- bind_rows(
  cluster_centroids |>
    transmute(
      umap_1,
      umap_2,
      label = as.character(Macro_Cluster),
      badge_fill = "white",
      text_color = "black"
    ),
  sample_centroids |>
    transmute(
      umap_1,
      umap_2,
      label = Sample_Clean,
      badge_fill = ifelse(
        tissue_type == "Adjacent_Brain_Tissue",
        "#B2182B",
        "#08519C"
      ),
      text_color = "white"
    )
)

# Panel A: Single Unified UMAP (Zero dropped labels, zero collisions!)
p_unified_umap <- ggplot(
  umap_df,
  aes(x = umap_1, y = umap_2, color = tissue_type)
) +
  geom_point(size = 0.35, alpha = 0.5) +
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
    data = all_badges,
    aes(x = umap_1, y = umap_2, label = label),
    fill = all_badges$badge_fill,
    color = all_badges$text_color,
    fontface = "bold",
    size = 3.5,
    box.padding = 0.55,
    point.padding = 0.3,
    force = 4,
    max.overlaps = Inf,
    min.segment.length = 0,
    segment.color = "grey25",
    segment.size = 0.35,
    inherit.aes = FALSE
  ) +
  guides(color = guide_legend(override.aes = list(size = 3.5, alpha = 1))) +
  theme_bw(base_size = 11) +
  labs(
    title = "A. Unified Cellular Atlas: Clusters (0–10) & Patient Localization",
    x = "UMAP 1",
    y = "UMAP 2"
  ) +
  theme(legend.position = "bottom", plot.title = element_text(face = "bold"))


# 5. Panel B: Tissue Composition (%) per Cluster
prop_df <- umap_df |>
  count(Macro_Cluster, tissue_type) |>
  group_by(Macro_Cluster) |>
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
    title = "B. Tissue Composition per Cluster",
    x = "Cluster ID (0–10)",
    y = "Percentage of Cells (%)",
    fill = NULL
  ) +
  theme(legend.position = "bottom", plot.title = element_text(face = "bold"))

# 6. Combine & Save Final 2-Panel Figure
final_umap_fig <- (p_unified_umap | p_prop) + plot_layout(widths = c(1.5, 0.85))
ggsave(
  file.path(DIR_OUT_FIG, "03_umap_unified_clusters_and_samples.png"),
  plot = final_umap_fig,
  width = 14,
  height = 6.2,
  dpi = 300,
  bg = "white"
)
cat(
  "[SUCCESS] Saved unified 2-panel UMAP to 'outputs/figures/03_umap_unified_clusters_and_samples.png'!\n"
)
