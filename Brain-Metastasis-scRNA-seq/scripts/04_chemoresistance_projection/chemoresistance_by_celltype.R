source("config.R")
suppressPackageStartupMessages({
  library(Seurat)
  library(tidyverse)
  library(patchwork)
  library(ggrepel)
})

bcbm <- readRDS(file.path(DIR_DATA_SER, "04_annotated_bcbm_seurat.rds"))

# 1. Compute Island Centroids for Direct Cell-Type Callouts on Panel A
umap_df <- cbind(Embeddings(bcbm, "umap"), bcbm@meta.data)
centroids <- umap_df |>
  group_by(cell_type) |>
  summarise(umap_1 = median(umap_1), umap_2 = median(umap_2), .groups = "drop")

# 2. Panel A: Chemoresistance UMAP + Direct Cell-Type Indicators
p_umap_z <- ggplot(
  umap_df |> arrange(Resistance_Z),
  aes(x = umap_1, y = umap_2, color = Resistance_Z)
) +
  geom_point(size = 0.4, alpha = 0.8) +
  scale_color_viridis_c(option = "inferno", name = "Resistance\n(Z-Score)") +
  geom_label_repel(
    data = centroids,
    aes(x = umap_1, y = umap_2, label = cell_type),
    fill = alpha("white", 0.88),
    color = "black",
    fontface = "bold",
    size = 3.0,
    box.padding = 0.45,
    point.padding = 0.25,
    force = 4,
    max.overlaps = Inf,
    min.segment.length = 0,
    segment.color = "grey30",
    inherit.aes = FALSE
  ) +
  theme_bw(base_size = 11) +
  labs(
    title = "A. Primary Chemoresistance Signature on UMAP (Bright = High)",
    x = "UMAP 1",
    y = "UMAP 2"
  ) +
  theme(plot.title = element_text(face = "bold"))

# 3. Panel B: Spacious Side-by-Side Ranked Boxplot by Annotated Cell Type
plot_df <- umap_df |>
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
  geom_boxplot(outlier.size = 0.25, width = 0.7) +
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
  theme(
    legend.position = "bottom",
    axis.text.y = element_text(face = "bold", size = 10, color = "black"),
    plot.title = element_text(face = "bold")
  )

# 4. Combine Side-by-Side (17 x 7.5 in) so Every Label Has Clean Breathing Room
ggsave(
  file.path(DIR_OUT_FIG, "04_chemoresistance_by_celltype.png"),
  (p_umap_z | p_box_z) + plot_layout(widths = c(1.15, 1.15)),
  width = 17,
  height = 7.5,
  dpi = 300,
  bg = "white"
)

cat(
  "[SUCCESS] Saved spacious side-by-side '04_chemoresistance_by_celltype.png' with direct cell-type indicators!\n"
)
