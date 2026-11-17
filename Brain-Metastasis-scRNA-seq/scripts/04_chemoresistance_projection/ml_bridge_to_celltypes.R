source("config.R")
suppressPackageStartupMessages({
  library(Seurat)
  library(tidyverse)
  library(patchwork)
  library(ggrepel)
})

bcbm <- readRDS(file.path(DIR_DATA_SER, "04_annotated_bcbm_seurat.rds"))

# 1. Assign Each Cell Type to Its Spatial Microenvironment Niche (Blue / Purple / Red)
niche_map <- c(
  "Malignant_Epithelial_1" = "Tumor-Dominated Niche",
  "Malignant_Epithelial_2" = "Tumor-Dominated Niche",
  "Malignant_Plastic_EMT" = "Tumor-Dominated Niche",
  "T_NK_Lymphocytes" = "Shared Tumor–Host Niche",
  "Microglia_TAMs" = "Shared Tumor–Host Niche",
  "Pericytes_Vascular" = "Shared Tumor–Host Niche",
  "Resident_Myeloid_Glia" = "Normal-Brain-Dominated Niche",
  "BBB_Endothelial_1" = "Normal-Brain-Dominated Niche",
  "BBB_Endothelial_2" = "Normal-Brain-Dominated Niche",
  "Oligodendrocytes" = "Normal-Brain-Dominated Niche"
)
bcbm$Spatial_Niche <- unname(niche_map[bcbm$cell_type])

# Order Y-axis cleanly from Normal Brain (bottom) -> Shared Host (middle) -> Tumor (top)
cell_order <- c(
  "Oligodendrocytes",
  "BBB_Endothelial_2",
  "BBB_Endothelial_1",
  "Resident_Myeloid_Glia",
  "Pericytes_Vascular",
  "Microglia_TAMs",
  "T_NK_Lymphocytes",
  "Malignant_Plastic_EMT",
  "Malignant_Epithelial_2",
  "Malignant_Epithelial_1"
)
bcbm$cell_type_ordered <- factor(bcbm$cell_type, levels = cell_order)

# 2. Panel A: Niche-Colored UMAP (Blue = Tumor, Purple = Shared, Red = Normal Brain)
umap_df <- cbind(Embeddings(bcbm, "umap"), bcbm@meta.data)
centroids <- umap_df |>
  group_by(cell_type_ordered, Spatial_Niche) |>
  summarise(umap_1 = median(umap_1), umap_2 = median(umap_2), .groups = "drop")

p_niche_umap <- ggplot(
  umap_df,
  aes(x = umap_1, y = umap_2, color = Spatial_Niche)
) +
  geom_point(size = 0.35, alpha = 0.6) +
  scale_color_manual(
    name = "Microenvironment Niche",
    values = c(
      "Tumor-Dominated Niche" = "#1F78B4", # Royal Blue
      "Shared Tumor–Host Niche" = "#7570B3", # Purple Mix
      "Normal-Brain-Dominated Niche" = "#E41A1C" # Coral Red
    )
  ) +
  geom_label_repel(
    data = centroids,
    aes(
      x = umap_1,
      y = umap_2,
      label = cell_type_ordered,
      fill = Spatial_Niche
    ),
    color = "white",
    fontface = "bold",
    size = 3.0,
    box.padding = 0.45,
    point.padding = 0.25,
    force = 4,
    max.overlaps = Inf,
    min.segment.length = 0,
    segment.color = "grey25",
    show.legend = FALSE
  ) +
  scale_fill_manual(
    values = c(
      "Tumor-Dominated Niche" = "#08519C",
      "Shared Tumor–Host Niche" = "#54278F",
      "Normal-Brain-Dominated Niche" = "#B2182B"
    )
  ) +
  guides(color = guide_legend(override.aes = list(size = 3.5, alpha = 1))) +
  theme_bw(base_size = 11) +
  labs(
    title = "A. Spatial Niche & Cell-Type Atlas (13,186 Cells)",
    x = "UMAP 1",
    y = "UMAP 2"
  ) +
  theme(legend.position = "bottom", plot.title = element_text(face = "bold"))

# 3. Panel B: Unified Split DotPlot Connecting Brain Lineage Markers + Breast ML Biomarkers
brain_markers <- c("EPCAM", "KRT18", "CD3D", "C1QA", "RGS5", "PECAM1", "MBP")
breast_ml_markers <- c(
  "TBC1D9",
  "DSC2",
  "GATA3",
  "TXNDC5",
  "SKP1",
  "PGK1",
  "TMSB15A",
  "TTLL4",
  "EPHB3"
)
all_genes <- c(brain_markers, breast_ml_markers)
all_genes <- all_genes[all_genes %in% rownames(bcbm)]

Idents(bcbm) <- "cell_type_ordered"
dot_data <- DotPlot(bcbm, features = all_genes)$data |>
  mutate(
    Gene_Panel = ifelse(
      features.plot %in% brain_markers,
      "1. Brain & Lineage Identity Markers",
      "2. Breast Cancer ML Biomarkers (GSE25066)"
    )
  )

p_unified_dot <- ggplot(dot_data, aes(x = features.plot, y = id)) +
  geom_point(aes(size = pct.exp, color = avg.exp.scaled)) +
  scale_color_gradient2(
    low = "grey90",
    mid = "#FDB863",
    high = "#B2182B",
    name = "Scaled\nExpr"
  ) +
  scale_size_continuous(range = c(0.5, 7), name = "% Cells\nExpressing") +
  facet_grid(~Gene_Panel, scales = "free_x", space = "free_x") +
  theme_bw(base_size = 11) +
  labs(
    title = "B. Integrated Lineage Proof & Breast ML Biomarker Bridge",
    x = NULL,
    y = NULL
  ) +
  theme(
    axis.text.x = element_text(
      angle = 45,
      hjust = 1,
      face = "bold",
      color = "black"
    ),
    axis.text.y = element_text(face = "bold", color = "black"),
    strip.background = element_rect(fill = "#2C3E50"),
    strip.text = element_text(color = "white", face = "bold", size = 10),
    plot.title = element_text(face = "bold")
  )

# 4. Combine & Save Master Bridge Figure
ggsave(
  file.path(DIR_OUT_FIG, "04_complete_master_bridge.png"),
  (p_niche_umap | p_unified_dot) + plot_layout(widths = c(1.15, 1.45)),
  width = 18,
  height = 7.2,
  dpi = 300,
  bg = "white"
)
cat(
  "[SUCCESS] Saved unified master bridge figure '04_complete_master_bridge.png'!\n"
)
