# ==============================================================================
# Synchronized Phase 4B (Clean Cell Annotation) + Phase 5 (CellChat Interactome)
# ==============================================================================
source("config.R")
suppressPackageStartupMessages({
  library(Seurat)
  library(CellChat)
  library(patchwork)
  library(tidyverse)
  library(future)
})

plan("sequential")
options(future.globals.maxSize = 2000 * 1024^2)

# 1. Load Scored Seurat Object from Phase 4 (Macro-Clusters 0-10)
bcbm <- readRDS(file.path(DIR_DATA_SER, "04_scored_bcbm_seurat.rds"))
Idents(bcbm) <- "seurat_clusters"

# 2. Map Verified Canonical Markers (image_p_1qBM.png) onto Clusters 0-10
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

bcbm$cell_type <- unname(cluster_annotations[as.character(Idents(bcbm))])
Idents(bcbm) <- "cell_type"

# Export Clean Annotated UMAP (Replaces old 17-cluster draft)
p_annot <- DimPlot(
  bcbm,
  reduction = "umap",
  label = TRUE,
  label.box = TRUE,
  repel = TRUE,
  pt.size = 0.35
) +
  ggtitle("GSE268662 Annotated Tumor & Brain Microenvironment (13,186 Cells)") +
  theme_bw(base_size = 11) +
  theme(legend.position = "right", plot.title = element_text(face = "bold"))

ggsave(
  file.path(DIR_OUT_FIG, "04B_annotated_celltypes_umap.png"),
  p_annot,
  width = 12,
  height = 6.5,
  dpi = 300,
  bg = "white"
)
saveRDS(
  bcbm,
  file = file.path(DIR_DATA_SER, "04_annotated_bcbm_seurat.rds"),
  compress = FALSE
)
cat("[+] Saved clean annotated UMAP and '04_annotated_bcbm_seurat.rds'!\n")

# 3. Subset to Core Tumor, Immune & BBB Populations for CellChat (8GB RAM Safe)
target_populations <- c(
  "Malignant_Epithelial_1",
  "Malignant_Epithelial_2",
  "Malignant_Plastic_EMT",
  "Microglia_TAMs",
  "Resident_Myeloid_Glia",
  "T_NK_Lymphocytes",
  "BBB_Endothelial_1",
  "Pericytes_Vascular"
)
target_populations <- target_populations[
  target_populations %in% unique(bcbm$cell_type)
]

bcbm_sub <- subset(bcbm, subset = cell_type %in% target_populations)
bcbm_sub$cell_type <- droplevels(factor(bcbm_sub$cell_type))
rm(bcbm)
gc()

# 4. Run CellChat Ligand-Receptor Inference
data_input <- GetAssayData(bcbm_sub, assay = "RNA", layer = "data")
meta_input <- bcbm_sub@meta.data %>% dplyr::select(cell_type)

cellchat <- createCellChat(
  object = data_input,
  meta = meta_input,
  group.by = "cell_type"
)
cellchat@DB <- subsetDB(
  CellChatDB.human,
  search = c("Secreted Signaling", "Cell-Cell Contact")
)
rm(bcbm_sub, data_input)
gc()

cellchat <- subsetData(cellchat)
cellchat <- identifyOverExpressedGenes(cellchat)
cellchat <- identifyOverExpressedInteractions(cellchat)
cellchat <- computeCommunProb(
  cellchat,
  type = "triMean",
  population.size = TRUE
)
cellchat <- filterCommunication(cellchat, min.cells = 10)
cellchat <- computeCommunProbPathway(cellchat)
cellchat <- aggregateNet(cellchat)

# 5. Export Significant Tumor -> Brain Immune/BBB Interactions (MIF-CD74, APP-CD74)
df_net <- subsetCommunication(cellchat)
write.csv(
  df_net,
  file.path(DIR_OUT_TAB, "05_cellchat_all_interactions.csv"),
  row.names = FALSE
)

tumor_sources <- c(
  "Malignant_Epithelial_1",
  "Malignant_Epithelial_2",
  "Malignant_Plastic_EMT"
)
immune_targets <- c(
  "Microglia_TAMs",
  "Resident_Myeloid_Glia",
  "T_NK_Lymphocytes",
  "BBB_Endothelial_1",
  "Pericytes_Vascular"
)

tumor_to_niche <- df_net %>%
  filter(source %in% tumor_sources, target %in% immune_targets) %>%
  arrange(pval, desc(prob))

write.csv(
  tumor_to_niche,
  file.path(DIR_OUT_TAB, "05_tumor_to_neuro_immune_signaling.csv"),
  row.names = FALSE
)

# 6. Export Publication-Grade CellChat Figures
png(
  file.path(DIR_OUT_FIG, "05_cellchat_tumor_circle.png"),
  width = 950,
  height = 850,
  res = 120,
  bg = "white"
)
groupSize <- as.numeric(table(cellchat@idents))
sources_idx <- which(levels(cellchat@idents) %in% tumor_sources)
netVisual_circle(
  cellchat@net$weight,
  vertex.weight = groupSize,
  weight.scale = TRUE,
  sources.use = sources_idx,
  label.edge = FALSE,
  title.name = "Signaling Sent by Metastatic Breast Tumor Clones to Brain Niche"
)
dev.off()

targets_idx <- which(levels(cellchat@idents) %in% immune_targets)
p_bubble <- netVisual_bubble(
  cellchat,
  sources.use = sources_idx,
  targets.use = targets_idx,
  remove.isolate = TRUE
) +
  theme_bw(base_size = 10) +
  theme(axis.text.x = element_text(angle = 45, hjust = 1, color = "black")) +
  ggtitle("Metastatic Breast Tumor -> Brain Immune & BBB Ligand-Receptor Axes")

ggsave(
  file.path(DIR_OUT_FIG, "05_cellchat_ligand_receptor_bubble.png"),
  p_bubble,
  width = 12,
  height = 8.5,
  dpi = 300,
  bg = "white"
)
saveRDS(
  cellchat,
  file = file.path(DIR_DATA_SER, "05_cellchat_bcbm.rds"),
  compress = FALSE
)

cat("\n[=== PHASE 4B & 5 COMPLETE ===]\n")
print(head(
  tumor_to_niche[, c(
    "source",
    "target",
    "ligand",
    "receptor",
    "pathway_name",
    "prob",
    "pval"
  )],
  15
))
