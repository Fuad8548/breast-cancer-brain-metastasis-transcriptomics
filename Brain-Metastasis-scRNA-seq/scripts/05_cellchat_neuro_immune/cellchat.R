# ==============================================================================
# Pipeline Phase 5: Tumor-to-Brain Niche Ligand-Receptor Modeling (CellChat)
# Script: scripts/05_cellchat_neuro_immune.R
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

# 1. Load Annotated Seurat Object Directly from Milestone 04
bcbm <- readRDS(file.path(DIR_DATA_SER, "04_annotated_bcbm_seurat.rds"))
Idents(bcbm) <- "cell_type"

# 2. Subset to Core Interacting Tumor, Immune & BBB Populations (8GB RAM Safe)
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
bcbm_sub <- subset(bcbm, subset = cell_type %in% target_populations)
bcbm_sub$cell_type <- droplevels(factor(bcbm_sub$cell_type))
rm(bcbm)
gc()

# 3. Initialize CellChat & Load Secreted / Cell-Cell Contact Databases
data_input <- GetAssayData(bcbm_sub, assay = "RNA", layer = "data")
meta_input <- bcbm_sub@meta.data |> dplyr::select(cell_type)

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

# 4. Run Ligand-Receptor Statistical Inference
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

# 5. Export Full Interactome & Ranked Tumor-to-Niche Signaling Tables
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

tumor_to_niche <- df_net |>
  filter(source %in% tumor_sources, target %in% immune_targets) |>
  arrange(pval, desc(prob))
write.csv(
  tumor_to_niche,
  file.path(DIR_OUT_TAB, "05_tumor_to_neuro_immune_signaling.csv"),
  row.names = FALSE
)

# 6. Export Publication Figures (with margin = 0.25 for clean circle labels)
png(
  file.path(DIR_OUT_FIG, "05_cellchat_tumor_circle.png"),
  width = 1000,
  height = 900,
  res = 120,
  bg = "white"
)
group_size <- as.numeric(table(cellchat@idents))
sources_idx <- which(levels(cellchat@idents) %in% tumor_sources)
netVisual_circle(
  cellchat@net$weight,
  vertex.weight = group_size,
  weight.scale = TRUE,
  sources.use = sources_idx,
  label.edge = FALSE,
  margin = 0.25,
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

cat("\n[=== PHASE 5 COMPLETE ===]\n")
