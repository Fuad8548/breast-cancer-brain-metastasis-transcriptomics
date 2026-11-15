#!/usr/bin/env Rscript

# ==============================================================================
# Pipeline Phase 5: Ligand-Receptor Crosstalk Modeling (CellChat)
# Script: scripts/05_cellchat_neuro_immune.R
# Focus: Chemoresistant Breast Tumor Clones <-> Brain-Resident Immune Niche (8GB RAM Safe)
# ==============================================================================

cat("[===] Initializing Phase 5: Neuro-Immune Ligand-Receptor Analysis (CellChat)...\n")

source("config.R")
suppressPackageStartupMessages({
  library(Seurat)
  library(CellChat)
  library(patchwork)
  library(tidyverse)
  library(future)
})

# Enforce sequential execution to prevent 8GB RAM duplication crash
plan("sequential")
options(future.globals.maxSize = 2000 * 1024^2)

# 1. Load Annotated Seurat Object from Phase 4B
input_rds <- file.path(DIR_DATA_SER, "04_annotated_bcbm_seurat.rds")
if (!file.exists(input_rds)) {
  stop("[!] Fatal Error: '04_annotated_bcbm_seurat.rds' missing. Run Phase 4B first!")
}
bcbm <- readRDS(input_rds)

# 2. Subset to Key Tumor, Immune & BBB Microenvironment Populations (Saves RAM & Sharpens Signal)
target_populations <- c(
  "Tumor_Chemoresistant_C11",
  "Tumor_Chemoresistant_C13",
  "Tumor_Epithelial_LowRes",
  "Tumor_Epithelial_KRT18",
  "Microglia_TAMs_C1QA",
  "Myeloid_AIF1",
  "T_Cells_CD3D",
  "NK_Cytotoxic_T_Cells",
  "BBB_Endothelial_1"
)
target_populations <- target_populations[target_populations %in% unique(bcbm$cell_type)]

cat("[-] Subsetting to core tumor-immune-BBB interacting populations...\n")
bcbm_sub <- subset(bcbm, subset = cell_type %in% target_populations)
bcbm_sub$cell_type <- droplevels(factor(bcbm_sub$cell_type))

rm(bcbm)
gc()

# 3. Initialize CellChat Object from Normalized RNA Data
cat("[-] Initializing CellChat object with human ligand-receptor database...\n")
data_input <- GetAssayData(bcbm_sub, assay = "RNA", layer = "data")
meta_input <- bcbm_sub@meta.data %>% select(cell_type)

cellchat <- createCellChat(object = data_input, meta = meta_input, group.by = "cell_type")
CellChatDB <- CellChatDB.human
# Focus on Secreted Signaling & Cell-Cell Contact (most relevant to neuro-immune crosstalk)
CellChatDB.use <- subsetDB(CellChatDB, search = c("Secreted Signaling", "Cell-Cell Contact"))
cellchat@DB <- CellChatDB.use

rm(bcbm_sub, data_input)
gc()

# 4. Identify Overexpressed Ligands, Receptors & Compute Communication Probabilities
cat("[-] Identifying overexpressed ligand-receptor pairs...\n")
cellchat <- subsetData(cellchat)
cellchat <- identifyOverExpressedGenes(cellchat)
cellchat <- identifyOverExpressedInteractions(cellchat)

cat("[-] Computing cell-cell communication probabilities (takes ~2-4 mins on 8GB RAM)...\n")
cellchat <- computeCommunProb(cellchat, type = "triMean", population.size = TRUE)
cellchat <- filterCommunication(cellchat, min.cells = 10)
cellchat <- computeCommunProbPathway(cellchat)
cellchat <- aggregateNet(cellchat)

# 5. Export All Significant Interactions Sent by Chemoresistant Tumor Clones to Immune Cells
df_net <- subsetCommunication(cellchat)
write.csv(df_net, file.path(DIR_OUT_TAB, "05_cellchat_all_interactions.csv"), row.names = FALSE)

chemo_to_immune <- df_net %>%
  filter(
    source %in% c("Tumor_Chemoresistant_C11", "Tumor_Chemoresistant_C13"),
    target %in% c("Microglia_TAMs_C1QA", "Myeloid_AIF1", "T_Cells_CD3D", "NK_Cytotoxic_T_Cells", "BBB_Endothelial_1")
  ) %>%
  arrange(pval, desc(prob))

write.csv(chemo_to_immune, file.path(DIR_OUT_TAB, "05_chemoresistant_tumor_to_immune_signaling.csv"), row.names = FALSE)

# 6. Export Publication-Ready CellChat Figures (White Background)
cat("[-] Generating CellChat interaction network and bubble plots...\n")

# Figure A: Circle Plot of Interaction Weights from Chemoresistant Clones
png(file.path(DIR_OUT_FIG, "05_cellchat_chemoresistant_circle.png"), width = 900, height = 800, res = 120, bg = "white")
groupSize <- as.numeric(table(cellchat@idents))
sources_idx <- which(levels(cellchat@idents) %in% c("Tumor_Chemoresistant_C11", "Tumor_Chemoresistant_C13"))
netVisual_circle(
  cellchat@net$weight,
  vertex.weight = groupSize,
  weight.scale = TRUE,
  sources.use = sources_idx,
  label.edge = FALSE,
  title.name = "Signaling Sent by Chemoresistant Breast Tumor Clones (C11 & C13)"
)
dev.off()

# Figure B: Bubble Plot of Top Ligand-Receptor Pairs to Microglia & T-Cells
targets_idx <- which(levels(cellchat@idents) %in% c("Microglia_TAMs_C1QA", "Myeloid_AIF1", "T_Cells_CD3D", "BBB_Endothelial_1"))
p_bubble <- netVisual_bubble(
  cellchat,
  sources.use = sources_idx,
  targets.use = targets_idx,
  remove.isolate = TRUE
) +
  theme_bw(base_size = 10) +
  theme(
    axis.text.x = element_text(angle = 45, hjust = 1, color = "black"),
    plot.background = element_rect(fill = "white", color = NA)
  ) +
  ggtitle("Chemoresistant Tumor -> Brain Immune & BBB Ligand-Receptor Axes")

ggsave(
  file.path(DIR_OUT_FIG, "05_cellchat_ligand_receptor_bubble.png"),
  p_bubble,
  width = 11, height = 8, dpi = 150, bg = "white"
)

# Save Final CellChat Object
saveRDS(cellchat, file = file.path(DIR_DATA_SER, "05_cellchat_bcbm.rds"), compress = FALSE)

cat("\n[=== PHASE 5 COMPLETE: CELLCHAT NEURO-IMMUNE SIGNALING FINISHED ===]\n")
print(head(chemo_to_immune[, c("source", "target", "ligand", "receptor", "pathway_name", "prob", "pval")], 12))
