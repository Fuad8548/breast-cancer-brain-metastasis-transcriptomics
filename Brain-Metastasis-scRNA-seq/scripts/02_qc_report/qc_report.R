source("config.R")
library(Seurat)
library(tidyverse)
library(patchwork)

bcbm <- readRDS(file.path(DIR_DATA_SER, "03_clustered_seurat.rds"))
meta <- bcbm@meta.data

p1 <- ggplot(meta, aes(x = sample_id, y = nFeature_RNA, fill = tissue_type)) +
    geom_violin(trim = TRUE) +
    geom_boxplot(width = 0.1, fill = "white", outlier.size = 0.2) +
    theme_minimal() +
    theme(axis.text.x = element_text(angle = 45, hjust = 1), legend.position = "none") +
    labs(title = "Genes Detected per Cell", x = "", y = "nFeature_RNA")

p2 <- ggplot(meta, aes(x = sample_id, y = nCount_RNA, fill = tissue_type)) +
    geom_violin(trim = TRUE) +
    geom_boxplot(width = 0.1, fill = "white", outlier.size = 0.2) +
    theme_minimal() +
    theme(axis.text.x = element_text(angle = 45, hjust = 1), legend.position = "none") +
    labs(title = "UMI Counts per Cell", x = "", y = "nCount_RNA")

p3 <- ggplot(meta, aes(x = sample_id, y = percent.mt, fill = tissue_type)) +
    geom_violin(trim = TRUE) +
    geom_boxplot(width = 0.1, fill = "white", outlier.size = 0.2) +
    theme_minimal() +
    theme(axis.text.x = element_text(angle = 45, hjust = 1)) +
    labs(title = "Mitochondrial %", x = "", y = "percent.mt")

ggsave(file.path(DIR_OUT_FIG, "02_qc_postfilter_violin.png"), p1 + p2 + p3, width = 14, height = 5, dpi = 150)
cat("[+] Fixed QC Violin Plot saved!\n")
