#!/usr/bin/env Rscript

# ==============================================================================
# Pipeline Phase 1: Automated GEO Data Ingestion, Per-Sample QC & Balancing
# Script: 01_data_ingestion.R
# Target Cohort: GSE268662 (Breast Cancer Brain Metastasis scRNA-seq - 8GB RAM)
# ==============================================================================

cat("[===] Initializing Phase 1: Ingesting GSE268662 single-cell data...\n")

if (!file.exists("config.R")) {
    stop("[!] Fatal Error: Root 'config.R' missing! Ensure working directory is project root.")
}
source("config.R")

suppressPackageStartupMessages({
    library(Seurat)
    library(Matrix)
    library(tidyverse)
})

# 1. Unpack .tar archive ONLY if .mtx.gz files are not already extracted
mtx_files <- list.files(DIR_DATA_RAW, pattern = "matrix\\.mtx(\\.gz)?$", full.names = TRUE, recursive = TRUE)

if (length(mtx_files) == 0) {
    tar_files <- list.files(DIR_DATA_RAW, pattern = "\\.tar$", full.names = TRUE)
    if (length(tar_files) == 0) {
        stop("[!] Fatal Error: Place GSE268662_RAW.tar inside 'data/raw/' first!")
    }
    for (tf in tar_files) {
        cat(sprintf("[-] Extracting archive: %s\n", basename(tf)))
        untar(tf, exdir = DIR_DATA_RAW)
    }
    mtx_files <- list.files(DIR_DATA_RAW, pattern = "matrix\\.mtx(\\.gz)?$", full.names = TRUE, recursive = TRUE)
}

if (length(mtx_files) == 0) {
    stop("[!] Fatal Error: No matrix.mtx.gz files found after extraction.")
}

# 2. Match Prefixed 10x Triplets via Exact File Paths (No Regex)
cat(sprintf("[+] Detected %d matrix file(s). Pairing with barcode and feature files...\n", length(mtx_files)))
seurat_list <- list()

for (i in seq_along(mtx_files)) {
    m_file <- mtx_files[i]
    m_dir <- dirname(m_file)

    # Extract exact prefix before 'matrix.mtx.gz' (e.g., "GSM8296168_DN_")
    prefix <- sub("matrix\\.mtx(\\.gz)?$", "", basename(m_file))
    sample_id <- ifelse(prefix == "", paste0("BCBM_Sample_", i), sub("_$", "", prefix))

    # Construct exact companion filenames
    b_file <- file.path(m_dir, paste0(prefix, "barcodes.tsv.gz"))
    f_file <- file.path(m_dir, paste0(prefix, "features.tsv.gz"))

    # Fallback if uncompressed or named 'genes.tsv.gz'
    if (!file.exists(b_file)) b_file <- file.path(m_dir, paste0(prefix, "barcodes.tsv"))
    if (!file.exists(f_file)) f_file <- file.path(m_dir, paste0(prefix, "genes.tsv.gz"))
    if (!file.exists(f_file)) f_file <- file.path(m_dir, paste0(prefix, "features.tsv"))

    if (!file.exists(b_file) || !file.exists(f_file)) {
        warning(sprintf("[!] Skipping %s: missing %s or %s", basename(m_file), basename(b_file), basename(f_file)))
        next
    }

    cat(sprintf("    -> Reading triplet for: %s\n", sample_id))
    counts <- ReadMtx(
        mtx = m_file,
        cells = b_file,
        features = f_file
    )

    obj <- CreateSeuratObject(counts = counts, project = sample_id, min.cells = MIN_CELLS, min.features = 100)
    obj$sample_id <- sample_id
    seurat_list[[sample_id]] <- obj
}

if (length(seurat_list) == 0) {
    cat("[!] Files inside data/raw/:\n")
    print(list.files(DIR_DATA_RAW))
    stop("[!] Fatal Error: No Seurat objects could be constructed.")
}

# 3. Per-Sample QC Filtering & Balanced Downsampling (8GB RAM Safe)
cat("[-] Filtering low-quality cells and balancing sample sizes for 8GB RAM...\n")
set.seed(42)

for (sample_id in names(seurat_list)) {
    obj <- seurat_list[[sample_id]]
    obj[["percent.mt"]] <- PercentageFeatureSet(obj, pattern = "^MT-")

    obj <- subset(
        obj,
        subset = nFeature_RNA > MIN_FEATURES &
            nFeature_RNA < MAX_FEATURES &
            percent.mt < MAX_MT_PCT
    )

    if (ncol(obj) > MAX_CELLS_PER_SAMPLE) {
        keep_cells <- sample(colnames(obj), size = MAX_CELLS_PER_SAMPLE)
        obj <- subset(obj, cells = keep_cells)
    }

    seurat_list[[sample_id]] <- obj
    cat(sprintf("    -> %s: %d viable cells ready\n", sample_id, ncol(obj)))
}
gc()

# 4. Memory-Safe Merge (Note double brackets [[1]] for Seurat object)
if (length(seurat_list) > 1) {
    cat(sprintf("[-] Merging %d balanced patient samples into unified BCBM cohort...\n", length(seurat_list)))
    bcbm_seurat <- merge(
        x = seurat_list[[1]],
        y = seurat_list[-1],
        add.cell.ids = names(seurat_list),
        project = PROJECT_NAME
    )
} else {
    bcbm_seurat <- seurat_list[[1]]
}

rm(seurat_list)
gc()

# 5. Compute Ribosomal QC & Save Uncompressed
cat("[-] Computing ribosomal (RPS/RPL) QC metrics and saving object...\n")
bcbm_seurat[["percent.ribo"]] <- PercentageFeatureSet(bcbm_seurat, pattern = "^RP[SL]")

output_rds_path <- file.path(DIR_DATA_SER, "01_filtered_merged_seurat.rds")
saveRDS(bcbm_seurat, file = output_rds_path, compress = FALSE)

cat("\n[=== PHASE 1 COMPLETE (8GB RAM MODE) ===]\n")
print(bcbm_seurat)
