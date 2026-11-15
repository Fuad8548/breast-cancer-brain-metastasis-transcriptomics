# ==============================================================================
# Global Configuration & Environment Anchors
# Project: Single-Cell Mapping of Chemoresistance in Breast Cancer Brain Metastasis
# ==============================================================================

# Directory Structure Anchors
DIR_DATA_RAW <- file.path("data", "raw")
DIR_DATA_SER <- file.path("data", "serialized")
DIR_OUT_FIG <- file.path("outputs", "figures")
DIR_OUT_TAB <- file.path("outputs", "tables")

# Ensure all required directories exist
dirs_to_create <- c(DIR_DATA_RAW, DIR_DATA_SER, DIR_OUT_FIG, DIR_OUT_TAB)
for (d in dirs_to_create) {
    if (!dir.exists(d)) dir.create(d, recursive = TRUE)
}

# Dataset Metadata
GEO_ACCESSION <- "GSE268662"
PROJECT_NAME <- "BCBM_Chemoresistance_Atlas"

# Quality Control Thresholds (Calibrated for 10x Genomics Human Brain Metastasis)
MIN_CELLS <- 3 # Keep genes detected in at least 3 cells
MIN_FEATURES <- 300 # Remove empty droplets / low-complexity cellular debris
MAX_FEATURES <- 7500 # Remove potential doublets / multiplets
MAX_MT_PCT <- 15.0 # Remove lysed/dying cells (>15% mitochondrial reads)

# Dimensionality & Clustering Parameters
TARGET_PCS <- 30 # Number of Principal Components to compute


# Add this line to config.R
MAX_CELLS_PER_SAMPLE <- 1500 # Caps each patient sample at 1,500 viable cells for 8GB RAM safety
