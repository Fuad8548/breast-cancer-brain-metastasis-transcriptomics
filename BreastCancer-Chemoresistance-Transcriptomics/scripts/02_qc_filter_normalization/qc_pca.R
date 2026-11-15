# Milestone 2: GSE25066 Cleaning & Resilient PCA Processing
library(GEOquery)
library(tidyverse)
library(ggplot2)

# 1. Ensure data is loaded
if (!exists("gse_data")) {
    gse <- getGEO("GSE25066", GSEMatrix = TRUE)
    gse_data <- gse[[1]]
}
metadata <- pData(gse_data)
expression_matrix <- exprs(gse_data)

cat("--- Processing Metadata and Cleaning Strings ---\n")

# 2. Extracting and formatting our target response column cleanly
meta_clean <- metadata %>%
    mutate(
        # Strip away the prefix text to isolate just "RD" or "pCR"
        Raw_Response = `characteristics_ch1.11`,
        Response = case_when(
            grepl("pCR", Raw_Response) ~ "Sensitive (pCR)",
            grepl("RD", Raw_Response) ~ "Resistant (RD)",
            TRUE ~ NA_character_
        ),
        # Clean the PAM50 molecular subtype labels for plotting
        Subtype = gsub("pam50_class: ", "", `characteristics_ch1.20`)
    ) %>%
    filter(!is.na(Response)) # Drop any missing samples

# 3. Synchronizing our expression matrix with the cleaned metadata
common_samples <- intersect(colnames(expression_matrix), rownames(meta_clean))
expression_aligned <- expression_matrix[, common_samples]
meta_aligned <- meta_clean[common_samples, ]

# 4. Enforcing Log2 normalization safety check
if (max(expression_aligned, na.rm = TRUE) > 50) {
    cat("Applying Log2 normalization transformation...\n")
    expression_aligned <- log2(expression_aligned + 1)
}

# 5. Transposing and clearing out absolute zero/constant variance columns
pca_input <- t(expression_aligned)
column_sds <- apply(pca_input, 2, sd)
pca_input_clean <- pca_input[, column_sds > 0 & !is.na(column_sds)]

cat("\n=== Cohort Distribution Metrics ===\n")
print(table(meta_aligned$Response))

# 6. Executing PCA safely
pca_result <- prcomp(pca_input_clean, scale. = TRUE, center = TRUE)
cat("\n[Success] PCA executed perfectly!\n")

# 7. Calculating Percent Variance Explained
var_explained <- (pca_result$sdev^2) / sum(pca_result$sdev^2) * 100

# 8. Creating Plotting Framework Dataframe
pca_df <- data.frame(
    PC1 = pca_result$x[, 1],
    PC2 = pca_result$x[, 2],
    Response = meta_aligned$Response,
    Subtype = meta_aligned$Subtype
)

# 9. Generating and Saving Publication Graphic
pca_plot <- ggplot(pca_df, aes(x = PC1, y = PC2, color = Response, shape = Subtype)) +
    geom_point(size = 3.5, alpha = 0.75) +
    scale_color_manual(values = c("Resistant (RD)" = "#D95F02", "Sensitive (pCR)" = "#1B9E77")) +
    labs(
        title = "PCA: Baseline Pre-treatment Transcriptomic Layout",
        subtitle = "NCBI GSE25066 Breast Cancer Validation Cohort",
        x = paste0("PC1 (", round(var_explained[1], 1), "% Variance)"),
        y = paste0("PC2 (", round(var_explained[2], 1), "% Variance)")
    ) +
    theme_bw() +
    theme(
        plot.title = element_text(face = "bold", size = 12, hjust = 0.5),
        plot.subtitle = element_text(size = 10, hjust = 0.5)
    )

ggsave("results/figures/Milestone2_GSE25066_PCA.png", plot = pca_plot, width = 7, height = 5, dpi = 300)
cat("QC Graph successfully generated and saved as 'results/figures/Milestone2_GSE25066_PCA.png'.\n")
