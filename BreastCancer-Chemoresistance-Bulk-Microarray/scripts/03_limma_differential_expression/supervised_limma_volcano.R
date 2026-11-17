# Milestone 3: Generating a Publication-Grade Volcano Plot
library(ggplot2)

# Prepare a plotting dataframe from the topTable results
volcano_df <- top_table %>%
    rownames_to_column(var = "Gene") %>%
    mutate(
        Significance = case_when(
            adj.P.Val < 0.05 & logFC > 0.5 ~ "Upregulated in Resistance",
            adj.P.Val < 0.05 & logFC < -0.5 ~ "Downregulated in Resistance",
            TRUE ~ "Not Significant"
        )
    )

# Build the Volcano Plot
volcano_plot <- ggplot(volcano_df, aes(x = logFC, y = -log10(P.Value), color = Significance)) +
    geom_point(alpha = 0.6, size = 1.8) +
    scale_color_manual(values = c(
        "Upregulated in Resistance" = "#E41A1C", # Red for active resistance drivers
        "Downregulated in Resistance" = "#377EB8", # Blue for suppressed targets
        "Not Significant" = "#999999" # Grey for genetic noise
    )) +
    # Add threshold lines
    geom_vline(xintercept = c(-0.5, 0.5), linetype = "dashed", color = "black", alpha = 0.5) +
    geom_hline(yintercept = -log10(0.005), linetype = "dashed", color = "black", alpha = 0.5) +
    labs(
        title = "Volcano Plot: Target Transcriptomic Shifts in Chemoresistance",
        subtitle = "Supervised limma Contrast (Resistant vs. Sensitive) | GSE25066",
        x = "Log2 Fold Change (Magnitude of Expression Shift)",
        y = "-Log10 Raw P-Value (Statistical Certainty)",
        caption = "Dashed lines represent exploratory biological thresholds."
    ) +
    theme_bw(base_size = 12) +
    theme(
        plot.title = element_text(face = "bold", size = 13, hjust = 0.5),
        plot.subtitle = element_text(size = 10, hjust = 0.5),
        legend.position = "bottom"
    )

ggsave("Biomarker_Volcano.png", plot = volcano_plot, width = 7, height = 6, dpi = 300)
cat("[SUCCESS] Volcano plot visual generated and saved as 'Biomarker_Volcano.png'.\n")
