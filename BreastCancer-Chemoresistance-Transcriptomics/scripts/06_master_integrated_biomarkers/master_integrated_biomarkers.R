library(tidyverse)
library(patchwork)

# 1. Load your already-generated Milestone 3, 4, and 5 CSV files
limma_df <- read.csv("results/tables/03_Biomarkers.csv", row.names = 1) %>%
    rownames_to_column("Probe_ID")
ml_top10 <- read.csv("results/tables/05_GSE25066_Top10_ML_Biomarkers_Annotated.csv")
go_bridge <- read.csv("results/tables/Milestone4_ML_to_GO_Pathway_Bridge.csv")
ego_df <- as.data.frame(go_enrich) # From your Milestone 4 script

# 2. Merge Milestone 5 (ML Gini) + Milestone 3 (limma logFC & adj.P.Val)
master_df <- ml_top10 %>%
    left_join(limma_df %>% dplyr::select(Probe_ID, logFC, AveExpr, adj.P.Val), by = "Probe_ID") %>%
    mutate(
        Direction = ifelse(logFC < 0, "Resistant-High (-logFC)", "Sensitive-High (+logFC)"),
        Top_GO_Pathway = sapply(Gene_Symbol, function(g) {
            hits <- go_bridge$Description[grepl(g, go_bridge$Matched_ML_Genes)]
            if (length(hits) == 0) "Orthogonal Single-Gene Driver (Non-Cluster)" else hits[1]
        })
    )

write.csv(master_df, "GSE25066_Master_Integrated_Biomarkers.csv", row.names = FALSE)
print(master_df)

# 3. Build ONE Combined 3-Panel Figure Telling the Whole Story
# Panel A: Milestone 3 (limma Direction & Fold-Change of Top 10 ML Biomarkers)
pA <- ggplot(master_df, aes(x = reorder(Gene_Symbol, Gini_Importance), y = logFC, fill = Direction)) +
    geom_col(width = 0.68) +
    geom_hline(yintercept = 0, linetype = "dashed", color = "grey40") +
    coord_flip() +
    scale_fill_manual(values = c("Resistant-High (-logFC)" = "#B2182B", "Sensitive-High (+logFC)" = "#2166AC")) +
    theme_bw(base_size = 10) +
    labs(title = "A. Milestone 3: Expression Direction (limma)", x = NULL, y = "Log2 Fold Change (pCR vs RD)", fill = NULL) +
    theme(legend.position = "bottom")

# Panel B: Milestone 4 (Top 8 Enriched GO Biological Programs)
top_go <- ego_df %>%
    arrange(p.adjust) %>%
    head(8)
pB <- ggplot(top_go, aes(x = reorder(str_wrap(Description, 32), -log10(p.adjust)), y = -log10(p.adjust), size = Count)) +
    geom_point(color = "#6A3D9A") +
    coord_flip() +
    theme_bw(base_size = 10) +
    labs(title = "B. Milestone 4: Core GO Pathways", x = NULL, y = "-log10(Adjusted P-value)", size = "Gene Count") +
    theme(legend.position = "bottom")

# Panel C: Milestone 5 (Random Forest Gini Diagnostic Ranking)
pC <- ggplot(master_df, aes(x = reorder(Gene_Symbol, Gini_Importance), y = Gini_Importance, fill = Direction)) +
    geom_col(width = 0.68) +
    coord_flip() +
    scale_fill_manual(values = c("Resistant-High (-logFC)" = "#B2182B", "Sensitive-High (+logFC)" = "#2166AC")) +
    theme_bw(base_size = 10) +
    labs(title = "C. Milestone 5: ML Gini Diagnostic Rank", x = NULL, y = "Gini Importance Score", fill = NULL) +
    theme(legend.position = "bottom")

# Combine all 3 panels into one single figure
combined_story <- (pA | pB | pC) +
    plot_annotation(
        title = "Integrated Multi-Stage Biomarker Discovery in Chemoresistant Breast Cancer (GSE25066, ROC-AUC = 0.90)",
        subtitle = "Combining Differential Expression (limma), Functional Pathways (GO), and Multivariate Random Forest Classification",
        theme = theme(plot.title = element_text(face = "bold", size = 13))
    )

ggsave("Milestone_3_4_5_Integrated_Story.png", combined_story, width = 16, height = 6, dpi = 300, bg = "white")
cat("[+] Saved 'GSE25066_Master_Integrated_Biomarkers.csv' and 'Milestone_3_4_5_Integrated_Story.png'!\n")
