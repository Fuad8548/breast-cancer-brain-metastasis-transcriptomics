library(clusterProfiler)
library(org.Hs.eg.db)
library(hgu133a.db) # Exact GPL96 chip database for GSE25066!
library(enrichplot)
library(ggplot2)
library(dplyr)

cat("--- Running Split Functional Enrichment (Resistant vs. Sensitive) ---\n")

# 1. Read Milestone 3 significant DEGs
deg_table <- read.csv("results/tables/03_Biomarkers.csv", row.names = 1)
deg_table$Probe_ID <- rownames(deg_table)

# 2. Map Affymetrix HG-U133A Probe IDs to Entrez IDs
gene_map <- bitr(deg_table$Probe_ID,
  fromType = "PROBEID",
  toType   = "ENTREZID",
  OrgDb    = hgu133a.db
)

deg_mapped <- deg_table %>%
  inner_join(gene_map, by = c("Probe_ID" = "PROBEID")) %>%
  mutate(Phenotype = ifelse(logFC < 0, "Chemoresistant (RD-High)", "Chemosensitive (pCR-High)"))

# 3. Split genes into two lists FIRST, then take unique Entrez IDs per group
gene_clusters <- lapply(
  split(as.character(deg_mapped$ENTREZID), as.character(deg_mapped$Phenotype)),
  unique
)

# Verify both groups look clean
str(gene_clusters)

# 4. Run Side-by-Side GO Enrichment (compareCluster)
comp_go <- compareCluster(
  geneCluster   = gene_clusters,
  fun           = "enrichGO",
  OrgDb         = org.Hs.eg.db,
  ont           = "BP",
  pAdjustMethod = "BH",
  pvalueCutoff  = 0.01,
  qvalueCutoff  = 0.05,
  readable      = TRUE
)

# 5. Save Publication-Grade Split Dotplot
p_split_go <- dotplot(comp_go, showCategory = 10) +
  ggtitle("Functional Pathways: Chemoresistant (RD) vs. Chemosensitive (pCR) Tumors") +
  theme(axis.text.x = element_text(angle = 15, hjust = 1, face = "bold"))

ggsave("results/figures/Milestone4_GO_Split_Pathways.png", plot = p_split_go, width = 10, height = 7.5, dpi = 300, bg = "white")
write.csv(as.data.frame(comp_go), "results/tables/Milestone4_GO_Split_Results.csv", row.names = FALSE)
cat("[SUCCESS] Saved split GO plot and table!\n")
