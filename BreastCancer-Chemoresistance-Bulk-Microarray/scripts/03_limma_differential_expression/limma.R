# Milestone 3: Supervised Statistical Modeling with limma
library(limma)
library(tidyverse)

cat("--- Setting up Supervised Statistical Design Matrix ---\n")

# 1. Converting the clinical Response variable into a strict computational factor
# We set up the comparison so that it calculates: Resistant relative to Sensitive
response_factor <- factor(meta_aligned$Response, levels = c("Sensitive (pCR)", "Resistant (RD)"))

# 2. Build the design matrix
# This tells limma to model expression as a function of drug response status
design <- model.matrix(~response_factor)
colnames(design) <- c("Intercept", "Resistant_vs_Sensitive")

# 3. Fit the linear model for all 19,000+ genes simultaneously
fit <- lmFit(expression_aligned, design)

# 4. Apply Empirical Bayes moderation
# This borrows information across genes to make the statistical variance estimates highly robust
fit_ebayes <- eBayes(fit)

# 5. Extract the complete statistical summary table
top_table <- topTable(fit_ebayes, coef = "Resistant_vs_Sensitive", number = Inf, adjust.method = "BH")

# 6. Apply strict publication-grade statistical filters
# Threshold: Adjusted p-value (False Discovery Rate) < 0.05
significant_genes <- top_table %>%
    filter(adj.P.Val < 0.05) %>%
    arrange(adj.P.Val)

cat("\n=== Differential Expression Results ===\n")
cat("Total transcripts scanned across genome: ", nrow(top_table), "\n")
cat("Statistically significant resistance biomarkers found: ", nrow(significant_genes), "\n")

# 7. Save the candidate biomarkers to disk for our Python Machine Learning phase
write.csv(significant_genes, "Significant_Biomarkers.csv", row.names = TRUE)
cat("\n[SUCCESS] Targeted biomarker signature exported to 'Significant_Biomarkers.csv'!\n")
