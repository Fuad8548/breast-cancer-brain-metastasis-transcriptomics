# Part 1: Line-by-Line Explanation of the Codebase

## 1. Loading Settings, Libraries, and Data:r
```R
meta <- bcbm@meta.data
```
- `meta <- bcbm@meta.data` extracts the metadata table—a spreadsheet where every row is one individual cell, and the columns store each cell's ID, sample name (sample_id), tissue group (tissue_type), gene count (nFeature_RNA), UMI molecule count (nCount_RNA), and mitochondrial percentage (percent.mt).

## 2. Building Each Panel (p1, p2, p3) with ggplot2:r
```R
p1 <- ggplot(meta, aes(x = sample_id, y = nFeature_RNA, fill = tissue_type)) +
    geom_violin(trim = TRUE) +
    geom_boxplot(width = 0.1, fill = "white", outlier.size = 0.2) + ...
```
- aes(x = sample_id, y = ..., fill = tissue_type): Tells R to put the 9 patient samples on the horizontal X-axis, the QC metric on the vertical Y-axis, and color-code each sample by Coral/Red (Adjacent_Brain_Tissue) vs. Teal/Cyan (Brain_Metastasis_Tumor).
- geom_violin(trim = TRUE): Draws the colored "violin" silhouette showing the density (shape) of cells at every value—where the violin is wide, thousands of cells share that value; where it is skinny, very few cells have that value.
- geom_boxplot(width = 0.1, fill = "white", ...): Draws the tiny white box inside the center of each violin. The black horizontal line inside the white box is the median (the exact middle cell), and the top/bottom edges of the white box show where the middle 50% of cells sit.


# Part 2: How to Read the 3-Panel Plot 

The X-axis lists the 9 individual patient tissue samples from GSE268662, color-coded into two groups:
- Coral / Red (4 samples — Normal Adjacent Brain): GSM8296165_BN1, GSM8296166_BN2, GSM8296168_DN, GSM8296170_EN (N = Normal brain).
- Teal / Cyan (5 samples — Brain Metastasis Tumor): GSM8296167_BT, GSM8296169_DT, GSM8296171_ET, GSM8296172_FT, GSM8296173_GT (T = Tumor).

## Panel 1 (Left): Genes Detected per Cell (`nFeature_RNA`)

- Y-Axis (0 to ~7,000): How many different genes were turned ON (detected) inside each individual cell.
- Most cells in the brain (especially resting immune cells and neurons) express around 500 to 2,000 genes — that's why the violin bulges wide near the bottom and the white median box sits low.
- Rapidly dividing cancer cells express 3,000 to 6,500 genes, creating the tall, thin upward needle.
- Why it's flat at the very bottom (~200–300): Because Phase 1 deleted empty droplets (< MIN_FEATURES), so no cells exist below the cutoff.

## Panel 2 (Middle): UMI Counts per Cell (`nCount_RNA`)

- Y-Axis (0 to ~55,000): Total number of individual mRNA molecules (UMIs) captured inside each cell.
- Why the violins look flat at the bottom with tall thin spikes:
	- Most cells contain 1,000 to 5,000 total mRNA molecules (the wide base at the bottom), while a small subset of giant, hyper-active tumor cells contain 20,000 to 50,000 mRNA molecules (the thin vertical spikes).
	- Because all 9 samples sit at roughly the same baseline height, this proves sequencing depth is consistent across all 9 patient samples.

## Panel 3 (Right): Mitochondrial % (`percent.mt`)

- Y-Axis (0% to 15%): Percentage of a cell's RNA coming from mitochondrial genes (MT-...).
- Why the violins are chopped flat right at 15%:
	- When a cell ruptures or dies during surgery, its outer cytoplasmic RNA leaks out while mitochondrial RNA stays trapped inside — causing percent.mt to spike to 30%–80%.
	- In Phase 1, we set a strict filter (percent.mt < 15%), which sliced off every dying cell above 15%—leaving a clean flat ceiling at 15% across all 9 samples.
- Why Teal (Tumor) violins bulge higher (10%–14%) than Red (Normal Brain, 2%–6%):
	- Metastatic breast cancer cells are metabolically hyper-active and packed with mitochondria compared to normal brain tissue, which is a classic biological signature of tumor cells.


### Summary
- Genes & UMI Counts (Bottom-Heavy Bulges): In the first two plots, the bulges sit near the bottom (around 500–1,000 genes and 1,000–3,000 UMIs), and the white boxes/medians are compressed low down. Because scRNA-seq has low capture efficiency and many cells are naturally small or transcriptionally quiet, the vast majority of cells have low-to-moderate counts, while a small minority of large or highly active cells create the long, skinny upward tails (a right-skewed distribution).
- Mitochondrial % (Top-Heavy Bulges & 15% Cutoff): In the third plot, almost all violins have a sharp flat top at 15%, and several tumor samples (teal) bulge near the top (10%–15%). This happens because your dataset has already been filtered using a standard QC threshold (percent.mt < 15), cutting off dead or dying cells that leak cytoplasmic RNA while retaining mitochondrial RNA. Tumor samples often experience higher metabolic stress or cell death during tissue dissociation, pushing their median and bulk density closer to that 15% ceiling.