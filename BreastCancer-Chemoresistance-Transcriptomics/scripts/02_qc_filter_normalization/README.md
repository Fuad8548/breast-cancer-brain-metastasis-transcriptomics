# The Clinical Theory: Why pCR vs. RD Dictates Response
- **Pathologic Complete Response (pCR)**: This is the gold standard endpoint in neoadjuvant (pre-surgery) chemotherapy clinical trials. It means that after the patient completes their course of chemotherapy, the surgeon removes the breast tissue and lymph nodes, and the pathologist finds zero invasive cancer cells remaining under the microscope. Biologically, the patient's tumor cells were completely sensitive to the drug.
- **Residual Disease (RD)**: This means that despite completing the full course of chemotherapy, a viable tumor mass still remains in the breast or lymph nodes at the time of surgery. Biologically, sub-populations of cancer cells used survival mechanisms (such as drug efflux pumps, DNA repair activation, or mutation changes) to withstand the drug treatment. These tumors are classified as clinically resistant.


# Understanding the PCA plot ([pca plot](results/figures/GSE25066_Faceted_PCA.png))

## Part 1: The Biology of the 5 Subtypes
Breast cancer is not a single disease. In 2000, scientists discovered that breast tumors can be classified into 5 distinct molecular subtypes (called the Perou-Sorlie Molecular Subtypes or **PAM50** classification) based entirely on which master-regulator genes they turn on.
- **Basal-like** (Top Left Panel): These tumors mimic the basal layer of the breast ducts. They are heavily associated with **Triple-Negative Breast Cancer (TNBC)** because they lack Estrogen Receptors (ER), Progesterone Receptors (PR), and HER2 amplification. They are highly aggressive, mutate rapidly, and have high proliferation rates.
- **Her2-enriched** (Top Middle Panel): These tumors feature high amplification and overexpression of the `ERBB2` (`HER2`) oncogene. Historically highly aggressive, they are now targeted directly with anti-HER2 therapies (like **Trastuzumab**).
- **Luminal A** (`LumA` - Top Right Panel): These tumors mimic the inner luminal lining of the breast ducts. They are ER-positive and/or PR-positive, grow slowly, and generally have the best clinical prognosis.
- **Luminal B** (`LumB` - Bottom Left Panel): Also ER-positive, but unlike Luminal A, they express high levels of cell-proliferation markers (like MKI67). They grow faster and are more aggressive than Luminal A.
- **Normal-like** (Bottom Middle Panel): These tumors resemble normal, healthy breast tissue profiles. They are often smaller tumors with favorable outcomes but are less common in clinical datasets.


## Part 2: Understanding the Math (PC1, PC2, and the 10.2% Variance)
Our plot shows that PC1 explains `10.2%` of the variance, and PC2 also explains `10.2%` of the variance.

**What is a Principal Component (PC)?**
A single patient has `19,502` different gene expression numbers in this dataset. We cannot look at a 19,502-dimensional graph.
- PC1 is a mathematical formula that combines all 19,502 genes into a single score. It finds the single biggest pattern or "trend" that makes these `306` patients look different from each other.
- PC2 is the second biggest pattern, completely independent of PC1.

### Why is the variance only 10.2%?
In simple datasets (like comparing a heart cell to a brain cell), PC1 might explain 70% of the variance because the groups are drastically different.
In our data, a value of `10.2%` means the data is **incredibly noisy and complex**. There is no single "super-gene" or dominant group of genes dictating the layout. Instead, thousands of small, independent biological pathways are fluctuating simultaneously across different patients.


## Part 3: How to Interpret the Interaction of the Boxes, Axes, and Colors

1. Subtypes are Driven by PC1 (The Horizontal Axis)
- Luminal A is shifted entirely to the right side of the graph (centered around +25 to +50 on the X-axis).
- Basal is shifted entirely to the left side of the graph (centered around -25 to -50 on the X-axis).
- Her2, LumB, and Normal sit right in the middle (around 0).

The Scientific Interpretation:
This tells us that PC1 captures the biological identity of the tumor. The single largest variance trend in our data (10.2%) is the difference between Hormone-Receptor Positive tumors (Luminal A) on the right and Hormone-Receptor Negative tumors (Basal) on the left.

2. The Relationship Between Red and Blue Dots
- Red dots (Resistant - Residual Disease) and Blue dots (Sensitive - pCR) are completely swirled together within every single panel.
- If drug resistance were a major driver of the overall transcriptome, all the red dots would gather at the top of the box and all the blue dots would gather at the bottom. But they do not. They sit right next to each other.

The Scientific Interpretation:
This proves that at baseline (before treatment), a resistant tumor looks almost exactly like a sensitive tumor when viewed globally. The genetic signals that allow a tumor to survive chemotherapy are hidden beneath the massive noise of the tumor's subtype identity.

## Core narrative
"Preliminary Exploratory Data Analysis (EDA) via Principal Component Analysis (PCA) on the baseline GSE25066 cohort (n = 306) demonstrated that the primary axes of transcriptomic variance (PC1 = 10.2%, PC2 = 10.2%) are driven by intrinsic molecular subtype architectures rather than therapeutic responses. As visualized in the faceted projection, Luminal A and Basal subtypes are strongly segregated along the PC1 horizontal axis, reflecting baseline hormonal receptor status. Concurrently, chemosensitive (pCR) and chemoresistant (RD) clinical phenotypes exhibit profound spatial overlap across all panels. This architectural confounding underscores a critical limitation of unsupervised analysis and directly justifies the supervised pipeline proposed herein: utilizing a localized limma linear model to regress out global subtype noise, followed by Random Forest feature selection to isolate the discrete, low-variance transcriptomic biomarkers driving clinical drug resistance."


