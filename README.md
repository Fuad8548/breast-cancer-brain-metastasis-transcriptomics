# Multi-Omics & Machine Learning Dissection of Breast Cancer Progression: From Primary Chemoresistance to Brain Metastasis

## 📌 Executive Summary
Treatment failure in locally advanced and metastatic breast cancer follows a two-step clinical trajectory: **intrinsic/acquired chemoresistance** in the primary tumor followed by **organotropic dissemination to the brain** (Breast Cancer Brain Metastasis, BCBM), where tumor cells evade neuro-immune surveillance.

This computational oncology portfolio integrates **bulk clinical transcriptomics (`GSE25066`, *n* = 306 patients)** with **single-cell RNA sequencing (`GSE268662`, human BCBM resections)** across a unified **R + Python** bioinformatics stack to map the molecular continuum from primary taxane-anthracycline resistance to brain microenvironmental immune suppression.

---

## 🧭 Unified Translational Research Continuum

```text
STAGE 1: PRIMARY BREAST TUMOR CHEMORESISTANCE (GSE25066 | Bulk Microarray, n = 306)
├── Cohort: 249 Chemoresistant (Residual Disease, RD) vs. 57 Chemosensitive (pCR)
├── Differential Expression (limma): 2,938 significant transcripts (adj.P.Val < 0.05)
├── Pathway Dichotomy (compareCluster GO-BP):
│    ├── RD-High: Chromosome Segregation (p.adj = 3.98e-22) & Mitotic Transition
│    └── pCR-High: Fatty Acid Metabolism & Microtubule Transport (p.adj = 7.61e-05)
└── Diagnostic ML Panel (Random Forest, ROC-AUC = 0.9033 | RD Precision = 0.94):
     └── Top Drivers: TMSB15A (#1), TTLL4 (#2), EPHB3 (#3) + Pathway Sentinels (SKP1, GATA3, PGK1)
                                          │
                                          ▼  [Dissemination & Blood-Brain Barrier Colonization]
STAGE 2: BREAST CANCER BRAIN METASTASIS NICHE (GSE268662 | Single-Cell RNA-seq)
├── Ecosystem Resolution (Seurat): Malignant Epithelial, Microglia/TAMs, T Cells, Astrocytes & Endothelium
└── Intercellular Interactome (CellChat):
     ├── MIF ──► CD74/(CXCR4/CD44): Polarizes brain microglia toward pro-tumorigenic M2 phenotype
     └── APP ──► CD74: Jams MHC Class II antigen presentation to evade cytotoxic T-cell clearance
```

---

## 📊 Comparative Overview of Both Studies

| Dimension                       | Study 1: Primary Chemoresistance (`GSE25066`)                                                                                                     | Study 2: Brain Metastasis Niche (`GSE268662`)                                                                        |
| :------------------------------ | :------------------------------------------------------------------------------------------------------------------------------------------------ | :------------------------------------------------------------------------------------------------------------------- |
| **Clinical Question**           | Which genes and pathways allow primary breast tumors to survive neoadjuvant taxane-anthracycline chemotherapy?                                    | How do disseminated breast cancer cells reprogram resident brain immune cells to colonize the brain parenchyma?      |
| **Data Modality & Scale**       | Clinical Bulk Microarray (`GPL96` / HG-U133A)<br>**306 patient tumors** (`249 RD` vs. `57 pCR`)                                                   | Single-Cell RNA-seq (`10x Genomics`)<br>**Human BCBM surgical resections** (multi-lineage)                           |
| **Core Tech Stack**             | **R:** `GEOquery`, `limma`, `clusterProfiler`, `hgu133a.db`<br>**Python:** `scikit-learn` (Random Forest), `pandas`, `seaborn`                    | **R:** `Seurat` (v5), `CellChat`, `SingleR`, `clusterProfiler`, `org.Hs.eg.db`                                       |
| **Primary Quantitative Output** | • **2,938 DEGs** (`adj.P.Val < 0.05`)<br>• **Holdout ROC-AUC = 0.9033** (`0.94` Resistant precision)                                              | • High-resolution UMAP of tumor, glial & immune states<br>• Ranked intercellular ligand-receptor probability network |
| **Key Biomarkers & Targets**    | **Resistance Drivers (`-logFC`):** `TMSB15A`, `TTLL4`, `EPHB3`, `DSC2`, `PGK1`<br>**Sensitivity Sentinels (`+logFC`):** `TBC1D9`, `SKP1`, `GATA3` | **Master Immune Checkpoints:** `MIF–CD74` and `APP–CD74` signaling axes converging on **Microglia / TAMs**           |
| **Therapeutic Strategy**        | Combine taxanes with **spindle/cell-cycle inhibitors (`AURKA/B`, `PLK1`, `CDK4/6`)** or `TTLL4`/`TMSB15A` targeting                               | **`CD74` receptor blockade** (e.g., milatuzumab) or **`MIF` inhibition** to restore microglial anti-tumor immunity   |

---

## 🔬 Integrated Biological Synthesis

1. **Cytoskeletal & Mitotic Escape Drives Primary Chemoresistance (`GSE25066`):**
   * Chemosensitive (`pCR-High`) tumors rely on intact **microtubule-based transport (`GO:0099111`, `MAPT`)** and differentiated **fatty-acid metabolism (`GO:0006631`, `GATA3`)**, making them directly susceptible to taxane-induced microtubule arrest.
   * Conversely, chemoresistant (`RD-High`) tumors hyper-activate **chromosome segregation (`GO:0007059`, `p.adjust = 3.98e-22`)**, **hypoxic glycolysis (`PGK1`)**, and orthogonal cytoskeletal/tubulin modifiers (**`TMSB15A`**, **`TTLL4`**) that remodel microtubules to evade taxane binding and promote invasive motility.
2. **Surviving Resistant Clones Hijack Brain Microglia via `CD74` (`GSE268662`):**
   * Once invasive, chemoresistant cells cross the blood-brain barrier, single-cell `CellChat` modeling reveals that malignant epithelial cells and reactive astrocytes establish an immunosuppressive sanctuary by targeting the **`CD74` receptor** on brain-resident microglia and tumor-associated macrophages (TAMs) via **`MIF–CD74`** (M2 polarization) and **`APP–CD74`** (MHC-II antigen-presentation blockade).

---

## 📂 Portfolio Repositories & Artifacts

* **Repository 1 — Breast Cancer Chemoresistance ML (`GSE25066`):**
  * Full R statistical/pathway pipeline (`limma` + `compareCluster`) and Python `RandomForestClassifier` (`ROC-AUC = 0.9033`).
  * Key artifacts: `GSE25066_Master_Integrated_Biomarkers.csv`, `Milestone4_GO_Split_Pathways.png`, `Milestone_3_4_5_Integrated_Story.png`.
* **Repository 2 — Breast Cancer Brain Metastasis scRNA-seq (`GSE268662`):**
  * Full R `Seurat` clustering, lineage annotation, and `CellChat` ligand-receptor interaction pipeline.
  * Key artifacts: `BCBM_CellChat_Significant_Interactions.csv`, `Milestone2_UMAP_CellTypes.png`, `Milestone4_CellChat_MIF_APP_CD74_Network.png`.

---

## 🛠️ Core Technical Competencies Demonstrated
* **Statistical Genomics & Differential Expression:** Empirical Bayes linear modeling (`limma`), multiple-testing FDR correction (`Benjamini-Hochberg`), platform-specific probe annotation (`hgu133a.db`).
* **Functional Pathway & Single-Cell Analysis:** Over-representation & split-cluster enrichment (`clusterProfiler`), single-cell QC/normalization/clustering (`Seurat`), and intercellular ligand-receptor inference (`CellChat`).
* **Machine Learning & Reproducible Engineering:** Stratified cross-validation, class-imbalance handling (`class_weight="balanced"`), ROC-AUC benchmarking, Gini feature selection, and modular R/Python pipelines.
