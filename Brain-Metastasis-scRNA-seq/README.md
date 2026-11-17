# Single-Cell Resolution of Primary Chemoresistance & Neuro-Immune Remodeling in Breast Cancer Brain Metastasis (`GSE268662`)

## 1. Executive Summary & Translational Bridge

Breast cancer brain metastasis (**BCBM**) represents a lethal clinical bottleneck driven by intrinsic therapy resistance and a specialized, immunosuppressive blood–brain barrier (**BBB**) microenvironment. This repository implements an end-to-end single-cell transcriptomic (**scRNA-seq**) and intercellular signaling pipeline analyzing **13,186 quality-controlled single cells** across **9 human clinical resections** (`5` Brain Metastasis Tumors vs. `4` Matched Adjacent Normal Brain tissues; **NCBI GEO: `GSE268662`**).

Directly bridging our upstream primary breast cancer machine-learning biomarker discovery (**`GSE25066`**, *n* = 508 neoadjuvant taxane–anthracycline patients; 2,938 significant `limma` DEGs & Random Forest Gini drivers), this study resolves:
1. **Intratumoral Phenotypic Plasticity:** How metastatic breast lesions split into differentiated luminal-like epithelial cores (`Malignant_Epithelial_1/2`, `EPCAM+ / KRT18+ / GATA3+ / TBC1D9+ / DSC2+`) and dedifferentiated, ER-stressed mesenchymal escape states (`Malignant_Plastic_EMT`, `TXNDC5+ / PGK1+`).
2. **Chemoresistance Program Persistence:** How the primary `GSE25066` chemoresistance module (`Resistance_Z`) and hypoxic glycolytic reprogramming (`PGK1`) concentrate within `Malignant_Epithelial_2` and the remodeling BBB vascular niche.
3. **Division of Labor in Neuro-Immune Reprogramming (`CellChat`):** How epithelial tumor cores suppress T/NK cells and microglia via **`APP–CD74`** and **`MIF–(CD74+CXCR4/CD44)`**, while `Malignant_Plastic_EMT` cells drive **`VEGFA–VEGFR1/2`** angiogenesis and recruit **`Microglia_TAMs`** (`CSF1–CSF1R`, `ANXA1–FPR1/2`).

---

## 2. Dataset Architecture & Cohort Design (`GSE268662`)

| Sample ID                                        | Tissue Compartment       | Clinical Condition                          | Post-QC Cells | Primary Biological Role                                                   |
| :----------------------------------------------- | :----------------------- | :------------------------------------------ | :-----------: | :------------------------------------------------------------------------ |
| **`BT`**, **`DT`**, **`ET`**, **`FT`**, **`GT`** | `Brain_Metastasis_Tumor` | Metastatic Breast Tumor Resection (*n* = 5) |     7,722     | Malignant epithelial subclones, EMT cells, infiltrating TAMs & T/NK cells |
| **`BN1`**, **`BN2`**, **`DN`**, **`EN`**         | `Adjacent_Brain_Tissue`  | Matched Adjacent Normal Brain (*n* = 4)     |     5,464     | Resting BBB endothelium, pericytes, oligodendrocytes & homeostatic glia   |
| **Total Cohort**                                 | **9 Patient Samples**    | **Paired Tumor vs. Normal Brain Atlas**     |  **13,186**   | **Harmony-integrated cross-patient single-cell ecosystem**                |

---

## 3. Five-Phase Reproducible Computational Pipeline

```text
BCBM-scRNAseq-Chemoresistance-Atlas/
├── config.R                                      # Centralized paths, QC thresholds & tissue color palettes
├── scripts/
│   ├── 01_data_ingestion_and_qc.R                # 10x H5 ingestion, MAD-based mitochondrial/feature QC & doublets
│   ├── 02_normalization_and_harmony.R            # LogNormalize, HVG selection, PCA & Harmony batch integration
│   ├── 03_clustering_and_marker_discovery.R      # Shared Nearest Neighbor (SNN) macro-clustering & UMAP embedding
│   ├── 04_chemoresistance_projection.R           # Affymetrix U133A -> HUGO symbol bridge, annotation & Z-scoring
│   └── 05_cellchat_neuro_immune.R                # Ligand-receptor probability inference & signaling network export
├── data/signatures/
│   ├── 03_Biomarkers.csv                         # Upstream GSE25066 limma differential expression table (2,938 probes)
│   └── 04_ML_Top_Biomarkers.csv                  # Top 10 Random Forest Gini importance drivers (TMSB15A -> PGK1)
└── results/
    ├── figures/                                  # High-resolution (300 DPI) publication figures
    └── tables/                                   # Annotated cell-type scores & ranked CellChat interactomes
```

### Pipeline Workflow Summary

* **Phase 01 — Multi-Sample Ingestion & Quality Control (`01_data_ingestion_and_qc.R`):** Ingests raw 10x Genomics HDF5 feature-barcode matrices across all 9 clinical specimens, filters low-quality droplets and high-mitochondrial stress artifacts, and retains **13,186 high-confidence cells**.
* **Phase 02 — Normalization & Harmony Integration (`02_normalization_and_harmony.R`):** Applies `LogNormalize`, selects the top 2,000 highly variable features, and runs **Harmony** across patient sample IDs to eliminate inter-patient batch effects while preserving biological separation between tumor and normal brain lineages.
* **Phase 03 — Macro-Clustering & Spatial Niche Structure (`03_clustering_and_marker_discovery.R`):** Resolves **11 macro-clusters (`0–10`)** at resolution `0.2`, partitioning the UMAP into three spatial compartments: **Tumor-Dominated Islands**, **Shared Tumor–Host Microenvironment**, and **Normal-Brain-Dominated Islands**.
* **Phase 04 — Cross-Platform Probe Bridge, Lineage Annotation & Chemoresistance Scoring (`04_chemoresistance_projection.R`):**
  * Translates **2,938 Affymetrix U133A (`GPL96`) probes** from `GSE25066` into official HUGO Gene Symbols via `hgu133a.db`.
  * Annotates all 13,186 cells into **10 canonical biological lineages** (`Malignant_Epithelial_1`, `Malignant_Epithelial_2`, `Malignant_Plastic_EMT`, `Microglia_TAMs`, `Resident_Myeloid_Glia`, `T_NK_Lymphocytes`, `BBB_Endothelial_1`, `BBB_Endothelial_2`, `Pericytes_Vascular`, `Oligodendrocytes`).
  * Computes standardized single-cell module scores (`Resistance_Z` & `Sensitivity_Z` via `AddModuleScore`) and maps the Top 10 Random Forest biomarkers across spatial compartments.
* **Phase 05 — Tumor-to-Brain Niche Interactome Modeling (`05_cellchat_neuro_immune.R`):** Uses **CellChat** (`triMean`, `population.size = TRUE`, `min.cells = 10`) across `Secreted Signaling` and `Cell-Cell Contact` databases to quantify statistically significant (`p < 0.01`) ligand–receptor axes sent from metastatic tumor clones to surrounding immune and BBB cells.

---

## 4. Key Biological & Quantitative Discoveries

### A. Intratumoral Heterogeneity & Primary ML Biomarker Mapping (Phase 04)
Projecting our `GSE25066` Random Forest biomarkers onto the annotated brain metastasis atlas (`04_complete_master_bridge.png` & `04_chemoresistance_by_celltype.png`) uncovers three distinct tumor states:

1. **Differentiated Luminal-Like Tumor Core (`Malignant_Epithelial_1`):** Expresses canonical epithelial markers (`EPCAM+`, `KRT18+`) alongside primary breast luminal and adhesion sentinels (**`GATA3`**, **`TBC1D9`**, **`DSC2`**).
2. **Hyper-Metabolic Chemoresistant Subclone (`Malignant_Epithelial_2`):** Exhibits the **highest standardized chemoresistance enrichment (`Resistance_Z` median `+0.8`, upper tail `+5.1`)** among malignant populations, pairing luminal retention (`KRT18+`, `GATA3+`) with intense hypoxic glycolysis (**`PGK1`**) and focal **`TMSB15A`** expression.
3. **Dedifferentiated EMT Escape State (`Malignant_Plastic_EMT`):** Silences epithelial junctions and chemosensitive luminal markers (`EPCAM-`, `DSC2-`, `GATA3-`, `TBC1D9-`) while upregulating the endoplasmic reticulum stress chaperone **`TXNDC5`** and glycolytic **`PGK1`** to survive therapeutic and microenvironmental stress.

### B. Functional Division of Labor in Neuro-Immune Crosstalk (Phase 05 — `CellChat`)
Filtering all significant (`pval = 0`, `p < 0.01`) tumor-to-niche interactions (`05_tumor_to_neuro_immune_signaling.csv`) reveals that distinct tumor subpopulations execute complementary niche-remodeling roles:

| Sender Tumor Clone           | Receiver Niche Population                  | Top Ligand → Receptor Axis                                       | Pathway / Class     | Comm. Prob (`p < 0.01`)  | Functional Impact on Brain Microenvironment                                               |
| :--------------------------- | :----------------------------------------- | :--------------------------------------------------------------- | :------------------ | :----------------------: | :---------------------------------------------------------------------------------------- |
| **`Malignant_Epithelial_1`** | `BBB_Endothelial_1`                        | **`APP` → `CD74`**                                               | `APP` (Contact)     |        `0.002615`        | Saturates `CD74` at the vascular interface, modulating endothelial antigen presentation   |
| **`Malignant_Epithelial_1`** | `T_NK_Lymphocytes`                         | **`APP` → `CD74`**                                               | `APP` (Contact)     |        `0.002278`        | Impairs MHC-II chaperone signaling and dampens local T/NK cytotoxic activation            |
| **`Malignant_Epithelial_1`** | `Microglia_TAMs`                           | **`APP` → `CD74`**                                               | `APP` (Contact)     |        `0.001535`        | Reprograms CNS microglia toward an antigen-presentation-defective, tumor-supportive state |
| **`Malignant_Epithelial_1`** | `T_NK_Lymphocytes`                         | **`MIF` → `CD74 + CXCR4`**<br>**`MIF` → `CD74 + CD44`**          | `MIF` (Secreted)    | `0.001297`<br>`0.000956` | Drives T-cell dysfunction and chemokine-axis arrest across infiltrating lymphocytes       |
| **`Malignant_Epithelial_1`** | `T_NK_Lymphocytes` & `BBB_Endothelial_1`   | **`MDK` → `NCL`**                                                | `MK` (Secreted)     | `0.000508`<br>`0.000489` | Promotes tumor survival, immune evasion, and endothelial activation via Midkine–Nucleolin |
| **`Malignant_Plastic_EMT`**  | `BBB_Endothelial_1`                        | **`VEGFA` → `FLT1` (`VEGFR1`)**<br>**`NAMPT` → `ITGA5 + ITGB1`** | `VEGF` / `VISFATIN` | `0.000214`<br>`0.000219` | Stimulates neo-angiogenesis and blood–tumor barrier remodeling                            |
| **`Malignant_Plastic_EMT`**  | `Resident_Myeloid_Glia` & `Microglia_TAMs` | **`ANXA1` → `FPR1`**<br>**`CSF1` → `CSF1R`**                     | `ANNEXIN` / `CSF`   | `0.000220`<br>`0.000082` | Recruits brain myeloid glia and sustains immunosuppressive `Microglia_TAMs` survival      |

