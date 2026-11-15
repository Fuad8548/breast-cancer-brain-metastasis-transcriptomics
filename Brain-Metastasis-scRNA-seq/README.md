# Single-Cell Transcriptomic & Cell-Cell Communication Profiling of the Breast Cancer Brain Metastasis Microenvironment (GSE268662)

## 📌 Project Overview
Breast Cancer Brain Metastasis (BCBM) carries a dismal prognosis because disseminated tumor cells cross the blood-brain barrier and reprogram resident neural and immune cells into an immunosuppressive niche. 

This end-to-end **R single-cell computational genomics pipeline** analyzes human BCBM tumor resections from NCBI GEO (**GSE268662**) using **Seurat** and **CellChat** to:
1. Perform single-cell quality control, normalization (`SCTransform` / LogNormalize), dimensionality reduction (PCA/UMAP), and unsupervised graph-based clustering to resolve intratumoral cellular heterogeneity.
2. Annotate malignant epithelial subclusters alongside brain-resident and infiltrating stromal/immune compartments (**Microglia/Macrophages, T Cells, Astrocytes, Oligodendrocytes, and Endothelial Cells**).
3. Reconstruct the intercellular signaling network (`CellChat`) to pinpoint dominant immunosuppressive ligand-receptor axes (**`MIF–CD74`** and **`APP–CD74`**) mediating tumor–microglia crosstalk and actionable immunotherapy targets.

---

## 🧬 End-to-End Single-Cell Architecture

```text
NCBI GEO (GSE268662: Human Breast Cancer Brain Metastasis scRNA-seq)
       │
       ▼
[Milestone 1: R / Seurat QC & Normalization] ──► Mitochondrial/Feature Filtering, Doublet Removal & Scaling
       │
       ▼
[Milestone 2: Clustering & Cell Annotation] ──► PCA, UMAP Embedding & Canonical Marker Identification
       │                                         ├─ Malignant Epithelial Tumor Subclusters
       │                                         ├─ Brain-Resident Microglia & Tumor-Associated Macrophages (TAMs)
       │                                         └─ Lymphoid (T/NK), Reactive Astrocytes & Vascular Endothelium
       │
       ├──► [Milestone 3: Subcluster Functional Profiling]
       │      Differential Marker Genes (FindAllMarkers) & GO/KEGG Pathway Enrichment
       │
       └──► [Milestone 4: R / CellChat Interactome Modeling]
              Quantitative Ligand-Receptor Inference & Signaling Network Topology:
              ├─ MIF ──► CD74 / (CD74+CXCR4 / CD74+CD44) Axis (Microglial M2 Reprogramming)
              └─ APP ──► CD74 Axis (MHC-II Antigen-Presentation Blockade)
```

---

## 🔬 Key Cellular & Microenvironmental Discoveries

### 1. Cellular Ecosystem Resolution (`Seurat` Clustering & Canonical Markers)
Unsupervised Leiden/Louvain clustering and UMAP manifold projection partitioned the BCBM tumor microenvironment into distinct functional compartments:
- **Malignant Epithelial Compartment:** Characterized by high expression of epithelial/tumor markers (`EPCAM`, `KRT8`, `KRT18`, `KRT19`, `MKI67`), exhibiting pronounced transcriptional plasticity and metabolic adaptation to the brain parenchyma.
- **Myeloid / Microglia-Macrophage Compartment:** Dominated by brain-resident microglia (`TMEM119`, `P2RY12`, `CX3CR1`) and infiltrating monocyte-derived macrophages (`CD68`, `CD163`, `APOE`, `C1QA`), forming the largest non-malignant immune population in the brain lesion.
- **Adaptive Immune & Neuro-Vascular Niche:** Comprises exhausted/effector T-cell subsets (`CD3D`, `CD8A`, `IL7R`), reactive astrocytes (`GFAP`, `AQP4`), and blood-tumor barrier endothelial cells (`PECAM1`, `CLDN5`).

### 2. Dominant Immunosuppressive Ligand-Receptor Axes (`CellChat` Interactome)
Systematic ligand-receptor probability modeling (`CellChat`) identified **Malignant Tumor Cells** and **Reactive Astrocytes** as primary signaling senders targeting **Microglia/Macrophages** via two master immunosuppressive checkpoints:

| Signaling Axis   | Ligand (Sender)                                      | Receptor Complex (Receiver)                | Primary Cell-Cell Direction                              | Biological Mechanism & Translational Implication                                                                                                                                                                 |
| :--------------- | :--------------------------------------------------- | :----------------------------------------- | :------------------------------------------------------- | :--------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------- |
| **`MIF – CD74`** | **`MIF`** *(Macrophage Migration Inhibitory Factor)* | **`CD74` + `CXCR4`** / **`CD74` + `CD44`** | Malignant Tumor Cells → Microglia & TAMs                 | **Microglial Immunosuppressive Polarization:** Tumor-secreted `MIF` binds `CD74` on brain microglia/macrophages, suppressing anti-tumor phagocytosis and driving pro-tumorigenic (`M2`-like) cytokine secretion. |
| **`APP – CD74`** | **`APP`** *(Amyloid Beta Precursor Protein)*         | **`CD74`**                                 | Malignant Cells & Reactive Astrocytes → Microglia & TAMs | **Antigen-Presentation Jamming:** Binding of `APP` to `CD74` blocks `CD74`-mediated `MHC Class II` trafficking to the cell surface, preventing microglia from presenting tumor antigens to cytotoxic T cells.    |

> **Therapeutic Takeaway:** Both dominant tumor-to-microglia communication highways converge on the **`CD74` receptor**, highlighting **`CD74`-targeted monoclonal antibodies (e.g., milatuzumab)** or **`MIF` inhibitors** as promising microenvironmental therapies to re-sensitize BCBM lesions to immune checkpoint blockade.

---

## 🚀 Translational Next Steps
* **Cross-Dataset Integration (`GSE25066` $\leftrightarrow$ `GSE268662`):** Map the bulk taxane-chemoresistance signature (`TMSB15A`, `TTLL4`, `EPHB3`, `PGK1`) onto malignant epithelial subclusters in `GSE268662` to identify pre-existing chemoresistant clones within brain metastases.
* **Spatial Transcriptomics Validation:** Project `MIF–CD74` and `APP–CD74` signaling niches onto spatial transcriptomic brain metastasis sections to map tumor-microglia border zones.


# Repository structure

```text
singlecell-brain-metastasis/
├── README.md                                      # Project documentation & key findings
├── config.R                                       # Global directory & 8GB RAM parameters
├── scripts/
│   ├── 01_data_ingestion.R                        # Per-sample QC & 8GB RAM balancing
│   ├── 02_qc_report.R                             # Violin plots & per-sample audit table
│   ├── 03_normalization_clustering.R              # LogNormalize, Harmony & 17-cluster UMAP
│   ├── 04_chemoresistance_projection.R            # hgu133a.db probe bridge & AddModuleScore
│   ├── 04b_annotate_and_fix_dotplot.R             # Biological cell-type labels & clean DotPlot
│   └── 05_cellchat_neuro_immune.R                 # CellChat ligand-receptor crosstalk
└── outputs/
    ├── figures/
    │   ├── 02_qc_postfilter_violin.png
    │   ├── 03_umap_clusters_and_tissue.png
    │   ├── 04_lineage_markers_dotplot_clean.png
    │   ├── 04_annotated_celltypes_umap.png
    │   ├── 04_chemoresistance_projection_umap.png
    │   ├── 04_top_ml_drivers_featureplot.png
    │   ├── 05_cellchat_chemoresistant_circle.png
    │   └── 05_cellchat_ligand_receptor_bubble.png
    └── tables/
        ├── 02_qc_sample_audit.csv
        ├── 04_GSE25066_translated_gene_bridge.csv
        ├── 04_cluster_chemoresistance_scores.csv
        ├── 05_cellchat_all_interactions.csv
        └── 05_chemoresistant_tumor_to_immune_signaling.csv
```











































