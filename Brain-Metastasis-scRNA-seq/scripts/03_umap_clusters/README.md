# Milestone 03: Low-Memory Normalization, Harmony Batch Integration & Macro-Cluster UMAP 

## 🎯 Objective
Transform the QC-filtered single-cell matrix (`01_filtered_merged_seurat.rds`, **13,186 cells across 9 patient samples**) into a batch-corrected, low-dimensional cellular atlas that clearly separates **malignant brain metastasis compartments** from **adjacent normal brain parenchyma** and **shared host immune/stromal niches**.

---

## ⚙️ Step-by-Step Codebase Architecture

| Step                                        | Function / Method                                                                 | What It Does & Why It Was Engineered This Way                                                                                                                                                                                   |
| :------------------------------------------ | :-------------------------------------------------------------------------------- | :------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------ |
| **1. Tissue Origin Tagging**                | `grepl("T$", sample_id)`                                                          | Classifies samples ending in `T` (`BT`, `DT`, `ET`, `FT`, `GT`) as **`Brain_Metastasis_Tumor`** (*n* = 5) and samples ending in `N`/`N1`/`N2` (`BN1`, `BN2`, `DN`, `EN`) as **`Adjacent_Brain_Tissue`** (*n* = 4).              |
| **2. Sparse Normalization & HVG Selection** | `NormalizeData(scale.factor = 10000)`<br>`FindVariableFeatures(nfeatures = 2000)` | Corrects for cell-to-cell sequencing depth differences and isolates the **top 2,000 highly variable genes (HVGs)** carrying real biological signal across cell types.                                                           |
| **3. RAM-Safe Feature Scaling**             | `ScaleData(features = VariableFeatures, vars.to.regress = "percent.mt")`          | Scales **only** the 2,000 variable features (preventing 8GB RAM overflow) while mathematically regressing out residual mitochondrial expression (`percent.mt`).                                                                 |
| **4. PCA & Harmony Batch Correction**       | `RunPCA(npcs = TARGET_PCS)`<br>`RunHarmony(group.by.vars = "sample_id")`          | Removes technical patient-to-patient batch effects across the 9 surgical resections so shared host cells co-cluster by biological lineage rather than donor ID.                                                                 |
| **5. Macro-Resolution Clustering**          | `FindNeighbors(dims = 1:20)`<br>`FindClusters(resolution = 0.2)`                  | Upgrades granularity from `0.5` (17 fragmented micro-clusters) to **`0.2` (11 clean macro-clusters, `0–10`)**, grouping cells into broad, interpretable anatomical compartments.                                                |
| **6. Largest-Island Centroid Anchoring**    | `group_by(Sample_Clean, Macro_Cluster) %>% slice_max(n_cells, n = 1)`             | Computes each patient's `(UMAP 1, UMAP 2)` badge coordinates using the median of its **single largest cluster** rather than a global average—preventing labels from floating in empty white space between disconnected islands. |
| **7. 3-Panel Publication Figure**           | `patchwork` (`p_macro`, `p_direct_labels`, `p_prop`)                              | Combines macro-cluster boundaries, binary tissue-colored patient localization, and exact per-cluster tissue percentages into **`03_umap_clusters.png`**.                                                                        |

---

## 📊 Visual & Biological Output Analysis 

### Panel A: Broad Cellular Clusters (`Resolution = 0.2`)
* Partitions all **13,186 cells** into **11 distinct macro-clusters (`Clusters 0–10`)** with non-overlapping boxed labels (`repel = TRUE`).
* Tight, well-separated islands confirm clean transcriptional separation between malignant epithelium, CNS glia/neurons, and infiltrating immune populations.

### Panel B: Direct Sample Mapping (`Red = Normal Brain | Blue = Brain Tumor`)
* Replaces hard-to-read 9-color rainbow palettes with a **high-contrast binary palette**:
  * **Coral Red (`#E41A1C`):** Adjacent Normal Brain cells (`BN1`, `BN2`, `DN`, `EN`).
  * **Royal Blue (`#1F78B4`):** Brain Metastasis Tumor cells (`BT`, `DT`, `ET`, `FT`, `GT`).
* Stamped patient badges (`ET`, `FT`, `GT`, `DT`, `DN`, `BN2`, etc.) anchor directly onto each sample's dominant cellular island, showing **patient-specific tumor subclusters** on the left/bottom (`Blue`) alongside **multi-patient integrated host islands** on the right/top (`Red` + `Blue`).

### Panel C: Tissue Composition per Cluster (`% Normal vs. % Tumor`)
Quantifies the exact cellular makeup of every macro-cluster, dividing the ecosystem into **three clear biological compartments**:

| Compartment Category                   | Macro-Clusters                       | Tissue Breakdown                                                                                                  | Biological Interpretation                                                                                                                      |
| :------------------------------------- | :----------------------------------- | :---------------------------------------------------------------------------------------------------------------- | :--------------------------------------------------------------------------------------------------------------------------------------------- |
| **Tumor-Dominated Islands**            | **Clusters `9`, `0`, `3`, `7`, `8`** | **~65% – 98% Brain Tumor (Blue)**<br>*(Cluster 9: ~98% Tumor; Clusters 0 & 3: ~86% Tumor)*                        | **Malignant Metastatic Epithelial Subclusters:** Disseminated breast cancer cells absent from normal brain tissue.                             |
| **Normal-Brain-Dominated Islands**     | **Clusters `10`, `6`, `2`, `1`**     | **~73% – 100% Normal Brain (Red)**<br>*(Cluster 10: 100% Normal; Cluster 6: ~95% Normal; Cluster 2: ~88% Normal)* | **Resident CNS Parenchyma:** Oligodendrocytes, astrocytes, and cortical neurons enriched in adjacent non-malignant brain.                      |
| **Shared Tumor–Host Microenvironment** | **Clusters `4`, `5`**                | **Blended Blue + Red**<br>*(Cluster 4: ~60% Tumor / 40% Normal; Cluster 5: ~33% Tumor / 67% Normal)*              | **Infiltrating Immune & Vascular Niche:** Microglia/macrophages, T cells, and endothelial cells co-integrated across normal and tumor tissues. |


