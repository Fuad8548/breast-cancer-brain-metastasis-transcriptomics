# Bulk Microarray Transcriptomic & Machine Learning Biomarker Discovery for Taxane-Anthracycline Chemoresistance in Breast Cancer (GSE25066)

## 📌 Project Overview
Neoadjuvant taxane-anthracycline chemotherapy is a cornerstone treatment for locally advanced breast cancer, yet **~80% of patients fail to achieve a pathological Complete Response (pCR)**, leaving Residual Disease (RD) with high relapse risk. 

This end-to-end **hybrid R + Python computational genomics pipeline** analyzes **306 clinical breast tumor biopsies** from NCBI GEO (**GSE25066**, Affymetrix Human Genome U133A / `GPL96`)—comprising **249 Chemoresistant (RD)** and **57 Chemosensitive (pCR)** tumors—to:
1. Identify statistically significant transcriptional drivers of treatment response using empirical Bayes linear modeling (`limma`).
2. Dissect opposing biological programs active in **Chemoresistant (`RD-High`)** vs. **Chemosensitive (`pCR-High`)** tumors using split Gene Ontology enrichment (`clusterProfiler::compareCluster`).
3. Distill a compact, non-redundant **10-gene diagnostic biomarker panel** using a class-balanced **Random Forest classifier (`ROC-AUC = 0.9033`)** that bridges core cell-cycle/metabolic pathway sentinels with orthogonal cytoskeletal and microtubule taxane-evasion targets.

---

## 🧬 End-to-End 5-Milestone Architecture

```text
NCBI GEO (GSE25066: 306 Patients | 249 RD vs. 57 pCR)
       │
       ▼
[Milestone 1 & 2: R / GEOquery] ──► Sample Alignment, Metadata Curation & QC
       │
       ▼
[Milestone 3: R / limma] ─────────► 2,938 Significant DEGs (adj.P.Val < 0.05)
       │                            ├─ 1,218 Chemoresistant-High (-logFC)
       │                            └─ 1,009 Chemosensitive-High (+logFC)
       │
       ├──► [Milestone 4: R / clusterProfiler + hgu133a.db]
       │      Split GO Enrichment (compareCluster):
       │      ├─ RD-High: Chromosome Segregation (p.adj = 3.98e-22) & Mitotic Transition
       │      └─ pCR-High: Fatty Acid Metabolism & Microtubule Transport (p.adj = 7.61e-05)
       │
       └──► [Milestone 5: Python / scikit-learn]
              Stratified 80/20 Random Forest (500 Trees, Class-Balanced):
              ├─ Holdout Test ROC-AUC = 0.9033 (Resistant Precision = 0.94)
              └─ Top 10 Gini Biomarkers + ML-to-GO Pathway Integration
```

---

## 🔬 Key Biological & Machine Learning Discoveries

### 1. Opposing Functional Programs (`Milestone 4: compareCluster GO-BP Pathway Enrichment`)
Splitting the 2,938 significant transcripts (`adj.P.Val < 0.05`) by fold-change direction revealed two sharply contrasting cellular states:
- **Chemoresistant (`RD-High`, `logFC < 0`, 1,218 genes):** Overwhelmingly enriched for **chromosome segregation (`GO:0007059`, `p.adjust = 3.98e-22`, Count = 92)**, **mitotic cell cycle phase transition (`GO:0044772`, `p.adjust = 4.17e-20`)**, **double-strand DNA repair (`GO:0006302`)**, and **response to hypoxia (`GO:0001666`, `p.adjust = 4.37e-06`)**. Resistant tumors survive taxane/anthracycline exposure by accelerating spindle checkpoints, repairing DNA breaks, and switching to hypoxic glycolysis.

- **Chemosensitive (`pCR-High`, `logFC > 0`, 1,009 genes):** Enriched for **fatty acid metabolic process (`GO:0006631`, `p.adjust = 7.61e-05`, Count = 49)**, **microtubule-based transport (`GO:0099111`, `p.adjust = 7.61e-05`, `MAPT`)**, and **protein localization to membrane (`GO:0072657`)**. Sensitive tumors retain well-differentiated luminal metabolism and intact microtubule tracks, making them directly vulnerable to taxane-mediated microtubule stabilization.

### 2. Top 10 Predictive Biomarkers (`Milestone 5: Random Forest Gini Importance`)
Using a stratified 80/20 train-test split with `class_weight="balanced"`, the Random Forest model achieved a **Holdout Test ROC-AUC of 0.9033** (**0.94 precision** on Resistant tumors). Integrating **Gini Importance** (multivariate diagnostic separability) with **`limma` `logFC`** (phenotype direction) and **GO enrichment** uncovered two complementary biomarker classes:

| Gini Rank | Gene Symbol   | Affymetrix Probe ID | Gini Importance | `limma` `logFC` (pCR vs RD) | `limma` `adj.P.Val` | Phenotype Status           | Linked Biological Program & Translational Role                                                                                        |
| :-------: | :------------ | :------------------ | :-------------: | :-------------------------: | :-----------------: | :------------------------- | :------------------------------------------------------------------------------------------------------------------------------------ |
|  **#1**   | **`TMSB15A`** | `205347_s_at`       |  **0.006574**   |           `-1.40`           |     `1.80e-08`      | **Resistant-High (`RD`)**  | *Orthogonal Cytoskeletal Switch:* Actin monomer sequestration driving tumor motility; **#1 diagnostic predictor & knockdown target**. |
|  **#2**   | **`TTLL4`**   | `203702_s_at`       |  **0.005051**   |           `-1.98`           |     `1.04e-08`      | **Resistant-High (`RD`)**  | *Orthogonal Spindle Switch:* Tubulin polyglutamylase remodeling microtubules to evade taxane binding (~3.9× upregulated in RD).       |
|  **#3**   | **`EPHB3`**   | `1438_at`           |  **0.004844**   |           `-0.94`           |     `1.36e-05`      | **Resistant-High (`RD`)**  | *Orthogonal Signaling Switch:* Ephrin receptor tyrosine kinase promoting cell-adhesion survival signaling.                            |
|  **#4**   | **`TBC1D9`**  | `212956_at`         |  **0.004739**   |           `+1.76`           |     `2.21e-07`      | **Sensitive-High (`pCR`)** | *Orthogonal Trafficking Switch:* Luminal Rab-GTPase regulator; loss serves as a strong marker of chemoresistance.                     |
|  **#5**   | **`DSC2`**    | `204750_s_at`       |  **0.004364**   |           `-1.46`           |     `2.85e-07`      | **Resistant-High (`RD`)**  | *Protein Localization (`GO:1903829`):* Desmosomal cadherin reinforcing intercellular survival junctions under drug stress.            |
|  **#6**   | **`SKP1`**    | `200718_s_at`       |  **0.004237**   |           `+0.39`           |     `6.65e-06`      | **Sensitive-High (`pCR`)** | *Mitotic Cell Cycle (`GO:0044772`):* SCF ubiquitin ligase core subunit; pathway sentinel for the 116-gene cell-cycle cluster.         |
|  **#7**   | **`GATA3`**   | `209604_s_at`       |  **0.004216**   |           `+1.83`           |     `4.33e-08`      | **Sensitive-High (`pCR`)** | *EMT & Differentiation (`GO:0001837`):* Master luminal transcription factor; loss triggers EMT-driven drug resistance.                |
|  **#8**   | **`KCNK5`**   | `219615_s_at`       |  **0.004069**   |           `-0.57`           |     `6.77e-04`      | **Resistant-High (`RD`)**  | *Orthogonal Ion-Channel Switch:* Two-pore potassium channel regulating membrane potential and apoptosis evasion.                      |
|  **#9**   | **`TXNDC5`**  | `221253_s_at`       |  **0.004044**   |           `-0.47`           |     `3.14e-05`      | **Resistant-High (`RD`)**  | *Protein Folding (`GO:0006457`):* Endoplasmic reticulum disulfide isomerase protecting tumors from chemotherapy-induced ER stress.    |
|  **#10**  | **`PGK1`**    | `200737_at`         |  **0.003943**   |           `-0.41`           |     `5.14e-04`      | **Resistant-High (`RD`)**  | *Response to Hypoxia (`GO:0001666`):* Glycolytic ATP-generating enzyme anchoring the 65-gene hypoxic Warburg survival cluster.        |

> **Why Hybrid R + Python Matters:** Hypergeometric GO enrichment (`clusterProfiler`) captures large multi-gene pathway crowds (`SKP1`, `GATA3`, `PGK1`), whereas Random Forest **Gini Importance** captures both pathway sentinels *and* orthogonal single-gene biochemical bottlenecks (`TMSB15A`, `TTLL4`, `EPHB3`) that single-pathway tests miss.

