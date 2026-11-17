Milestone 04 is the biological bridge between our two projects. It takes the `2,938` chemoresistance genes we discovered in bulk primary breast tumors (GSE25066) and projects them onto every single cell of our brain metastasis atlas (GSE268662) to prove which exact cell populations inside the brain carry the primary chemoresistance program.

# Step-by-Step Codebase Walkthrough 

## 1. Cross-Platform Probe Translation (GSE25066 → GSE268662)

- `read.csv(deg_files[1], ...)`: Loads our limma DEG table (03_Biomarkers.csv) from the primary breast cancer study (GSE25066).
- `mapIds(hgu133a.db, ...)`: Because GSE25066 used an older Affymetrix microarray (205347_s_at), while GSE268662 uses modern single-cell RNA-seq gene symbols (TMSB15A), `mapIds()` translates every Affymetrix Probe ID into its official HUGO Gene Symbol and keeps only genes actually detected in the brain metastasis cells (`Gene_Symbol %in% rownames(bcbm)`).
- resistant_genes vs. sensitive_genes: Splits the mapped genes into two opposing signatures:
	- Chemoresistant Signature (logFC < 0): Genes upregulated in RD tumors (including our top Random Forest drivers TMSB15A, TTLL4, EPHB3, PGK1).
	- Chemosensitive Signature (logFC > 0): Genes upregulated in pCR tumors (including TBC1D9, SKP1, GATA3).

## 2. Canonical Lineage Marker DotPlot 

- `DotPlot(bcbm, features = lineage_markers)`: Plots 18 classic "identity badge" genes across Clusters 0–10:
	- Epithelial Tumor: EPCAM, KRT8, KRT18, KRT19
	- Microglia / Macrophages: PTPRC, C1QA, AIF1, P2RY12
	- T/ NK Lymphocytes: CD3D, CD3E, NKG7
	- Blood-Brain Barrier (BBB) Vessels: PECAM1, VWF, RGS5
	- Normal Brain Glia: MBP, PLP1, GFAP, AQP4
- How to read the DotPlot (04_lineage_markers_dotplot.png): Dot size = percentage of cells in that cluster expressing the gene; Dot color (grey → dark red) = average expression level. This is how we proved Cluster 3 and 9 are epithelial tumors, Cluster 4 is microglia/TAMs, and Cluster 0 is T/NK cells.

## 3. Single-Cell Chemoresistance Scoring (AddModuleScore)

- `AddModuleScore(...)`: For each of the 13,186 individual cells, `Seurat` calculates the average expression of all resistant_genes minus randomly selected background control genes with matching expression levels—producing a normalized GSE25066_Resistant_Score1 for every single cell.
	- Positive score (> 0) = This individual cell actively expresses the primary breast cancer chemoresistance program.
	- Negative/zero score (<= 0) = This cell lacks the chemoresistance signature.


# Interpretation of the output

## 1. 04_annotated_umap_with_markers.png
Eliminates the need to squint at DotPlot percentages—every island on the UMAP carries its biological name + exact marker proof (EPCAM+, KRT18+, C1QA+, CD3D+, PECAM1+, RGS5+, MBP+) right inside its label box.
## 2. 04_chemoresistance_by_celltype.png
Shows at a glance which named cell types sit above the dashed 0 baseline (Malignant_Epithelial_1 and Malignant_Epithelial_2 in blue) versus normal brain cells below 0.
## 3. 04_top_ml_drivers_featureplot.png
Uses a high-contrast grey → orange → dark red scale with descriptive titles (#1 RF Driver: TMSB15A, #2 RF Driver: TTLL4, #3 RF Driver: EPHB3, #10 RF Sentinel: PGK1) so we immediately see each primary chemoresistance gene lighting up inside the metastatic brain tumor islands.

- **Why PGK1 was included**
	Alongside our Top 3 Random Forest Gini drivers (TMSB15A #1, TTLL4 #2, EPHB3 #3), PGK1 (#10 in our Random Forest ranking) is our hypoxic glycolysis metabolic sentinel.
  - Why PGK1 lights up heavily (Expr up to 4) while the other 3 appear as scattered red/orange dots (Expr 0–2.5):
    	1. mRNA Copy-Number Biology (Enzyme vs. Signaling Regulator): PGK1 (Phosphoglycerate Kinase 1) is a core glycolytic enzyme transcribed at thousands of mRNA copies per cell. Inside the hypoxic brain metastasis tumor islands (bottom-left and bottom-center islands), cancer cells crank up glycolysis (the Warburg effect), so PGK1 glows dark orange/red across almost every tumor cell.
    	2. 10x scRNA-seq "Dropout" on Low-Copy Regulators: TMSB15A, TTLL4, and EPHB3 are low-copy cytoskeletal/receptor signaling genes (only a few mRNA molecules per cell). In 10x droplet single-cell sequencing, low-abundance transcripts are captured in ~5–15% of cells (showing up as distinct red/orange dots rather than blanket coverage).
    	3. Why This Proves the Value of Your Multi-Gene Module Score: Because individual signaling drivers have single-cell dropout, combining all GSE25066 resistance genes into one composite module score (Resistance_Z_Score) overcomes single-gene dropout and reveals the full chemoresistant tumor population cleanly.

## 4. Breast Cancer ML Biomarkers, Spatial Cell Types and Brain Lineage Identity Markers Bridge integrated plot (04_complete_master_bridge.png)

Looking strictly at Milestone 04 (connecting our `GSE25066` primary breast cancer Random Forest biomarkers to the `GSE268662` brain metastasis cell types), here is the complete biological story told by the data in this plot:

### 1. Validation of Cellular Lineage Identity (Panel B, Left Block: Brain & Lineage Proof)

Before testing breast cancer chemoresistance genes, the left panel proves that our clustering cleanly separated metastatic tumor cells from native brain structures:
- Epithelial Breast Carcinoma (`Malignant_Epithelial_1 & 2`): Strongly express `EPCAM` and `KRT18` (cytokeratin 18). Because normal brain tissue never expresses epithelial keratins or EPCAM, this confirms these two Blue islands are true **breast-derived carcinoma** cells.
- Dedifferentiated Tumor State (`Malignant_Plastic_EMT`): Sits inside the Tumor-Dominated Niche (Blue), yet shows near-complete loss of `EPCAM` and `KRT18`. This is the classic signature of Epithelial-to-Mesenchymal Transition (EMT), where metastatic cells shed their rigid epithelial identity to become **migratory** and **therapy-tolerant**.
- Native Brain & Immune Compartments (Purple & Red Niches): Each host lineage is locked onto its own canonical identity gene with zero epithelial contamination: `CD3D` (`T_NK_Lymphocytes`), `C1QA` (`Microglia_TAMs`), `RGS5` (`Pericytes_Vascular`), `PECAM1` (`BBB_Endothelial_1/2`), and `MBP` (`Oligodendrocytes`).

### 2. How Primary Breast Cancer ML Biomarkers Behave Inside the Brain (Panel B, Right Block)

Once cell identities are established, scanning across the 9 GSE25066 Random Forest biomarkers reveals four distinct biological programs operating inside the brain lesion:

#### Discovery A: Retention vs. Loss of Luminal Breast Differentiation (GATA3, TBC1D9, DSC2)

- Where they appear: High in `Malignant_Epithelial_1` and `Malignant_Epithelial_2`, completely absent in `Malignant_Plastic_EMT`, and completely absent in all host brain cells.
- Biological meaning:
	- GATA3 is the master transcription factor that maintains luminal breast epithelial identity.
	- TBC1D9 is a calcium-binding/Rab-GAP regulator tightly co-expressed with luminal differentiation and favorable chemotherapy response (pCR) in primary breast tumors.
	- DSC2 (Desmocollin-2) encodes a desmosomal cell-cell adhesion junction protein that holds epithelial tumor sheets together.
- The Takeaway: Brain metastases do not exist as a single uniform tumor mass. Instead, they split into two co-existing states:
	1. Differentiated Luminal-Like Clones (`Malignant_Epithelial_1/2`): Retain primary breast architecture (`DSC2+`, `GATA3+`, `TBC1D9+`).
	2. Dedifferentiated EMT Clones (`Malignant_Plastic_EMT`): Silence `GATA3`, `TBC1D9`, and `DSC2` — dissolving cell-cell junctions (DSC2 loss) and abandoning the chemosensitive luminal program (TBC1D9/GATA3 loss) to survive in the brain.

#### Discovery B: Endoplasmic Reticulum (ER) Stress Adaptation in EMT Cells (`TXNDC5`)

- Where it appears: Peaks sharply in `Malignant_Plastic_EMT` (large orange-red dot).
- Biological meaning: `TXNDC5` (Thioredoxin Domain Containing 5) is an endoplasmic reticulum chaperone induced by hypoxia and oxidative stress that prevents protein misfolding and blocks apoptosis.
- The Takeaway: When `Malignant_Plastic_EMT` cells shed their epithelial protection (EPCAM-, DSC2-, GATA3-), they survive the harsh brain microenvironment by activating TXNDC5-mediated ER-stress buffering.

#### Discovery C: Hypoxic Glycolytic Reprogramming (`PGK1` as Metabolic Baseline & Tumor Amplifier)

- Where it appears: Expressed across all 10 cell types (proving healthy mRNA capture across Blue, Purple, and Red niches), but reaches its highest intensity (Scaled Expr dark red) in `Malignant_Epithelial_1/2`.
- Biological meaning: `PGK1` catalyzes ATP generation in glycolysis. Rapidly dividing brain metastases outgrow local blood vessels and face severe hypoxia, forcing Malignant_Epithelial_1/2 to hyper-activate PGK1 (the Warburg effect) while normal oligodendrocytes (MBP+) and resting vessels (PECAM1+) maintain low baseline expression.

#### Discovery D: Single-Cell Dropout vs. Composite Module Enrichment (`TMSB15A`, `TTLL4`, `EPHB3`)

- Where they appear: Show up as small, focal dots (`TMSB15A` dark red in ~5% of `Malignant_Epithelial_2`; `TTLL4` and `EPHB3` sparse across tumor/stroma).
- Biological meaning: Unlike high-abundance structural keratins (`KRT18`) or metabolic enzymes (`PGK1`), signaling receptors (`EPHB3`) and cytoskeletal modifiers (`TMSB15A`, `TTLL4`) operate at low mRNA copy numbers per cell. In 10x droplet scRNA-seq, low-copy transcripts experience technical dropout at the single-gene level—which is exactly why aggregating all ~1,400 GSE25066 resistance genes via `AddModuleScore` (`Resistance_Z`) was essential to capture the full chemoresistance signal across the tumor islands.


#### Discovery E: Broadly Detected Housekeeping/Cell-Cycle Machinery (`SKP1`)
- Like PGK1, SKP1 shows up across almost every row (25%–50% of cells) because every viable cell needs basal ubiquitin-proteasome turnover to stay alive.
- Peaks in Proliferative Tumor (`Malignant_Epithelial_1/2`) & Vascular/Glial Stroma (`BBB_Endothelial_1`, `Pericytes_Vascular`, `Oligodendrocytes`): Notice its bright orange-red intensity in Malignant_Epithelial_1 (the actively dividing tumor subclone) and remodeling blood-brain barrier vessels (BBB_Endothelial_1 / Pericytes_Vascular), whereas Malignant_Plastic_EMT shows slightly lower SKP1 expression — consistent with EMT tumor cells slowing down proliferation (entering a slow-cycling, TXNDC5-buffered survival state) to resist chemotherapy!

# How Raw_Resistant_Score, Resistance_Z_Score, and Sensitivity_Z_Score Work

| Metric Column in CSV | Formula / How R Computes It                                                                                                              | Why It Matters & How to Read It                                                                                                                                                              |
| -------------------- | ---------------------------------------------------------------------------------------------------------------------------------------- | -------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------- |
| Raw_Resistant_Score  | Seurat AddModuleScore(): Average expression of your GSE25066 chemoresistance genes minus 100 expression-matched random background genes. | Raw score (~0.03 to 0.14). Its baseline is not centered at zero, making raw numbers hard to compare visually.                                                                                |
| Resistance_Z_Score   | scale(Raw_Resistant_Score) = \\(\\frac{\\text{Cell\\ Score}-\\text{Cohort\\ Mean}}{\\text{Cohort\\ SD}}\\)                               | Standardized score centered at 0.00:<br>• 0.00 = Exact average across all 13,186 cells<br>• > 0 (Positive) = Enriched for chemoresistance<br>• < 0 (Negative) = Depleted (normal background) |
| Sensitivity_Z_Score  | scale(Raw_Sensitive_Score) (Z-score of the opposite pCR-enriched genes from GSE25066).                                                   | Built-in negative control: Shows whether a cell population expresses genes associated with chemosensitivity (pCR) vs. chemoresistance (RD).                                                  |
















