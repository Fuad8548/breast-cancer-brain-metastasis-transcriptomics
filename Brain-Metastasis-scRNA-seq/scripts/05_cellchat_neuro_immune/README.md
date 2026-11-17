# Step-by-Step Walkthrough of 05_cellchat_neuro_immune.R

## Step 1: Subsetting to Active Interacting Lineages

- `subset(bcbm, subset = cell_type %in% target_populations)`: Keeps the 8 interacting tumor, immune, and blood–brain barrier (BBB) cell populations while excluding static myelin (`Oligodendrocytes`) and normal capillary controls (`BBB_Endothelial_2`) so signaling probabilities reflect active tumor–niche crosstalk.
- `rm(bcbm)`; `gc()`: Frees the full Seurat object from RAM immediately before building the CellChat matrix.

## Step 2: Initializing CellChat & Filtering the Ligand–Receptor Database

- `GetAssayData(bcbm_sub, assay = "RNA", layer = "data")`: Extracts the log-normalized single-cell expression matrix from Seurat v5.
- `subsetDB(CellChatDB.human, search = c("Secreted Signaling", "Cell-Cell Contact"))`: Restricts the human ligand–receptor database to soluble cytokines/growth factors (MIF, VEGFA, MDK) and direct membrane contact pairs (APP–CD74), excluding structural ECM-only pairs.

## Step 3: Statistical Inference of Signaling Probabilities

- `identifyOverExpressedGenes()` & `identifyOverExpressedInteractions()`: Finds ligands and receptors significantly overexpressed in each of the 8 cell types.
- `computeCommunProb(cellchat, type = "triMean", population.size = TRUE)`: Calculates communication probability using the triMean (a robust, outlier-resistant average) while weighting by cell population size (`population.size = TRUE`) so abundant cell populations accurately reflect total ligand output.
- `filterCommunication(min.cells = 10)`: Removes noisy interactions supported by fewer than 10 cells in either the sender or receiver group.

## Step 4: Filtering Tumor-to-Niche Axes & Plotting

- `tumor_sources → immune_targets`: Filters df_net to isolate unidirectional signals sent from the 3 malignant breast tumor populations to the 5 brain immune and BBB populations, sorted by `pval` (0) and descending `prob`.
- `netVisual_circle()` & `netVisual_bubble()`: Draws the weighted network circle plot (05_cellchat_tumor_circle.png) and the dot-matrix ligand–receptor bubble plot (05_cellchat_ligand_receptor_bubble.png).


# Interpretation of the output

## How to read `05_tumor_to_neuro_immune_signaling.csv` and `05_cellchat_all_interactions.csv`.

### 1. How to Read the Columns in Both CSV Tables

- source & target: The sender cell population (source) releasing or displaying the signal, and the receiver cell population (target) expressing the matching receptor.
- ligand & receptor: The exact interacting molecule pair (e.g., MIF binding to the CD74_CXCR4 receptor complex).
- prob (Communication Probability): Quantitative interaction strength computed from normalized ligand/receptor expression and cell population size (triMean method); higher values indicate stronger predicted signaling.
- pval: Permutation test significance (0 means p < 0.01, the highest statistical confidence in CellChat).
- annotation: Distinguishes Secreted Signaling (diffusible cytokines/growth factors acting over short distances) from Cell-Cell Contact (membrane-bound proteins requiring direct physical cell contact).


### 2. Table 1: `05_tumor_to_neuro_immune_signaling.csv` — Tumor-Driven Niche Reprogramming

This table isolates 48 significant interactions (pval = 0) sent directly from our three breast tumor subpopulations (Malignant_Epithelial_1, Malignant_Epithelial_2, Malignant_Plastic_EMT) to the surrounding brain immune and vascular niche. It reveals four major biological mechanisms:

#### A. Dominant CD74-Centered Immune Suppression (Malignant_Epithelial_1 → Immune & BBB Cells)

Seven of the top nine strongest signals in the entire table converge on CD74 (pval = 0), establishing `Malignant_Epithelial_1` as the primary **immunosuppressive** hub:
- APP → CD74 (Ranks #1, #2, #3 overall): Sent from `Malignant_Epithelial_1` to `BBB_Endothelial_1` (prob = 0.002615), `T_NK_Lymphocytes` (prob = 0.002278), and `Microglia_TAMs` (prob = 0.001535) via direct Cell-Cell Contact. Because `CD74` is the invariant chain chaperone required for **MHC Class II antigen presentation**, saturating CD74 blocks antigen presentation on microglia/macrophages and dampens T/NK cell activation.
- MIF → CD74 + CXCR4 & CD74 + CD44 (Ranks #4, #5, #8, #9): `Malignant_Epithelial_1` secretes Macrophage Migration Inhibitory Factor (`MIF`) onto `T_NK_Lymphocytes` (prob = 0.001297 and 0.000956) and `Microglia_TAMs` (prob = 0.000384 and 0.000346). This **reprograms brain-resident microglia into tumor-supportive (M2-like) macrophages while inducing T-lymphocyte dysfunction**.
- `Malignant_Epithelial_2` Reinforces MIF–CD74 (Ranks #13, #17): Although smaller in cell number, the hyper-metabolic `Malignant_Epithelial_2` subclone sends the exact same MIF → CD74+CXCR4 (prob = 0.000163) and MIF → CD74+CD44 (prob = 0.000120) signals to `T_NK_Lymphocytes`.

#### B. Midkine Survival & Vascular Remodeling (MDK → NCL/ ITGA6_ITGB1)

- MDK → NCL (Ranks #6, #7, #18, #19): `Malignant_Epithelial_1` secretes Midkine (MDK) to Nucleolin (NCL) on `T_NK_Lymphocytes` (prob = 0.000508), `BBB_Endothelial_1` (prob = 0.000489), `Microglia_TAMs` (prob = 0.000116), and `Pericytes_Vascular` (prob = 0.000111). MDK–NCL signaling promotes **tumor cell survival, suppresses cytotoxic lymphocytes, and stimulates endothelial proliferation**.

#### C. Angiogenesis & Blood–Brain Barrier Remodeling by EMT Cells (Malignant_Plastic_EMT → BBB_Endothelial_1)

Unlike `Malignant_Epithelial_1`, dedifferentiated `Malignant_Plastic_EMT` cells specialize in vascular and myeloid remodeling:
- VEGFA → FLT1 (VEGFR1), FLT1_KDR (VEGFR1/2), & KDR (VEGFR2): `Malignant_Plastic_EMT` secretes VEGFA directly onto `BBB_Endothelial_1` (prob = 0.000214, 0.000127, and 0.000066), driving neo-angiogenesis and blood–tumor barrier permeability.
- NAMPT → ITGA5_ITGB1 (Rank #11, prob = 0.000219): Secretes extracellular NAMPT (visfatin) to integrin α₅β₁ on BBB_Endothelial_1, further stimulating endothelial migration and vessel sprouting (or sprouting angiogenesis, the process by which **a tumor builds its own new network of blood vessels from pre-existing capillaries to feed its growth**).

#### D. Myeloid Recruitment & Immune Checkpoint Contact (Malignant_Plastic_EMT → Microglia_TAMs / Resident_Myeloid_Glia)

- ANXA1 → FPR1 / FPR2 / FPR3: `Malignant_Plastic_EMT` secretes Annexin A1 (ANXA1) to formyl peptide receptors on `Resident_Myeloid_Glia` (prob = 0.000220 via FPR1) and `Microglia_TAMs` (prob = 0.000037 via FPR3), recruiting and polarizing myeloid cells toward an anti-inflammatory phenotype.
- CSF1 → CSF1R (prob = 0.000082): `Malignant_Plastic_EMT` secretes Macrophage Colony-Stimulating Factor (CSF1) to Microglia_TAMs, **sustaining tumor-associated macrophage survival**.
- HLA-DRB1 / HLA-DPA1 / HLA-DRA → CD4 & PTPRC → MRC1: Direct cell-cell contact between `Malignant_Plastic_EMT` and `Microglia_TAMs` via MHC-II–CD4 and CD45–MRC1 (CD206), **reinforcing M2 macrophage polarization**.


### 3. Table 2: 05_cellchat_all_interactions.csv — Full Ecosystem Crosstalk (Host-to-Host & Feedback Loops)


| Ecosystem Axis (source → target)                                                                           | Top Ligand → Receptor (prob, pval = 0)                                                                                 | Biological Meaning in the Brain Metastasis Ecosystem                                                                                                                                                                                                            |
| ---------------------------------------------------------------------------------------------------------- | ---------------------------------------------------------------------------------------------------------------------- | --------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------- |
| Vascular Autocrine & Endothelial–Immune Loop<br>(BBB_Endothelial_1 ↔ BBB_Endothelial_1 / T_NK / Microglia) | APP → CD74 (0.01015 & 0.00958)<br>SELE → CD44 (0.00689)<br>PECAM1 → PECAM1 (0.00561)<br>ESAM → ESAM (0.00541)          | Activated BBB_Endothelial_1 cells amplify APP–CD74 signaling onto T_NK_Lymphocytes and Microglia_TAMs, while expressing E-selectin (SELE–CD44) and tight-junction adhesion molecules (PECAM1, ESAM) to trap circulating immune cells at the vessel wall.        |
| Chemokine Trapping & T-Cell Sequestration<br>(BBB_Endothelial_1 / Pericytes → T_NK_Lymphocytes)            | CCL5 → ACKR1 (0.00654)<br>CCL14 / CCL2 → ACKR1 (0.00509, 0.00500)<br>CXCL12 → CXCR4 (0.00260, 0.00115)                 | BBB_Endothelial_1 uses the decoy receptor ACKR1 (DARC) to scavenge inflammatory chemokines (CCL2, CCL5, CCL14, CXCL8) and uses CXCL12–CXCR4 to sequester T_NK_Lymphocytes in the perivascular space—preventing them from infiltrating deep into the tumor core. |
| Host-Amplified Angiogenesis<br>(Microglia_TAMs & Pericytes_Vascular → BBB_Endothelial_1)                   | VEGFA → FLT1 (VEGFR1) (0.000356, 0.000240)<br>NAMPT → ITGA5_ITGB1 (0.00162, 0.00112)<br>ANGPT2 → ITGA5_ITGB1 (0.00130) | Once co-opted by the tumor, Microglia_TAMs actually secrete more VEGFA (0.000356) onto BBB_Endothelial_1 than the tumor cells themselves (0.000214), creating a self-sustaining angiogenic feed-forward loop.                                                   |
| Microglial Checkpoint Suppression of T/NK Cells<br>(Microglia_TAMs → T_NK_Lymphocytes)                     | LGALS9 → PTPRC (CD45) (0.001000)<br>LGALS9 → CD44 (0.000723)<br>LGALS9 → HAVCR2 (TIM-3) (0.000020)                     | Reprogrammed Microglia_TAMs express Galectin-9 (LGALS9) to engage TIM-3 (HAVCR2) and CD45 on T_NK_Lymphocytes, directly triggering T-cell exhaustion and apoptosis.                                                                                             |
| Pericyte–Endothelial BBB Stabilization<br>(BBB_Endothelial_1 ↔ Pericytes_Vascular)                         | JAG1 → NOTCH3 (0.000324)<br>PDGFA → PDGFRB (0.000032)<br>GJC1 → GJA1 (0.000112)                                        | Canonical JAG1–NOTCH3 and PDGFA–PDGFRB mural cell signaling maintains pericyte coverage around remodeled tumor blood vessels.                                                                                                                                   |



### 4. Unified Synthesis: How Milestone 04 + Milestone 05 Fit Together

1. Milestone 04 proved that metastatic breast cancer cells inside the brain split into differentiated epithelial clones (Malignant_Epithelial_1/2, GATA3+/TBC1D9+/DSC2+) and dedifferentiated ER-stressed EMT clones (Malignant_Plastic_EMT, TXNDC5+/PGK1+).
2. Milestone 05 (CellChat) proves that these two tumor states divide the workload of building a therapy-resistant sanctuary:
	- Malignant_Epithelial_1 suppresses T/NK cells and microglia via APP–CD74 and MIF–CD74 ().
	- Malignant_Plastic_EMT drives VEGFA–VEGFR1/2 angiogenesis and recruits Microglia_TAMs (CSF1–CSF1R, ANXA1–FPR1) ().
	- Co-opted Microglia_TAMs and BBB_Endothelial_1 then amplify VEGFA, trap chemokines via ACKR1, sequester T cells via CXCL12–CXCR4, and exhaust T cells via LGALS9–TIM-3 (HAVCR2) ().


## Interpretation of the plots

### 1. 05_cellchat_tumor_circle.png: Who Is the Master Communicator?

- What the visual elements represent:
	- Colored Circles (Nodes): Each represents one cell population, with circle size proportional to cell count (`T_NK_Lymphocytes` in brown and `BBB_Endothelial_1` in red are the largest receiver hubs).
	- Colored Arrows/Edges: Show directional signals sent from the three tumor populations (`Malignant_Epithelial_1` in blue, `Malignant_Epithelial_2` in green, and `Malignant_Plastic_EMT` in purple) to the surrounding brain microenvironment.
	- Line Thickness: Encodes total interaction weight (cumulative signaling strength).
- Key Takeaways:
	- `Malignant_Epithelial_1` (Blue Node, Upper Right) Is the Primary Signaling Hub: Three massive blue cables radiate directly from `Malignant_Epithelial_1` onto `T_NK_Lymphocytes`, `BBB_Endothelial_1`, and `Microglia_TAMs` — proving that differentiated epithelial tumor cores dominate neuro-immune reprogramming.
	- `Malignant_Plastic_EMT` (Purple Node, Upper Left) & `Malignant_Epithelial_2` (Green Node, Top) Act as Targeted Niche Remodelers: Their thinner purple and green edges project selectively onto `BBB_Endothelial_1`, `Microglia_TAMs`, and `Resident_Myeloid_Glia`, specializing in localized vascular and myeloid remodeling rather than broad T-cell engagement.


### 2. 05_cellchat_ligand_receptor_bubble.png: The Exact Molecular Axes

- How to read the matrix:
	- X-axis (source.target): Every sender–receiver pair grouped by tumor sender (Malignant_Epithelial_1 on the left, Malignant_Epithelial_2 in the middle, Malignant_Plastic_EMT on the right).
	- Y-axis (interaction_name_2): The exact Ligand - Receptor pair.
	- Dot Color (min blue → max red): Communication probability (prob), with every solid dot meeting p < 0.01.
- Four Distinct Molecular Modules Visible in the Columns:
	1. APP - CD74 Red Hotspot (Bottom-Left Columns 1, 2, 5): Three bright red/crimson dots sit at APP - CD74 from Malignant_Epithelial_1 to BBB_Endothelial_1, Microglia_TAMs, and T_NK_Lymphocytes, **blocking MHC-II antigen presentation and inducing local immune tolerance**.
	2. MIF - (CD74+CXCR4) & MIF - (CD74+CD44) Spanning All Three Tumor Clones (Upper-Middle Rows): Orange-red dots from `Malignant_Epithelial_1` onto `T_NK_Lymphocytes` (Column 5) and yellow-green dots extending across `Malignant_Epithelial_2` (Columns 6–8) prove that MIF–CD74 signaling is a **conserved pan-tumor mechanism suppressing T/NK cells and polarizing microglia**.
	3. MDK - NCL Survival & Vascular Axis (Middle Rows): Yellow-green dots from `Malignant_Epithelial_1` onto `BBB_Endothelial_1` and `T_NK_Lymphocytes` show Midkine (MDK) engaging Nucleolin (NCL) **to promote endothelial remodeling and immune evasion**.
	4. `Malignant_Plastic_EMT` Specialized Vascular & Myeloid Columns (Rightmost Columns 9–12):
		- Column 9 (Malignant_Plastic_EMT -> BBB_Endothelial_1): Lights up exclusively at the top left with VEGFA - VEGFR1, VEGFA - VEGFR1R2, VEGFA - VEGFR2, and NAMPT - (ITGA5+ITGB1), proving that **EMT-transitioned tumor cells drive neo-angiogenesis at the blood–brain barrier**.
		- Column 10 (Malignant_Plastic_EMT -> Microglia_TAMs): Forms a vertical column of green/blue dots across HLA-DR/DP/DQ - CD4, CSF1 - CSF1R, PTPRC - MRC1, and ICAM1 - ITGAM/ITGAX+ITGB2, demonstrating **direct physical contact and macrophage polarization**.
		- Column 11 (Malignant_Plastic_EMT -> Resident_Myeloid_Glia): Highlights ANXA1 - FPR1 and ANXA1 - FPR2 at the bottom right, **recruiting brain-resident myeloid glia**.









# Ligand Receptor Pair in Brain

|                    Ligand–Receptor Pair                    |                  Sender Cell (Ligand)                   |                                       Receiver Cell (Receptor)                                       |                                                                                                                   What This "Molecular Handshake" Actually Does in the Brain                                                                                                                   |
| :--------------------------------------------------------: | :-----------------------------------------------------: | :--------------------------------------------------------------------------------------------------: | :--------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------: |
| MIF – CD74 (Macrophage Migration Inhibitory Factor → CD74) |        Malignant Breast Cancer Cells secrete MIF        | Brain Microglia & Tumor-Associated Macrophages (TAMs) express CD74 (often paired with CXCR4 or CD44) |               Immune Brainwashing: Normally, brain microglia attack and eat invading cancer cells. When tumor cells flood the brain with MIF and bind CD74 on microglia, it switches the microglia from tumor-killing (M1) into tumor-protecting (M2 immunosuppressive) helpers.               |
|       APP – CD74 (Amyloid Precursor Protein → CD74)        | Malignant Tumor Cells & Reactive Astrocytes express APP |                                 Microglia & Macrophages express CD74                                 | Antigen-Presentation Blockade: CD74 is the master chaperone protein required for immune cells to present tumor antigens (MHC Class II) to killer T cells. When tumor-derived APP binds CD74, it jams antigen presentation, hiding the breast cancer cells from T-cell attack inside the brain. |




## 1. What our New CellChat Results Prove

Looking at the Top 9 strongest tumor-to-niche signaling axes across the entire brain metastasis microenvironment—7 out of the top 9 target CD74 (pval = 0):

| Rank   | Source Cell (Sender)   | Target Cell (Receiver) | Ligand → Receptor                            | Pathway          | Probability (prob)   | P-Value | Biological & Clinical Meaning in Brain Metastasis                                                                                                                                |
| ------ | ---------------------- | ---------------------- | -------------------------------------------- | ---------------- | -------------------- | ------- | -------------------------------------------------------------------------------------------------------------------------------------------------------------------------------- |
| #1     | Malignant_Epithelial_1 | BBB_Endothelial_1      | APP → CD74                                   | APP              | 0.002615             | 0       | Tumor cells bind CD74 on blood-brain barrier endothelium to remodel vascular permeability & adhesion.                                                                            |
| #2     | Malignant_Epithelial_1 | T_NK_Lymphocytes       | APP → CD74                                   | APP              | 0.002278             | 0       | Suppresses effector T/NK lymphocyte activation inside the tumor bed.                                                                                                             |
| #3     | Malignant_Epithelial_1 | Microglia_TAMs         | APP → CD74                                   | APP              | 0.001535             | 0       | Blocks MHC-II antigen presentation on brain microglia/macrophages so T cells cannot recognize tumor antigens.                                                                    |
| #4     | Malignant_Epithelial_1 | T_NK_Lymphocytes       | MIF → CD74 + CXCR4                           | MIF              | 0.001297             | 0       | Drives immunosuppressive chemokine signaling and T-cell dysfunction.                                                                                                             |
| #5     | Malignant_Epithelial_1 | T_NK_Lymphocytes       | MIF → CD74 + CD44                            | MIF              | 0.000956             | 0       | Reinforces T-cell exhaustion/anergy via CD74/CD44 co-receptor signaling.                                                                                                         |
| #6     | Malignant_Epithelial_1 | T_NK_Lymphocytes       | MDK → NCL                                    | MK               | 0.000508             | 0       | Midkine (MDK)-Nucleolin (NCL) signaling promoting tumor survival and immune suppression.                                                                                         |
| #7     | Malignant_Epithelial_1 | BBB_Endothelial_1      | MDK → NCL                                    | MK               | 0.000489             | 0       | Midkine-driven endothelial remodeling & angiogenesis at the blood-tumor barrier.                                                                                                 |
| #8     | Malignant_Epithelial_1 | Microglia_TAMs         | MIF → CD74 + CXCR4                           | MIF              | 0.000384             | 0       | Reprograms brain microglia into pro-tumorigenic M2 macrophages.                                                                                                                  |
| #9     | Malignant_Epithelial_1 | Microglia_TAMs         | MIF → CD74 + CD44                            | MIF              | 0.000346             | 0       | Sustains M2 microglial survival and anti-inflammatory cytokine release.                                                                                                          |
| #11–12 | Malignant_Plastic_EMT  | BBB_Endothelial_1      | NAMPT → ITGA5+ITGB1<br>VEGFA → FLT1 (VEGFR1) | VISFATIN<br>VEGF | 0.000219<br>0.000214 | 0       | Division of labor: While Malignant_Epithelial_1 suppresses immunity via APP/MIF–CD74, dedifferentiated Malignant_Plastic_EMT cells drive angiogenesis (VEGFA–VEGFR1) at the BBB! |