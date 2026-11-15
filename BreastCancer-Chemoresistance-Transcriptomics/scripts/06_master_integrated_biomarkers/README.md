# 1. The Complete Biological Story (Milestones 3 + 4 + 5 Together)

- Milestone 3 (limma, 2,938 DEGs): Scans all 22,283 probes and separates the tumor's genes into Resistance Accelerators (-logFC) that turn ON in RD tumors vs. Sensitivity Brakes (+logFC) that are lost in RD tumors.
- Milestone 4 (clusterProfiler GO Pathways): Explains HOW the tumor survives chemotherapy as a system—revealing three core survival programs (Count = 43–123 genes working together):
	1. Uncontrolled Cell-Cycle Division (GO:0044772, SKP1, CCND1)
	2. Epithelial-to-Mesenchymal Transition / EMT (GO:0001837, GATA3)
	3. Hypoxic Glycolytic Energy Survival (GO:0001666, PGK1)
- Milestone 5 (Random Forest Gini ML, AUC = 0.9033): Proves that you do not need to test 2,938 genes in the clinic. Just 10 non-redundant biomarkers capture both the direct drug-blocking bottlenecks (TMSB15A, TTLL4, EPHB3) AND the three GO survival programs (SKP1, GATA3, PGK1, DSC2, TXNDC5).

# 2. Master integrated summary table

| Rank (ML) | Gene Symbol (Probe_ID) | Milestone 5: ML Gini Score | Milestone 3: limma logFC & adj.P.Val |    Phenotype Status     |   Milestone 4: Linked GO Biological Program (p.adjust)   |                 Clinical & Translational Role                 |
| :-------: | :--------------------: | :------------------------: | :----------------------------------: | :---------------------: | :------------------------------------------------------: | :-----------------------------------------------------------: |
|    #1     | TMSB15A (205347_s_at)  |          0.006574          |           -1.40 (1.80e-08)           | High in Resistant (RD)  |   Orthogonal Cytoskeletal Switch (Actin sequestration)   |      #1 Diagnostic Predictor & Primary Knockdown Target       |
|    #2     |  TTLL4 (203702_s_at)   |          0.005051          |           -1.98 (1.04e-08)           | High in Resistant (RD)  |  Orthogonal Spindle Switch (Tubulin polyglutamylation)   |     #2 Diagnostic Predictor & Direct Taxane-Shield Target     |
|    #3     |    EPHB3 (1438_at)     |          0.004844          |           -0.94 (1.36e-05)           | High in Resistant (RD)  |  Orthogonal Signaling Switch (Ephrin adhesion receptor)  |     #3 Diagnostic Predictor (Invasion/Survival signaling)     |
|    #4     |   TBC1D9 (212956_at)   |          0.004739          |           +1.76 (2.21e-07)           | High in Sensitive (pCR) | Orthogonal Trafficking Switch (Lost in resistant tumors) |      Chemosensitivity Sentinel (Loss signals resistance)      |
|    #5     |   DSC2 (204750_s_at)   |          0.004364          |           -1.46 (2.85e-07)           | High in Resistant (RD)  |  Positive regulation of protein localization (1.32e-05)  |             Junctional Adhesion Resistance Driver             |
|    #6     |   SKP1 (200718_s_at)   |          0.004237          |           +0.39 (6.65e-06)           | High in Sensitive (pCR) |      Mitotic cell cycle phase transition (3.39e-14)      | Cell-Cycle Pathway Sentinel (116-gene cluster representative) |
|    #7     |  GATA3 (209604_s_at)   |          0.004216          |           +1.83 (4.33e-08)           | High in Sensitive (pCR) |     Epithelial to mesenchymal transition (1.15e-04)      |           Luminal Differentiation & EMT Gatekeeper            |
|    #8     |  KCNK5 (219615_s_at)   |          0.004069          |           -0.57 (6.77e-04)           | High in Resistant (RD)  |  Orthogonal Ion-Channel Switch (Potassium homeostasis)   |         Membrane Potential & Apoptosis-Evasion Marker         |
|    #9     |  TXNDC5 (221253_s_at)  |          0.004044          |           -0.47 (3.14e-05)           | High in Resistant (RD)  |          Protein folding / ER stress (5.99e-03)          |               ER Unfolded-Protein Stress Shield               |
|    #10    |    PGK1 (200737_at)    |          0.003943          |           -0.41 (5.14e-04)           | High in Resistant (RD)  |       Response to hypoxia & glycolysis (3.59e-06)        |     Metabolic Warburg Sentinel (65-gene hypoxia cluster)      |










