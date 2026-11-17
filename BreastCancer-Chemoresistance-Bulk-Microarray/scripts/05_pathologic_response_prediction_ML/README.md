Using Python (`scikit-learn`), we will train an advanced classifier to evaluate whether the transcriptomic profiles of our 2,938 differentially expressed genes hold real clinical utility. The goal is to predict if a newly diagnosed patient will respond completely to treatment (Sensitive / pCR) or fail treatment (Resistant / RD).

From [Export Metrics](export_metrices.R) we get two files:
- [Features](GSE25066_ML_Features.csv)
- [Labels](GSE25066_ML_Labels.csv)



## Ranking of the genes based on Gini Importance, linked with limma logFC

| Rank (Gini) | Gene Symbol | Gini Importance () |  limma logFC ()  | Which Group Is Higher? |                           Biological Meaning                           |
| :---------: | :---------: | :----------------: | :--------------: | :--------------------: | :--------------------------------------------------------------------: |
|     #1      |   TMSB15A   |      0.006574      | -1.40 (Negative) |  Chemoresistant (RD)   |            Active Resistance Driver (cytoskeleton/invasion)            |
|     #2      |    TTLL4    |      0.005051      | -1.98 (Negative) |  Chemoresistant (RD)   |          Active Resistance Driver (microtubule taxane shield)          |
|     #3      |    EPHB3    |      0.004844      | -0.94 (Negative) |  Chemoresistant (RD)   |           Active Resistance Driver (cell-adhesion survival)            |
|     #4      |   TBC1D9    |      0.004739      | +1.76 (Positive) |  Chemosensitive (pCR)  |          Sensitivity Marker (high in responders; lost in RD)           |
|     #5      |    DSC2     |      0.004364      | -1.46 (Negative) |  Chemoresistant (RD)   |             Active Resistance Driver (desmosomal survival)             |
|     #6      |    SKP1     |      0.004237      | +0.39 (Positive) |  Chemosensitive (pCR)  |             Sensitivity Marker (ubiquitin ligase complex)              |
|     #7      |    GATA3    |      0.004216      | +1.83 (Positive) |  Chemosensitive (pCR)  | Sensitivity Marker (luminal differentiation marker; loss = resistance) |


### Comparison based on logFC
Noticing the minus sign:
- **TTLL4** (203702_s_at) has `logFC = -1.98` (almost 4x higher in Resistant tumors).
- **TMSB15A** (205347_s_at) has `logFC = -1.40` (about 2.6x higher in Resistant tumors).
So by logFC magnitude alone, TTLL4 (1.98) actually has a LARGER average jump than TMSB15A (1.40)! Yet in our Random Forest(), TMSB15A (0.006574) beat TTLL4 (0.005051) to become the #1 most important predictor. Why?   

All 2,938 genes in our Random Forest share 100% (1.0000) of the decision-making credit across 500 trees.
- If all 2,938 genes were equally useful, each gene would get $\(1 / 2938 = \mathbf{0.00034}\)$.
- TMSB15A got 0.006574 (nearly 20 times more predictive power than average!) because across almost every single patient biopsy, TMSB15A sits at a clean threshold (> 9.0) that separates Resistant patients from Sensitive patients with very few mistakes.
- TTLL4 got 0.005051 (about 15 times more than average): even though TTLL4 jumps higher (-1.98) in some Resistant patients, it has wider patient-to-patient overlap with Sensitive patients in the middle range, so it makes slightly more classification errors than TMSB15A.


## Gene Profile Comparison (Took 2 patients at 5 genes for example)
```text
Gene Feature	GSM615096 (label 1)	GSM615099 (label 0)	Relative Difference
205225_at	    13.93	            8.19	            Much higher in Class 1
219051_x_at	    9.64	            7.16	            Higher in Class 1
209644_x_at	    9.11	            11.56	            Higher in Class 0
206373_at	    5.81	            9.04	            Higher in Class 0
203702_s_at	    5.10	            9.89	            Much higher in Class 0
```

In this setup, Label 1 represents **Chemo-Resistant** patients (tumor genes active/turned on), while Label 0 represents **Chemo-Sensitive** patients (tumor genes inactive/turned off). With an AUC (Area Under Curve) of `0.90`, this model is highly effective at identifying which patients are likely to resist standard chemotherapy before treatment even begins.














GO  pathway relation =====================================


## Biological Function Profiles for our First 5 ML Probes
These 5 probes are part of our machine learning dataset. Here is their official mapping and how they relate to the mitotic and cell-cycle pathways discovered in our GO enrichment results:
- `219051_x_at` $\(\rightarrow \)$ Gene Symbol: **EIF4G1** (Eukaryotic Translation Initiation Factor 4 Gamma 1)
	- Function: This gene initiates the translation of mRNAs into operational proteins.
	- Tumor Role: Highly aggressive breast cancer cells upregulate EIF4G1 to selectively translate pro-survival and cell-cycle mRNAs, helping the tumor continue growing even when under chemical stress from chemotherapy.
- `209644_x_at` $\(\rightarrow \)$ Gene Symbol: **CCND1** (Cyclin D1)
	- Function: A critical cell cycle regulator. It binds to CDK4/6 to drive the cell past the G1/S restriction point, triggering cellular division.
	- Tumor Role: This connects directly to our top GO pathways (mitotic cell cycle phase transition). Overexpression of Cyclin D1 is a classic marker of chemotherapy resistance, pushing cells to divide regardless of drug blocks.
- `205225_at` $\(\rightarrow \)$ Gene Symbol: **ESR1** (Estrogen Receptor 1)
	- Function: Encodes the Estrogen Receptor alpha (ER$\(\alpha \)$) protein, the key nuclear receptor driving hormone-responsive breast cancer growth.
	- Tumor Role: This gene is the primary biological driver of the Luminal A and Luminal B subtypes we analyzed in our PCA plot. Its baseline expression levels dictate whether a tumor responds to standard chemotherapy or requires hormone-targeted alternatives.
- `206373_at` $\(\rightarrow \)$ Gene Symbol: **PGR** (Progesterone Receptor)
	- Function: Encodes the steroid hormone receptor for progesterone, which interacts closely with ESR1 cascades inside luminal breast tissue cells.
	- Tumor Role: Serves as a key indicator of tumor differentiation. High expression typically signals standard hormone sensitivity, while a sudden loss of PGR expression can indicate a tumor is evolving into a more aggressive, chemoresistant state.
- `203702_s_at` $\(\rightarrow \)$ Gene Symbol: **TTLL4** (Tubulin Tyrosine Ligase Like 4)
	- Function: An enzyme responsible for polyglutamylation—a chemical modification made directly to the cell's internal microtubule skeleton.
	- Tumor Role: This is our direct connection to the sister chromatid segregation and chromosome segregation GO pathways. Taxane chemotherapies work by freezing microtubules to halt mitosis. By turning TTLL4 ON, resistant tumors alter their microtubules to escape the drug's physical block.

## Clinical Breakdown of the Two Patients
- Patient GSM615096: Chemo-Resistant (Label 1)
  - Biological State: The model classified this patient as 1 (**Resistant**) because their tumor exhibits a clear "gene active" signature.
  - Key Drivers: The exceptionally high expression of 205225_at (13.93) and 219051_x_at (9.64) likely represents up-regulated survival pathways, drug efflux pumps, or rapid cell proliferation mechanisms. This genetic activity shields the tumor from the impending chemotherapy.

- Patient GSM615099: Chemo-Sensitive (Label 0)
  - Biological State: The model classified this patient as 0 (**Sensitive**). Their high-risk tumor genes are effectively "turned off."
  - Key Drivers: Their expression of 205225_at drops to a low 8.19. Without this protective genetic shield, and with higher expression in baseline genes like 203702_s_at (9.89), the tumor remains highly vulnerable. This profile suggests the chemotherapy will likely be highly effective.


## Clinical & Translational Value

By catching these signatures early, an oncologist looking at these results could potentially bypass standard chemotherapy for patient GSM615096 in favor of targeted alternative therapies, while confidently proceeding with standard chemotherapy for patient GSM615099.


## Interpretation of the plot
An ROC Curve (Receiver Operating Characteristic) measures how good the model is at separating the two classes—in our case, telling the difference between a Chemo-Resistant patient (1) and a Chemo-Sensitive patient (0).


1. The Core Metrics (The Axes)

- Vertical Axis (True Positive Rate / Sensitivity): This is our "Catch Rate". It answers: Out of all the truly chemo-resistant patients, what percentage did the model correctly identify? You want this to be as close to 1.0 (100%) as possible.
- Horizontal Axis (False Positive Rate): This is our "False Alarm Rate". It answers: Out of all the chemo-sensitive patients, what percentage did the model mistakenly label as resistant? You want this to be as close to 0.0 (0%) as possible.

2. The Baseline vs. our Model

- The Diagonal Dashed Line: This represents pure random guessing (like flipping a coin). A useless model runs straight down this diagonal line.
- The Red Solid Line (our Model): This shows how our model actually performs. Notice how it shoots straight up dramatically at the beginning? That means the model can catch a huge chunk of resistant patients (up to 56% immediately) with zero false alarms.

3. Understanding the AUC score (0.90)

AUC stands for "Area Under the Curve." It measures the entire two-dimensional area underneath the red line.
- A perfect model has an AUC of 1.0 (the red line would go straight up to the top left corner, then flat across).
- A random guessing model has an AUC of 0.50 (the area under the dashed line).
- our score of 0.90 is excellent. Biologically, it means that if you randomly choose one Chemo-Resistant patient and one Chemo-Sensitive patient from our study, there is a 90% chance the model will correctly point out which one is resistant based purely on their gene data.

### Summary

The closer the red curve bows toward the top-left corner, the better the machine learning model is. Because your curve pushes high up into that corner, your model is highly capable of looking at the 2,938 genes in a new patient and accurately telling the doctor whether their tumor will resist or yield to chemotherapy.


# The Translational Story: Precision Medicine in Action

This is the ultimate clinical goal of our project, known as Precision Medicine. If a doctor can predict treatment failure before starting therapy, they can completely change the patient's clinical outcome.

## Part 1: The Clinical Decision Loop

Let's imagine a patient walks into a clinic in Malaysia, and a doctor takes a small needle biopsy of their breast tumor before any treatment begins. The lab runs a transcriptomic panel on that biopsy tissue.

### Patient 1 (Profiles like GSM615096)

- The ML Prediction: Our Random Forest classifier scans their genes and predicts a Label of 1 (Resistant / RD).
- The Doctor's Decision: Halt standard therapy. If the doctor gives them standard taxane chemotherapy, the tumor will resist it, meaning the patient will suffer severe toxic side effects for months while the tumor potentially continues to grow or spread. Instead, the doctor routes this patient to an alternative treatment strategy (e.g., matching them with a targeted immunotherapy trial, adding a cell-cycle inhibitor, or scheduling immediate surgery instead of waiting).

### Patient 2 (Profiles like GSM615099)

- The ML Prediction: Our Random Forest classifier scans their 2,938 genes and predicts a Label of 0 (Sensitive / pCR).
- The Doctor's Decision: Go ahead with the Standard Therapy (Taxane-based neoadjuvant chemotherapy). The model predicts the patient's tumor lacks structural defense mechanisms, meaning the drugs will likely destroy the tumor completely, leading to a Pathologic Complete Response (pCR).


## Part 2: Connecting the ML Predictions to our GO Pathways

Our machine learning code identified who is resistant, but our Gene Ontology (GO) enrichment results explain why they are resistant.
Let's look at how the specific values for Patient 1 (GSM615096) map directly to the molecular machinery uncovered in Milestone 4:

1. **Overcoming the Drug's Mechanism (The Taxane Connection)**

- The Biology: Standard chemotherapy for these breast cancer patients relies heavily on Taxanes (like Paclitaxel). Taxanes act like structural handcuffs inside a cell — they bind to microtubules to freeze the cell's skeleton, preventing it from dividing during mitosis. If a cell cannot divide, it triggers a self-destruct sequence and dies.
- Our GO Connection: Our top enrichment pathways were "mitotic cell cycle phase transition" and "sister chromatid segregation".
- The Story: In Patient 1 (Resistant), genes like `205347_s_at` (which controls actin/cytoskeleton modeling) and `203702_s_at` (which chemically modifies tubulin/microtubules) are heavily turned ON. Because these specific genes are hyperactive, the resistant tumor has essentially remodeled its cell-division machinery. They are actively pushing past the drug's structural handcuffs, forcing "mitotic nuclear division" and "chromosome segregation" to occur despite the chemotherapy being present.

2. **The Proliferation Shield**

- Because Patient 1 has turned these cell cycle transition pathways completely ON, their tumor cells divide rapidly and ignore the chemical signals that should be forcing them to stop. The tumor builds a transcriptomic "shield."
- Conversely, in Patient 2, these pathways are turned OFF or quiet. Their cell cycle machinery is vulnerable, so when the taxanes enter the cell, they easily break the structural handcuffs into place, halting mitosis completely and clearing the disease.