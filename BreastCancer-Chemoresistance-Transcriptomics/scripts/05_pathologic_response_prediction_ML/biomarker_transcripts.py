from pathlib import Path
import pandas as pd
import numpy as np
import matplotlib.pyplot as plt
import seaborn as sns
from sklearn.model_selection import train_test_split
from sklearn.ensemble import RandomForestClassifier
from sklearn.metrics import roc_auc_score, roc_curve, classification_report, confusion_matrix

print("--- Initializing Genomic Machine Learning Pipeline ---")

# 1. Dynamically find the absolute path of the directory containing this script
SCRIPT_DIR = Path(__file__).resolve().parent

# Build absolute paths to the data folder
features_path = SCRIPT_DIR / "05_GSE25066_ML_Features.csv"
labels_path = SCRIPT_DIR / "05_GSE25066_ML_Labels.csv"

# Load using the fully resolved absolute paths
X = pd.read_csv(features_path, index_col=0)
y_df = pd.read_csv(labels_path)

# Ensure samples line up perfectly
y = y_df.set_index("Sample_ID").loc[X.index]["Label"].values

print(f"Dataset Shape: {X.shape[0]} patient samples x {X.shape[1]} biological gene features")
print(f"Class Balance: {np.sum(y==1)} Resistant (RD) vs {np.sum(y==0)} Sensitive (pCR)")

# 2. Perform a Stratified Train/Test Split (80% Train, 20% Evaluation Test)
# Stratification ensures both slices maintain identical balances of resistant cases
X_train, X_test, y_train, y_test = train_test_split(
    X, y, test_size=0.20, random_state=42, stratify=y
)

# 3. Instantiate and Train the Random Forest Classifier
rf_model = RandomForestClassifier(
    n_estimators=500,        # Number of decision trees in the forest
    max_depth=6,             # Limit depth to control overfitting on low sample size
    random_state=42,
    class_weight="balanced", # Handles the dataset skew (249 RD vs 57 pCR)
    n_jobs=-1                # Use all available CPU cores for speed
)
rf_model.fit(X_train, y_train)
print("[Success] Random Forest model training complete.")

# 4. Generate Performance Evaluations on the Test Set
y_pred = rf_model.predict(X_test)
y_probs = rf_model.predict_proba(X_test)[:, 1] # Probability estimates for being resistant

test_auc = roc_auc_score(y_test, y_probs)
print(f"\nHoldout Test Set ROC-AUC Score: {test_auc:.4f}")
print("\n--- Detailed Classification Metrics ---")
print(classification_report(y_test, y_pred, target_names=["Sensitive (pCR)", "Resistant (RD)"]))

# 5. Dynamically Extract, Annotate & Export All Biomarkers by Gini Importance
importances = rf_model.feature_importances_
importance_df = pd.DataFrame({
    "Probe_ID": X.columns,
    "Gini_Importance": importances
}).sort_values(by="Gini_Importance", ascending=False)

# Merge official HUGO Gene Symbols if map exists
map_path = SCRIPT_DIR / "GSE25066_Probe_to_Gene_Map.csv"
if map_path.exists():
    probe_map = pd.read_csv(map_path)
    importance_df = importance_df.merge(probe_map, on="Probe_ID", how="left")
    importance_df["Label"] = importance_df["Gene_Symbol"] + " (" + importance_df["Probe_ID"] + ")"
else:
    importance_df["Gene_Symbol"] = importance_df["Probe_ID"]
    importance_df["Label"] = importance_df["Probe_ID"]

# Save complete dynamically ranked biomarker table to CSV
# importance_df.to_csv(SCRIPT_DIR / "GSE25066_RF_Gini_Importance_All.csv", index=False)
top_10_biomarkers = importance_df.head(10)

print("\n=== Top 10 Diagnostic Biomarker Transcripts Dynamically Extracted ===")
print(top_10_biomarkers[["Probe_ID", "Gene_Symbol", "Gini_Importance"]].to_string(index=False))

# 6. Plot Combined 2-Panel Publication Figure (Panel A: ROC Curve | Panel B: Top 10 Gini Biomarkers)
fpr, tpr, _ = roc_curve(y_test, y_probs)
fig, axes = plt.subplots(1, 2, figsize=(14, 5.5))

# Panel A: ROC-AUC Validation Curve
axes[0].plot(fpr, tpr, color="#E41A1C", lw=2.5, label=f"Random Forest (AUC = {test_auc:.2f})")
axes[0].plot([0, 1], [0, 1], color="#999999", linestyle="--", lw=1.2)
axes[0].set_xlim([0.0, 1.0])
axes[0].set_ylim([0.0, 1.05])
axes[0].set_xlabel("False Positive Rate (1 - Specificity)", fontsize=11)
axes[0].set_ylabel("True Positive Rate (Sensitivity / Recall)", fontsize=11)
axes[0].set_title("A. Predictive Classification Performance (GSE25066)", fontsize=12, fontweight="bold")
axes[0].legend(loc="lower right", fontsize=10)
axes[0].grid(alpha=0.3)

# Panel B: Dynamically Ranked Top 10 Biomarkers (Gini Importance)
top10_plot = top_10_biomarkers.iloc[::-1]
axes[1].barh(top10_plot["Label"], top10_plot["Gini_Importance"], color="#B2182B", height=0.68)
axes[1].set_xlabel("Random Forest Gini Importance Score", fontsize=11)
axes[1].set_title("B. Top 10 Predictive Chemoresistance Biomarkers", fontsize=12, fontweight="bold")
axes[1].grid(axis="x", alpha=0.3)

plt.tight_layout()
plt.savefig(SCRIPT_DIR / "Milestone5_Model_ROC_and_Biomarkers.png", dpi=300)
print("\n[SUCCESS] Combined ROC + Top 10 Biomarkers plot saved as 'Milestone5_Model_ROC_and_Biomarkers.png'.")
