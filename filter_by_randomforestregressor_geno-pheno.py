import pandas as pd
import numpy as np
import matplotlib.pyplot as plt
import seaborn as sns

# Load importance matrix (genes x conditions)
imp = pd.read_csv("rf_importances.csv", index_col=0)
X = pd.read_csv("genotype_matrix.csv", index_col=0).T  # shape: 23 samples x n_genes

# ── Thresholds ─────────────────────────────────────────────────────────────
# Rule 1: mean + 1SD across ALL values (global baseline)
global_mean = imp.values.mean()
global_sd   = imp.values.std()
threshold_mean = global_mean + global_sd

# Rule 2: max importance in any single condition
threshold_max = imp.values.max() * 0.20  # at least 10% of the best score seen

# Rule 3: must exceed global mean in at least N conditions
min_conditions = 1  # lower to 1 if mutations are condition-specific

# ── Apply filters ──────────────────────────────────────────────────────────
mean_imp      = imp.mean(axis=1)
max_imp       = imp.max(axis=1)
n_above_mean  = (imp.gt(global_mean)).sum(axis=1)

mask = (
    (mean_imp >= threshold_mean) |          # broadly important
    (max_imp  >= threshold_max)             # important in at least one condition
) & (n_above_mean >= min_conditions)        # not a one-hit wonder

filtered_genes = imp[mask]

print(f"Total genes:    {len(imp)}")
print(f"Filtered genes: {len(filtered_genes)}")
print(f"\nThresholds used:")
print(f"  Mean importance threshold : {threshold_mean:.4f}  (global mean+1SD)")
print(f"  Max importance threshold  : {threshold_max:.4f}  (10% of best score)")
print(f"  Min conditions above mean : {min_conditions}")
print(f"\nSelected genes:")
print(filtered_genes.mean(axis=1).sort_values(ascending=False).to_string())

# ── Heatmap of filtered genes ──────────────────────────────────────────────
fig, ax = plt.subplots(figsize=(12, max(4, len(filtered_genes) * 0.4)))
sns.heatmap(
    filtered_genes,
    cmap="YlOrRd",
    annot=True, fmt=".3f",
    linewidths=0.5,
    ax=ax
)
ax.set_title("RF Feature Importances — Filtered Genes", fontsize=13)
ax.set_xlabel("Condition")
ax.set_ylabel("Gene")
plt.tight_layout()
plt.savefig("rf_importances_filtered_heatmap.png", dpi=150)
plt.show()

# ── Save filtered gene list and matrix ────────────────────────────────────
filtered_genes.to_csv("rf_importances_filtered.csv")
filtered_gene_list = filtered_genes.index.tolist()
print(f"\nFinal gene list saved: {filtered_gene_list}")

genotype_filtered = X[filtered_gene_list]
genotype_filtered.to_csv("genotype_filtered_matrix.csv")
print(f"\nFiltered genotype matrix saved: {genotype_filtered.shape[0]} samples x {genotype_filtered.shape[1]} genes")