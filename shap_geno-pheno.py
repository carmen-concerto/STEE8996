import pandas as pd
import numpy as np
import matplotlib.pyplot as plt
import shap
from sklearn.ensemble import RandomForestRegressor

# ── Load data ──────────────────────────────────────────────────────────────
X = pd.read_csv("genotype_filtered_matrix.csv", index_col=0)  # samples x genes
Y = pd.read_csv("phenotype_matrix.csv", index_col=0).T        # samples x conditions

gene_names      = X.columns.tolist()
condition_names = Y.columns.tolist()
pos_ctrl        = 17519

print(f"Genotype matrix : {X.shape[0]} samples x {X.shape[1]} genes")
print(f"Phenotype matrix: {Y.shape[0]} samples x {Y.shape[1]} conditions")

# ── Storage ────────────────────────────────────────────────────────────────
shap_mean_abs   = pd.DataFrame(index=gene_names, columns=condition_names, dtype=float)
shap_mean_dir   = pd.DataFrame(index=gene_names, columns=condition_names, dtype=float)  # signed

# ── Per-condition SHAP ─────────────────────────────────────────────────────
for condition in condition_names:
    print(f"\nProcessing: {condition}")
    y = Y[condition].values

    rf = RandomForestRegressor(n_estimators=500, random_state=42)
    rf.fit(X, y)

    explainer   = shap.TreeExplainer(rf)
    shap_values = explainer.shap_values(X)  # shape: n_samples x n_genes

    # Mean absolute SHAP (magnitude of effect)
    shap_mean_abs[condition] = np.abs(shap_values).mean(axis=0)

    # Mean signed SHAP (direction of effect: positive = benefits fitness)
    shap_mean_dir[condition] = shap_values.mean(axis=0)

    # ── Per-condition summary plot ─────────────────────────────────────────
    fig, ax = plt.subplots(figsize=(8, max(4, len(gene_names) * 0.35)))
    shap.summary_plot(
        shap_values, X,
        feature_names=gene_names,
        show=False,
        plot_size=None
    )
    plt.title(f"SHAP Summary — {condition}", fontsize=12)
    plt.tight_layout()
    plt.savefig(f"shap_summary_{condition}.png", dpi=150, bbox_inches="tight")
    plt.close()

    # ── Dependence plot for positive control ──────────────────────────────
    if pos_ctrl in gene_names:
        fig, ax = plt.subplots(figsize=(6, 4))
        shap.dependence_plot(
            pos_ctrl, shap_values, X,
            feature_names=gene_names,
            ax=ax, show=False
        )
        ax.set_title(f"SHAP Dependence ({pos_ctrl}) — {condition}", fontsize=11)
        plt.tight_layout()
        plt.savefig(f"shap_dependence_{pos_ctrl}_{condition}.png", dpi=150, bbox_inches="tight")
        plt.close()

    # ── Save raw SHAP values per condition ────────────────────────────────
    shap_df = pd.DataFrame(shap_values, index=X.index, columns=gene_names)
    shap_df.to_csv(f"shap_raw_{condition}.csv")

# ── Summary matrices ───────────────────────────────────────────────────────
shap_mean_abs.to_csv("shap_mean_abs.csv")
shap_mean_dir.to_csv("shap_mean_dir.csv")
print("\nSaved: shap_mean_abs.csv, shap_mean_dir.csv")

# ── Heatmap: mean absolute SHAP (importance) ──────────────────────────────
import seaborn as sns

fig, ax = plt.subplots(figsize=(12, max(4, len(gene_names) * 0.4)))
sns.heatmap(
    shap_mean_abs.astype(float),
    cmap="YlOrRd",
    annot=True, fmt=".3f",
    linewidths=0.5,
    ax=ax
)
ax.set_title("Mean |SHAP| per Gene per Condition", fontsize=13)
ax.set_xlabel("Condition")
ax.set_ylabel("Gene")
plt.tight_layout()
plt.savefig("shap_heatmap_importance.png", dpi=150, bbox_inches="tight")
plt.show()
print("Saved: shap_heatmap_importance.png")

# ── Heatmap: mean signed SHAP (direction) ─────────────────────────────────
fig, ax = plt.subplots(figsize=(12, max(4, len(gene_names) * 0.4)))
sns.heatmap(
    shap_mean_dir.astype(float),
    cmap="RdBu",       # blue = negative (fitness cost), red = positive (fitness benefit)
    center=0,
    annot=True, fmt=".3f",
    linewidths=0.5,
    ax=ax
)
ax.set_title("Mean Signed SHAP per Gene per Condition\n(red = fitness benefit, blue = fitness cost)", fontsize=13)
ax.set_xlabel("Condition")
ax.set_ylabel("Gene")
plt.tight_layout()
plt.savefig("shap_heatmap_direction.png", dpi=150, bbox_inches="tight")
plt.show()
print("Saved: shap_heatmap_direction.png")

# ── Positive control summary ───────────────────────────────────────────────
if pos_ctrl in gene_names:
    print(f"\nPositive control ({pos_ctrl}) SHAP summary:")
    ctrl_abs = shap_mean_abs.loc[pos_ctrl].astype(float)
    ctrl_dir = shap_mean_dir.loc[pos_ctrl].astype(float)
    summary  = pd.DataFrame({
        "mean_abs_shap" : ctrl_abs,
        "mean_signed_shap" : ctrl_dir,
        "direction" : ctrl_dir.apply(lambda v: "beneficial" if v > 0 else "deleterious")
    })
    print(summary.to_string())
    summary.to_csv(f"shap_summary_{pos_ctrl}.csv")