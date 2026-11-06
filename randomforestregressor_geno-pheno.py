import pandas as pd
import numpy as np
from sklearn.ensemble import RandomForestRegressor
from sklearn.model_selection import LeaveOneOut
from sklearn.metrics import r2_score

# Load data
# X: samples x genes (allele frequencies)
# Y: samples x conditions (log2 fitness)
X = pd.read_csv("genotype_matrix.csv", index_col=0).T  # shape: 23 samples x n_genes
Y = pd.read_csv("phenotype_matrix.csv", index_col=0).T  # shape: 23 samples x 14 conditions

gene_names = X.columns
condition_names = Y.columns

importance_df = pd.DataFrame(index=gene_names, columns=condition_names)

for condition in condition_names:
    y = Y[condition].values
    
    # Leave-one-out CV (appropriate for N=23)
    loo = LeaveOneOut()
    rf = RandomForestRegressor(n_estimators=500, random_state=42)
    
    preds = []
    for train_idx, test_idx in loo.split(X):
        rf.fit(X.iloc[train_idx], y[train_idx])
        preds.append(rf.predict(X.iloc[test_idx])[0])
    
    print(f"{condition} LOO R²: {r2_score(y, preds):.3f}")
    
    # Fit on full data for importances
    rf.fit(X, y)
    importance_df[condition] = rf.feature_importances_

importance_df.to_csv("rf_importances.csv")
print(importance_df)

# Check where veA ranks in each condition
pos_ctrl = 17519

for condition in condition_names:
    col = importance_df[condition].sort_values(ascending=False)
    rank = col.index.get_loc(pos_ctrl) + 1
    score = col[pos_ctrl]
    print(f"{condition}: {pos_ctrl} ranked #{rank} with score {score:.3f}")