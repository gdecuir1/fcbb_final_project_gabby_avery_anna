"""
step2_fisher_pairwise_heatmap_tp53_binary.py


Phase 2 (core stats MVP): Run Fisher's exact test for all pairs of genes,
stratified by strict binary TP53 status (WT vs mutated).

Inputs
------
- `config/gene_lists/mvp_driver_genes.csv` (genes to test; ~30-50)
- `data/processed/gene_mutation_binarized_matrix.parquet` (or tidy equivalent)
- `data/processed/tp53_binary_status.parquet`

Outputs
--------
1. Pairwise p-value tables:
   - `outputs/tables/fisher_tp53_binary_pairwise_pvalues.csv`
   - (optionally) `outputs/tables/fisher_tp53_binary_pairwise_effects.csv`
2. Heatmap(s):
   - `outputs/figures/fisher_tp53_binary_pvalues_heatmap.png`

Functionality (to implement)
-------------------------------
- For each TP53 stratum:
  - Build a 2x2 contingency table for each gene pair:
    - Presence/absence of mutation in gene A
    - Presence/absence of mutation in gene B
  - Compute Fisher's exact test p-values.
- Combine/plot results:
  - At minimum: p-value heatmap for each stratum, or a single combined result.

Notes
-----
- Decide whether to:
  - Use sample-level binary mutation presence (recommended for co-occurrence).
  - Exclude gene pairs with zero counts in any contingency cell.
- Apply multiple-testing correction if required by your class rubric (optional step).
"""
import pandas as pd
import numpy as np
import seaborn as sns
import matplotlib.pyplot as plt
from scipy.stats import fisher_exact
from itertools import combinations

from step1_digestion_and_processing import preprocess

def plot_fisher_heatmap(df, title):
    # drop status column (only want binary data)
    df_genes = df.drop(columns=["TP53_status"], errors='ignore')
    genes = df_genes.columns.tolist()
    n_genes = len(genes)

    # initialize p-val matrix
    p_val_matrix = pd.DataFrame(1.0, index=genes, columns=genes)
    for g1, g2 in combinations(genes, 2):
        # build 2x2 contingency table
        a = ((df_genes[g1] == 1) & (df_genes[g2] == 1)).sum()
        b = ((df_genes[g1] == 1) & (df_genes[g2] == 0)).sum()
        c = ((df_genes[g1] == 0) & (df_genes[g2] == 1)).sum()
        d = ((df_genes[g1] == 0) & (df_genes[g2] == 0)).sum() 

        table = [[a,b],[c,d]]

        # calculate fisher's exact test
        _, pval = fisher_exact(table)
        # populate matrix
        p_val_matrix.loc[g1, g2] = pval 
        p_val_matrix.loc[g2, g1] = pval

    # convert to -log10(p-value) for visualization
    log_p_matrix = -np.log10(p_val_matrix.astype(float) + 1e-10)
    # hide upper triangle and diagonal
    mask = np.triu(np.ones_like(log_p_matrix, dtype=bool))

    # plot the heatmap
    plt.figure(figsize=(10,8))
    ax = sns.heatmap(
        log_p_matrix,
        mask=mask,
        cmap="Reds",
        annot=False,
        linewidths=0.5,
        cbar_kws={'label': '-log10(p-value)'}
    )
    # add asterisks for statistically significant pairs
    for i in range(n_genes):
        for j in range(i+1, n_genes):
            if p_val_matrix.iloc[j, i] <= 0.05:
                plt.text(i + 0.5, j + 0.5, '*', ha='center', va='center', color='black', fontsize=20)
    plt.title(title)
    plt.tight_layout()
    plt.show()

    return p_val_matrix

if __name__ == "__main__":
    tp53_mut, tp53_wt = preprocess()
    print("Generating Heatmap for TP53 Mutated Cohort...")
    p_vals_mut = plot_fisher_heatmap(tp53_mut, title="Gene Co-occurrence (TP53 Mutated)")

    print("Generating Heatmap for TP53 Wild-Type Cohort...")
    p_vals_wt = plot_fisher_heatmap(tp53_wt, title="Gene Co-occurrence (TP53 Wild-Type)")