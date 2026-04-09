"""
step2_fisher_pairwise_heatmap_tp53_binary.py


Phase 2 (core stats MVP): Run Fisher's exact test for all pairs of genes,
stratified by strict binary TP53 status (WT vs mutated).

Inputs
------
- In-memory mutation matrix from `step1_digestion_and_processing.preprocess()`

Outputs
-------
- Heatmaps saved under `outputs/figures/`:
  - `fisher_tp53_binary_pairwise_heatmap_TP53_mut.png`
  - `fisher_tp53_binary_pairwise_heatmap_TP53_WT.png`

Functionality
-------------
- For each TP53 stratum:
  - Build a 2x2 contingency table for each gene pair:
    - Presence/absence of mutation in gene A
    - Presence/absence of mutation in gene B
  - Compute Fisher's exact test p-values.
  - Plot lower-triangle -log10(p) heatmap (same logic as before).
"""
from __future__ import annotations

import argparse
import sys
from itertools import combinations
from pathlib import Path

import matplotlib.pyplot as plt
import numpy as np
import pandas as pd
import seaborn as sns
from scipy.stats import fisher_exact

# Allow `python scripts/step2_...py` from repo root: import step1 from `scripts/`.
_SCRIPT_DIR = Path(__file__).resolve().parent
REPO_ROOT = _SCRIPT_DIR.parent
if str(_SCRIPT_DIR) not in sys.path:
    sys.path.insert(0, str(_SCRIPT_DIR))

from step1_digestion_and_processing import preprocess


def plot_fisher_heatmap(
    df: pd.DataFrame,
    title: str,
    save_path: Path | None = None,
    show: bool = True,
):
    # Fisher is run only on gene columns (0/1 mutation presence); drop cohort label.
    # drop status column (only want binary data)
    df_genes = df.drop(columns=["TP53_status"], errors="ignore")
    genes = df_genes.columns.tolist()
    n_genes = len(genes)

    # Symmetric p-value matrix (default 1.0 = not updated; pairs get filled in the loop).
    # initialize p-val matrix
    p_val_matrix = pd.DataFrame(1.0, index=genes, columns=genes)
    for g1, g2 in combinations(genes, 2):
        # 2x2 table for co-occurrence of mutations in gene1 vs gene2 within this stratum:
        #              gene2=1   gene2=0
        #   gene1=1       a        b
        #   gene1=0       c        d
        # build 2x2 contingency table
        a = ((df_genes[g1] == 1) & (df_genes[g2] == 1)).sum()
        b = ((df_genes[g1] == 1) & (df_genes[g2] == 0)).sum()
        c = ((df_genes[g1] == 0) & (df_genes[g2] == 1)).sum()
        d = ((df_genes[g1] == 0) & (df_genes[g2] == 0)).sum()

        table = [[a, b], [c, d]]

        # calculate fisher's exact test
        _, pval = fisher_exact(table)
        # populate matrix
        p_val_matrix.loc[g1, g2] = pval
        p_val_matrix.loc[g2, g1] = pval

    # Brighter cells = smaller p-values; tiny epsilon avoids log(0).
    # convert to -log10(p-value) for visualization
    log_p_matrix = -np.log10(p_val_matrix.astype(float) + 1e-10)
    # hide upper triangle and diagonal
    mask = np.triu(np.ones_like(log_p_matrix, dtype=bool))

    # Lower triangle only: each unordered gene pair appears once in the plot.
    # plot the heatmap
    plt.figure(figsize=(10, 8))
    sns.heatmap(
        log_p_matrix,
        mask=mask,
        cmap="Reds",
        annot=False,
        linewidths=0.5,
        cbar_kws={"label": "-log10(p-value)"},
    )
    # Mark pairs with p <= 0.05 (uncorrected) in the lower triangle.
    # add asterisks for statistically significant pairs
    for i in range(n_genes):
        for j in range(i + 1, n_genes):
            if p_val_matrix.iloc[j, i] <= 0.05:
                plt.text(
                    i + 0.5,
                    j + 0.5,
                    "*",
                    ha="center",
                    va="center",
                    color="black",
                    fontsize=20,
                )
    plt.title(title)
    plt.tight_layout()
    # Persist figure next to step 4 outputs; optional interactive display.
    if save_path is not None:
        save_path = Path(save_path)
        save_path.parent.mkdir(parents=True, exist_ok=True)
        plt.savefig(save_path, dpi=150, bbox_inches="tight")
    if show:
        plt.show()
    else:
        plt.close()

    return p_val_matrix


def main() -> None:
    parser = argparse.ArgumentParser(description="TP53 binary stratified Fisher heatmaps.")
    parser.add_argument(
        "--out-dir",
        default=str(REPO_ROOT / "outputs" / "figures"),
        help="Directory for saved PNGs.",
    )
    parser.add_argument(
        "--no-show",
        action="store_true",
        help="Do not open interactive plot windows (save PNG only).",
    )
    args = parser.parse_args()

    # Resolve output directory relative to repo root when a relative path is passed.
    out_dir = Path(args.out_dir)
    if not out_dir.is_absolute():
        out_dir = REPO_ROOT / out_dir

    # Same mutation matrices as step 1, split by TP53_mut vs TP53_WT.
    tp53_mut, tp53_wt = preprocess()
    show = not args.no_show

    print("Generating Heatmap for TP53 Mutated Cohort...")
    plot_fisher_heatmap(
        tp53_mut,
        title="Gene Co-occurrence (TP53 Mutated)",
        save_path=out_dir / "fisher_tp53_binary_pairwise_heatmap_TP53_mut.png",
        show=show,
    )

    print("Generating Heatmap for TP53 Wild-Type Cohort...")
    plot_fisher_heatmap(
        tp53_wt,
        title="Gene Co-occurrence (TP53 Wild-Type)",
        save_path=out_dir / "fisher_tp53_binary_pairwise_heatmap_TP53_WT.png",
        show=show,
    )

    print(f"Saved heatmaps under {out_dir}")


if __name__ == "__main__":
    main()
