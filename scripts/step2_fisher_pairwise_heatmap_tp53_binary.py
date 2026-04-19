"""
step2_fisher_pairwise_heatmap_tp53_binary.py

Pipeline — **Step 2** (LUAD; TP53 binary stratification)

**Phase 2 — core co-occurrence screen:** pairwise Fisher exact tests on the curated gene list,
stratified by **binary TP53** (mutated vs wild-type), mirroring the logic used later in step 4
for LoF/GoF strata.

Inputs
------
- In-memory matrices from ``step1_digestion_and_processing.preprocess()`` (TP53 Mut vs WT splits).

Outputs
-------
- ``outputs/figures/fisher_tp53_binary_pairwise_heatmap_TP53_mut.png``
- ``outputs/figures/fisher_tp53_binary_pairwise_heatmap_TP53_WT.png``

Method (per stratum)
--------------------
For each unordered gene pair, build a 2×2 table of mutation presence/absence, run Fisher's exact
test, visualize **−log₁₀(p)** on a heatmap (lower triangle), and mark pairs with uncorrected *p* ≤ 0.05.
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
    # Gene columns only (0/1); drop TP53_status cohort label used for stratification.
    df_genes = df.drop(columns=["TP53_status"], errors="ignore")
    genes = df_genes.columns.tolist()
    n_genes = len(genes)

    # Symmetric p-value and odds-ratio matrices (1.0 default until each pair is filled).
    p_val_matrix = pd.DataFrame(1.0, index=genes, columns=genes)
    or_matrix = pd.DataFrame(1.0, index=genes, columns=genes)
    for g1, g2 in combinations(genes, 2):
        # 2x2 table for co-occurrence of mutations in gene1 vs gene2 within this stratum:
        #              gene2=1   gene2=0
        #   gene1=1       a        b
        #   gene1=0       c        d
        a = ((df_genes[g1] == 1) & (df_genes[g2] == 1)).sum()
        b = ((df_genes[g1] == 1) & (df_genes[g2] == 0)).sum()
        c = ((df_genes[g1] == 0) & (df_genes[g2] == 1)).sum()
        d = ((df_genes[g1] == 0) & (df_genes[g2] == 0)).sum()

        table = [[a, b], [c, d]]

        odds_ratio, pval = fisher_exact(table)
        p_val_matrix.loc[g1, g2] = pval
        p_val_matrix.loc[g2, g1] = pval
        or_matrix.loc[g1, g2] = odds_ratio
        or_matrix.loc[g2, g1] = odds_ratio

    # −log₁₀(p); small epsilon avoids log(0). Mask upper triangle (each pair once).
    log_p_matrix = -np.log10(p_val_matrix.astype(float) + 1e-10)
    mask = np.triu(np.ones_like(log_p_matrix, dtype=bool))

    plt.figure(figsize=(10, 8))
    sns.heatmap(
        log_p_matrix,
        mask=mask,
        cmap="Reds",
        annot=False,
        linewidths=0.5,
        cbar_kws={"label": "-log10(p-value)"},
    )
    # Uncorrected p ≤ 0.05: print summary and mark cell with '*'.
    for i in range(n_genes):
        for j in range(i + 1, n_genes):
            if p_val_matrix.iloc[j, i] <= 0.05:
                g1 = genes[i]
                g2 = genes[j]
                pval = p_val_matrix.iloc[j, i]
                odds = or_matrix.iloc[j, i]
                if odds > 1:
                    relationship = "Co-occurring"
                elif odds < 1:
                    relationship = "Mutually exclusive"
                else:
                    relationship = "Independent"
                print(f"{g1} & {g2}: p={pval:.4f}, OR={odds:.2f} ({relationship})")
                
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
    parser = argparse.ArgumentParser(
        description="Step 2: TP53 binary (Mut vs WT) stratified pairwise Fisher heatmaps (LUAD).",
    )
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
