"""
07_alternative_hypothesis_pathway_redundancy.py

Phase 4 backup plan (Alternative Hypothesis #1):
"Pathway redundancy" co-occurrence.

Hypothesis
----------
TP53 LoF tumors will show significantly higher rates of co-occurring mutations
across parallel signaling pathways (e.g., RTK/RAS/PI3K + Wnt) compared to TP53 WT
tumors, which may show more mutual exclusivity between pathways.

Inputs
------
- `config/pathways/pathway_to_genes.csv`
- `data/processed/gene_mutation_binarized_matrix.parquet`
- `data/processed/lof_gof/tp53_functional_status.csv` (from step 3)

Outputs
--------
- Pathway-level co-occurrence results (tables):
  - `outputs/tables/pathway_redundancy_tp53_lof_vs_wt.tsv`
- Optional visualization(s):
  - `outputs/figures/pathway_redundancy_heatmap.png`

Functionality (to implement)
-------------------------------
- Aggregate gene-level mutation binaries into pathway-level binaries per sample:
  - pathway_mutated = 1 if ANY gene in pathway is mutated in the sample
  - pathway_mutated = 0 otherwise
- Compare pathway pairs (parallel pathways):
  - For each pathway pair (A, B), compute enrichment/co-occurrence differences
    between `TP53_LoF` vs `TP53_WT` groups.
- Run Fisher's exact test (or chi-square) per pathway pair.

Notes
-----
- Decide whether to allow overlapping genes across pathways.
- Decide how many pathway pairs to test (all vs curated parallel pairs).
"""

from __future__ import annotations
import argparse
from itertools import combinations
from pathlib import Path

import matplotlib.pyplot as plt
import numpy as np
import pandas as pd
import seaborn as sns
from scipy.stats import fisher_exact

_REPO = Path(__file__).resolve().parent.parent
_DEFAULT_TP53_FUNCTIONAL = str(
    _REPO / "data" / "processed" / "lof_gof" / "tp53_functional_status.csv"
)

def benjamini_hochberg_fdr(p_values: np.ndarray) -> np.ndarray:
    p = np.asarray(p_values, dtype=float)
    n = len(p)
    if n == 0:
        return p
    order = np.argsort(p)
    sorted_p = p[order]
    ranks = np.arange(1, n + 1, dtype=float)
    q_sorted = sorted_p * n / ranks
    for i in range(n - 2, -1, -1):
        q_sorted[i] = min(q_sorted[i], q_sorted[i + 1])
    q = np.empty(n)
    q[order] = np.clip(q_sorted, 0.0, 1.0)
    return q

def aggregate_to_pathways(mut_matrix: pd.DataFrame, pathway_map: pd.DataFrame) -> pd.DataFrame:
    pathway_dict = pathway_map.groupby('pathway')['gene'].apply(list).to_dict()
    pathway_matrix = pd.DataFrame(index=mut_matrix.index)

    for pathway, genes in pathway_dict.items():
        valid_genes = [g for g in genes if g in mut_matrix.columns]
        if valid_genes:
            # mutated if any valid gene in pathway is mutated
            pathway_matrix[pathway] = mut_matrix[valid_genes].max(axis=1)
    return pathway_matrix.astype(np.int8)

def compute_pairwise_fisher(matrix: pd.DataFrame) -> pd.DataFrame:
    pathways = matrix.columns.tolist()
    rows = []

    for p1, p2 in combinations(pathways, 2):
        a = int(((matrix[p1] == 1) & (matrix[p2] == 1)).sum())
        b = int(((matrix[p1] == 1) & (matrix[p2] == 0)).sum())
        c = int(((matrix[p1] == 0) & (matrix[p2] == 1)).sum())
        d = int(((matrix[p1] == 0) & (matrix[p2] == 0)).sum())
        table = [[a, b], [c, d]]
        oddsratio, pval = fisher_exact(table)

        rows.append({
            "pathway_a": p1,
            "pathway_b": p2,
            "n_both": a,
            "p1_only": b,
            "p2_only": c,
            "neither": d,
            "odds_ratio": oddsratio,
            "p_value": pval
        })
    out = pd.DataFrame(rows)
    if len(out) > 0:
        out["q_value_bh"] = benjamini_hochberg_fdr(out["p_value"].values)
    return out

def plot_comparative_heatmap(lof_df: pd.DataFrame, wt_df: pd.DataFrame, out_path: Path):
    out_path.parent.mkdir(parents=True, exist_ok=True)
    pathways = list(set(lof_df['pathway_a']).union(set(lof_df['pathway_b'])))
    pathways.sort()

    def build_matrix(df):
        mat = pd.DataFrame(np.nan, index=pathways, columns=pathways)
        for _, r in df.iterrows():
            p1, p2 = r["pathway_a"], r["pathway_b"]
            sign = 1 if r["odds_ratio"] > 1 else -1
            val = sign * -np.log10(max(float(r["p_value"]), 1e-300))
            mat.loc[p1, p2] = val
            mat.loc[p2, p1] = val
        return mat
        
    lof_mat = build_matrix(lof_df)
    wt_mat = build_matrix(wt_df)

    fig, axes = plt.subplots(1, 2, figsize=(20, 8), sharey=True)
    mask = np.triu(np.ones_like(lof_mat.values, dtype=bool))

    # Calculate vmax, but enforce a minimum of 1.0 so the colorbar doesn't collapse to 0
    vmax_calc = max(np.nanmax(np.abs(lof_mat.values)), np.nanmax(np.abs(wt_mat.values)))
    vmax = max(vmax_calc, 1.0) 

    sns.heatmap(lof_mat, mask=mask, cmap="RdBu_r", center=0, vmin=-vmax, vmax=vmax,
                ax=axes[0], cbar_kws={"label": "Signed -log10(p)"})
    axes[0].set_title("TP53 LoF: Pathway Interactions")
    
    sns.heatmap(wt_mat, mask=mask, cmap="RdBu_r", center=0, vmin=-vmax, vmax=vmax,
                ax=axes[1], cbar_kws={"label": "Signed -log10(p)"})
    axes[1].set_title("TP53 WT: Pathway Interactions") # FIXED: This was mistakenly axes[0]

    plt.tight_layout()
    plt.savefig(out_path, dpi=150, bbox_inches="tight")
    plt.close()

def main() -> None:
    parser = argparse.ArgumentParser(description="Alternative hypothesis: pathway redundancy (TP53 LoF vs WT).")
    parser.add_argument("--pathway-map", default="config/pathways/pathway_to_genes.csv", help="Pathway to genes CSV.")
    parser.add_argument(
        "--tp53-status",
        default=_DEFAULT_TP53_FUNCTIONAL,
        help="TP53 functional groups (step 3 CSV under processed/lof_gof/).",
    )
    parser.add_argument(
        "--mutation-matrix",
        default="data/processed/gene_mutation_binarized_matrix.parquet",
        help="Gene-level binarized mutations per sample.",
    )
    parser.add_argument("--out-dir", default="outputs", help="Base output directory.")
    args = parser.parse_args()

    out_tables = Path(args.out_dir) / "tables"
    out_figures = Path(args.out_dir) / "figures"
    out_tables.mkdir(parents=True, exist_ok=True)
    out_figures.mkdir(parents=True, exist_ok=True)

    print("Loading data...")
    mut_matrix = pd.read_parquet(args.mutation_matrix)
    pathway_map = pd.read_csv(args.pathway_map)
    tp53_status = pd.read_csv(args.tp53_status)

    sample_col = "Tumor_Sample_Barcode" if "Tumor_Sample_Barcode" in tp53_status.columns else tp53_status.columns[0]
    tp53_status = tp53_status.set_index(sample_col)

    print("Aggregating genes to pathway level...")
    pathway_matrix = aggregate_to_pathways(mut_matrix, pathway_map)
    print(f"Collapsed {mut_matrix.shape[1]} genes into {pathway_matrix.shape[1]} pathways.")

    print("Splitting cohorts by TP53 status...")
    lof_samples = tp53_status[tp53_status['tp53_group'] == 'TP53_LoF'].index.intersection(pathway_matrix.index)
    wt_samples = tp53_status[tp53_status['tp53_group'] == 'TP53_WT'].index.intersection(pathway_matrix.index)
    
    pathway_matrix_lof = pathway_matrix.loc[lof_samples]
    pathway_matrix_wt = pathway_matrix.loc[wt_samples]
    
    print(f"TP53 LoF samples: {len(pathway_matrix_lof)}")
    print(f"TP53 WT samples: {len(pathway_matrix_wt)}")

    print("Running pairwise Fisher exact tests...")
    lof_results = compute_pairwise_fisher(pathway_matrix_lof)
    wt_results = compute_pairwise_fisher(pathway_matrix_wt)
    
    # Merge for easy comparison
    merged = pd.merge(
        lof_results, wt_results, 
        on=["pathway_a", "pathway_b"], 
        suffixes=("_lof", "_wt")
    )
    
    # Calculate the difference in odds ratios
    merged["odds_ratio_diff"] = merged["odds_ratio_lof"] - merged["odds_ratio_wt"]
    merged = merged.sort_values("odds_ratio_diff", ascending=False)

    out_file = out_tables / "pathway_redundancy_tp53_lof_vs_wt.tsv"
    merged.to_csv(out_file, sep="\t", index=False)
    print(f"Saved tabular results to {out_file}")

    print("Generating comparative heatmap...")
    heatmap_file = out_figures / "pathway_redundancy_heatmap.png"
    plot_comparative_heatmap(lof_results, wt_results, heatmap_file)
    print(f"Saved figure to {heatmap_file}")

if __name__ == "__main__":
    main()

