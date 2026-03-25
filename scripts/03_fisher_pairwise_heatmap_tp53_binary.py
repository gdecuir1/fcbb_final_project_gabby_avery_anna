"""
03_fisher_pairwise_heatmap_tp53_binary.py

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

from __future__ import annotations

from pathlib import Path
import argparse


def main() -> None:
    parser = argparse.ArgumentParser(description="Run Fisher tests and plot heatmap (TP53 binary).")
    parser.add_argument("--genes", default="config/gene_lists/mvp_driver_genes.csv", help="Path to driver genes CSV.")
    parser.add_argument("--tp53-status", default="data/processed/tp53_binary_status.parquet", help="TP53 binary labels file.")
    parser.add_argument(
        "--mutation-matrix",
        default="data/processed/gene_mutation_binarized_matrix.parquet",
        help="Binarized mutation calls table.",
    )
    parser.add_argument("--out-dir", default="outputs", help="Base output directory.")
    args = parser.parse_args()

    # TODO: load inputs, compute all gene-pair Fisher tests, write outputs, generate heatmaps.
    raise NotImplementedError("Fill in Fisher test + heatmap generation (TP53 binary).")


if __name__ == "__main__":
    main()

