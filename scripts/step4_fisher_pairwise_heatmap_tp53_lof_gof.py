"""
05_fisher_pairwise_heatmap_tp53_lof_gof.py

Phase 3 (core stats refinement): Re-run co-occurrence Fisher tests using
three TP53 groups:
1) TP53_WT
2) TP53_LoF
3) TP53_GoF_or_missense

Inputs
------
- `config/gene_lists/mvp_driver_genes.csv`
- `data/processed/gene_mutation_binarized_matrix.parquet`
- `data/processed/lof_gof/tp53_functional_status.csv` (from step 3)

Outputs
--------
1. Pairwise p-value tables:
   - `outputs/tables/fisher_tp53_lof_gof_pairwise_pvalues.csv`
2. Heatmap(s):
   - `outputs/figures/fisher_tp53_lof_gof_pvalues_heatmap_TP53_WT.png`
   - `outputs/figures/fisher_tp53_lof_gof_pvalues_heatmap_TP53_LoF.png`
   - `outputs/figures/fisher_tp53_lof_gof_pvalues_heatmap_TP53_GoF_missense.png`
   (or a combined multi-panel figure)

Functionality (to implement)
-------------------------------
- For each TP53 group independently:
  - compute gene-pair Fisher exact p-values using sample-level binary mutation presence.
- Optionally compare distributions between groups (e.g., delta p-values or effect sizes).

Notes
-----
- If group sample sizes are small, decide whether to:
  - restrict the tested gene set further,
  - use an alternative method,
  - or apply robust multiple-testing correction.
"""

from __future__ import annotations

from pathlib import Path
import argparse

_REPO = Path(__file__).resolve().parent.parent
_DEFAULT_TP53_FUNCTIONAL = str(
    _REPO / "data" / "processed" / "lof_gof" / "tp53_functional_status.csv"
)


def main() -> None:
    parser = argparse.ArgumentParser(description="Run Fisher tests + heatmaps (TP53 LoF vs GoF).")
    parser.add_argument("--genes", default="config/gene_lists/mvp_driver_genes.csv", help="Driver genes CSV.")
    parser.add_argument(
        "--tp53-status",
        default=_DEFAULT_TP53_FUNCTIONAL,
        help="TP53 functional group labels per sample (step 3 CSV under processed/lof_gof/).",
    )
    parser.add_argument(
        "--mutation-matrix",
        default="data/processed/gene_mutation_binarized_matrix.parquet",
        help="Binarized mutation calls table.",
    )
    parser.add_argument("--out-dir", default="outputs", help="Base output directory.")
    args = parser.parse_args()

    # TODO: load inputs, compute gene-pair Fisher tests per TP53 group, write tables and plots.
    raise NotImplementedError("Fill in Fisher test + heatmap generation (TP53 LoF/GoF).")


if __name__ == "__main__":
    main()

