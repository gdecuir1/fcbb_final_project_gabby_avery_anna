"""
08_alternative_hypothesis_structural_vs_point.py

Phase 4 backup plan (Alternative Hypothesis #2):
"Structural vs point" mechanism shift.

Hypothesis
----------
TP53-mutant tumors will favor copy-number alterations (CNA amplifications in PI3K/RTK)
due to genomic instability, whereas TP53-WT tumors will favor point mutations.

Inputs
------
- `data/processed/tp53_functional_status.parquet`
- `data/processed/gene_mutation_binarized_matrix.parquet` (point mutations)
- A CNA-derived table (you will create / preprocess):
  - e.g., `data/processed/cna_amplification_binarized_matrix.parquet`
  - sample-level indicators for gene/pathway amplifications
- Optional mapping/config for which CNA genes count as "PI3K/RTK amplifications":
  - `config/gene_lists/rtk_ras_pi3k_cna_genes.csv` (you can add later)

Outputs
--------
- Summary tables comparing TP53 LoF vs WT:
  - `outputs/tables/structural_vs_point_tp53_lof_vs_wt.tsv`
- Optional stats:
  - per-gene odds ratios and p-values

Functionality (to implement)
-------------------------------
- For each sample, generate:
  - point_mutated = 1 if a point mutation exists in RTK/PI3K genes
  - cna_amplified = 1 if a CNA amplification exists in RTK/PI3K genes
- Compare rates between:
  - TP53_LoF vs TP53_WT (and optionally GoF vs WT)
- Statistical tests:
  - Fisher's exact test per gene (or aggregated pathway level)

Notes
-----
- CNA datasets may come in different formats depending on MAF/CNA source.
- Decide on amplification thresholds later (e.g., log2 copy ratio > X).
"""

from __future__ import annotations

from pathlib import Path
import argparse


def main() -> None:
    parser = argparse.ArgumentParser(description="Alternative hypothesis: structural vs point (CNA vs SNV).")
    parser.add_argument("--tp53-status", default="data/processed/tp53_functional_status.parquet", help="TP53 functional groups.")
    parser.add_argument("--mutation-matrix", default="data/processed/gene_mutation_binarized_matrix.parquet", help="Point mutations matrix.")
    parser.add_argument(
        "--cna-matrix",
        default="data/processed/cna_amplification_binarized_matrix.parquet",
        help="CNA amplification matrix (to be created).",
    )
    parser.add_argument("--out-dir", default="outputs", help="Base output directory.")
    args = parser.parse_args()

    # TODO: implement CNA vs SNV comparison by TP53 group.
    raise NotImplementedError("Fill in structural vs point analysis logic.")


if __name__ == "__main__":
    main()

