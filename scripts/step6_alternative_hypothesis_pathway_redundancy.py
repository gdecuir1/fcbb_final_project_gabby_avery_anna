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
- `data/processed/tp53_functional_status.parquet`

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

from pathlib import Path
import argparse


def main() -> None:
    parser = argparse.ArgumentParser(description="Alternative hypothesis: pathway redundancy (TP53 LoF vs WT).")
    parser.add_argument("--pathway-map", default="config/pathways/pathway_to_genes.csv", help="Pathway to genes CSV.")
    parser.add_argument("--tp53-status", default="data/processed/tp53_functional_status.parquet", help="TP53 functional groups.")
    parser.add_argument(
        "--mutation-matrix",
        default="data/processed/gene_mutation_binarized_matrix.parquet",
        help="Gene-level binarized mutations per sample.",
    )
    parser.add_argument("--out-dir", default="outputs", help="Base output directory.")
    args = parser.parse_args()

    # TODO: implement pathway aggregation + Fisher tests for TP53 LoF vs WT.
    raise NotImplementedError("Fill in pathway redundancy analysis logic.")


if __name__ == "__main__":
    main()

