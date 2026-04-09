"""
02_preprocess_mvp_primary_coding_binary_tp53.py

Phase 2 (setup): Preprocess MVP MAF and binarize mutation calls.

Inputs
------
- `config/cohort_mvp.yaml`
- Raw mutation files from:
  - `data/raw/<cohort>_mutations.maf(.gz)?`
  - `data/raw/<cohort>_sample_metadata.(tsv/csv)(.gz)?` (if you downloaded it)

Outputs (fill in exact filenames after inspecting your MAF schema)
---------------------------------------------------------------
1. `data/processed/<cohort>_primary_coding_mutations.parquet`
   - Filtered to primary solid tumors and coding (protein-changing) variants only.
2. `data/processed/tp53_binary_status.parquet`
   - One row per sample:
     - `sample_id`
     - `tp53_binary` (1 = mutated, 0 = WT; per `tp53_binary_definition`)
3. `data/processed/gene_mutation_binarized_matrix.parquet` (or long/tidy format)
   - For downstream pairwise Fisher tests:
     - Either wide matrix: rows=samples, cols=genes, values ∈ {0,1}
     - Or long table: (`sample_id`, `gene`, `mutated_flag`)

Functionality (to implement)
-------------------------------
- Filter for primary solid tumors only.
- Filter for coding regions only.
- Binarize mutation data:
  - For each gene X and sample:
    - mutated = 1
    - WT = 0
- Build strict TP53 binary labels for Phase 2 MVP:
  - TP53 WT vs any non-silent coding mutation in TP53
  - Ignore non-coding/promoter regions for now (as requested).

Notes
-----
- The hardest part here is mapping your MAF columns to:
  - sample identifier key used across all tables
  - variant class/type to decide "coding vs non-coding"
  - TP53 variant classification for the binary definition.
"""

from __future__ import annotations

from pathlib import Path
import argparse


def main() -> None:
    parser = argparse.ArgumentParser(description="Preprocess MVP MAF and create TP53 binary labels.")
    parser.add_argument("--config", default="config/cohort_mvp.yaml", help="Path to cohort config YAML.")
    parser.add_argument("--cohort", default="TCGA-LUAD", help="Cohort study id (e.g. TCGA-LUAD).")
    args = parser.parse_args()

    config_path = Path(args.config)
    if not config_path.exists():
        raise FileNotFoundError(f"Config not found: {config_path}")

    # TODO: Implement reading raw MAF/metadata, filtering, and producing binarized outputs.
    raise NotImplementedError("Fill in preprocessing logic (filtering + coding-only + binarization).")


if __name__ == "__main__":
    main()

