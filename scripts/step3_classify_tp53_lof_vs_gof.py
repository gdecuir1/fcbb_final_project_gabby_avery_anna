"""
04_classify_tp53_lof_vs_gof.py

Phase 3: Refine TP53 definition into three functional categories.

Inputs
------
- `config/tp53_classification.yaml`
- Processed mutation data from Phase 2 (must retain TP53 variant-level info needed
  to decide LoF vs GoF/missense categories):
  - ideally: an intermediate TP53-annotated table derived from your MAF
  - e.g. `data/processed/tp53_variants.parquet` (you will create in Phase 2)

Outputs
--------
`data/processed/tp53_functional_status.parquet`
One row per sample:
- `sample_id`
- `tp53_group` (labels per `tp53_classification.yaml`)
  - `TP53_WT`
  - `TP53_LoF`
  - `TP53_GoF_missense`

Functionality (to implement)
-------------------------------
- WT: no non-silent coding TP53 mutation (per Phase 2 definition).
- LoF: nonsense / frameshift / splice-site (per your chosen variant mapping).
- GoF / Missense:
  - Start with missense as a proxy for GoF/DN (as requested).
  - Later: optionally refine classification using curated annotations or predictive algorithms
    (SIFT/PolyPhen/REVEL) if you have them in the dataset.

Notes
-----
- Decide a rule for samples with multiple TP53 alterations:
  - Example placeholder rule: prioritize LoF over missense, or mark as "mixed" then decide later.
- Ensure you are classifying at the sample level (not variant level) for downstream co-occurrence tests.
"""

from __future__ import annotations

from pathlib import Path
import argparse


def main() -> None:
    parser = argparse.ArgumentParser(description="Classify TP53 into WT / LoF / GoF+missense groups.")
    parser.add_argument("--config", default="config/tp53_classification.yaml", help="Path to TP53 classification YAML.")
    parser.add_argument(
        "--tp53-variants",
        default="data/processed/tp53_variants.parquet",
        help="TP53 variant-level table produced during preprocessing.",
    )
    parser.add_argument("--out", default="data/processed/tp53_functional_status.parquet", help="Output path.")
    args = parser.parse_args()

    # TODO: implement LoF vs GoF/missense classification logic.
    raise NotImplementedError("Fill in TP53 functional grouping logic.")


if __name__ == "__main__":
    main()

