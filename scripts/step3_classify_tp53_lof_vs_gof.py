"""
step3_classify_tp53_lof_vs_gof.py

Phase 3: Classify TP53 into three sample-level groups for stratified co-occurrence analysis.

Uses the same MAF filters and LUAD cohort as step1 (coding variants only in the curated gene list).

Groups
------
- TP53_WT: no coding TP53 variant in the merged mutation table (TP53 column == 0).
- TP53_LoF: at least one TP53 variant classified as loss-of-function:
  Nonsense_Mutation, Frame_Shift_Ins, Frame_Shift_Del, Splice_Site.
- TP53_GoF_missense: has TP53 coding mutation but no LoF class above; missense is treated as
  the GoF / dominant-negative proxy (refine later with SIFT/PolyPhen/REVEL or TP53 databases).

Multi-mutation rule
--------------------
If a sample carries both LoF and missense TP53 hits, it is labeled TP53_LoF (LoF prioritized).

Inputs
------
- Same paths as step1 (`build_mutation_matrix`), or override via CLI.

Outputs (default)
-----------------
All step-3 artifacts go under `data/processed/lof_gof/` so step 1–2 behavior and any other
`data/processed/` files are untouched.

- data/processed/lof_gof/tp53_functional_status.csv — sample_id, tp53_group
- data/processed/lof_gof/mutation_matrix_with_tp53_group.csv — full binary matrix + TP53_status + tp53_group

Downstream
----------
Import `classify_tp53_functional_groups` and `attach_tp53_group_to_matrix` from this module, or
load the CSVs under `lof_gof/` in step 4+ (drop `tp53_group` and `TP53_status` like step2).
"""

from __future__ import annotations

import argparse
import sys
from pathlib import Path

import pandas as pd

# Repo layout + import path (same pattern as step 2).
_SCRIPT_DIR = Path(__file__).resolve().parent
REPO_ROOT = _SCRIPT_DIR.parent
# Step-3-only outputs live here so other `data/processed/` paths stay untouched.
LOF_GOF_DIR = REPO_ROOT / "data" / "processed" / "lof_gof"
if str(_SCRIPT_DIR) not in sys.path:
    sys.path.insert(0, str(_SCRIPT_DIR))

from step1_digestion_and_processing import build_mutation_matrix

# Subset of MAF `Variant_Classification` strings: treat as loss-of-function TP53.
# MAF Variant_Classification values used in step1
LOF_VARIANT_CLASSES = frozenset(
    {
        "Nonsense_Mutation",
        "Frame_Shift_Ins",
        "Frame_Shift_Del",
        "Splice_Site",
    }
)
# Missense here proxies GoF / dominant-negative until you add external scores (REVEL, etc.).
MISSENSE_VARIANT_CLASSES = frozenset({"Missense_Mutation"})

TP53_GENE = "TP53"

# String labels written to CSV and used by step 4 for filtering.
LABEL_WT = "TP53_WT"
LABEL_LOF = "TP53_LoF"
LABEL_GOF = "TP53_GoF_missense"


def _sample_tp53_group(subdf: pd.DataFrame) -> str:
    """Assign TP53 functional group from all TP53 rows for one sample."""
    # Union of variant classes observed for this sample at TP53.
    classes = set(subdf["Variant_Classification"].dropna().astype(str).unique())
    # If both LoF and missense exist, LoF wins (explicit priority rule).
    if classes & LOF_VARIANT_CLASSES:
        return LABEL_LOF
    if classes & MISSENSE_VARIANT_CLASSES:
        return LABEL_GOF
    if classes:
        return "TP53_other_coding"
    return LABEL_WT


def classify_tp53_functional_groups(df: pd.DataFrame, all_luad_sample_ids) -> pd.DataFrame:
    """
    Build one row per sample with tp53_group from long-format merged mutation table.
    Samples with no TP53 rows are TP53_WT.
    """
    # Long-format merged table from step 1: only TP53 gene rows carry functional class info.
    tp53_rows = df[df["Hugo_Symbol"] == TP53_GENE]
    labels: dict[str, str] = {}
    for sample_id, grp in tp53_rows.groupby("Tumor_Sample_Barcode"):
        labels[str(sample_id)] = _sample_tp53_group(grp)

    # One row per LUAD patient; missing key => never had a TP53 hit in filtered data => WT.
    out = pd.DataFrame({"sample_id": pd.Index(all_luad_sample_ids).astype(str)})
    out["tp53_group"] = out["sample_id"].map(lambda s: labels.get(s, LABEL_WT))
    return out


def attach_tp53_group_to_matrix(
    mutation_matrix: pd.DataFrame,
    tp53_status: pd.DataFrame,
) -> pd.DataFrame:
    """Join tp53_group; align by matrix index (tumor sample barcode)."""
    id_to_group = dict(zip(tp53_status["sample_id"], tp53_status["tp53_group"]))
    mm = mutation_matrix.copy()
    # mutation_matrix index matches `Tumor_Sample_Barcode` / clinical `Sample ID`.
    mm["tp53_group"] = mm.index.astype(str).map(lambda s: id_to_group.get(s, LABEL_WT))
    return mm


def _validate_binary_vs_group(mutation_matrix: pd.DataFrame, tp53_status: pd.DataFrame) -> None:
    """TP53 column 0 should match TP53_WT; TP53 1 should not be WT."""
    mm = mutation_matrix.copy()
    # Quick invariant: binary TP53 gene column must agree with three-way label.
    mp = dict(zip(tp53_status["sample_id"].astype(str), tp53_status["tp53_group"]))
    for sid in mm.index.astype(str):
        g = mp.get(sid, LABEL_WT)
        tp53_bin = int(mm.loc[sid, TP53_GENE])
        if tp53_bin == 0 and g != LABEL_WT:
            raise ValueError(f"Sample {sid}: TP53=0 but tp53_group={g}")
        if tp53_bin == 1 and g == LABEL_WT:
            raise ValueError(f"Sample {sid}: TP53=1 but tp53_group=WT (contradiction)")


def split_by_tp53_group(mutation_matrix_with_group: pd.DataFrame):
    """
    Return three dataframes (WT / LoF / GoF+missense), each including gene columns
    plus tp53_group (and TP53_status if present). Use for Fisher heatmaps like step2.
    """
    # Convenience helper if you call step 3 functions from a notebook instead of step 4.
    mm = mutation_matrix_with_group
    wt = mm[mm["tp53_group"] == LABEL_WT].copy()
    lof = mm[mm["tp53_group"] == LABEL_LOF].copy()
    gof = mm[mm["tp53_group"] == LABEL_GOF].copy()
    return wt, lof, gof


def main() -> None:
    parser = argparse.ArgumentParser(description="Classify TP53 LoF vs GoF/missense (sample level).")
    parser.add_argument("--maf", default=None, help="Path to MAF-style mutations TSV.")
    parser.add_argument("--clinical", default=None, help="Path to clinical TSV.")
    parser.add_argument(
        "--out-status",
        default=str(LOF_GOF_DIR / "tp53_functional_status.csv"),
        help="Output table: sample_id, tp53_group.",
    )
    parser.add_argument(
        "--out-matrix",
        default=str(LOF_GOF_DIR / "mutation_matrix_with_tp53_group.csv"),
        help="Full matrix + TP53_status + tp53_group.",
    )
    args = parser.parse_args()

    # Rebuild the same objects as step 1 (long table + binary matrix + clinical LUAD list).
    maf_path = args.maf or None
    clin_path = args.clinical or None
    if maf_path and clin_path:
        df, mutation_matrix, clinical_la = build_mutation_matrix(maf_path, clin_path)
    elif maf_path or clin_path:
        raise ValueError("Pass both --maf and --clinical, or neither for defaults.")
    else:
        df, mutation_matrix, clinical_la = build_mutation_matrix()

    # Functional labels derived from TP53 rows' Variant_Classification in `df`.
    all_ids = clinical_la["Sample ID"].unique()
    tp53_status = classify_tp53_functional_groups(df, all_ids)

    # Match step 1 naming: Mut/WT column for readability in exported matrix.
    mutation_matrix = mutation_matrix.copy()
    mutation_matrix["TP53_status"] = mutation_matrix[TP53_GENE].apply(
        lambda x: "Mut" if int(x) == 1 else "WT"
    )
    mm_out = attach_tp53_group_to_matrix(mutation_matrix, tp53_status)
    _validate_binary_vs_group(mutation_matrix.drop(columns=["TP53_status"], errors="ignore"), tp53_status)

    # CSVs for step 4 (and notebooks): status table + full matrix with group column.
    out_status = Path(args.out_status)
    out_matrix = Path(args.out_matrix)
    out_status.parent.mkdir(parents=True, exist_ok=True)
    out_matrix.parent.mkdir(parents=True, exist_ok=True)

    tp53_status.to_csv(out_status, index=False)
    mm_out.to_csv(out_matrix, index_label="sample_id")

    print("TP53 functional group counts:")
    print(tp53_status["tp53_group"].value_counts())
    print(f"Wrote {out_status}")
    print(f"Wrote {out_matrix}")
    out_parquet = LOF_GOF_DIR.parent / "gene_mutation_binarized_matrix.parquet"
    mutation_matrix.drop(columns=["TP53_status"], errors="ignore").to_parquet(out_parquet)
    print(f"Wrote {out_parquet}")

if __name__ == "__main__":
    main()
