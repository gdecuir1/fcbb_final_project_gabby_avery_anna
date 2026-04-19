"""
step1_digestion_and_processing.py

Pipeline — **Step 1** (LUAD cohort, MSK-style tables)

Build a reproducible **sample × gene binary mutation matrix** for lung adenocarcinoma (LUAD)
from a MAF-like mutation table plus clinical metadata (e.g. cBioPortal exports).

What this step does
--------------------
- Loads mutation calls and clinical metadata (tab-separated).
- Keeps **coding / protein-altering** variant classes only (missense, nonsense, frameshift, splice).
- Restricts to a **curated driver gene list** (``GENE_LIST``; shared with steps 3–5).
- Restricts to **LUAD** patients via clinical ``Cancer Type`` / ``Cancer Type Detailed``.
- Merges mutations with LUAD samples and pivots to a **0/1 matrix** (all LUAD samples as rows;
  missing gene entries = 0).

Downstream
----------
- **Step 2** calls ``preprocess()`` for TP53-mut vs TP53-WT stratified Fisher heatmaps.
- **Step 3** calls ``build_mutation_matrix()`` for LoF / GoF / WT TP53 labels and exports.

Optional: MC3 (pan-cancer) data download
-----------------------------------------
To fetch NCI GDC MC3 publication supplements into ``data/raw/mc3/`` (used by **step 5**, not by
the LUAD ``build_mutation_matrix`` path here), run::

  python scripts/step1_digestion_and_processing.py --download-mc3

Use ``--mc3-open-only`` for open-access files only. Controlled files need ``GDC_TOKEN`` or
``--gdc-token``. See: https://gdc.cancer.gov/about-data/publications/mc3-2017

The downloader lives in ``mc3_gdc_download.py`` in the same directory as this script.
"""

from __future__ import annotations

import argparse
import os
import sys
from pathlib import Path

import pandas as pd

# Project root (parent of `scripts/`) so data paths resolve from repo root or from `scripts/`.
REPO_ROOT = Path(__file__).resolve().parent.parent

# ---------------------------------------------------------------------------
# Analysis parameters: which genes and which MAF variant classes to keep
# ---------------------------------------------------------------------------
GENE_LIST = [
    "TP53", "RB1", "ATM",
    "KRAS", "EGFR", "PIK3CA", "AKT3", "PTEN", "ERBB4",
    "APC", "CTNNB1",
    "STK11", "SMARCA4", "PBRM1", "CREBBP",
    "MGA", "PTPRD",
]

CODING_VARIANT_CLASSES = [
    "Missense_Mutation", "Nonsense_Mutation",
    "Frame_Shift_Ins", "Frame_Shift_Del",
    "Splice_Site",
]

DEFAULT_MAF_PATH = str(REPO_ROOT / "data" / "data_mutations.txt")
DEFAULT_CLINICAL_PATH = str(
    REPO_ROOT / "data" / "lung_msk_mind_2020_clinical_data (1).tsv"
)

# Default folder for GDC MC3 publication file downloads (ingestion only; not used by LUAD preprocess).
DEFAULT_MC3_RAW_DIR = REPO_ROOT / "data" / "raw" / "mc3"


def build_mutation_matrix(
    maf_path: str = DEFAULT_MAF_PATH,
    clinical_path: str = DEFAULT_CLINICAL_PATH,
) -> tuple[pd.DataFrame, pd.DataFrame, pd.DataFrame]:
    """
    Load MAF + clinical data, filter to LUAD and curated genes, return long merged
    table and sample × gene binary matrix (all LUAD samples, fill 0 for missing).

    Returns
    -------
    df : DataFrame
        Merged mutation rows (one row per variant event) after filters.
    mutation_matrix : DataFrame
        Rows = tumor sample barcode, columns = gene symbols, values 0/1.
    clinical_data_LA : DataFrame
        LUAD-only clinical subset used for the merge.
    """
    # Load full MAF-style mutation table (one row per variant call).
    data_all = pd.read_csv(Path(maf_path).expanduser(), sep="\t")

    # Keep only coding / protein-altering classes and the curated gene list.
    data_filtered = data_all[data_all["Variant_Classification"].isin(CODING_VARIANT_CLASSES)]
    data_filtered = data_filtered[data_filtered["Hugo_Symbol"].isin(GENE_LIST)]

    # Clinical metadata: restrict to lung adenocarcinoma (LUAD) rows.
    clinical_data_all = pd.read_csv(Path(clinical_path).expanduser(), sep="\t")
    clinical_data_LA = clinical_data_all[
        (clinical_data_all["Cancer Type"] == "Lung Adenocarcinoma")
        & (clinical_data_all["Cancer Type Detailed"] == "Lung Adenocarcinoma")
    ].copy()

    # Inner join: only samples that appear in both mutation and clinical tables.
    df = pd.merge(
        data_filtered,
        clinical_data_LA,
        left_on="Tumor_Sample_Barcode",
        right_on="Sample ID",
        how="inner",
    )

    # Wide matrix: index = tumor sample, columns = genes, values = 1 if any qualifying mutation.
    mutation_matrix = (
        df.assign(mut=1)
        .pivot_table(
            index="Tumor_Sample_Barcode",
            columns="Hugo_Symbol",
            values="mut",
            aggfunc="max",
            fill_value=0,
        )
    )

    # Reindex so every LUAD patient is a row (0 = no mutation in listed genes for that patient).
    all_luad_samples = clinical_data_LA["Sample ID"].unique()
    mutation_matrix = mutation_matrix.reindex(all_luad_samples, fill_value=0)

    return df, mutation_matrix, clinical_data_LA


def preprocess():
    """
    Interactive / demo path for **step 2**: build matrix, print sanity checks, split by TP53 Mut vs WT.

    Returns two matrices (TP53 mutated vs wild-type) each with a ``TP53_status`` column for stratified plots.
    """
    df, mutation_matrix, clinical_data_LA = build_mutation_matrix()

    # Size and ID checks after the merge.
    print(df.shape)
    print(df["Tumor_Sample_Barcode"].nunique())
    print(df["Sample ID"].nunique())

    print(mutation_matrix.shape)
    print(mutation_matrix.head())

    print(mutation_matrix.index.is_unique)
    print(sorted(pd.unique(mutation_matrix.values.ravel())))

    print((mutation_matrix.sum(axis=1) > 0).sum())
    print(mutation_matrix.columns.tolist())

    # Consistency: every (sample, gene) in long `df` must be 1 in the pivot.
    check_pairs = df[["Tumor_Sample_Barcode", "Hugo_Symbol"]].drop_duplicates()
    all_checks_pass = True
    for _, row in check_pairs.iterrows():
        sample = row["Tumor_Sample_Barcode"]
        gene = row["Hugo_Symbol"]
        if mutation_matrix.loc[sample, gene] != 1:
            all_checks_pass = False
            print(f"Mismatch found: sample={sample}, gene={gene}")

    print(all_checks_pass)

    # Binary TP53 stratification for step 2 (mutated vs wild-type column).
    mutation_matrix = mutation_matrix.copy()
    mutation_matrix["TP53_status"] = mutation_matrix["TP53"].apply(
        lambda x: "Mut" if x == 1 else "WT"
    )

    print("TP53 status counts:")
    print(mutation_matrix["TP53_status"].value_counts())

    tp53_mut = mutation_matrix[mutation_matrix["TP53_status"] == "Mut"].copy()
    tp53_wt = mutation_matrix[mutation_matrix["TP53_status"] == "WT"].copy()

    return tp53_mut, tp53_wt


def _run_mc3_gdc_download(
    out_dir: Path,
    *,
    gdc_token: str | None,
    open_only: bool,
    verify_md5: bool,
    overwrite: bool,
) -> None:
    """Delegate to ``mc3_gdc_download`` (same directory as this script)."""
    script_dir = Path(__file__).resolve().parent
    if str(script_dir) not in sys.path:
        sys.path.insert(0, str(script_dir))
    from mc3_gdc_download import download_mc3_publication_supplements

    download_mc3_publication_supplements(
        out_dir,
        gdc_token=gdc_token,
        skip_controlled=open_only,
        verify_md5=verify_md5,
        overwrite=overwrite,
    )


if __name__ == "__main__":
    parser = argparse.ArgumentParser(
        description="Step 1: LUAD mutation preprocessing, or MC3 GDC raw data download.",
    )
    parser.add_argument(
        "--download-mc3",
        action="store_true",
        help="Download MC3 publication files from GDC into data/raw/mc3 (see mc3_gdc_download).",
    )
    parser.add_argument(
        "--mc3-out",
        default=str(DEFAULT_MC3_RAW_DIR),
        help="Output directory for --download-mc3 (default: <repo>/data/raw/mc3).",
    )
    parser.add_argument(
        "--gdc-token",
        default=os.environ.get("GDC_TOKEN"),
        help="GDC download token for controlled MC3 files (or set env GDC_TOKEN).",
    )
    parser.add_argument(
        "--mc3-open-only",
        action="store_true",
        help="With --download-mc3: only download open-access targets (skip controlled).",
    )
    parser.add_argument(
        "--mc3-verify-md5",
        action="store_true",
        help="With --download-mc3: verify MD5 after each file (slow for large archives).",
    )
    parser.add_argument(
        "--mc3-overwrite",
        action="store_true",
        help="With --download-mc3: re-download even if file exists with expected size.",
    )
    args = parser.parse_args()

    if args.download_mc3:
        out = Path(args.mc3_out).expanduser()
        if not out.is_absolute():
            out = REPO_ROOT / out
        _run_mc3_gdc_download(
            out,
            gdc_token=args.gdc_token,
            open_only=args.mc3_open_only,
            verify_md5=args.mc3_verify_md5,
            overwrite=args.mc3_overwrite,
        )
    else:
        tp53_mut, tp53_wt = preprocess()
