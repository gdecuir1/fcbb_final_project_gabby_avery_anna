"""
step5_scale_discover_pan_cancer.py

Pan-cancer scaling using the MC3 open-access public MAF (after step 1 ``--download-mc3``).

The Bioconductor DISCOVER method (R) is not bundled here; this step implements a
transparent Python analogue used in many papers: sample-level binary mutation
matrices, pairwise Fisher exact tests for co-occurrence / mutual exclusivity, and
Benjamini–Hochberg FDR on the resulting p-values. Results are written to
``outputs/tables/`` and summarized in ``outputs/reports/``.

Inputs
------
- ``data/raw/mc3/mc3.v0.2.8.PUBLIC.maf.gz`` (preferred), or pass ``--mc3-maf``.
- If that file is missing and you omit ``--mc3-maf``, the script falls back to
  ``tests/fixtures/minimal_mc3_public.maf`` with a **warning** (for local runs without the full download).

Outputs
-------
- ``outputs/tables/step5_pancancer_pairwise_fisher.tsv``
- ``outputs/tables/step5_mutation_frequency_by_gene.tsv``
- ``outputs/tables/step5_mutation_frequency_by_study_gene.tsv`` (if a study column exists)
- ``outputs/figures/step5_pancancer_fisher_logp_heatmap.png``
- ``outputs/reports/step5_pan_cancer_analytics_report.md``
"""

from __future__ import annotations

import argparse
import gzip
import sys
from collections import defaultdict
from datetime import datetime, timezone
from itertools import combinations
from pathlib import Path

import matplotlib.pyplot as plt
import numpy as np
import pandas as pd
import seaborn as sns
from scipy.stats import fisher_exact

_SCRIPT_DIR = Path(__file__).resolve().parent
REPO_ROOT = _SCRIPT_DIR.parent
if str(_SCRIPT_DIR) not in sys.path:
    sys.path.insert(0, str(_SCRIPT_DIR))

from step1_digestion_and_processing import CODING_VARIANT_CLASSES, GENE_LIST
from step3_classify_tp53_lof_vs_gof import (
    LABEL_GOF,
    LABEL_LOF,
    LABEL_WT,
    LOF_VARIANT_CLASSES,
    MISSENSE_VARIANT_CLASSES,
)

DEFAULT_MC3_MAF = REPO_ROOT / "data" / "raw" / "mc3" / "mc3.v0.2.8.PUBLIC.maf.gz"
DEFAULT_MC3_MAF_UNCOMPRESSED = REPO_ROOT / "data" / "raw" / "mc3" / "mc3.v0.2.8.PUBLIC.maf"
MINIMAL_MC3_FIXTURE = REPO_ROOT / "tests" / "fixtures" / "minimal_mc3_public.maf"


def resolve_mc3_maf_path(explicit: str | None) -> Path:
    """
    If ``explicit`` is set, use that path (must exist).
    Otherwise try: public .maf.gz → public .maf → bundled minimal fixture.
    """
    if explicit:
        p = Path(explicit).expanduser()
        if not p.is_absolute():
            p = REPO_ROOT / p
        if not p.exists():
            raise FileNotFoundError(
                f"MC3 MAF not found: {p}\n"
                "Download open-access data with:\n"
                "  python scripts/step1_digestion_and_processing.py --download-mc3 --mc3-open-only\n"
                "or pass a valid file:  --mc3-maf path/to/your.maf.gz"
            )
        return p

    candidates = [
        DEFAULT_MC3_MAF,
        DEFAULT_MC3_MAF_UNCOMPRESSED,
        MINIMAL_MC3_FIXTURE,
    ]
    for c in candidates:
        if c.exists():
            if c == MINIMAL_MC3_FIXTURE:
                print(
                    "WARNING: Full MC3 public MAF not found; using tests/fixtures/minimal_mc3_public.maf "
                    "for a smoke run. Download the real file for pan-cancer results:\n"
                    "  python scripts/step1_digestion_and_processing.py --download-mc3 --mc3-open-only",
                    file=sys.stderr,
                )
            return c

    raise FileNotFoundError(
        "No MAF found. Tried:\n"
        f"  - {DEFAULT_MC3_MAF}\n"
        f"  - {DEFAULT_MC3_MAF_UNCOMPRESSED}\n"
        f"  - {MINIMAL_MC3_FIXTURE}\n"
        "Download open-access MC3 with:\n"
        "  python scripts/step1_digestion_and_processing.py --download-mc3 --mc3-open-only\n"
        "or pass an existing file:  --mc3-maf /path/to/mc3.v0.2.8.PUBLIC.maf.gz"
    )

REQUIRED_MAF_COLS = ("Hugo_Symbol", "Tumor_Sample_Barcode", "Variant_Classification")

# Prefer explicit MC3 / TCGA-style study labels if present (column names vary by MAF version).
STUDY_COLUMN_CANDIDATES = (
    "CODE",
    "Code",
    "Study",
    "study",
    "ONCOTREE_CODE",
    "project_id",
    "Project",
    "Disease",
    "disease_type",
    "primary_site",
)


def _df_to_markdown(df: pd.DataFrame) -> str:
    try:
        return df.to_markdown(index=False)
    except Exception:
        return "```\n" + df.to_string(index=False) + "\n```"


def _series_to_markdown(s: pd.Series, *, key: str, val: str) -> str:
    t = s.reset_index()
    t.columns = [key, val]
    return _df_to_markdown(t)


def benjamini_hochberg_fdr(p_values: np.ndarray) -> np.ndarray:
    """Benjamini–Hochberg FDR; return q-values in original order."""
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


def _detect_study_column(columns: list[str]) -> str | None:
    colset = set(columns)
    for cand in STUDY_COLUMN_CANDIDATES:
        if cand in colset:
            return cand
    for c in columns:
        cl = c.lower()
        if "project" in cl and "uuid" not in cl:
            return c
        if cl in ("code", "study", "cancer", "tumor_type", "primary_diagnosis"):
            return c
    return None


def _tp53_group_from_variant_classes(classes: set[str]) -> str:
    if classes & LOF_VARIANT_CLASSES:
        return LABEL_LOF
    if classes & MISSENSE_VARIANT_CLASSES:
        return LABEL_GOF
    if classes:
        return "TP53_other_coding"
    return LABEL_WT


def load_mc3_pairs_and_tp53(
    maf_path: Path,
    genes: list[str],
    coding_classes: list[str],
    chunksize: int,
    max_rows: int | None,
) -> tuple[pd.DataFrame, pd.Series, list[str], pd.Series]:
    """
    Stream MAF in chunks; return long (sample, gene) pairs, study per sample, column list, TP53 group.
    """
    gene_set = set(genes)
    coding_set = set(coding_classes)
    open_fn = gzip.open if str(maf_path).endswith(".gz") else open

    with open_fn(maf_path, "rt", encoding="utf-8", errors="replace") as f:
        header = f.readline()
        if not header:
            raise ValueError("Empty MAF")
        all_columns = header.strip().split("\t")

    missing = [c for c in REQUIRED_MAF_COLS if c not in all_columns]
    if missing:
        raise ValueError(f"MAF missing required columns {missing}; found: {all_columns[:40]}…")

    study_col = _detect_study_column(all_columns)
    usecols = list(REQUIRED_MAF_COLS)
    if study_col and study_col not in usecols:
        usecols.append(study_col)

    pair_parts: list[pd.DataFrame] = []
    study_first: dict[str, str] = {}
    tp53_classes: dict[str, set[str]] = defaultdict(set)

    rows_read = 0
    compression = "gzip" if str(maf_path).endswith(".gz") else None

    for chunk in pd.read_csv(
        maf_path,
        sep="\t",
        compression=compression,
        usecols=lambda c: c in usecols,
        chunksize=chunksize,
        low_memory=False,
    ):
        if max_rows is not None and rows_read >= max_rows:
            break
        if max_rows is not None:
            remain = max_rows - rows_read
            if len(chunk) > remain:
                chunk = chunk.iloc[:remain].copy()
        rows_read += len(chunk)

        chunk = chunk[
            chunk["Hugo_Symbol"].isin(gene_set)
            & chunk["Variant_Classification"].isin(coding_set)
        ]
        if chunk.empty:
            continue

        if study_col:
            for barcode, val in zip(chunk["Tumor_Sample_Barcode"], chunk[study_col]):
                if pd.isna(val):
                    continue
                b = str(barcode)
                if b not in study_first:
                    study_first[b] = str(val)

        tp53_sub = chunk[chunk["Hugo_Symbol"] == "TP53"]
        for barcode, vc in zip(tp53_sub["Tumor_Sample_Barcode"], tp53_sub["Variant_Classification"]):
            tp53_classes[str(barcode)].add(str(vc))

        pairs = chunk[["Tumor_Sample_Barcode", "Hugo_Symbol"]].drop_duplicates()
        pair_parts.append(pairs)

    if not pair_parts:
        raise ValueError("No rows passed filters (genes + coding classes). Check MAF path and filters.")

    all_pairs = pd.concat(pair_parts, ignore_index=True).drop_duplicates()
    all_pairs["Tumor_Sample_Barcode"] = all_pairs["Tumor_Sample_Barcode"].astype(str)

    tp53_group = pd.Series(
        {
            s: _tp53_group_from_variant_classes(tp53_classes.get(s, set()))
            for s in all_pairs["Tumor_Sample_Barcode"].unique()
        },
        dtype=object,
    )

    study_series = pd.Series(study_first, dtype=object) if study_first else pd.Series(dtype=object)

    return all_pairs, study_series, all_columns, tp53_group


def pairs_to_matrix(pairs: pd.DataFrame, genes: list[str]) -> pd.DataFrame:
    samples = pairs["Tumor_Sample_Barcode"].unique()
    wide = (
        pairs.assign(v=1)
        .pivot_table(
            index="Tumor_Sample_Barcode",
            columns="Hugo_Symbol",
            values="v",
            aggfunc="max",
            fill_value=0,
        )
        .reindex(columns=genes, fill_value=0)
        .reindex(samples, fill_value=0)
    )
    return wide.astype(np.int8)


def pairwise_fisher_results(matrix: pd.DataFrame, genes: list[str]) -> pd.DataFrame:
    rows = []
    for g1, g2 in combinations(genes, 2):
        a = int(((matrix[g1] == 1) & (matrix[g2] == 1)).sum())
        b = int(((matrix[g1] == 1) & (matrix[g2] == 0)).sum())
        c = int(((matrix[g1] == 0) & (matrix[g2] == 1)).sum())
        d = int(((matrix[g1] == 0) & (matrix[g2] == 0)).sum())
        table = [[a, b], [c, d]]
        oddsratio, pval = fisher_exact(table)
        rows.append(
            {
                "gene_a": g1,
                "gene_b": g2,
                "n_both": a,
                "g1_only": b,
                "g2_only": c,
                "neither": d,
                "odds_ratio": oddsratio,
                "p_value": pval,
            }
        )
    out = pd.DataFrame(rows)
    if len(out):
        out["q_value_bh"] = benjamini_hochberg_fdr(out["p_value"].values)
    return out


def plot_logp_heatmap(pivot_logp: pd.DataFrame, out_path: Path, title: str) -> None:
    out_path.parent.mkdir(parents=True, exist_ok=True)
    mask = np.triu(np.ones_like(pivot_logp.values, dtype=bool))
    plt.figure(figsize=(12, 10))
    sns.heatmap(
        pivot_logp,
        mask=mask,
        cmap="RdYlBu_r",
        center=np.nanmedian(pivot_logp.values),
        linewidths=0.3,
        cbar_kws={"label": "-log10(p)"},
    )
    plt.title(title)
    plt.tight_layout()
    plt.savefig(out_path, dpi=150, bbox_inches="tight")
    plt.close()


def write_report(
    report_path: Path,
    *,
    maf_path: Path,
    n_samples: int,
    n_pairs_rows: int,
    n_fisher_rows: int,
    study_col: str | None,
    n_studies: int | None,
    top_exclusive: pd.DataFrame,
    top_cooccur: pd.DataFrame,
    gene_prev: pd.Series,
    tp53_group_counts: pd.Series,
) -> None:
    report_path.parent.mkdir(parents=True, exist_ok=True)
    lines = [
        "# Step 5 — Pan-cancer analytics (MC3 public MAF)",
        "",
        f"- **Generated (UTC):** {datetime.now(timezone.utc).strftime('%Y-%m-%d %H:%M:%SZ')}",
        f"- **MAF:** `{maf_path}`",
        f"- **Samples (after filters):** {n_samples:,}",
        f"- **Curated genes:** {len(GENE_LIST)}",
        f"- **Study column used:** {study_col or '*not detected — pan-cancer only*'}",
    ]
    if n_studies is not None:
        lines.append(f"- **Distinct study labels:** {n_studies}")
    lines.extend(
        [
            "",
            "## Methods (Python stand-in for DISCOVER screening)",
            "",
            "1. Stream the MC3 open-access MAF in chunks; keep **coding** variant classes matching step 1.",
            "2. Restrict to the same **driver gene list** as step 1.",
            "3. Build a **sample × gene** binary mutation matrix.",
            "4. **TP53** is classified per sample (LoF vs missense proxy vs WT) using the same rules as step 3.",
            "5. Run **pairwise Fisher exact tests** on all unordered gene pairs (co-occurrence table).",
            "6. Apply **Benjamini–Hochberg FDR** (`q_value_bh`) across tests.",
            "",
            "The R/Bioconductor **DISCOVER** package uses a dedicated statistical model; this pipeline is a",
            "lighter-weight **Fisher-based screen** suitable for coursework and for spotting pairs to follow up.",
            "",
            "## Pan-cancer mutation prevalence (curated genes)",
            "",
            _series_to_markdown(
                gene_prev.sort_values(ascending=False),
                key="gene",
                val="fraction_samples_mutated",
            ),
            "",
            "## TP53 functional groups (same rules as step 3)",
            "",
            _series_to_markdown(tp53_group_counts, key="tp53_group", val="n_samples"),
            "",
            "## Strongest mutual exclusivity signals (low odds ratio, small q)",
            "",
            "*(Fisher odds ratio below 1 suggests fewer co-mutations than independence would predict.)*",
            "",
        ]
    )
    if top_exclusive.empty:
        lines.append("*No pairs met reporting thresholds.*")
    else:
        lines.append(_df_to_markdown(top_exclusive))
    lines.extend(
        [
            "",
            "## Strongest co-occurrence signals (high odds ratio, small q)",
            "",
        ]
    )
    if top_cooccur.empty:
        lines.append("*No pairs met reporting thresholds.*")
    else:
        lines.append(_df_to_markdown(top_cooccur))
    lines.extend(
        [
            "",
            "## Outputs",
            "",
            "- `outputs/tables/step5_pancancer_pairwise_fisher.tsv`",
            "- `outputs/tables/step5_mutation_frequency_by_gene.tsv`",
            "- `outputs/tables/step5_mutation_frequency_by_study_gene.tsv` (if study labels available)",
            "- `outputs/tables/step5_tp53_group_counts.tsv`",
            "- `outputs/figures/step5_pancancer_fisher_logp_heatmap.png`",
            "",
        ]
    )
    report_path.write_text("\n".join(lines), encoding="utf-8")


def main() -> None:
    parser = argparse.ArgumentParser(
        description="Pan-cancer MC3 analytics + report (Fisher screen; DISCOVER-style summary).",
    )
    parser.add_argument(
        "--mc3-maf",
        default=None,
        help="Path to MAF (.maf or .maf.gz). If omitted, auto-detect: data/raw/mc3/*.PUBLIC.maf* then tests/fixtures/minimal_mc3_public.maf",
    )
    parser.add_argument(
        "--chunksize",
        type=int,
        default=150_000,
        help="Rows per chunk when streaming the MAF.",
    )
    parser.add_argument(
        "--max-rows",
        type=int,
        default=None,
        help="Optional cap on MAF rows read (debug / laptops).",
    )
    parser.add_argument(
        "--out-reports",
        default=str(REPO_ROOT / "outputs" / "reports"),
        help="Directory for Markdown report.",
    )
    parser.add_argument(
        "--out-tables",
        default=str(REPO_ROOT / "outputs" / "tables"),
        help="Directory for TSV tables.",
    )
    parser.add_argument(
        "--out-figures",
        default=str(REPO_ROOT / "outputs" / "figures"),
        help="Directory for figures.",
    )
    parser.add_argument("--no-show", action="store_true", help="Do not display plots interactively.")
    args = parser.parse_args()

    maf_path = resolve_mc3_maf_path(args.mc3_maf)

    out_reports = Path(args.out_reports)
    out_tables = Path(args.out_tables)
    out_figures = Path(args.out_figures)
    if not out_reports.is_absolute():
        out_reports = REPO_ROOT / out_reports
    if not out_tables.is_absolute():
        out_tables = REPO_ROOT / out_tables
    if not out_figures.is_absolute():
        out_figures = REPO_ROOT / out_figures

    out_tables.mkdir(parents=True, exist_ok=True)
    out_figures.mkdir(parents=True, exist_ok=True)
    out_reports.mkdir(parents=True, exist_ok=True)

    print(f"Streaming MAF: {maf_path}")
    pairs, study_series, all_columns, tp53_group = load_mc3_pairs_and_tp53(
        maf_path,
        GENE_LIST,
        CODING_VARIANT_CLASSES,
        chunksize=args.chunksize,
        max_rows=args.max_rows,
    )

    matrix = pairs_to_matrix(pairs, GENE_LIST)
    n_samples = len(matrix)
    print(f"Samples: {n_samples}; genes: {len(GENE_LIST)}")

    tp53_aligned = tp53_group.reindex(matrix.index).fillna(LABEL_WT)
    tp53_group_counts = tp53_aligned.value_counts()
    tp53_group_counts.to_frame("n_samples").to_csv(
        out_tables / "step5_tp53_group_counts.tsv",
        sep="\t",
        index_label="tp53_group",
    )

    # Mutation frequencies
    gene_prev = matrix.mean(axis=0).sort_values(ascending=False)
    prev_df = gene_prev.reset_index()
    prev_df.columns = ["gene", "fraction_samples_mutated"]
    prev_df["n_samples_mutated"] = (matrix[prev_df["gene"]] == 1).sum().values
    prev_df.to_csv(out_tables / "step5_mutation_frequency_by_gene.tsv", sep="\t", index=False)

    study_col_name = _detect_study_column(list(all_columns))
    if study_col_name and len(study_series):
        study_series = study_series.reindex(matrix.index)
        by_study = []
        for study in study_series.dropna().unique():
            idx = matrix.index[study_series == study]
            sub = matrix.loc[idx]
            for g in GENE_LIST:
                by_study.append(
                    {
                        "study": study,
                        "gene": g,
                        "n_samples": len(sub),
                        "n_mutated": int((sub[g] == 1).sum()),
                        "fraction": float((sub[g] == 1).mean()) if len(sub) else 0.0,
                    }
                )
        pd.DataFrame(by_study).to_csv(
            out_tables / "step5_mutation_frequency_by_study_gene.tsv",
            sep="\t",
            index=False,
        )
        n_studies = study_series.nunique(dropna=True)
    else:
        n_studies = None
        study_col_name = None

    fisher_df = pairwise_fisher_results(matrix, GENE_LIST)
    fisher_df.to_csv(out_tables / "step5_pancancer_pairwise_fisher.tsv", sep="\t", index=False)

    # Heatmap of -log10 p
    logp = pd.DataFrame(np.nan, index=GENE_LIST, columns=GENE_LIST, dtype=float)
    for _, r in fisher_df.iterrows():
        g1, g2 = r["gene_a"], r["gene_b"]
        lp = -np.log10(max(float(r["p_value"]), 1e-300))
        logp.loc[g1, g2] = lp
        logp.loc[g2, g1] = lp
    fig_path = out_figures / "step5_pancancer_fisher_logp_heatmap.png"
    plot_logp_heatmap(
        logp,
        fig_path,
        title=f"Pan-cancer Fisher pairwise -log10(p), n={n_samples} samples (MC3 public)",
    )

    # Top tables for report
    sig = fisher_df[(fisher_df["q_value_bh"] <= 0.25) & (fisher_df["n_both"] + fisher_df["g1_only"] + fisher_df["g2_only"] > 0)].copy()
    top_ex = sig.nsmallest(15, "odds_ratio")[
        ["gene_a", "gene_b", "odds_ratio", "p_value", "q_value_bh", "n_both", "g1_only", "g2_only"]
    ]
    top_co = sig.nlargest(15, "odds_ratio")[
        ["gene_a", "gene_b", "odds_ratio", "p_value", "q_value_bh", "n_both", "g1_only", "g2_only"]
    ]

    report_path = out_reports / "step5_pan_cancer_analytics_report.md"
    write_report(
        report_path,
        maf_path=maf_path,
        n_samples=n_samples,
        n_pairs_rows=len(pairs),
        n_fisher_rows=len(fisher_df),
        study_col=study_col_name,
        n_studies=int(n_studies) if n_studies is not None else None,
        top_exclusive=top_ex,
        top_cooccur=top_co,
        gene_prev=gene_prev,
        tp53_group_counts=tp53_group_counts,
    )

    print(f"Wrote tables under {out_tables}")
    print(f"Wrote figure {fig_path}")
    print(f"Wrote report {report_path}")

    if not args.no_show:
        plt.figure(figsize=(10, 8))
        sns.heatmap(logp, mask=np.triu(np.ones_like(logp, dtype=bool)), cmap="RdYlBu_r", cbar_kws={"label": "-log10(p)"})
        plt.title("Pan-cancer Fisher (preview)")
        plt.tight_layout()
        plt.show()


if __name__ == "__main__":
    main()
