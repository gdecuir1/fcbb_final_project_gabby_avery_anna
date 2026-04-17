"""
step7_alternative_hypothesis_structural_vs_point.py

Hypothesis
----------
TP53-mutant tumors will favor copy-number alterations (CNA amplifications in PI3K/RTK), 
whereas TP53-WT tumors will favor point mutations.

Inputs
------
- `data/processed/lof_gof/tp53_functional_status.csv` (from step 3)
- `data/processed/gene_mutation_binarized_matrix.parquet` (point mutations from step 3)
- A CNA-derived table (created here using data_cna.txt):

Outputs
--------
- Summary tables comparing TP53 LoF vs WT:
  - `outputs/tables/structural_vs_point_tp53_lof_vs_wt.tsv`

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
- 
"""

from __future__ import annotations

from pathlib import Path
import argparse

import matplotlib.pyplot as plt
import numpy as np
import pandas as pd
from scipy.stats import fisher_exact

_REPO = Path(__file__).resolve().parent.parent

_DEFAULT_TP53_FUNCTIONAL = _REPO / "data" / "processed" / "lof_gof" / "tp53_functional_status.csv" # from step 3
_DEFAULT_MAF = _REPO / "data" / "data_mutations.txt" # downloaded from cbioportal 
_DEFAULT_CLINICAL = _REPO / "data" / "lung_msk_mind_2020_clinical_data (1).tsv" # downloaded from cbioportal
_DEFAULT_POINT_MUT = _REPO / "data" / "processed" / "expanded_point_mutation_matrix.parquet" # used the same logic from step 1
_DEFAULT_RAW_CNA = _REPO / "data" / "data_cna.txt" # downloaded from cbioportal
_DEFAULT_CNA_OUT = _REPO / "data" / "processed" / "cna_amplification_binarized_matrix.parquet" # created in this script

LABEL_LOF = "TP53_LoF"
LABEL_WT = "TP53_WT"

CODING_VARIANT_CLASSES = [
    "Missense_Mutation", "Nonsense_Mutation",
    "Frame_Shift_Ins", "Frame_Shift_Del",
    "Splice_Site",
]

STRUCTURAL_POINT_GENE_LIST = [
    "EGFR", "ERBB2", "ERBB3", "ERBB4", "MET", "ALK", "ROS1", "RET", "AXL", "KIT", "IGF1R",
    "FGFR1", "FGFR2", "FGFR3", "FGFR4",
    "KRAS", "NRAS", "HRAS", "BRAF", "RAF1", "ARAF", "MAPK1", "MAPK3", "MAP2K1", "MAP2K2",
    "SOS1", "PTPN11",
    "PIK3CA", "PIK3CB", "PIK3CD", "AKT1", "AKT2", "AKT3", "PDPK1", "MTOR", "RPTOR", "RICTOR", "RHEB",
]

# converts an input path (string) into a Path 
def _resolve_repo_path(p: str | Path) -> Path:
    path = Path(p).expanduser()
    return path if path.is_absolute() else _REPO / path

# add script for compatibility 
from script.step5_scale_discover_pan_cancer import benjamini_hochberg_fdr, _df_to_markdown

# tests whether a binary feature is more frequent/enriched in one group vs. another using Fisher's exact test 
def run_fisher_binary_feature(
    feature: pd.Series,
    groups: pd.Series,
    group_a: str = LABEL_LOF,
    group_b: str = LABEL_WT,
) -> dict[str, object]:
    # keep only two groups to compare (drop samples in other TP53 class)
    mask = groups.isin([group_a, group_b])
    feature = feature.loc[mask].astype(int)
    groups = groups.loc[mask].astype(str)

    # count positives/negatives in each group 
    a_pos = int(((groups == group_a) & (feature == 1)).sum())
    a_neg = int(((groups == group_a) & (feature == 0)).sum())
    b_pos = int(((groups == group_b) & (feature == 1)).sum())
    b_neg = int(((groups == group_b) & (feature == 0)).sum())

    # build a 2x2 contingency table
    table = [[a_pos, a_neg], [b_pos, b_neg]]

    # run fisher's exact test
    odds_ratio, p_value = fisher_exact(table)

    return {
        "group_a": group_a,
        "group_b": group_b,
        "n_group_a_positive": a_pos,
        "n_group_a_negative": a_neg,
        "n_group_b_positive": b_pos,
        "n_group_b_negative": b_neg,
        "odds_ratio": odds_ratio, # odds_ratio > 1 = feature is more common in group_a than group_b
        "p_value": p_value, 
    }

# loads tp53_functional_status.csv and returns table mapping each sample to its TP53 group
def load_tp53_status(path: str | Path) -> pd.DataFrame:
    p = _resolve_repo_path(path)
    df = pd.read_csv(p, usecols=["sample_id", "tp53_group"], dtype=str)

    # basic cleanup
    df["sample_id"] = df["sample_id"].str.strip()
    df["tp53_group"] = df["tp53_group"].str.strip()

    # drop missing / duplicates
    df = df.dropna(subset=["sample_id", "tp53_group"]).drop_duplicates(subset=["sample_id"])

    return df


# recreate a sample × gene binary mutation matrix from step 1 with STRUCTURAL_POINT_GENE_LIST
def build_point_mutation_matrix(
    maf_path: str | Path,
    clinical_path: str | Path,
    out_path: str | Path
) -> pd.DataFrame:
    maf_path = Path(maf_path).expanduser()
    clinical_path = Path(clinical_path).expanduser()
    out_path = Path(out_path).expanduser()

    # load raw files
    maf = pd.read_csv(maf_path, sep="\t", low_memory=False) # data_mutations.txt
    clinical = pd.read_csv(clinical_path, sep="\t", low_memory=False) # lung_msk_mind_2020_clinical_data (1).tsv

    # restrict to only include coding mutations 
    maf = maf[
        maf["Variant_Classification"].isin(CODING_VARIANT_CLASSES)
    ].copy()

    # filter for only genes in STRUCTURAL_POINT_GENE_LIST
    maf = maf[
        maf["Hugo_Symbol"].isin(STRUCTURAL_POINT_GENE_LIST)
    ].copy()

    # filter for only LUAD patients 
    clinical_luad = clinical[
        (clinical["Cancer Type"] == "Lung Adenocarcinoma") &
        (clinical["Cancer Type Detailed"] == "Lung Adenocarcinoma")
    ].copy()

    # keep only LUAD tumor samples 
    merged = pd.merge(
        maf,
        clinical_luad,
        left_on="Tumor_Sample_Barcode",
        right_on="Sample ID",
        how="inner",
    )

    # build binary matrix
    mutation_matrix = (
        merged.assign(mut=1)
        .pivot_table(
            index="Tumor_Sample_Barcode",
            columns="Hugo_Symbol",
            values="mut",
            aggfunc="max",
            fill_value=0,
        )
    )

    # reindex so all LUAD samples appear, even if they have no selected mutations
    all_luad_samples = clinical_luad["Sample ID"].astype(str).unique()
    mutation_matrix = mutation_matrix.reindex(
        all_luad_samples,
        fill_value=0,
    )

    # ensure all genes in STRUCTURAL_POINT_GENE_LIST are present as columns
    mutation_matrix = mutation_matrix.reindex(
        columns=STRUCTURAL_POINT_GENE_LIST,
        fill_value=0,
    )

    mutation_matrix.index = mutation_matrix.index.astype(str)
    mutation_matrix = mutation_matrix.astype("int8")

    # save parquet
    out_path.parent.mkdir(parents=True, exist_ok=True)
    mutation_matrix.to_parquet(out_path)

    print(f"Wrote point mutation matrix: {out_path}")
    print(f"Shape: {mutation_matrix.shape}")
    print(f"Samples: {mutation_matrix.shape[0]}")
    print(f"Genes: {mutation_matrix.shape[1]}")
    print("Mutated LUAD samples per gene:")
    print(mutation_matrix.sum(axis=0).sort_values(ascending=False))

    return mutation_matrix

# make sure we've restricted to out 17 patients 
def load_raw_cna_matrix(
        path: str | Path, 
        clinical_path: str | Path
    ) -> pd.DataFrame:

    # load clinical + filter LUAD
    clinical_p = _resolve_repo_path(clinical_path)
    clinical = pd.read_csv(clinical_p, sep="\t", low_memory=False)

    # filter for only LUAD patients 
    clinical_luad = clinical[
        (clinical["Cancer Type"] == "Lung Adenocarcinoma") &
        (clinical["Cancer Type Detailed"] == "Lung Adenocarcinoma")
    ].copy()

    luad_sample_ids = set(clinical_luad["Sample ID"].astype(str).str.strip())

    # load CNA matrix
    p = _resolve_repo_path(path)
    df = pd.read_csv(p, sep="\t", low_memory=False)

    # ensure/clean gene symbol column
    if "Hugo_Symbol" not in df.columns:
        df = df.rename(columns={df.columns[0]: "Hugo_Symbol"})
    df["Hugo_Symbol"] = df["Hugo_Symbol"].astype(str).str.strip()

    # restrict to LUAD sample columns
    keep_sample_cols = [c for c in df.columns[1:] if str(c).strip() in luad_sample_ids]
    out = df[["Hugo_Symbol"] + keep_sample_cols].copy()

    return out

# converts raw CNA table into a sample x gene matrix 
def cna_gene_by_sample_to_sample_by_gene(cna_raw: pd.DataFrame) -> pd.DataFrame:
    cna = cna_raw.copy()
    cna = cna.drop_duplicates(subset=["Hugo_Symbol"], keep="first")
    cna = cna.set_index("Hugo_Symbol")
    cna = cna.apply(pd.to_numeric, errors="coerce")
    cna = cna.transpose()
    cna.index = cna.index.astype(str)
    cna.columns = cna.columns.astype(str)
    return cna


def build_cna_amplification_matrix(cna_sample_by_gene: pd.DataFrame, threshold: int = 2) -> pd.DataFrame:
    amp = (cna_sample_by_gene >= threshold).fillna(False).astype(np.int8)
    amp.index = amp.index.astype(str)
    amp.columns = amp.columns.astype(str)
    return amp


def save_cna_amplification_matrix(cna_amp: pd.DataFrame, out_path: str | Path) -> Path:
    out = _resolve_repo_path(out_path)
    out.parent.mkdir(parents=True, exist_ok=True)
    cna_amp.to_parquet(out)
    return out

def restrict_to_gene_set(df: pd.DataFrame, gene_list: list[str]) -> pd.DataFrame:
    keep = [g for g in gene_list if g in df.columns]
    return df.loc[:, keep].copy()


def harmonize_inputs(
    tp53_status: pd.DataFrame,
    point_mut: pd.DataFrame,
    cna_amp: pd.DataFrame,
    gene_list: list[str],
) -> tuple[pd.DataFrame, pd.DataFrame, pd.DataFrame, list[str]]:
    point_mut = restrict_to_gene_set(point_mut, gene_list)
    cna_amp = restrict_to_gene_set(cna_amp, gene_list)

    final_genes = sorted(set(point_mut.columns).intersection(cna_amp.columns))
    if not final_genes:
        raise ValueError("No overlapping genes between point mutation and CNA matrices in the selected gene set.")

    point_mut = point_mut.loc[:, final_genes]
    cna_amp = cna_amp.loc[:, final_genes]

    tp53_status = tp53_status[tp53_status["tp53_group"].isin([LABEL_LOF, LABEL_WT])].copy()

    shared_samples = sorted(set(tp53_status["sample_id"]).intersection(point_mut.index).intersection(cna_amp.index))
    if not shared_samples:
        raise ValueError("No overlapping samples across TP53 status, point mutation matrix, and CNA matrix.")

    tp53_status = tp53_status.set_index("sample_id").loc[shared_samples].reset_index()
    point_mut = point_mut.loc[shared_samples]
    cna_amp = cna_amp.loc[shared_samples]

    return tp53_status, point_mut, cna_amp, final_genes


def build_sample_level_mechanism_table(
    tp53_status: pd.DataFrame,
    point_mut: pd.DataFrame,
    cna_amp: pd.DataFrame,
) -> pd.DataFrame:
    sample_ids = point_mut.index.astype(str)
    group_map = dict(zip(tp53_status["sample_id"].astype(str), tp53_status["tp53_group"].astype(str)))

    any_point = (point_mut.max(axis=1) > 0).astype(np.int8)
    any_cna = (cna_amp.max(axis=1) > 0).astype(np.int8)

    out = pd.DataFrame(index=sample_ids)
    out["sample_id"] = sample_ids
    out["tp53_group"] = [group_map[s] for s in sample_ids]
    out["any_point_mutation"] = any_point.values
    out["any_cna_amplification"] = any_cna.values
    out["point_only"] = ((out["any_point_mutation"] == 1) & (out["any_cna_amplification"] == 0)).astype(np.int8)
    out["cna_only"] = ((out["any_point_mutation"] == 0) & (out["any_cna_amplification"] == 1)).astype(np.int8)
    out["both"] = ((out["any_point_mutation"] == 1) & (out["any_cna_amplification"] == 1)).astype(np.int8)
    out["neither"] = ((out["any_point_mutation"] == 0) & (out["any_cna_amplification"] == 0)).astype(np.int8)
    return out.reset_index(drop=True)


def run_primary_summary_tests(sample_level_df: pd.DataFrame) -> pd.DataFrame:
    features = [
        "any_point_mutation",
        "any_cna_amplification",
        "point_only",
        "cna_only",
        "both",
    ]
    rows: list[dict[str, object]] = []
    groups = sample_level_df["tp53_group"]

    for feature in features:
        res = run_fisher_binary_feature(sample_level_df[feature], groups, LABEL_LOF, LABEL_WT)
        res["feature"] = feature
        rows.append(res)

    cols = [
        "feature",
        "group_a",
        "group_b",
        "n_group_a_positive",
        "n_group_a_negative",
        "n_group_b_positive",
        "n_group_b_negative",
        "odds_ratio",
        "p_value",
    ]
    return pd.DataFrame(rows)[cols]


def run_per_gene_tests(
    point_mut: pd.DataFrame,
    cna_amp: pd.DataFrame,
    tp53_status: pd.DataFrame,
) -> pd.DataFrame:
    group_series = tp53_status.set_index("sample_id")["tp53_group"].reindex(point_mut.index)
    rows: list[dict[str, object]] = []

    for gene in point_mut.columns:
        for mechanism, matrix in [("point_mutation", point_mut), ("cna_amplification", cna_amp)]:
            feature = matrix[gene].astype(int)
            res = run_fisher_binary_feature(feature, group_series, LABEL_LOF, LABEL_WT)
            rows.append(
                {
                    "gene": gene,
                    "mechanism": mechanism,
                    "n_lof_positive": res["n_group_a_positive"],
                    "n_lof_negative": res["n_group_a_negative"],
                    "n_wt_positive": res["n_group_b_positive"],
                    "n_wt_negative": res["n_group_b_negative"],
                    "odds_ratio": res["odds_ratio"],
                    "p_value": res["p_value"],
                }
            )

    out = pd.DataFrame(rows)
    if len(out):
        out["q_value_bh"] = benjamini_hochberg_fdr(out["p_value"].values)
    return out.sort_values(["q_value_bh", "p_value", "gene", "mechanism"]).reset_index(drop=True)


def plot_mechanism_stacked_bar(sample_level_df: pd.DataFrame, out_path: str | Path) -> Path:
    out = _resolve_repo_path(out_path)
    out.parent.mkdir(parents=True, exist_ok=True)

    categories = ["point_only", "cna_only", "both", "neither"]
    plot_df = sample_level_df.groupby("tp53_group")[categories].mean().reindex([LABEL_LOF, LABEL_WT])

    fig, ax = plt.subplots(figsize=(8, 6))
    bottom = np.zeros(len(plot_df), dtype=float)

    for cat in categories:
        ax.bar(plot_df.index, plot_df[cat].values, bottom=bottom, label=cat)
        bottom += plot_df[cat].values

    ax.set_ylabel("Fraction of samples")
    ax.set_xlabel("TP53 group")
    ax.set_title("Structural vs point alteration patterns by TP53 group")
    ax.set_ylim(0, 1.0)
    ax.legend(title="Mechanism class")
    plt.tight_layout()
    plt.savefig(out, dpi=150, bbox_inches="tight")
    plt.close(fig)
    return out


def plot_per_gene_heatmap(per_gene_df: pd.DataFrame, out_path: str | Path) -> Path:
    out = _resolve_repo_path(out_path)
    out.parent.mkdir(parents=True, exist_ok=True)

    plot_df = per_gene_df.copy()
    signed_logp = []
    for _, r in plot_df.iterrows():
        p = max(float(r["p_value"]), 1e-300)
        orr = float(r["odds_ratio"])
        sign = 1.0 if orr > 1 else -1.0
        signed_logp.append(sign * (-np.log10(p)))
    plot_df["signed_log10_p"] = signed_logp

    pivot = plot_df.pivot(index="gene", columns="mechanism", values="signed_log10_p").fillna(0.0)

    fig, ax = plt.subplots(figsize=(6, max(8, len(pivot) * 0.35)))
    im = ax.imshow(pivot.values, aspect="auto")
    ax.set_xticks(np.arange(len(pivot.columns)))
    ax.set_xticklabels(pivot.columns, rotation=45, ha="right")
    ax.set_yticks(np.arange(len(pivot.index)))
    ax.set_yticklabels(pivot.index)
    ax.set_title("Per-gene mechanism enrichment in TP53 LoF vs WT\n(signed -log10 p)")
    cbar = fig.colorbar(im, ax=ax)
    cbar.set_label("signed -log10(p)")
    plt.tight_layout()
    plt.savefig(out, dpi=150, bbox_inches="tight")
    plt.close(fig)
    return out


def write_report(
    report_path: str | Path,
    *,
    tp53_path: Path,
    point_path: Path,
    raw_cna_path: Path,
    processed_cna_path: Path,
    final_genes: list[str],
    tp53_counts: pd.Series,
    n_samples: int,
    summary_df: pd.DataFrame,
    per_gene_df: pd.DataFrame,
    stacked_bar_path: Path,
    heatmap_path: Path,
) -> Path:
    out = _resolve_repo_path(report_path)
    out.parent.mkdir(parents=True, exist_ok=True)

    top_gene = per_gene_df.head(12)[
        ["gene", "mechanism", "n_lof_positive", "n_lof_negative", "n_wt_positive", "n_wt_negative", "odds_ratio", "p_value", "q_value_bh"]
    ]

    lines = [
        "# Step 8: Structural vs point mechanism shift",
        "",
        "## Inputs",
        f"- TP53 status: `{tp53_path}`",
        f"- Point mutation matrix: `{point_path}`",
        f"- Raw CNA matrix: `{raw_cna_path}`",
        f"- Processed CNA amplification matrix: `{processed_cna_path}`",
        "",
        "## Analysis setup",
        f"- Primary comparison: `{LABEL_LOF}` vs `{LABEL_WT}`",
        "- CNA threshold: amplification = 1 if discrete CNA call >= 2",
        f"- Number of overlapping samples: {n_samples}",
        f"- Number of final genes analyzed: {len(final_genes)}",
        f"- Final genes: {', '.join(final_genes)}",
        "",
        "## TP53 group counts",
        "",
        _df_to_markdown(tp53_counts.reset_index().rename(columns={"index": "tp53_group", "tp53_group": "n_samples"})),
        "",
        "## Primary sample-level summary tests",
        "",
        _df_to_markdown(summary_df),
        "",
        "## Top per-gene results",
        "",
        _df_to_markdown(top_gene),
        "",
        "## Outputs",
        f"- Summary table: `outputs/tables/structural_vs_point_tp53_lof_vs_wt_summary.tsv`",
        f"- Per-gene table: `outputs/tables/structural_vs_point_tp53_lof_vs_wt_per_gene.tsv`",
        f"- Sample-level mechanisms: `outputs/tables/structural_vs_point_sample_level_mechanisms.tsv`",
        f"- Stacked bar figure: `{stacked_bar_path}`",
        f"- Per-gene heatmap: `{heatmap_path}`",
        "",
    ]
    out.write_text("\n".join(lines), encoding="utf-8")
    return out


def main() -> None:
    parser = argparse.ArgumentParser(description="Alternative hypothesis: structural vs point (CNA amplification vs point mutation).")
    parser.add_argument("--maf", default=str(_DEFAULT_MAF), help="Raw mutation MAF/TSV.")
    parser.add_argument("--clinical", default=str(_DEFAULT_CLINICAL), help="Clinical TSV.")
    parser.add_argument("--tp53-status", default=str(_DEFAULT_TP53_FUNCTIONAL), help="TP53 functional groups CSV.")
    parser.add_argument("--mutation-matrix", default=str(_DEFAULT_POINT_MUT), help="Point mutation matrix parquet.")
    parser.add_argument("--raw-cna", default=str(_DEFAULT_RAW_CNA), help="Raw gene-level CNA matrix.")
    parser.add_argument("--cna-threshold", type=int, default=2, help="Amplification threshold; default = 2.")
    parser.add_argument("--out-dir", default="outputs", help="Base output directory.")
    parser.add_argument(
        "--out-cna-matrix",
        default=str(_DEFAULT_CNA_OUT),
        help="Processed binary CNA amplification matrix parquet.",
    )
    args = parser.parse_args()

    out_base = _resolve_repo_path(args.out_dir)
    out_tables = out_base / "tables"
    out_figures = out_base / "figures"
    out_reports = out_base / "reports"
    out_tables.mkdir(parents=True, exist_ok=True)
    out_figures.mkdir(parents=True, exist_ok=True)
    out_reports.mkdir(parents=True, exist_ok=True)

    tp53_path = _resolve_repo_path(args.tp53_status)
    point_path = _resolve_repo_path(args.mutation_matrix)
    raw_cna_path = _resolve_repo_path(args.raw_cna)
    clinical_path = _resolve_repo_path(args.clinical)

    print("Loading inputs...")
    tp53_status = load_tp53_status(tp53_path)
    point_mut = build_point_mutation_matrix(
        maf_path=args.maf,
        clinical_path=args.clinical,
        out_path=args.mutation_matrix,
    )
    raw_cna = load_raw_cna_matrix(
        path=raw_cna_path,
        clinical_path=clinical_path,
    )
    point_path = _resolve_repo_path(args.mutation_matrix)

    print("Processing CNA matrix...")
    cna_sample_by_gene = cna_gene_by_sample_to_sample_by_gene(raw_cna)
    cna_amp = build_cna_amplification_matrix(cna_sample_by_gene, threshold=args.cna_threshold)
    processed_cna_path = save_cna_amplification_matrix(cna_amp, args.out_cna_matrix)
    print(f"Wrote processed CNA matrix: {processed_cna_path}")

    print("Harmonizing samples and genes...")
    tp53_status, point_mut, cna_amp, final_genes = harmonize_inputs(
        tp53_status=tp53_status,
        point_mut=point_mut,
        cna_amp=cna_amp,
        gene_list=STRUCTURAL_POINT_GENE_LIST,
    )

    print(f"Overlapping samples retained: {len(point_mut)}")
    print(f"Final genes retained ({len(final_genes)}): {', '.join(final_genes)}")

    print("Building sample-level mechanism table...")
    sample_level_df = build_sample_level_mechanism_table(tp53_status, point_mut, cna_amp)

    print("Running primary summary tests...")
    summary_df = run_primary_summary_tests(sample_level_df)

    print("Running per-gene tests...")
    per_gene_df = run_per_gene_tests(point_mut, cna_amp, tp53_status)

    summary_path = out_tables / "structural_vs_point_tp53_lof_vs_wt_summary.tsv"
    per_gene_path = out_tables / "structural_vs_point_tp53_lof_vs_wt_per_gene.tsv"
    sample_level_path = out_tables / "structural_vs_point_sample_level_mechanisms.tsv"

    summary_df.to_csv(summary_path, sep="\t", index=False)
    per_gene_df.to_csv(per_gene_path, sep="\t", index=False)
    sample_level_df.to_csv(sample_level_path, sep="\t", index=False)

    print(f"Wrote {summary_path}")
    print(f"Wrote {per_gene_path}")
    print(f"Wrote {sample_level_path}")

    print("Making figures...")
    stacked_bar_path = plot_mechanism_stacked_bar(
        sample_level_df,
        out_figures / "structural_vs_point_mechanism_stacked_bar.png",
    )
    heatmap_path = plot_per_gene_heatmap(
        per_gene_df,
        out_figures / "structural_vs_point_per_gene_heatmap.png",
    )

    print(f"Wrote {stacked_bar_path}")
    print(f"Wrote {heatmap_path}")

    tp53_counts = sample_level_df["tp53_group"].value_counts()

    report_path = write_report(
        out_reports / "step8_structural_vs_point_report.md",
        tp53_path=tp53_path,
        point_path=point_path,
        raw_cna_path=raw_cna_path,
        processed_cna_path=processed_cna_path,
        final_genes=final_genes,
        tp53_counts=tp53_counts,
        n_samples=len(sample_level_df),
        summary_df=summary_df,
        per_gene_df=per_gene_df,
        stacked_bar_path=stacked_bar_path,
        heatmap_path=heatmap_path,
    )
    print(f"Wrote {report_path}")

if __name__ == "__main__":
    main()

