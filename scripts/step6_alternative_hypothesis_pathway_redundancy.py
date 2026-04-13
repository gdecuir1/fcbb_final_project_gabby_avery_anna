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
- Gene × sample matrix from earlier steps, resolved in this order if ``--mutation-matrix`` is omitted:
  1. ``data/processed/lof_gof/mutation_matrix_with_tp53_group.csv`` (step 3)
  2. ``data/processed/mutation_matrix_with_tp53_group.csv``
  3. ``outputs/tables/step5_pancancer_sample_gene_matrix.tsv`` (step 5 pan-cancer)
  4. ``data/processed/gene_mutation_binarized_matrix.parquet`` (step 3, needs pyarrow)
- ``data/processed/lof_gof/tp53_functional_status.csv`` (step 3; default for ``--tp53-status``)

Outputs
--------
- `outputs/tables/pathway_redundancy_tp53_lof_vs_wt.tsv` — merged Fisher results (LoF vs WT cohorts)
- `outputs/figures/pathway_redundancy_heatmap.png` — comparative heatmaps
- `outputs/reports/step6_pathway_redundancy_report.md` — short run summary

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
import argparse
from itertools import combinations
from pathlib import Path

import matplotlib.pyplot as plt
import numpy as np
import pandas as pd
import seaborn as sns
from scipy.stats import fisher_exact

_REPO = Path(__file__).resolve().parent.parent
_DEFAULT_TP53_FUNCTIONAL = str(
    _REPO / "data" / "processed" / "lof_gof" / "tp53_functional_status.csv"
)


def _resolve_repo_path(p: str | Path) -> Path:
    path = Path(p).expanduser()
    return path if path.is_absolute() else _REPO / path


def _df_to_markdown_table(df: pd.DataFrame) -> str:
    cols = list(df.columns)
    lines = ["| " + " | ".join(cols) + " |", "| " + " | ".join("---" for _ in cols) + " |"]
    for _, row in df.iterrows():
        cells = [f"{row[c]:.4g}" if isinstance(row[c], (float, np.floating)) else str(row[c]) for c in cols]
        lines.append("| " + " | ".join(cells) + " |")
    return "\n".join(lines)

def benjamini_hochberg_fdr(p_values: np.ndarray) -> np.ndarray:
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

def load_mutation_matrix(path: str | Path) -> pd.DataFrame:
    p = Path(path).expanduser()
    suf = p.suffix.lower()
    if suf == ".parquet":
        return pd.read_parquet(p)
    if suf in (".tsv", ".txt") or p.name.lower().endswith(".tsv.gz"):
        return pd.read_csv(p, sep="\t", index_col=0, low_memory=False)
    if suf == ".csv":
        return pd.read_csv(p, index_col=0, low_memory=False)
    if suf == ".gz" and "tsv" in p.name.lower():
        return pd.read_csv(p, sep="\t", index_col=0, low_memory=False)
    raise ValueError(
        f"Unsupported mutation matrix format: {p} "
        "(use .tsv, .csv, or .parquet; .parquet requires pyarrow or fastparquet)."
    )


_METADATA_COLS = frozenset({"TP53_status", "tp53_group"})


def gene_only_mutation_matrix(df: pd.DataFrame) -> pd.DataFrame:
    """Drop step-3 annotation columns; keep binary gene columns only."""
    drop = [c for c in df.columns if c in _METADATA_COLS]
    out = df.drop(columns=drop, errors="ignore")
    out.index = out.index.astype(str)
    return out


_MUTATION_MATRIX_CANDIDATES: tuple[Path, ...] = (
    _REPO / "data" / "processed" / "lof_gof" / "mutation_matrix_with_tp53_group.csv",
    _REPO / "data" / "processed" / "mutation_matrix_with_tp53_group.csv",
    _REPO / "outputs" / "tables" / "step5_pancancer_sample_gene_matrix.tsv",
    _REPO / "data" / "processed" / "gene_mutation_binarized_matrix.parquet",
)


def resolve_mutation_matrix_path(explicit: str | None) -> Path:
    """Prefer explicit path; otherwise first existing candidate under the repo."""
    if explicit:
        p = _resolve_repo_path(explicit)
        if not p.is_file():
            raise FileNotFoundError(
                f"Mutation matrix not found: {p}\n"
                "Omit --mutation-matrix to auto-detect from data/processed (step 3) or outputs (step 5)."
            )
        return p
    for c in _MUTATION_MATRIX_CANDIDATES:
        if c.is_file():
            print(f"Using mutation matrix: {c.relative_to(_REPO)}")
            return c
    tried = "\n".join(f"  - {c.relative_to(_REPO)}" for c in _MUTATION_MATRIX_CANDIDATES)
    raise FileNotFoundError(
        "No mutation matrix found. Run step 1–3 (LUAD pipeline) or step 5 (pan-cancer), or pass --mutation-matrix.\n"
        "Searched:\n"
        f"{tried}"
    )


def aggregate_to_pathways(mut_matrix: pd.DataFrame, pathway_map: pd.DataFrame) -> pd.DataFrame:
    pathway_dict = pathway_map.groupby('pathway')['gene'].apply(list).to_dict()
    pathway_matrix = pd.DataFrame(index=mut_matrix.index)

    for pathway, genes in pathway_dict.items():
        valid_genes = [g for g in genes if g in mut_matrix.columns]
        if valid_genes:
            # mutated if any valid gene in pathway is mutated
            pathway_matrix[pathway] = mut_matrix[valid_genes].max(axis=1)
    return pathway_matrix.astype(np.int8)

def compute_pairwise_fisher(matrix: pd.DataFrame) -> pd.DataFrame:
    pathways = matrix.columns.tolist()
    rows = []

    for p1, p2 in combinations(pathways, 2):
        a = int(((matrix[p1] == 1) & (matrix[p2] == 1)).sum())
        b = int(((matrix[p1] == 1) & (matrix[p2] == 0)).sum())
        c = int(((matrix[p1] == 0) & (matrix[p2] == 1)).sum())
        d = int(((matrix[p1] == 0) & (matrix[p2] == 0)).sum())
        table = [[a, b], [c, d]]
        oddsratio, pval = fisher_exact(table)

        rows.append({
            "pathway_a": p1,
            "pathway_b": p2,
            "n_both": a,
            "p1_only": b,
            "p2_only": c,
            "neither": d,
            "odds_ratio": oddsratio,
            "p_value": pval
        })
    out = pd.DataFrame(rows)
    if len(out) > 0:
        out["q_value_bh"] = benjamini_hochberg_fdr(out["p_value"].values)
    return out

def plot_comparative_heatmap(lof_df: pd.DataFrame, wt_df: pd.DataFrame, out_path: Path):
    out_path.parent.mkdir(parents=True, exist_ok=True)
    pathways = list(set(lof_df['pathway_a']).union(set(lof_df['pathway_b'])))
    pathways.sort()

    def build_matrix(df):
        mat = pd.DataFrame(np.nan, index=pathways, columns=pathways)
        for _, r in df.iterrows():
            p1, p2 = r["pathway_a"], r["pathway_b"]
            sign = 1 if r["odds_ratio"] > 1 else -1
            val = sign * -np.log10(max(float(r["p_value"]), 1e-300))
            mat.loc[p1, p2] = val
            mat.loc[p2, p1] = val
        return mat
        
    lof_mat = build_matrix(lof_df)
    wt_mat = build_matrix(wt_df)

    fig, axes = plt.subplots(1, 2, figsize=(20, 8), sharey=True)
    mask = np.triu(np.ones_like(lof_mat.values, dtype=bool))

    # Calculate vmax, but enforce a minimum of 1.0 so the colorbar doesn't collapse to 0
    vmax_calc = max(np.nanmax(np.abs(lof_mat.values)), np.nanmax(np.abs(wt_mat.values)))
    vmax = max(vmax_calc, 1.0) 

    sns.heatmap(lof_mat, mask=mask, cmap="RdBu_r", center=0, vmin=-vmax, vmax=vmax,
                ax=axes[0], cbar_kws={"label": "Signed -log10(p)"})
    axes[0].set_title("TP53 LoF: Pathway Interactions")
    
    sns.heatmap(wt_mat, mask=mask, cmap="RdBu_r", center=0, vmin=-vmax, vmax=vmax,
                ax=axes[1], cbar_kws={"label": "Signed -log10(p)"})
    axes[1].set_title("TP53 WT: Pathway Interactions") # FIXED: This was mistakenly axes[0]

    plt.tight_layout()
    plt.savefig(out_path, dpi=150, bbox_inches="tight")
    plt.close()

def main() -> None:
    parser = argparse.ArgumentParser(description="Alternative hypothesis: pathway redundancy (TP53 LoF vs WT).")
    parser.add_argument("--pathway-map", default="config/pathways/pathway_to_genes.csv", help="Pathway to genes CSV.")
    parser.add_argument(
        "--tp53-status",
        default=_DEFAULT_TP53_FUNCTIONAL,
        help="TP53 functional groups (step 3 CSV under processed/lof_gof/).",
    )
    parser.add_argument(
        "--mutation-matrix",
        default=None,
        help=(
            "Gene-level binarized mutations (samples × genes). If omitted, uses the first file found: "
            "data/processed/lof_gof/mutation_matrix_with_tp53_group.csv (step 3), "
            "data/processed/mutation_matrix_with_tp53_group.csv, "
            "outputs/tables/step5_pancancer_sample_gene_matrix.tsv, "
            "or data/processed/gene_mutation_binarized_matrix.parquet."
        ),
    )
    parser.add_argument(
        "--out-dir",
        default="outputs",
        help="Repo-relative base directory; writes under <out-dir>/tables, /figures, /reports.",
    )
    args = parser.parse_args()

    out_base = _resolve_repo_path(args.out_dir)
    out_tables = out_base / "tables"
    out_figures = out_base / "figures"
    out_reports = out_base / "reports"
    out_tables.mkdir(parents=True, exist_ok=True)
    out_figures.mkdir(parents=True, exist_ok=True)
    out_reports.mkdir(parents=True, exist_ok=True)

    mut_path = resolve_mutation_matrix_path(args.mutation_matrix)
    pathway_path = _resolve_repo_path(args.pathway_map)
    tp53_path = _resolve_repo_path(args.tp53_status)

    print("Loading data...")
    mut_matrix = gene_only_mutation_matrix(load_mutation_matrix(mut_path))
    pathway_map = pd.read_csv(pathway_path, comment="#")
    pathway_map = pathway_map.dropna(subset=["pathway", "gene"]).astype({"pathway": str, "gene": str})
    tp53_status = pd.read_csv(tp53_path)

    sample_col = "Tumor_Sample_Barcode" if "Tumor_Sample_Barcode" in tp53_status.columns else tp53_status.columns[0]
    tp53_status = tp53_status.set_index(sample_col)
    tp53_status.index = tp53_status.index.astype(str)

    print("Aggregating genes to pathway level...")
    pathway_matrix = aggregate_to_pathways(mut_matrix, pathway_map)
    print(f"Collapsed {mut_matrix.shape[1]} genes into {pathway_matrix.shape[1]} pathways.")

    print("Splitting cohorts by TP53 status...")
    lof_samples = tp53_status[tp53_status['tp53_group'] == 'TP53_LoF'].index.intersection(pathway_matrix.index)
    wt_samples = tp53_status[tp53_status['tp53_group'] == 'TP53_WT'].index.intersection(pathway_matrix.index)
    
    pathway_matrix_lof = pathway_matrix.loc[lof_samples]
    pathway_matrix_wt = pathway_matrix.loc[wt_samples]
    
    print(f"TP53 LoF samples: {len(pathway_matrix_lof)}")
    print(f"TP53 WT samples: {len(pathway_matrix_wt)}")

    print("Running pairwise Fisher exact tests...")
    lof_results = compute_pairwise_fisher(pathway_matrix_lof)
    wt_results = compute_pairwise_fisher(pathway_matrix_wt)
    
    # Merge for easy comparison
    merged = pd.merge(
        lof_results, wt_results, 
        on=["pathway_a", "pathway_b"], 
        suffixes=("_lof", "_wt")
    )
    
    # Calculate the difference in odds ratios
    merged["odds_ratio_diff"] = merged["odds_ratio_lof"] - merged["odds_ratio_wt"]
    merged = merged.sort_values("odds_ratio_diff", ascending=False)

    out_file = out_tables / "pathway_redundancy_tp53_lof_vs_wt.tsv"
    merged.to_csv(out_file, sep="\t", index=False)
    print(f"Saved tabular results to {out_file}")

    print("Generating comparative heatmap...")
    heatmap_file = out_figures / "pathway_redundancy_heatmap.png"
    plot_comparative_heatmap(lof_results, wt_results, heatmap_file)
    print(f"Saved figure to {heatmap_file}")

    report_path = out_reports / "step6_pathway_redundancy_report.md"
    top = merged.head(10)[
        ["pathway_a", "pathway_b", "odds_ratio_diff", "odds_ratio_lof", "odds_ratio_wt", "p_value_lof", "p_value_wt"]
    ]
    report_lines = [
        "# Step 6: Pathway redundancy (TP53 LoF vs WT)",
        "",
        "## Inputs",
        f"- Mutation matrix: `{mut_path}`",
        f"- Pathway map: `{pathway_path}`",
        f"- TP53 status: `{tp53_path}`",
        "",
        "## Cohort sizes",
        f"- TP53 LoF samples (with pathway data): {len(pathway_matrix_lof)}",
        f"- TP53 WT samples (with pathway data): {len(pathway_matrix_wt)}",
        "",
        "## Outputs",
        f"- Table: `{out_file}`",
        f"- Figure: `{heatmap_file}`",
        "",
        "## Top pathway pairs by `odds_ratio_diff` (LoF OR − WT OR)",
        "",
        _df_to_markdown_table(top),
        "",
    ]
    report_path.write_text("\n".join(report_lines), encoding="utf-8")
    print(f"Saved report to {report_path}")

if __name__ == "__main__":
    main()

