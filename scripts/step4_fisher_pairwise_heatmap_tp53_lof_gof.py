"""
step4_fisher_pairwise_heatmap_tp53_lof_gof.py

Pipeline — **Step 4** (LUAD; TP53 LoF vs GoF strata)

**Phase 4 — functional TP53**

1. **Pairwise Fisher** (same core as step 2): gene–gene co-occurrence **within** each stratum
   (**TP53_LoF** only, **TP53_GoF_missense** only). **TP53_WT** is excluded from these heatmaps.

2. **TP53_LoF vs TP53_GoF_missense** (TP53 **loss-of-function** vs **gain-of-function / missense proxy**)
   on the **RTK/RAS/PI3K** gene list: per-gene **prevalence in each TP53 stratum** plus **Fisher** panels
   for **coding MAF** and **CNA amplification** (binary amp, score ≥ threshold). Figures under
   ``outputs/figures/lof_vs_gof/`` — columns compare **TP53 cohorts**, not “CNA vs point” as the x-axis.

Reads the annotated matrix from **step 3**; MAF/clinical paths default to **step 1** defaults for
rebuilding the wide RTK panel (genes outside step-1 ``GENE_LIST``).

Inputs
------
- ``data/processed/lof_gof/mutation_matrix_with_tp53_group.csv`` (default)
- Optional: MAF + clinical (defaults from step 1), ``data/data_cna.txt`` (or ``--raw-cna``)

Outputs
-------
- ``outputs/figures/fisher_tp53_lof_pairwise_heatmap.png``
- ``outputs/figures/fisher_tp53_gof_missense_pairwise_heatmap.png``
- ``outputs/figures/lof_vs_gof/rtk_TP53_LoF_vs_GoF_missense_RTKpathway_enrichment.png``
- ``outputs/figures/lof_vs_gof/rtk_TP53_LoF_vs_GoF_missense_coding_Fisher_neglog10p_magnitude.png``
- ``outputs/figures/lof_vs_gof/rtk_structural_vs_point_patterns_TP53_LoF_vs_GoF_missense.png`` (stacked bar)
- ``outputs/figures/lof_vs_gof/{GENE}_TP53_LoF_vs_GoF_missense_exploration.png`` (default gene **PDPK1**; ``--focus-gene``)
- ``outputs/tables/step4_rtk_*_lof_vs_gof.tsv`` and ``step4_{GENE}_TP53_LoF_vs_GoF_contingency.tsv``
"""

from __future__ import annotations

import argparse
import sys
from itertools import combinations
from pathlib import Path

import matplotlib

matplotlib.use("Agg")
import matplotlib.pyplot as plt
import numpy as np
import pandas as pd
import seaborn as sns
from scipy.stats import fisher_exact

_SCRIPT_DIR = Path(__file__).resolve().parent
REPO_ROOT = _SCRIPT_DIR.parent
if str(_SCRIPT_DIR) not in sys.path:
    sys.path.insert(0, str(_SCRIPT_DIR))

from step1_digestion_and_processing import DEFAULT_CLINICAL_PATH, DEFAULT_MAF_PATH
from step3_classify_tp53_lof_vs_gof import (
    LABEL_GOF,
    LABEL_LOF,
    _align_cna_to_samples,
    _cna_amp_submatrix,
    _read_rtk_gene_panel,
    _write_rtk_pathway_enrichment_long_table,
    cna_gene_by_sample_to_sample_by_gene,
    fisher_amp_per_gene_pair,
    fisher_per_gene_lof_vs_gof,
    load_luad_point_matrix_for_genes,
    load_raw_cna_matrix,
    plot_rtk_pathway_cohort_enrichment_figure,
    split_by_tp53_group,
)

META_COLUMNS = ("TP53_status", "tp53_group")

_DEFAULT_MATRIX = str(
    REPO_ROOT / "data" / "processed" / "lof_gof" / "mutation_matrix_with_tp53_group.csv"
)


def plot_fisher_heatmap(
    df: pd.DataFrame,
    title: str,
    save_path: Path | None = None,
    show: bool = True,
):
    """Pairwise Fisher over gene columns (0/1) within one TP53 stratum."""
    df_genes = df.drop(columns=list(META_COLUMNS), errors="ignore")
    genes = df_genes.columns.tolist()
    n_genes = len(genes)

    if n_genes < 2:
        print(f"Skipping heatmap — need ≥2 genes, got {n_genes}. ({title})")
        return None

    p_val_matrix = pd.DataFrame(1.0, index=genes, columns=genes)
    for g1, g2 in combinations(genes, 2):
        a = ((df_genes[g1] == 1) & (df_genes[g2] == 1)).sum()
        b = ((df_genes[g1] == 1) & (df_genes[g2] == 0)).sum()
        c = ((df_genes[g1] == 0) & (df_genes[g2] == 1)).sum()
        d = ((df_genes[g1] == 0) & (df_genes[g2] == 0)).sum()

        table = [[a, b], [c, d]]
        _, pval = fisher_exact(table)
        p_val_matrix.loc[g1, g2] = pval
        p_val_matrix.loc[g2, g1] = pval

    log_p_matrix = -np.log10(p_val_matrix.astype(float) + 1e-10)
    mask = np.triu(np.ones_like(log_p_matrix, dtype=bool))

    plt.figure(figsize=(10, 8))
    sns.heatmap(
        log_p_matrix,
        mask=mask,
        cmap="Reds",
        annot=False,
        linewidths=0.5,
        cbar_kws={"label": "-log10(p-value)"},
    )
    for i in range(n_genes):
        for j in range(i + 1, n_genes):
            if p_val_matrix.iloc[j, i] <= 0.05:
                plt.text(
                    i + 0.5,
                    j + 0.5,
                    "*",
                    ha="center",
                    va="center",
                    color="black",
                    fontsize=20,
                )
    plt.title(title)
    plt.tight_layout()
    if save_path is not None:
        save_path.parent.mkdir(parents=True, exist_ok=True)
        plt.savefig(save_path, dpi=150, bbox_inches="tight")
    if show:
        plt.show()
    else:
        plt.close()

    return p_val_matrix


MECHANISM_COLORS = {
    "point_only": "#5B9BD5",
    "cna_only": "#ED7D31",
    "both": "#70AD47",
    "neither": "#C0504D",
}


def _mechanism_class_series(
    sample_ids: pd.Index,
    point_df: pd.DataFrame,
    cna_amp_df: pd.DataFrame | None,
    genes: list[str],
) -> pd.Series:
    """
    Per sample: ``point_only`` | ``cna_only`` | ``both`` | ``neither`` on RTK/RAS/PI3K ``genes``.
    Point = any coding MAF hit in panel; structural = any CNA amplification in panel.
    """
    genes_p = [g for g in genes if g in point_df.columns]
    genes_c = [g for g in genes if cna_amp_df is not None and g in cna_amp_df.columns]
    labels: list[str] = []
    for sid in sample_ids.astype(str):
        sid = str(sid)
        pt = 0
        if sid in point_df.index and genes_p:
            pt = int(point_df.loc[sid, genes_p].astype(int).clip(0, 1).max())
        cn = 0
        if cna_amp_df is not None and sid in cna_amp_df.index and genes_c:
            cn = int(cna_amp_df.loc[sid, genes_c].astype(int).clip(0, 1).max())
        if pt >= 1 and cn >= 1:
            labels.append("both")
        elif pt >= 1:
            labels.append("point_only")
        elif cn >= 1:
            labels.append("cna_only")
        else:
            labels.append("neither")
    return pd.Series(labels, index=sample_ids.astype(str), dtype=str)


def plot_structural_vs_point_stacked_lof_vs_gof(
    mech_lof: pd.Series,
    mech_gof: pd.Series,
    out_path: Path,
    show_plots: bool,
) -> None:
    """
    Stacked bar chart: fraction of samples in each mechanism class, **TP53_LoF** vs **TP53_GoF_missense**.
    Matches reference layout (point_only / cna_only / both / neither on RTK/RAS/PI3K genes).
    """
    from matplotlib.patches import Patch

    # Bottom of stack first (point_only at base, as in reference).
    stack_order = ["point_only", "cna_only", "both", "neither"]

    def fracs(series: pd.Series) -> dict[str, float]:
        n = len(series)
        if n == 0:
            return {c: 0.0 for c in stack_order}
        vc = series.value_counts()
        return {c: float(vc.get(c, 0)) / n for c in stack_order}

    f_lof = fracs(mech_lof)
    f_gof = fracs(mech_gof)

    x_labels = ["TP53_LoF", "TP53_GoF_missense"]
    x = np.arange(len(x_labels))
    width = 0.55

    fig, ax = plt.subplots(figsize=(6.5, 5.2))
    bottoms = np.zeros(2)
    for cat in stack_order:
        heights = np.array([f_lof[cat], f_gof[cat]])
        ax.bar(
            x,
            heights,
            width,
            bottom=bottoms,
            color=MECHANISM_COLORS[cat],
            edgecolor="white",
            linewidth=0.6,
        )
        bottoms += heights

    ax.set_xticks(x)
    ax.set_xticklabels(x_labels)
    ax.set_ylabel("Fraction of samples")
    ax.set_xlabel("TP53 group")
    ax.set_ylim(0.0, 1.0)
    ax.set_title("Structural vs point alteration patterns by TP53 group")
    ax.legend(
        handles=[
            Patch(facecolor=MECHANISM_COLORS[c], edgecolor="white", label=c.replace("_", " "))
            for c in stack_order
        ],
        title="Mechanism class",
        loc="upper left",
        bbox_to_anchor=(1.02, 1.0),
        borderaxespad=0.0,
    )
    fig.tight_layout()
    out_path.parent.mkdir(parents=True, exist_ok=True)
    fig.savefig(out_path, dpi=150, bbox_inches="tight")
    if show_plots:
        plt.show()
    plt.close(fig)


def plot_focus_gene_tp53_lof_vs_gof_exploration(
    gene: str,
    lof_idx: pd.Index,
    gof_idx: pd.Index,
    point_gene_df: pd.DataFrame,
    cna_amp_full: pd.DataFrame | None,
    maf_path: str,
    clinical_path: str,
    cohort_index: pd.Index,
    cna_amp_threshold: int,
    out_fig_path: Path,
    tables_dir: Path,
    show_plots: bool,
) -> None:
    """
    **TP53_LoF** vs **TP53_GoF_missense**: prevalence + Fisher + per-sample coding status for one gene
    (default **PDPK1**), plus CNA amp bars when data exist.
    """
    g = gene.strip().upper()
    if g in point_gene_df.columns:
        s_point = point_gene_df[g].reindex(cohort_index.astype(str)).fillna(0).astype(int)
    else:
        one = load_luad_point_matrix_for_genes(
            maf_path, clinical_path, [g], cohort_index
        )
        s_point = (
            one[g].reindex(cohort_index.astype(str)).fillna(0).astype(int)
            if g in one.columns
            else pd.Series(0, index=cohort_index.astype(str), dtype=int)
        )

    lof_ix = lof_idx.astype(str)
    gof_ix = gof_idx.astype(str)
    n_lof, n_gof = len(lof_ix), len(gof_ix)
    mut_lof = int(s_point.reindex(lof_ix).fillna(0).astype(int).sum())
    mut_gof = int(s_point.reindex(gof_ix).fillna(0).astype(int).sum())
    a, b = mut_lof, n_lof - mut_lof
    c, d = mut_gof, n_gof - mut_gof
    or_coding, p_coding = fisher_exact([[a, b], [c, d]])

    amp_lof = amp_gof = 0
    or_amp = float("nan")
    p_amp = float("nan")
    has_cna_gene = (
        cna_amp_full is not None
        and g in cna_amp_full.columns
    )
    if has_cna_gene and cna_amp_full is not None:
        s_amp = cna_amp_full[g].reindex(cohort_index.astype(str)).fillna(0).astype(int)
        amp_lof = int(s_amp.reindex(lof_ix).fillna(0).astype(int).sum())
        amp_gof = int(s_amp.reindex(gof_ix).fillna(0).astype(int).sum())
        or_amp, p_amp = fisher_exact(
            [[amp_lof, n_lof - amp_lof], [amp_gof, n_gof - amp_gof]]
        )

    pd.DataFrame(
        [
            {
                "gene": g,
                "layer": "coding_MAF",
                "n_lof": n_lof,
                "n_gof": n_gof,
                "mut_lof": mut_lof,
                "mut_gof": mut_gof,
                "odds_ratio": float(or_coding),
                "p_value": float(p_coding),
            },
            {
                "gene": g,
                "layer": "CNA_amp",
                "n_lof": n_lof,
                "n_gof": n_gof,
                "mut_lof": amp_lof,
                "mut_gof": amp_gof,
                "odds_ratio": float(or_amp) if has_cna_gene else np.nan,
                "p_value": float(p_amp) if has_cna_gene else np.nan,
            },
        ]
    ).to_csv(
        tables_dir / f"step4_{g}_TP53_LoF_vs_GoF_contingency.tsv",
        sep="\t",
        index=False,
    )

    fig, axes = plt.subplots(1, 3, figsize=(12.5, 4.2), gridspec_kw={"width_ratios": [1, 1, 1.15]})

    prev_df = pd.DataFrame(
        {
            "TP53 group": ["TP53_LoF", "TP53_GoF_missense"],
            "fraction": [mut_lof / n_lof if n_lof else 0.0, mut_gof / n_gof if n_gof else 0.0],
            "n_mut": [mut_lof, mut_gof],
            "n_tot": [n_lof, n_gof],
        }
    )
    sns.barplot(
        data=prev_df,
        x="TP53 group",
        y="fraction",
        hue="TP53 group",
        palette={"TP53_LoF": "#c0392b", "TP53_GoF_missense": "#27ae60"},
        ax=axes[0],
        dodge=False,
        legend=False,
    )
    axes[0].set_ylim(0, 1.05)
    axes[0].set_ylabel("Fraction with coding mutation")
    axes[0].set_title(f"{g} — coding (MAF)")
    for i, row in prev_df.iterrows():
        axes[0].text(
            i,
            row["fraction"] + 0.04,
            f'{int(row["n_mut"])}/{int(row["n_tot"])}',
            ha="center",
            fontsize=9,
        )
    axes[0].text(
        0.5,
        -0.22,
        f"Fisher OR={or_coding:.3g}, p={p_coding:.3g}",
        transform=axes[0].transAxes,
        ha="center",
        fontsize=9,
    )

    if has_cna_gene and cna_amp_full is not None:
        prev_amp = pd.DataFrame(
            {
                "TP53 group": ["TP53_LoF", "TP53_GoF_missense"],
                "fraction": [
                    amp_lof / n_lof if n_lof else 0.0,
                    amp_gof / n_gof if n_gof else 0.0,
                ],
            }
        )
        sns.barplot(
            data=prev_amp,
            x="TP53 group",
            y="fraction",
            hue="TP53 group",
            palette={"TP53_LoF": "#c0392b", "TP53_GoF_missense": "#27ae60"},
            ax=axes[1],
            dodge=False,
            legend=False,
        )
        axes[1].set_ylim(0, 1.05)
        axes[1].set_ylabel(f"Fraction amplified (CNA≥{cna_amp_threshold})")
        axes[1].set_title(f"{g} — CNA amplification")
        axes[1].text(
            0.5,
            -0.22,
            f"Fisher OR={or_amp:.3g}, p={p_amp:.3g}",
            transform=axes[1].transAxes,
            ha="center",
            fontsize=9,
        )
    else:
        axes[1].axis("off")
        axes[1].text(
            0.5,
            0.5,
            "CNA not available or\ngene absent from CNA matrix",
            ha="center",
            va="center",
            transform=axes[1].transAxes,
            fontsize=10,
        )

    ycol = f"{g} coding (0/1)"
    spf_lof = s_point.reindex(lof_ix).fillna(0).astype(int)
    spf_gof = s_point.reindex(gof_ix).fillna(0).astype(int)
    strip_df = pd.concat(
        [
            pd.DataFrame({"TP53 group": "TP53_LoF", ycol: spf_lof.values}),
            pd.DataFrame({"TP53 group": "TP53_GoF_missense", ycol: spf_gof.values}),
        ],
        ignore_index=True,
    )
    sns.stripplot(
        data=strip_df,
        x="TP53 group",
        y=ycol,
        hue="TP53 group",
        palette={"TP53_LoF": "#c0392b", "TP53_GoF_missense": "#27ae60"},
        jitter=0.22,
        dodge=False,
        ax=axes[2],
        legend=False,
        size=9,
    )
    axes[2].set_ylim(-0.15, 1.15)
    axes[2].set_ylabel(f"{g} coding mutation (binary)")
    axes[2].set_title("Per-sample calls (LUAD)")
    axes[2].set_yticks([0, 1])

    fig.suptitle(
        f"{g}: TP53 loss-of-function vs gain-of-function (missense proxy)\n"
        "LUAD — TP53_LoF vs TP53_GoF_missense (TP53_WT excluded)",
        fontsize=11,
        y=1.08,
    )
    fig.tight_layout()
    out_fig_path.parent.mkdir(parents=True, exist_ok=True)
    fig.savefig(out_fig_path, dpi=150, bbox_inches="tight")
    if show_plots:
        plt.show()
    plt.close(fig)


def _neglog10p_magnitude(p_value: float, cap: float = 3.0) -> float:
    if p_value is None or (isinstance(p_value, float) and (np.isnan(p_value) or p_value >= 1.0)):
        return 0.0
    return float(min(-np.log10(max(float(p_value), 1e-300)), cap))


def plot_rtk_coding_fisher_neglog10_magnitude_column(
    genes_ordered: list[str],
    tab_point: pd.DataFrame,
    title: str,
    out_path: Path,
    show_plots: bool,
    nl_cap: float = 3.0,
) -> None:
    """Single column: min(−log10(p), cap) for **coding MAF** Fisher TP53_LoF vs TP53_GoF_missense."""
    pp = (
        tab_point.set_index("gene")
        if len(tab_point) and "gene" in tab_point.columns
        else pd.DataFrame()
    )
    n = len(genes_ordered)
    mat = np.zeros((n, 1))
    for i, g in enumerate(genes_ordered):
        if not pp.empty and g in pp.index:
            mat[i, 0] = _neglog10p_magnitude(float(pp.loc[g, "p_value"]), nl_cap)
    vmax = max(0.15, float(np.nanmax(mat)) * 1.05) if mat.size else nl_cap
    vmax = min(nl_cap, vmax)
    fig_h = max(10.0, 0.32 * len(genes_ordered))
    fig, ax = plt.subplots(figsize=(5.0, fig_h))
    sns.heatmap(
        mat,
        yticklabels=genes_ordered,
        xticklabels=["Coding MAF\nFisher\nTP53_LoF vs\nTP53_GoF_missense"],
        cmap="viridis",
        vmin=0.0,
        vmax=vmax,
        cbar_kws={"label": "min(−log10(p), cap)"},
        ax=ax,
        linewidths=0.15,
        linecolor="lightgray",
    )
    ax.set_title(title + f"\n(uncorrected p; cap={nl_cap})")
    plt.setp(ax.get_xticklabels(), rotation=25, ha="right")
    fig.tight_layout()
    fig.savefig(out_path, dpi=150, bbox_inches="tight")
    if show_plots:
        plt.show()
    plt.close(fig)


def run_rtk_lof_vs_gof_pathway_enrichment(
    mm: pd.DataFrame,
    maf_path: str,
    clinical_path: str,
    raw_cna_path: Path | None,
    cna_path_attempted: Path | None,
    cna_amp_threshold: int,
    rtk_gene_list_csv: Path | None,
    out_fig_dir: Path,
    tables_dir: Path,
    show_plots: bool,
    focus_gene: str = "PDPK1",
) -> None:
    """
    TP53_LoF vs TP53_GoF_missense within **LUAD** (matrix rows only); **TP53_WT excluded**.
    RTK/RAS/PI3K genes; MAF-derived binary point mutations + CNA binary amplification Fisher tests.
    """
    out_fig_dir.mkdir(parents=True, exist_ok=True)
    tables_dir.mkdir(parents=True, exist_ok=True)

    genes_rtk = _read_rtk_gene_panel(rtk_gene_list_csv)
    if not genes_rtk:
        print("[step4] RTK gene list is empty; skipping LoF vs GoF pathway heatmaps.", file=sys.stderr)
        return

    point_panel = load_luad_point_matrix_for_genes(
        maf_path, clinical_path, genes_rtk, mm.index
    )
    point_panel = point_panel.copy()
    point_panel["tp53_group"] = mm["tp53_group"].values

    wt, lof, gof = split_by_tp53_group(point_panel)
    lof_g = lof.drop(columns=["tp53_group"], errors="ignore")
    gof_g = gof.drop(columns=["tp53_group"], errors="ignore")

    if len(lof_g) == 0 or len(gof_g) == 0:
        print(
            "[step4] Need both TP53_LoF and TP53_GoF_missense LUAD samples for RTK LoF vs GoF "
            f"(n_LoF={len(lof_g)}, n_GoF={len(gof_g)}); skipping pathway heatmaps.",
            file=sys.stderr,
        )
        return

    print(
        "[step4] RTK/RAS/PI3K: TP53_LoF vs TP53_GoF_missense within LUAD — "
        f"n_LoF={len(lof_g)}, n_GoF_missense={len(gof_g)}; {len(genes_rtk)} genes; TP53_WT excluded.",
        file=sys.stderr,
    )

    tab_point = fisher_per_gene_lof_vs_gof(lof_g, gof_g, genes_rtk)
    tab_point = tab_point.copy()
    tab_point["cohort"] = "LUAD"
    tab_point["comparison"] = "TP53_LoF_vs_TP53_GoF_missense"
    tab_point["pathway_scope"] = "RTK_RAS_PI3K_gene_list"
    tab_point.to_csv(
        tables_dir / "step4_rtk_point_fisher_lof_vs_gof.tsv",
        sep="\t",
        index=False,
    )

    tab_cna: pd.DataFrame | None = None
    amp_full_all: pd.DataFrame | None = None
    if raw_cna_path is not None and raw_cna_path.is_file():
        try:
            cna_raw = load_raw_cna_matrix(raw_cna_path)
            cna_sg = cna_gene_by_sample_to_sample_by_gene(cna_raw)
            amp_full_all = _cna_amp_submatrix(cna_sg, genes_rtk, cna_amp_threshold)
            amp_full_all.index = amp_full_all.index.astype(str)
            amp_full_all = (
                amp_full_all.reindex(mm.index.astype(str)).fillna(0).astype(int)
            )
            cna_lof = amp_full_all.reindex(lof.index.astype(str)).fillna(0).astype(int)
            cna_gof = amp_full_all.reindex(gof.index.astype(str)).fillna(0).astype(int)
            tab_cna = fisher_amp_per_gene_pair(cna_lof, cna_gof, genes_rtk)
            tab_cna = tab_cna.rename(
                columns={
                    "n_left": "n_lof",
                    "n_right": "n_gof",
                    "alt_left": "amp_lof",
                    "alt_right": "amp_gof",
                    "prev_left": "prev_lof",
                    "prev_right": "prev_gof",
                }
            )
            tab_cna = tab_cna.copy()
            tab_cna["cohort"] = "LUAD"
            tab_cna["comparison"] = "TP53_LoF_vs_TP53_GoF_missense"
            tab_cna["pathway_scope"] = "RTK_RAS_PI3K_gene_list"
            tab_cna.to_csv(
                tables_dir / "step4_rtk_cna_amp_fisher_lof_vs_gof.tsv",
                sep="\t",
                index=False,
            )
        except OSError as e:
            print(f"[step4] Could not read CNA from {raw_cna_path}: {e}", file=sys.stderr)
            amp_full_all = None
        except ValueError as e:
            print(f"[step4] CNA parse error ({raw_cna_path}): {e}", file=sys.stderr)
            amp_full_all = None
    else:
        hint = f" Path tried: {cna_path_attempted}" if cna_path_attempted else ""
        print(
            "[step4] CNA file not found or not readable; enrichment figure omits CNA panels."
            + hint,
            file=sys.stderr,
        )

    n_lof, n_gof = len(lof_g), len(gof_g)
    sub = (
        f"LUAD — RTK/RAS/PI3K (n_TP53_LoF={n_lof}, n_TP53_GoF_missense={n_gof}; TP53_WT excluded)\n"
        "TP53 loss-of-function vs TP53 gain-of-function (missense proxy)\n"
        "Per RTK gene: cohort prevalence + Fisher (coding MAF; CNA if available)"
    )
    plot_rtk_pathway_cohort_enrichment_figure(
        genes_rtk,
        tab_point,
        tab_cna,
        title=sub,
        out_path=out_fig_dir / "rtk_TP53_LoF_vs_GoF_missense_RTKpathway_enrichment.png",
        show_plots=show_plots,
        prev_a_col="prev_lof",
        prev_b_col="prev_gof",
        xlabel_a="TP53_LoF",
        xlabel_b="TP53_GoF_missense",
    )
    plot_rtk_coding_fisher_neglog10_magnitude_column(
        genes_rtk,
        tab_point,
        title=sub,
        out_path=out_fig_dir
        / "rtk_TP53_LoF_vs_GoF_missense_coding_Fisher_neglog10p_magnitude.png",
        show_plots=show_plots,
    )
    _write_rtk_pathway_enrichment_long_table(
        genes_rtk,
        tab_point,
        tab_cna,
        tables_dir / "step4_rtk_pathway_enrichment_scores_lof_vs_gof.tsv",
        "TP53_LoF_vs_TP53_GoF_missense",
        "TP53_LoF",
        "TP53_GoF_missense",
        prev_a_col="prev_lof",
        prev_b_col="prev_gof",
    )

    pt_genes = point_panel.drop(columns=["tp53_group"], errors="ignore")
    mech_lof = _mechanism_class_series(lof.index, pt_genes, amp_full_all, genes_rtk)
    mech_gof = _mechanism_class_series(gof.index, pt_genes, amp_full_all, genes_rtk)
    mech_rows: list[dict] = []
    stack_order = ["point_only", "cna_only", "both", "neither"]
    for grp, ser in (
        ("TP53_LoF", mech_lof),
        ("TP53_GoF_missense", mech_gof),
    ):
        n = len(ser)
        for m in stack_order:
            c = int((ser == m).sum())
            mech_rows.append(
                {
                    "tp53_group": grp,
                    "mechanism_class": m,
                    "n": c,
                    "fraction": (c / n) if n else 0.0,
                }
            )
    pd.DataFrame(mech_rows).to_csv(
        tables_dir / "step4_rtk_mechanism_class_counts_lof_vs_gof.tsv",
        sep="\t",
        index=False,
    )
    plot_structural_vs_point_stacked_lof_vs_gof(
        mech_lof,
        mech_gof,
        out_fig_dir / "rtk_structural_vs_point_patterns_TP53_LoF_vs_GoF_missense.png",
        show_plots,
    )

    fg = (focus_gene or "").strip()
    if fg:
        plot_focus_gene_tp53_lof_vs_gof_exploration(
            gene=fg,
            lof_idx=lof.index,
            gof_idx=gof.index,
            point_gene_df=pt_genes,
            cna_amp_full=amp_full_all,
            maf_path=maf_path,
            clinical_path=clinical_path,
            cohort_index=mm.index,
            cna_amp_threshold=cna_amp_threshold,
            out_fig_path=out_fig_dir / f"{fg.upper()}_TP53_LoF_vs_GoF_missense_exploration.png",
            tables_dir=tables_dir,
            show_plots=show_plots,
        )

    print(f"[step4] LoF vs GoF RTK pathway heatmaps written under {out_fig_dir}")


def main() -> None:
    parser = argparse.ArgumentParser(
        description="Step 4: TP53 LoF / GoF pairwise heatmaps + RTK LoF vs GoF per-gene enrichment (LUAD).",
    )
    parser.add_argument(
        "--matrix",
        default=_DEFAULT_MATRIX,
        help="Step 3 CSV: mutation matrix + tp53_group.",
    )
    parser.add_argument(
        "--out-dir",
        default=str(REPO_ROOT / "outputs" / "figures"),
        help="Directory for stratum pairwise Fisher PNGs.",
    )
    parser.add_argument(
        "--lof-vs-gof-rtk-dir",
        default=str(REPO_ROOT / "outputs" / "figures" / "lof_vs_gof"),
        help="Directory for TP53_LoF vs TP53_GoF RTK/CNA enrichment heatmaps.",
    )
    parser.add_argument(
        "--maf",
        default=None,
        help="MAF path for wide RTK gene panel (default: step 1 default MAF).",
    )
    parser.add_argument(
        "--clinical",
        default=None,
        help="Clinical TSV path (default: step 1 default).",
    )
    parser.add_argument(
        "--raw-cna",
        default=str(REPO_ROOT / "data" / "data_cna.txt"),
        help="CNA matrix (genes × samples). If missing, CNA panels are omitted from the enrichment figure.",
    )
    parser.add_argument(
        "--cna-amp-threshold",
        type=int,
        default=2,
        help="Binary amplification if discrete CNA score >= this value.",
    )
    parser.add_argument(
        "--rtk-gene-list",
        default=str(REPO_ROOT / "config" / "pathways" / "rtk_ras_pi3k_genes.csv"),
        help="CSV with ``gene`` column for RTK/RAS/PI3K panel.",
    )
    parser.add_argument(
        "--skip-lof-vs-gof-rtk",
        action="store_true",
        help="Skip TP53_LoF vs TP53_GoF RTK/CNA enrichment heatmaps.",
    )
    parser.add_argument(
        "--focus-gene",
        default="PDPK1",
        help=(
            "Gene symbol for the single-gene TP53_LoF vs GoF exploration figure "
            "(coding prevalence, CNA amp if available, per-sample strip). "
            "Empty string skips."
        ),
    )
    parser.add_argument(
        "--no-show",
        action="store_true",
        help="Do not open interactive plot windows (save PNG only).",
    )
    args = parser.parse_args()

    path = Path(args.matrix).expanduser()
    if not path.is_absolute():
        path = REPO_ROOT / path
    if not path.exists():
        raise FileNotFoundError(
            f"Matrix not found: {path}. Run step3_classify_tp53_lof_vs_gof.py first.",
        )

    mm = pd.read_csv(path, index_col="sample_id")
    if "tp53_group" not in mm.columns:
        raise ValueError(f"Expected column 'tp53_group' in {path}")

    tp53_lof = mm[mm["tp53_group"] == LABEL_LOF].copy()
    tp53_gof = mm[mm["tp53_group"] == LABEL_GOF].copy()

    out_dir = Path(args.out_dir)
    if not out_dir.is_absolute():
        out_dir = REPO_ROOT / out_dir

    print(f"TP53 LoF samples: {len(tp53_lof)} | TP53 GoF/missense samples: {len(tp53_gof)}")

    if len(tp53_lof) == 0:
        print("No TP53_LoF samples; skipping LoF pairwise heatmap.")
    else:
        print("Generating heatmap for TP53 LoF cohort...")
        plot_fisher_heatmap(
            tp53_lof,
            title="Gene co-occurrence (TP53 loss-of-function)",
            save_path=out_dir / "fisher_tp53_lof_pairwise_heatmap.png",
            show=not args.no_show,
        )

    if len(tp53_gof) == 0:
        print("No TP53_GoF_missense samples; skipping GoF pairwise heatmap.")
    else:
        print("Generating heatmap for TP53 GoF / missense cohort...")
        plot_fisher_heatmap(
            tp53_gof,
            title="Gene co-occurrence (TP53 GoF / missense)",
            save_path=out_dir / "fisher_tp53_gof_missense_pairwise_heatmap.png",
            show=not args.no_show,
        )

    if not args.skip_lof_vs_gof_rtk:
        rtk_csv = Path(args.rtk_gene_list)
        rtk_csv = rtk_csv if rtk_csv.is_file() else None
        maf_res = str(Path(args.maf).expanduser()) if args.maf else DEFAULT_MAF_PATH
        clin_res = (
            str(Path(args.clinical).expanduser())
            if args.clinical
            else DEFAULT_CLINICAL_PATH
        )
        raw_cna_arg = (args.raw_cna or "").strip()
        cna_attempted = Path(raw_cna_arg).expanduser() if raw_cna_arg else None
        if cna_attempted is not None and not cna_attempted.is_absolute():
            cna_attempted = REPO_ROOT / cna_attempted
        raw_cna_path = cna_attempted if cna_attempted and cna_attempted.is_file() else None
        lof_gof_dir = Path(args.lof_vs_gof_rtk_dir)
        if not lof_gof_dir.is_absolute():
            lof_gof_dir = REPO_ROOT / lof_gof_dir
        run_rtk_lof_vs_gof_pathway_enrichment(
            mm,
            maf_res,
            clin_res,
            raw_cna_path=raw_cna_path,
            cna_path_attempted=cna_attempted,
            cna_amp_threshold=int(args.cna_amp_threshold),
            rtk_gene_list_csv=rtk_csv,
            out_fig_dir=lof_gof_dir,
            tables_dir=REPO_ROOT / "outputs" / "tables",
            show_plots=not args.no_show,
            focus_gene=str(args.focus_gene),
        )


if __name__ == "__main__":
    main()
