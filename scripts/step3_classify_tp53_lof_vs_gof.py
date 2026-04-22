"""
step3_classify_tp53_lof_vs_gof.py

Pipeline — **Step 3** (LUAD; TP53 functional stratification)

**Phase 3:** assign each sample a **TP53 functional group** for stratified analyses (steps 4, 6, 7).

Uses the same MAF filters and LUAD cohort as **step 1** (coding variants in ``GENE_LIST`` only).

Groups
------
- **TP53_WT:** no coding TP53 row in the merged long table (binary TP53 column == 0).
- **TP53_LoF (loss-of-function):** at least one TP53 variant in ``LOF_VARIANT_CLASSES``:
  ``Nonsense_Mutation``, ``Frame_Shift_Ins`` / ``Frame_Shift_Del``, ``Splice_Site``.
- **TP53_GoF_missense (gain-of-function, missense proxy):** at least one ``Missense_Mutation`` at TP53
  and **no** LoF-class variants on that sample (if both LoF and missense are present, **TP53_LoF** wins).
- **TP53_other_coding:** other coding TP53 classes not mapped to LoF or missense.

**Rule:** if a sample has both LoF and missense TP53 calls, label **TP53_LoF** (LoF wins).

Inputs
------
- Same MAF + clinical paths as ``step1_digestion_and_processing.build_mutation_matrix``, or CLI overrides.

Outputs (default: ``data/processed/lof_gof/``)
----------------------------------------------
Keeps step 1–2 paths unchanged by isolating new artifacts here.

- ``tp53_functional_status.csv`` — ``sample_id``, ``tp53_group``
- ``mutation_matrix_with_tp53_group.csv`` — gene binary matrix + ``TP53_status`` + ``tp53_group``
- ``../gene_mutation_binarized_matrix.parquet`` — gene-only matrix (parquet) for optional tools / step 7
- ``outputs/figures/lof_vs_gof/rtk_TP53_LoF_vs_WT_RTKpathway_enrichment.png`` — RTK panel by **TP53 LoF vs WT**
  cohorts (prevalence + Fisher; see folder ``README.md``). **TP53_LoF vs GoF_missense** figure from **step 4**.
- ``outputs/tables/step3_rtk_*.tsv`` — Fisher summaries backing the figures.

Downstream
----------
**Step 4** reads the matrix CSV for LoF vs GoF heatmaps. **Steps 6–7** use ``tp53_functional_status.csv``
(and matrices with aligned sample IDs). For Fisher on genes only, drop metadata columns like step 2.
"""

from __future__ import annotations

import argparse
import sys
from pathlib import Path

import matplotlib

matplotlib.use("Agg")
import matplotlib.pyplot as plt
import numpy as np
import pandas as pd
import seaborn as sns
from scipy.stats import fisher_exact

# Repo layout + import path (same pattern as step 2).
_SCRIPT_DIR = Path(__file__).resolve().parent
REPO_ROOT = _SCRIPT_DIR.parent
# Step-3-only outputs live here so other `data/processed/` paths stay untouched.
LOF_GOF_DIR = REPO_ROOT / "data" / "processed" / "lof_gof"
if str(_SCRIPT_DIR) not in sys.path:
    sys.path.insert(0, str(_SCRIPT_DIR))

from step1_digestion_and_processing import (
    CODING_VARIANT_CLASSES,
    DEFAULT_CLINICAL_PATH,
    DEFAULT_MAF_PATH,
    build_mutation_matrix,
)


def benjamini_hochberg_fdr(p_values: np.ndarray) -> np.ndarray:
    """Benjamini–Hochberg FDR; return q-values in original order (same logic as step 5)."""
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

# MAF Variant_Classification strings treated as TP53 loss-of-function (aligned with step 1 filters).
LOF_VARIANT_CLASSES = frozenset(
    {
        "Nonsense_Mutation",
        "Frame_Shift_Ins",
        "Frame_Shift_Del",
        "Splice_Site",
    }
)
# Missense used as GoF / dominant-negative proxy for this pipeline (see module docstring).
MISSENSE_VARIANT_CLASSES = frozenset({"Missense_Mutation"})

TP53_GENE = "TP53"

# String labels written to CSV and used by step 4 for filtering.
LABEL_WT = "TP53_WT"
LABEL_LOF = "TP53_LoF"
LABEL_GOF = "TP53_GoF_missense"


def _sample_tp53_group(subdf: pd.DataFrame) -> str:
    """
    Assign TP53 functional group from all TP53 rows for one sample.

    **LoF:** any LoF class (nonsense / frameshift / splice). **GoF (missense proxy):** missense present
    and no LoF. If both, **LoF** is returned.
    """
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


# Default RTK/RAS/PI3K panel (Hugo symbols; alphabetical). Used if gene list CSV is missing/empty.
RTK_RAS_PI3K_FULL_GENES: tuple[str, ...] = tuple(
    sorted(
        {
            "AKT1",
            "AKT2",
            "AKT3",
            "ALK",
            "ARAF",
            "AXL",
            "BRAF",
            "EGFR",
            "ERBB2",
            "ERBB3",
            "ERBB4",
            "FGFR1",
            "FGFR2",
            "FGFR3",
            "FGFR4",
            "HRAS",
            "IGF1R",
            "KIT",
            "KRAS",
            "MAP2K1",
            "MAP2K2",
            "MAPK1",
            "MAPK3",
            "MET",
            "MTOR",
            "NRAS",
            "PDPK1",
            "PIK3CA",
            "PIK3CB",
            "PIK3CD",
            "PTPN11",
            "RAF1",
            "RET",
            "RHEB",
            "RICTOR",
            "ROS1",
            "RPTOR",
            "SOS1",
        }
    )
)


def _read_rtk_gene_panel(rtk_gene_list_csv: Path | None) -> list[str]:
    """
    Ordered gene list for RTK/RAS/PI3K figures: non-comment rows from CSV ``gene`` column,
    else ``RTK_RAS_PI3K_FULL_GENES``. Symbols uppercased; order preserved from file or alphabetical.
    """
    if rtk_gene_list_csv is not None and rtk_gene_list_csv.is_file():
        df = pd.read_csv(rtk_gene_list_csv, comment="#")
        if "gene" in df.columns:
            genes = [
                str(x).strip().upper()
                for x in df["gene"].dropna().tolist()
                if str(x).strip() and not str(x).strip().startswith("#")
            ]
            if genes:
                return genes
    return list(RTK_RAS_PI3K_FULL_GENES)


def load_luad_point_matrix_for_genes(
    maf_path: str,
    clinical_path: str,
    genes: list[str],
    cohort_index: pd.Index,
) -> pd.DataFrame:
    """
    Sample × gene binary matrix (coding MAF hits only) for ``genes`` on **LUAD** only
    (same lung adenocarcinoma clinical filter as step 1), reindexed to ``cohort_index``.
    Genes absent from MAF → all 0.
    """
    data_all = pd.read_csv(Path(maf_path).expanduser(), sep="\t")
    data_filtered = data_all[
        data_all["Variant_Classification"].isin(CODING_VARIANT_CLASSES)
    ]
    gene_set = set(genes)
    data_filtered = data_filtered[data_filtered["Hugo_Symbol"].isin(gene_set)]
    clinical_data_all = pd.read_csv(Path(clinical_path).expanduser(), sep="\t")
    clinical_la = clinical_data_all[
        (clinical_data_all["Cancer Type"] == "Lung Adenocarcinoma")
        & (clinical_data_all["Cancer Type Detailed"] == "Lung Adenocarcinoma")
    ].copy()
    if data_filtered.empty:
        out = pd.DataFrame(0, index=cohort_index, columns=genes, dtype=int)
        return out
    df = pd.merge(
        data_filtered,
        clinical_la,
        left_on="Tumor_Sample_Barcode",
        right_on="Sample ID",
        how="inner",
    )
    if df.empty:
        return pd.DataFrame(0, index=cohort_index, columns=genes, dtype=int)
    mut = df.assign(mut=1).pivot_table(
        index="Tumor_Sample_Barcode",
        columns="Hugo_Symbol",
        values="mut",
        aggfunc="max",
        fill_value=0,
    )
    idx = cohort_index.astype(str)
    mut.index = mut.index.astype(str)
    out = mut.reindex(idx).fillna(0).astype(int)
    for g in genes:
        if g not in out.columns:
            out[g] = 0
    return out[genes]


def _series_01(df: pd.DataFrame, g: str) -> pd.Series:
    if g not in df.columns:
        return pd.Series(0, index=df.index, dtype=int)
    return df[g].astype(int).clip(0, 1)


def _meta_cols() -> list[str]:
    return ["TP53_status", "tp53_group"]


def _fisher_2x2(a: int, b: int, c: int, d: int) -> tuple[float, float]:
    oddsr, p = fisher_exact([[a, b], [c, d]])
    return float(oddsr), float(p)


def signed_neglog10p(odds_ratio: float, p_value: float, cap: float = 8.0) -> float:
    """
    Directional significance from Fisher 2×2 with **LoF as the first row** (mutated vs not).

    Positive value ⇒ higher alteration rate in **TP53_LoF** than in the comparator (**TP53_WT** or
    **TP53_GoF_missense**). Negative ⇒ higher in the comparator. ``TP53_WT`` tumors are never used
    when the comparison is LoF vs GoF (those tests use only the two TP53-mutant subsets).
    """
    if p_value is None or (isinstance(p_value, float) and np.isnan(p_value)):
        return 0.0
    p = float(p_value)
    if p >= 1.0:
        return 0.0
    nl = min(-np.log10(max(p, 1e-300)), cap)
    if odds_ratio is None or (isinstance(odds_ratio, float) and np.isnan(odds_ratio)):
        return 0.0
    orf = float(odds_ratio)
    if orf > 1.0:
        return nl
    if orf < 1.0:
        return -nl
    return 0.0


def fisher_amp_per_gene_pair(
    left_amp: pd.DataFrame,
    right_amp: pd.DataFrame,
    genes: list[str],
) -> pd.DataFrame:
    """Fisher on binary amplification: first cohort vs second (same gene symbols)."""
    n_l, n_r = len(left_amp), len(right_amp)
    rows: list[dict] = []
    for g in genes:
        mut_l = int(_series_01(left_amp, g).sum())
        mut_r = int(_series_01(right_amp, g).sum())
        a, b = mut_l, n_l - mut_l
        c, d = mut_r, n_r - mut_r
        if n_l == 0 or n_r == 0:
            oratio, p = np.nan, np.nan
        else:
            oratio, p = _fisher_2x2(a, b, c, d)
        rows.append(
            {
                "gene": g,
                "n_left": n_l,
                "n_right": n_r,
                "alt_left": mut_l,
                "alt_right": mut_r,
                "prev_left": mut_l / n_l if n_l else np.nan,
                "prev_right": mut_r / n_r if n_r else np.nan,
                "odds_ratio": oratio,
                "p_value": p,
            }
        )
    out = pd.DataFrame(rows)
    if not out.empty and out["p_value"].notna().any():
        valid = out["p_value"].notna()
        out.loc[valid, "fdr_bh"] = benjamini_hochberg_fdr(
            out.loc[valid, "p_value"].astype(float).values
        )
    return out


def plot_rtk_pathway_cohort_enrichment_figure(
    genes_ordered: list[str],
    tab_point: pd.DataFrame,
    tab_cna: pd.DataFrame | None,
    title: str,
    out_path: Path,
    show_plots: bool,
    *,
    prev_a_col: str,
    prev_b_col: str,
    xlabel_a: str,
    xlabel_b: str,
    cap: float = 8.0,
) -> None:
    """
    RTK/RAS/PI3K pathway figure: **TP53 strata** (not “point vs CNA” as the main columns).

    For each RTK gene: (1) fraction with **coding mutation** in stratum A vs stratum B;
    (2) Fisher signed −log10(p) for **coding** LoF-vs-other; optionally (3–4) same for **CNA amp**
    if ``tab_cna`` is provided with matching ``prev_*`` columns.
    """
    pp = (
        tab_point.set_index("gene")
        if len(tab_point) and "gene" in tab_point.columns
        else pd.DataFrame()
    )
    has_cna = (
        tab_cna is not None
        and len(tab_cna) > 0
        and "gene" in tab_cna.columns
        and prev_a_col in tab_cna.columns
        and prev_b_col in tab_cna.columns
    )
    cc = tab_cna.set_index("gene") if has_cna else None

    n = len(genes_ordered)
    mat_prev = np.zeros((n, 2))
    mat_fisher_pt = np.zeros((n, 1))
    for i, g in enumerate(genes_ordered):
        if not pp.empty and g in pp.index:
            if prev_a_col in pp.columns:
                mat_prev[i, 0] = float(np.nan_to_num(pp.loc[g, prev_a_col], nan=0.0))
            if prev_b_col in pp.columns:
                mat_prev[i, 1] = float(np.nan_to_num(pp.loc[g, prev_b_col], nan=0.0))
            mat_fisher_pt[i, 0] = signed_neglog10p(
                float(pp.loc[g, "odds_ratio"]),
                float(pp.loc[g, "p_value"]),
                cap=cap,
            )

    mat_cna_prev: np.ndarray | None = None
    mat_fisher_cna: np.ndarray | None = None
    if has_cna and cc is not None:
        mat_cna_prev = np.zeros((n, 2))
        mat_fisher_cna = np.zeros((n, 1))
        for i, g in enumerate(genes_ordered):
            if g in cc.index:
                mat_cna_prev[i, 0] = float(np.nan_to_num(cc.loc[g, prev_a_col], nan=0.0))
                mat_cna_prev[i, 1] = float(np.nan_to_num(cc.loc[g, prev_b_col], nan=0.0))
                mat_fisher_cna[i, 0] = signed_neglog10p(
                    float(cc.loc[g, "odds_ratio"]),
                    float(cc.loc[g, "p_value"]),
                    cap=cap,
                )

    abs_f = max(
        float(np.nanmax(np.abs(mat_fisher_pt))) if mat_fisher_pt.size else 0.0,
        float(np.nanmax(np.abs(mat_fisher_cna))) if mat_fisher_cna is not None else 0.0,
    )
    if not np.isfinite(abs_f):
        abs_f = 0.0
    flim = min(float(cap), max(abs_f * 1.2, 0.35)) if abs_f > 0 else min(float(cap), 0.5)

    n_panels = 4 if has_cna else 2
    fig_h = max(12.0, 0.26 * n * n_panels + 4)
    fig, axes = plt.subplots(n_panels, 1, figsize=(7.5, fig_h))

    def _one_ax(ax, data, vmin, vmax, cmap, center, cbar_label, xlabs, row_title):
        hm_kw: dict = dict(
            data=data,
            yticklabels=genes_ordered,
            xticklabels=xlabs,
            cmap=cmap,
            vmin=vmin,
            vmax=vmax,
            cbar_kws={"label": cbar_label},
            ax=ax,
            linewidths=0.12,
            linecolor="lightgray",
        )
        if center is not None:
            hm_kw["center"] = center
        sns.heatmap(**hm_kw)
        ax.set_ylabel("")
        ax.set_title(row_title, fontsize=10)

    _one_ax(
        axes[0],
        mat_prev,
        0.0,
        1.0,
        "Blues",
        None,
        "fraction samples",
        [xlabel_a, xlabel_b],
        "Coding mutation (MAF) — prevalence by TP53 stratum",
    )
    _one_ax(
        axes[1],
        mat_fisher_pt,
        -flim,
        flim,
        "RdBu_r",
        0.0,
        "signed −log10(p)",
        [f"Fisher exact\ncoding MAF\n{xlabel_a} vs {xlabel_b}"],
        "Coding mutation — Fisher (same 2×2 as prevalence)",
    )
    if has_cna and mat_cna_prev is not None and mat_fisher_cna is not None:
        _one_ax(
            axes[2],
            mat_cna_prev,
            0.0,
            1.0,
            "Purples",
            None,
            "fraction samples",
            [xlabel_a, xlabel_b],
            "CNA amplification — prevalence by TP53 stratum",
        )
        _one_ax(
            axes[3],
            mat_fisher_cna,
            -flim,
            flim,
            "RdBu_r",
            0.0,
            "signed −log10(p)",
            [f"Fisher exact\nCNA amp\n{xlabel_a} vs {xlabel_b}"],
            "CNA amplification — Fisher",
        )

    for ax in axes:
        plt.setp(ax.get_xticklabels(), rotation=25, ha="right")
    fig.suptitle(title, fontsize=11, y=1.002)
    fig.tight_layout()
    fig.savefig(out_path, dpi=150, bbox_inches="tight")
    if show_plots:
        plt.show()
    plt.close(fig)


def fisher_per_gene_lof_vs_wt(
    lof: pd.DataFrame,
    wt: pd.DataFrame,
    genes: list[str],
) -> pd.DataFrame:
    """Per-gene Fisher: mutated (1) vs not in LoF vs WT (missing gene column ⇒ all wild-type)."""
    n_lof, n_wt = len(lof), len(wt)
    rows: list[dict] = []
    for g in genes:
        mut_lof = int(_series_01(lof, g).sum())
        mut_wt = int(_series_01(wt, g).sum())
        a, b = mut_lof, n_lof - mut_lof
        c, d = mut_wt, n_wt - mut_wt
        if n_lof == 0 or n_wt == 0:
            oratio, p = np.nan, np.nan
        else:
            oratio, p = _fisher_2x2(a, b, c, d)
        rows.append(
            {
                "gene": g,
                "n_lof": n_lof,
                "n_wt": n_wt,
                "mut_lof": mut_lof,
                "mut_wt": mut_wt,
                "prev_lof": mut_lof / n_lof if n_lof else np.nan,
                "prev_wt": mut_wt / n_wt if n_wt else np.nan,
                "odds_ratio": oratio,
                "p_value": p,
            }
        )
    out = pd.DataFrame(rows)
    if not out.empty and out["p_value"].notna().any():
        valid = out["p_value"].notna()
        fdr = benjamini_hochberg_fdr(
            out.loc[valid, "p_value"].astype(float).values
        )
        out.loc[valid, "fdr_bh"] = fdr
    return out


def fisher_per_gene_lof_vs_gof(
    lof: pd.DataFrame,
    gof: pd.DataFrame,
    genes: list[str],
) -> pd.DataFrame:
    """
    Per-gene Fisher exact on **TP53_LoF** vs **TP53_GoF_missense** tumors only (LUAD rows in ``lof`` / ``gof``).

    **TP53_WT** samples must not appear in either dataframe. Contingency per gene:
    ``[[mut_lof, non_mut_lof], [mut_gof, non_mut_gof]]`` for RTK/RAS/PI3K pathway genes (``genes``).
    """
    n_lof, n_gof = len(lof), len(gof)
    rows: list[dict] = []
    for g in genes:
        mut_lof = int(_series_01(lof, g).sum())
        mut_gof = int(_series_01(gof, g).sum())
        a, b = mut_lof, n_lof - mut_lof
        c, d = mut_gof, n_gof - mut_gof
        if n_lof == 0 or n_gof == 0:
            oratio, p = np.nan, np.nan
        else:
            oratio, p = _fisher_2x2(a, b, c, d)
        rows.append(
            {
                "gene": g,
                "n_lof": n_lof,
                "n_gof": n_gof,
                "mut_lof": mut_lof,
                "mut_gof": mut_gof,
                "prev_lof": mut_lof / n_lof if n_lof else np.nan,
                "prev_gof": mut_gof / n_gof if n_gof else np.nan,
                "odds_ratio": oratio,
                "p_value": p,
            }
        )
    out = pd.DataFrame(rows)
    if not out.empty and out["p_value"].notna().any():
        valid = out["p_value"].notna()
        fdr = benjamini_hochberg_fdr(
            out.loc[valid, "p_value"].astype(float).values
        )
        out.loc[valid, "fdr_bh"] = fdr
    return out


def _pathway_any_hit(df: pd.DataFrame, genes: list[str]) -> pd.Series:
    if not genes:
        return pd.Series(0, index=df.index, dtype=int)
    stacked = pd.concat([_series_01(df, g) for g in genes], axis=1)
    return (stacked.max(axis=1) >= 1).astype(int)


def load_raw_cna_matrix(path: Path) -> pd.DataFrame:
    """Read cBioPortal-style CNA matrix: rows = genes, columns = sample IDs."""
    return pd.read_csv(path, sep="\t", index_col=0)


def cna_gene_by_sample_to_sample_by_gene(cna: pd.DataFrame) -> pd.DataFrame:
    return cna.T


def build_cna_amplification_matrix(
    cna_sample_gene: pd.DataFrame,
    genes: list[str],
    threshold: int = 2,
) -> pd.DataFrame:
    present = [g for g in genes if g in cna_sample_gene.columns]
    if not present:
        return pd.DataFrame(index=cna_sample_gene.index)
    sub = cna_sample_gene[present].apply(pd.to_numeric, errors="coerce").fillna(0)
    return (sub >= threshold).astype(int)


def _align_cna_to_samples(
    cna_sample_gene: pd.DataFrame,
    sample_ids: pd.Index | list[str],
) -> pd.DataFrame:
    ids = [str(s) for s in sample_ids]
    idx = [s for s in ids if s in cna_sample_gene.index]
    return cna_sample_gene.reindex(idx).fillna(0)


def _cna_amp_submatrix(
    cna_sg: pd.DataFrame,
    genes: list[str],
    threshold: int,
) -> pd.DataFrame:
    """Sample × gene binary amplification for ``genes`` (missing genes → 0)."""
    amp = build_cna_amplification_matrix(cna_sg, genes, threshold=threshold)
    for g in genes:
        if g not in amp.columns:
            amp[g] = 0
    return amp.reindex(columns=genes, fill_value=0).astype(int)


def run_rtk_lof_vs_wt_gof_visualizations(
    mm_out: pd.DataFrame,
    maf_path: str,
    clinical_path: str,
    out_fig_dir: Path,
    raw_cna_path: Path | None,
    cna_amp_threshold: int,
    rtk_gene_list_csv: Path | None,
    show_plots: bool,
) -> None:
    """
    **LUAD-only** cohort (same clinical filter as step 1). **RTK/RAS/PI3K** gene panel from CSV / defaults.

    **TP53_LoF vs TP53_WT** only here (point + optional CNA signed heatmap, pathway burden, bars if ≤12 genes).
    **TP53_LoF vs TP53_GoF_missense** RTK enrichment is implemented in **step 4**.
    """
    out_fig_dir.mkdir(parents=True, exist_ok=True)
    tables_dir = REPO_ROOT / "outputs" / "tables"
    tables_dir.mkdir(parents=True, exist_ok=True)

    genes_rtk = _read_rtk_gene_panel(rtk_gene_list_csv)
    if not genes_rtk:
        print("[step3] RTK gene list is empty; skipping lof_vs_gof figures.", file=sys.stderr)
        return

    point_panel = load_luad_point_matrix_for_genes(
        maf_path, clinical_path, genes_rtk, mm_out.index
    )
    point_panel = point_panel.copy()
    point_panel["tp53_group"] = mm_out["tp53_group"].values

    wt, lof, gof = split_by_tp53_group(point_panel)
    # LoF vs GoF branch below uses only ``lof`` and ``gof`` rows (WT never enters those tests).
    lof_g = lof.drop(columns=["tp53_group"], errors="ignore")
    wt_g = wt.drop(columns=["tp53_group"], errors="ignore")
    gof_g = gof.drop(columns=["tp53_group"], errors="ignore")

    tab_pw = fisher_per_gene_lof_vs_wt(lof_g, wt_g, genes_rtk)
    tab_pw.to_csv(
        tables_dir / "step3_rtk_point_fisher_lof_vs_wt.tsv", sep="\t", index=False
    )

    hit_lof = _pathway_any_hit(lof_g, genes_rtk)
    hit_wt = _pathway_any_hit(wt_g, genes_rtk)
    n_lof, n_wt = len(lof_g), len(wt_g)
    a, b = int(hit_lof.sum()), n_lof - int(hit_lof.sum())
    c, d = int(hit_wt.sum()), n_wt - int(hit_wt.sum())
    or_path, p_path = (
        _fisher_2x2(a, b, c, d) if n_lof and n_wt else (np.nan, np.nan)
    )
    pd.DataFrame(
        [
            {
                "panel": "RTK_RAS_PI3K_any_gene_mutated",
                "n_lof": n_lof,
                "n_wt": n_wt,
                "mut_lof": a,
                "mut_wt": c,
                "odds_ratio": or_path,
                "p_value": p_path,
            }
        ]
    ).to_csv(
        tables_dir / "step3_rtk_pathway_any_mut_fisher_lof_vs_wt.tsv",
        sep="\t",
        index=False,
    )
    prev_path = pd.DataFrame(
        [
            {"group": "TP53_LoF", "prevalence": a / n_lof if n_lof else 0.0},
            {"group": "TP53_WT", "prevalence": c / n_wt if n_wt else 0.0},
        ]
    )
    fig, ax = plt.subplots(figsize=(4.2, 4.5))
    sns.barplot(
        data=prev_path,
        x="group",
        y="prevalence",
        hue="group",
        ax=ax,
        palette={"TP53_LoF": "#c0392b", "TP53_WT": "#2980b9"},
        legend=False,
    )
    ax.set_ylabel("Fraction with ≥1 coding mutation in panel")
    ax.set_xlabel("")
    ax.set_title("LUAD RTK/RAS/PI3K pathway burden\n(any listed gene mutated, MAF)")
    fig.tight_layout()
    fig.savefig(
        out_fig_dir / "rtk_pathway_any_mut_prevalence_lof_vs_wt.png",
        dpi=150,
        bbox_inches="tight",
    )
    if show_plots:
        plt.show()
    plt.close(fig)

    compact_bars = len(genes_rtk) <= 12
    if not compact_bars:
        for fn in (
            "rtk_point_prevalence_lof_vs_wt.png",
            "rtk_point_fisher_neglog10p_lof_vs_wt.png",
            "rtk_cna_amp_prevalence_lof_vs_wt.png",
        ):
            p = out_fig_dir / fn
            if p.is_file():
                p.unlink()

    if compact_bars:
        long_rows: list[dict] = []
        for _, r in tab_pw.iterrows():
            long_rows.append(
                {"gene": r["gene"], "group": "TP53_LoF", "prevalence": r["prev_lof"]}
            )
            long_rows.append(
                {"gene": r["gene"], "group": "TP53_WT", "prevalence": r["prev_wt"]}
            )
        fig, ax = plt.subplots(figsize=(max(6, 0.45 * len(tab_pw) + 3), 4.5))
        sns.barplot(
            data=pd.DataFrame(long_rows),
            x="gene",
            y="prevalence",
            hue="group",
            ax=ax,
            palette={"TP53_LoF": "#c0392b", "TP53_WT": "#2980b9"},
        )
        ax.set_ylabel("Fraction with coding mutation (MAF)")
        ax.set_xlabel("Gene (RTK / RAS / PI3K)")
        ax.set_title("LUAD: point mutations — TP53_LoF vs TP53_WT")
        ax.legend(title="")
        plt.xticks(rotation=25, ha="right")
        fig.tight_layout()
        fig.savefig(
            out_fig_dir / "rtk_point_prevalence_lof_vs_wt.png",
            dpi=150,
            bbox_inches="tight",
        )
        if show_plots:
            plt.show()
        plt.close(fig)

        tab_plot = tab_pw.copy()
        tab_plot["neglog10p"] = -np.log10(
            np.clip(tab_plot["p_value"].astype(float), 1e-300, None)
        )
        fig, ax = plt.subplots(figsize=(6, max(3, 0.35 * len(tab_plot) + 1)))
        y = np.arange(len(tab_plot))
        colors = np.where(tab_plot["odds_ratio"].astype(float) > 1, "#c0392b", "#2980b9")
        ax.barh(y, tab_plot["neglog10p"], color=colors, alpha=0.85)
        ax.set_yticks(y)
        ax.set_yticklabels(tab_plot["gene"].astype(str))
        ax.set_xlabel(r"$-\log_{10}(p)$ Fisher exact (LoF vs WT)")
        ax.set_title("LUAD RTK/RAS/PI3K: coding mutation enrichment LoF vs WT")
        ax.axvline(-np.log10(0.05), color="gray", linestyle="--", linewidth=0.8, alpha=0.7)
        fig.tight_layout()
        fig.savefig(
            out_fig_dir / "rtk_point_fisher_neglog10p_lof_vs_wt.png",
            dpi=150,
            bbox_inches="tight",
        )
        if show_plots:
            plt.show()
        plt.close(fig)

    tab_cna_lw: pd.DataFrame | None = None

    if raw_cna_path is not None and raw_cna_path.is_file():
        try:
            cna_raw = load_raw_cna_matrix(raw_cna_path)
            cna_sg = cna_gene_by_sample_to_sample_by_gene(cna_raw)
            if not any(g in cna_sg.columns for g in genes_rtk):
                print(
                    "[step3] CNA matrix has no RTK panel gene columns; CNA heatmap cells will be 0.",
                    file=sys.stderr,
                )
            amp_all = _cna_amp_submatrix(cna_sg, genes_rtk, cna_amp_threshold)
            cna_lof = _align_cna_to_samples(amp_all, lof.index)
            cna_wt = _align_cna_to_samples(amp_all, wt.index)
            tab_cna_lw = fisher_amp_per_gene_pair(cna_lof, cna_wt, genes_rtk)
            tab_cna_lw = tab_cna_lw.rename(
                columns={
                    "n_left": "n_lof",
                    "n_right": "n_wt",
                    "alt_left": "amp_lof",
                    "alt_right": "amp_wt",
                    "prev_left": "prev_lof",
                    "prev_right": "prev_wt",
                }
            )
            tab_cna_lw.to_csv(
                tables_dir / "step3_rtk_cna_amp_fisher_lof_vs_wt.tsv",
                sep="\t",
                index=False,
            )

            if compact_bars:
                long_c: list[dict] = []
                for _, r in tab_cna_lw.iterrows():
                    long_c.append(
                        {
                            "gene": r["gene"],
                            "group": "TP53_LoF",
                            "prevalence": r["prev_lof"],
                        }
                    )
                    long_c.append(
                        {
                            "gene": r["gene"],
                            "group": "TP53_WT",
                            "prevalence": r["prev_wt"],
                        }
                    )
                fig, ax = plt.subplots(
                    figsize=(max(6, 0.45 * len(tab_cna_lw) + 3), 4.5)
                )
                sns.barplot(
                    data=pd.DataFrame(long_c),
                    x="gene",
                    y="prevalence",
                    hue="group",
                    ax=ax,
                    palette={"TP53_LoF": "#c0392b", "TP53_WT": "#2980b9"},
                )
                ax.set_ylabel(
                    f"Fraction amplified (CNA ≥ {cna_amp_threshold})"
                )
                ax.set_xlabel("Gene (RTK / RAS / PI3K)")
                ax.set_title("LUAD: CNA amplification — TP53_LoF vs TP53_WT")
                ax.legend(title="")
                plt.xticks(rotation=25, ha="right")
                fig.tight_layout()
                fig.savefig(
                    out_fig_dir / "rtk_cna_amp_prevalence_lof_vs_wt.png",
                    dpi=150,
                    bbox_inches="tight",
                )
                if show_plots:
                    plt.show()
                plt.close(fig)

            tab_pw_idx = tab_pw.set_index("gene")
            tab_cna_idx = tab_cna_lw.set_index("gene")
            common = [g for g in tab_pw_idx.index if g in tab_cna_idx.index]
            if common:
                dx = (
                    tab_pw_idx.loc[common, "prev_lof"]
                    - tab_pw_idx.loc[common, "prev_wt"]
                ).astype(float)
                dy = (
                    tab_cna_idx.loc[common, "prev_lof"]
                    - tab_cna_idx.loc[common, "prev_wt"]
                ).astype(float)
                fig, ax = plt.subplots(figsize=(6.5, 5.5))
                ax.axhline(0, color="gray", lw=0.6)
                ax.axvline(0, color="gray", lw=0.6)
                ax.scatter(
                    dx,
                    dy,
                    s=80,
                    alpha=0.85,
                    c="#8e44ad",
                    edgecolors="k",
                    linewidths=0.4,
                )
                for g, x, y in zip(common, dx, dy):
                    ax.annotate(
                        str(g),
                        (float(x), float(y)),
                        textcoords="offset points",
                        xytext=(4, 4),
                        fontsize=9,
                    )
                ax.set_xlabel("Δ coding mutation prevalence (LoF − WT)")
                ax.set_ylabel(
                    f"Δ CNA amp prevalence (LoF − WT), CNA≥{cna_amp_threshold}"
                )
                ax.set_title(
                    "LUAD RTK/RAS/PI3K: point vs structural (amp) contrast"
                )
                fig.tight_layout()
                fig.savefig(
                    out_fig_dir / "rtk_point_vs_cna_delta_scatter.png",
                    dpi=150,
                    bbox_inches="tight",
                )
                if show_plots:
                    plt.show()
                plt.close(fig)
        except OSError as e:
            print(
                f"[step3] Could not read CNA from {raw_cna_path}: {e}",
                file=sys.stderr,
            )
        except ValueError as e:
            print(
                f"[step3] CNA parse/plot error ({raw_cna_path}): {e}",
                file=sys.stderr,
            )
    else:
        print(
            "[step3] CNA path missing or file not found; "
            "RTK enrichment figure will omit CNA panels (zeros).",
            file=sys.stderr,
        )

    plot_rtk_pathway_cohort_enrichment_figure(
        genes_rtk,
        tab_pw,
        tab_cna_lw,
        title=(
            "LUAD — RTK/RAS/PI3K genes\n"
            "TP53 loss-of-function vs TP53 wild-type (coding MAF + optional CNA amp)"
        ),
        out_path=out_fig_dir / "rtk_TP53_LoF_vs_WT_RTKpathway_enrichment.png",
        show_plots=show_plots,
        prev_a_col="prev_lof",
        prev_b_col="prev_wt",
        xlabel_a="TP53_LoF",
        xlabel_b="TP53_WT",
    )

    _write_rtk_pathway_enrichment_long_table(
        genes_rtk,
        tab_pw,
        tab_cna_lw,
        tables_dir / "step3_rtk_pathway_enrichment_scores_lof_vs_wt.tsv",
        "TP53_LoF_vs_TP53_WT",
        "TP53_LoF",
        "TP53_WT",
        prev_a_col="prev_lof",
        prev_b_col="prev_wt",
    )

    print(f"[step3] RTK LoF vs WT figures written under {out_fig_dir}")


def _write_rtk_pathway_enrichment_long_table(
    genes: list[str],
    tab_point: pd.DataFrame,
    tab_cna: pd.DataFrame | None,
    out_path: Path,
    comparison: str,
    label_a: str,
    label_b: str,
    prev_a_col: str,
    prev_b_col: str,
) -> None:
    """Long TSV: per gene, cohort prevalences + Fisher signed scores (coding and CNA)."""
    rows: list[dict] = []
    pp = tab_point.set_index("gene") if len(tab_point) else pd.DataFrame()
    cc = (
        tab_cna.set_index("gene")
        if tab_cna is not None and len(tab_cna)
        else None
    )
    sa = label_a.replace(" ", "_")
    sb = label_b.replace(" ", "_")
    for g in genes:
        if not pp.empty and g in pp.index:
            va = float(pp.loc[g, prev_a_col]) if prev_a_col in pp.columns else np.nan
            vb = float(pp.loc[g, prev_b_col]) if prev_b_col in pp.columns else np.nan
            sp = signed_neglog10p(
                float(pp.loc[g, "odds_ratio"]),
                float(pp.loc[g, "p_value"]),
            )
            rows.append(
                {
                    "comparison": comparison,
                    "gene": g,
                    "layer": "coding_MAF",
                    "metric": f"fraction_{sa}",
                    "value": va,
                }
            )
            rows.append(
                {
                    "comparison": comparison,
                    "gene": g,
                    "layer": "coding_MAF",
                    "metric": f"fraction_{sb}",
                    "value": vb,
                }
            )
            rows.append(
                {
                    "comparison": comparison,
                    "gene": g,
                    "layer": "coding_MAF",
                    "metric": "Fisher_signed_neglog10p",
                    "value": sp,
                    "odds_ratio": pp.loc[g, "odds_ratio"],
                    "p_value": pp.loc[g, "p_value"],
                }
            )
        if cc is not None and not cc.empty and g in cc.index:
            if prev_a_col in cc.columns and prev_b_col in cc.columns:
                va = float(cc.loc[g, prev_a_col])
                vb = float(cc.loc[g, prev_b_col])
                sc = signed_neglog10p(
                    float(cc.loc[g, "odds_ratio"]),
                    float(cc.loc[g, "p_value"]),
                )
                rows.append(
                    {
                        "comparison": comparison,
                        "gene": g,
                        "layer": "CNA_amp",
                        "metric": f"fraction_{sa}",
                        "value": va,
                    }
                )
                rows.append(
                    {
                        "comparison": comparison,
                        "gene": g,
                        "layer": "CNA_amp",
                        "metric": f"fraction_{sb}",
                        "value": vb,
                    }
                )
                rows.append(
                    {
                        "comparison": comparison,
                        "gene": g,
                        "layer": "CNA_amp",
                        "metric": "Fisher_signed_neglog10p",
                        "value": sc,
                        "odds_ratio": cc.loc[g, "odds_ratio"],
                        "p_value": cc.loc[g, "p_value"],
                    }
                )
    pd.DataFrame(rows).to_csv(out_path, sep="\t", index=False)


def main() -> None:
    parser = argparse.ArgumentParser(
        description="Step 3: TP53 LoF vs GoF/missense vs WT labels and annotated mutation matrix (LUAD).",
    )
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
    parser.add_argument(
        "--skip-rtk-figures",
        action="store_true",
        help="Skip RTK/RAS/PI3K LoF vs WT (and LoF vs GoF) figures under outputs/figures/lof_vs_gof/.",
    )
    parser.add_argument(
        "--raw-cna",
        default=str(REPO_ROOT / "data" / "data_cna.txt"),
        help="cBioPortal-style CNA matrix (genes × samples). Empty string skips CNA plots.",
    )
    parser.add_argument(
        "--cna-amp-threshold",
        type=int,
        default=2,
        help="Discrete amplification if CNA integer score >= this value (typical cBioPortal: 2).",
    )
    parser.add_argument(
        "--rtk-gene-list",
        default=str(REPO_ROOT / "config" / "pathways" / "rtk_ras_pi3k_genes.csv"),
        help="CSV with a ``gene`` column (# comments OK). If missing/empty, uses built-in RTK/RAS/PI3K list (see module).",
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

    if not args.skip_rtk_figures:
        rtk_csv = Path(args.rtk_gene_list)
        rtk_csv = rtk_csv if rtk_csv.is_file() else None
        raw_cna_arg = (args.raw_cna or "").strip()
        raw_cna_path = Path(raw_cna_arg) if raw_cna_arg else None
        if raw_cna_path is not None and not raw_cna_path.is_file():
            raw_cna_path = None
        maf_res = str(Path(args.maf).expanduser()) if args.maf else DEFAULT_MAF_PATH
        clin_res = (
            str(Path(args.clinical).expanduser())
            if args.clinical
            else DEFAULT_CLINICAL_PATH
        )
        run_rtk_lof_vs_wt_gof_visualizations(
            mm_out,
            maf_res,
            clin_res,
            REPO_ROOT / "outputs" / "figures" / "lof_vs_gof",
            raw_cna_path=raw_cna_path,
            cna_amp_threshold=int(args.cna_amp_threshold),
            rtk_gene_list_csv=rtk_csv,
            show_plots=False,
        )
        print(
            "[step3] TP53_LoF vs TP53_GoF_missense RTK/CNA heatmaps: "
            "run scripts/step4_fisher_pairwise_heatmap_tp53_lof_gof.py",
            file=sys.stderr,
        )


if __name__ == "__main__":
    main()
