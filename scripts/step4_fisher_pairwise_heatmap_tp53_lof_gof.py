"""
step4_fisher_pairwise_heatmap_tp53_lof_gof.py

Pairwise Fisher exact tests and heatmaps, stratified by TP53 **LoF** vs **GoF/missense**
(two separate cohorts; WT samples are excluded).

Based on step2 logic; reads the step 3 matrix under `data/processed/lof_gof/`.

Inputs
------
- `data/processed/lof_gof/mutation_matrix_with_tp53_group.csv` (from step 3)

Outputs
-------
- Figures under `outputs/figures/` (optional): LoF and GoF heatmaps
"""

from __future__ import annotations

import argparse
import sys
from itertools import combinations
from pathlib import Path

import matplotlib.pyplot as plt
import numpy as np
import pandas as pd
import seaborn as sns
from scipy.stats import fisher_exact

# Allow imports from `scripts/` when invoked as `python scripts/step4_...py`.
_SCRIPT_DIR = Path(__file__).resolve().parent
REPO_ROOT = _SCRIPT_DIR.parent
if str(_SCRIPT_DIR) not in sys.path:
    sys.path.insert(0, str(_SCRIPT_DIR))

from step3_classify_tp53_lof_vs_gof import LABEL_GOF, LABEL_LOF

# Columns to drop before Fisher: cohort labels, not mutation indicators.
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
    """Same core logic as step2: pairwise Fisher over gene columns (0/1)."""
    # Gene-only 0/1 matrix for this stratum (one row = one patient in LoF or GoF cohort).
    df_genes = df.drop(columns=list(META_COLUMNS), errors="ignore")
    genes = df_genes.columns.tolist()
    n_genes = len(genes)

    if n_genes < 2:
        print(f"Skipping heatmap — need ≥2 genes, got {n_genes}. ({title})")
        return None

    # Mirror step 2: fill symmetric p-values for each unordered pair.
    p_val_matrix = pd.DataFrame(1.0, index=genes, columns=genes)
    for g1, g2 in combinations(genes, 2):
        # 2x2 co-occurrence within this TP53 functional stratum only.
        a = ((df_genes[g1] == 1) & (df_genes[g2] == 1)).sum()
        b = ((df_genes[g1] == 1) & (df_genes[g2] == 0)).sum()
        c = ((df_genes[g1] == 0) & (df_genes[g2] == 1)).sum()
        d = ((df_genes[g1] == 0) & (df_genes[g2] == 0)).sum()

        table = [[a, b], [c, d]]
        _, pval = fisher_exact(table)
        p_val_matrix.loc[g1, g2] = pval
        p_val_matrix.loc[g2, g1] = pval

    # Visualize strength of association via -log10(p); mask duplicates upper triangle.
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
    # Uncorrected p <= 0.05 markers (same convention as step 2).
    for i in range(n_genes):
        for j in range(i + 1, n_genes):
            if p_val_matrix.iloc[j, i] <= 0.05:
                plt.text(
                    i + 0.5, j + 0.5, "*",
                    ha="center", va="center", color="black", fontsize=20,
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


def main() -> None:
    parser = argparse.ArgumentParser(
        description="Fisher pairwise heatmaps: TP53 LoF vs GoF/missense strata.",
    )
    parser.add_argument(
        "--matrix",
        default=_DEFAULT_MATRIX,
        help="Step 3 CSV: mutation matrix + tp53_group.",
    )
    parser.add_argument(
        "--out-dir",
        default=str(REPO_ROOT / "outputs" / "figures"),
        help="Directory for saved PNGs.",
    )
    parser.add_argument(
        "--no-show",
        action="store_true",
        help="Do not open interactive plot windows (save PNG only).",
    )
    args = parser.parse_args()

    # Default matrix is step 3 export under `data/processed/lof_gof/`.
    path = Path(args.matrix).expanduser()
    if not path.is_absolute():
        path = REPO_ROOT / path
    if not path.exists():
        raise FileNotFoundError(
            f"Matrix not found: {path}. Run step3_classify_tp53_lof_vs_gof.py first.",
        )

    # Rows = patients; columns = gene binaries + `TP53_status` + `tp53_group`.
    mm = pd.read_csv(path, index_col="sample_id")
    if "tp53_group" not in mm.columns:
        raise ValueError(f"Expected column 'tp53_group' in {path}")

    # WT patients are excluded: only compare LoF vs GoF/missense strata side by side.
    tp53_lof = mm[mm["tp53_group"] == LABEL_LOF].copy()
    tp53_gof = mm[mm["tp53_group"] == LABEL_GOF].copy()

    out_dir = Path(args.out_dir)
    if not out_dir.is_absolute():
        out_dir = REPO_ROOT / out_dir

    print(f"TP53 LoF samples: {len(tp53_lof)} | TP53 GoF/missense samples: {len(tp53_gof)}")

    # Each stratum is optional: tiny cohorts may still run but yield noisy p-values.
    if len(tp53_lof) == 0:
        print("No TP53_LoF samples; skipping LoF heatmap.")
    else:
        print("Generating heatmap for TP53 LoF cohort...")
        plot_fisher_heatmap(
            tp53_lof,
            title="Gene co-occurrence (TP53 loss-of-function)",
            save_path=out_dir / "fisher_tp53_lof_pairwise_heatmap.png",
            show=not args.no_show,
        )

    if len(tp53_gof) == 0:
        print("No TP53_GoF_missense samples; skipping GoF heatmap.")
    else:
        print("Generating heatmap for TP53 GoF / missense cohort...")
        plot_fisher_heatmap(
            tp53_gof,
            title="Gene co-occurrence (TP53 GoF / missense)",
            save_path=out_dir / "fisher_tp53_gof_missense_pairwise_heatmap.png",
            show=not args.no_show,
        )


if __name__ == "__main__":
    main()
