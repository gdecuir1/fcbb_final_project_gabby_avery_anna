# FCBB Final Project: TP53 modulation of RTK/RAS/PI3K and Wnt / cell-cycle co-occurrence

This repo contains a **bare-bones pipeline scaffold** for your class project. The numbered scripts (`scripts/01_...` onward) are intended to be filled in incrementally.

## Execution order (MVP)
1. `python scripts/01_acquire_mvp_data.py`
2. `python scripts/02_preprocess_mvp_primary_coding_binary_tp53.py`
3. `python scripts/03_fisher_pairwise_heatmap_tp53_binary.py`
4. `python scripts/04_classify_tp53_lof_vs_gof.py`
5. `python scripts/05_fisher_pairwise_heatmap_tp53_lof_gof.py`

## Optional scaling / backups
- DISCOVER pan-cancer (may be compute-heavy):
  - `python scripts/06_scale_discover_pan_cancer.py`
- Alternative hypotheses if pan-cancer is too noisy / slow:
  - `python scripts/07_alternative_hypothesis_pathway_redundancy.py`
  - `python scripts/08_alternative_hypothesis_structural_vs_point.py`

## Key files / folders to fill in
File(s) you will edit:
- `config/cohort_mvp.yaml`
- `config/tp53_classification.yaml`
- `config/gene_lists/mvp_driver_genes.csv` (start with ~30–50 RTK/RAS/PI3K, Wnt, and core cell-cycle drivers)
- `config/pathways/pathway_to_genes.csv` (pathway grouping for redundancy hypothesis)

Generated artifacts (by scripts):
- `data/raw/` (downloaded MAF/CNA sources)
- `data/processed/` (cleaned and binarized mutation tables, TP53 group assignments)
- `outputs/tables/` (pairwise stats tables)
- `outputs/figures/` (heatmaps, sanity-check plots)

