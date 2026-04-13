# Step 6: Pathway redundancy (TP53 LoF vs WT)

## Inputs
- Mutation matrix: `/Users/gabrielledecuir/Desktop/Comp Bio/fcbb_final_project_gabby_avery_anna/data/processed/lof_gof/mutation_matrix_with_tp53_group.csv`
- Pathway map: `/Users/gabrielledecuir/Desktop/Comp Bio/fcbb_final_project_gabby_avery_anna/config/pathways/pathway_to_genes.csv`
- TP53 status: `/Users/gabrielledecuir/Desktop/Comp Bio/fcbb_final_project_gabby_avery_anna/data/processed/lof_gof/tp53_functional_status.csv`

## Cohort sizes
- TP53 LoF samples (with pathway data): 4
- TP53 WT samples (with pathway data): 9

## Outputs
- Table: `/Users/gabrielledecuir/Desktop/Comp Bio/fcbb_final_project_gabby_avery_anna/outputs/tables/pathway_redundancy_tp53_lof_vs_wt.tsv`
- Figure: `/Users/gabrielledecuir/Desktop/Comp Bio/fcbb_final_project_gabby_avery_anna/outputs/figures/pathway_redundancy_heatmap.png`

## Top pathway pairs by `odds_ratio_diff` (LoF OR − WT OR)

| pathway_a | pathway_b | odds_ratio_diff | odds_ratio_lof | odds_ratio_wt | p_value_lof | p_value_wt |
| --- | --- | --- | --- | --- | --- | --- |
| Chromatin_remodeling | DNA_damage_cell_cycle | nan | nan | nan | 1 | 1 |
| Chromatin_remodeling | Other_tumor_suppressor | nan | inf | nan | 1 | 1 |
| Chromatin_remodeling | RTK_RAS_PI3K | nan | inf | inf | 1 | 0.007937 |
| Chromatin_remodeling | Wnt_signaling | nan | inf | inf | 1 | 0.4444 |
| DNA_damage_cell_cycle | Other_tumor_suppressor | nan | nan | nan | 1 | 1 |
| DNA_damage_cell_cycle | RTK_RAS_PI3K | nan | nan | nan | 1 | 1 |
| DNA_damage_cell_cycle | Wnt_signaling | nan | nan | nan | 1 | 1 |
| Other_tumor_suppressor | RTK_RAS_PI3K | nan | inf | nan | 1 | 1 |
| Other_tumor_suppressor | Wnt_signaling | nan | inf | nan | 1 | 1 |
| RTK_RAS_PI3K | Wnt_signaling | nan | inf | inf | 0.25 | 0.4444 |
