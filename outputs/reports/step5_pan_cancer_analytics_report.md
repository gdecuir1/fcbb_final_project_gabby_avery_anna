# Step 5 — Pan-cancer analytics (MC3 public MAF)

- **Generated (UTC):** 2026-04-12 19:30:38Z
- **MAF:** `/Users/gabrielledecuir/Desktop/Comp Bio/fcbb_final_project_gabby_avery_anna/data/raw/mc3/mc3.v0.2.8.PUBLIC.maf.gz`
- **Samples (after filters):** 6,562
- **Curated genes:** 17
- **Study column used:** *not detected — pan-cancer only*

## Methods (Python stand-in for DISCOVER screening)

1. Stream the MC3 open-access MAF in chunks; keep **coding** variant classes matching step 1.
2. Restrict to the same **driver gene list** as step 1.
3. Build a **sample × gene** binary mutation matrix.
4. **TP53** is classified per sample (LoF vs missense proxy vs WT) using the same rules as step 3.
5. Run **pairwise Fisher exact tests** on all unordered gene pairs (co-occurrence table).
6. Apply **Benjamini–Hochberg FDR** (`q_value_bh`) across tests.

The R/Bioconductor **DISCOVER** package uses a dedicated statistical model; this pipeline is a
lighter-weight **Fisher-based screen** suitable for coursework and for spotting pairs to follow up.

## Pan-cancer mutation prevalence (curated genes)

```
   gene  fraction_samples_mutated
   TP53                  0.579549
 PIK3CA                  0.206644
   PTEN                  0.131210
    APC                  0.125876
   KRAS                  0.119781
    ATM                  0.085187
  PTPRD                  0.079397
 CREBBP                  0.071625
  ERBB4                  0.070253
    RB1                  0.068272
    MGA                  0.065986
  PBRM1                  0.063395
 CTNNB1                  0.062786
SMARCA4                  0.058976
   EGFR                  0.055319
  STK11                  0.023011
   AKT3                  0.017525
```

## TP53 functional groups (same rules as step 3)

```
       tp53_group  n_samples
          TP53_WT       2759
TP53_GoF_missense       2460
         TP53_LoF       1343
```

## Strongest mutual exclusivity signals (low odds ratio, small q)

*(Fisher odds ratio below 1 suggests fewer co-mutations than independence would predict.)*

```
gene_a gene_b  odds_ratio      p_value   q_value_bh  n_both  g1_only  g2_only
  TP53 CTNNB1    0.226780 2.717179e-43 9.238410e-42     105     3698      307
  TP53  PBRM1    0.254637 3.680545e-38 1.001108e-36     115     3688      301
  TP53   PTEN    0.258586 2.436650e-71 3.313844e-69     257     3546      604
  TP53 PIK3CA    0.339231 5.313718e-68 3.613329e-66     502     3301      854
  TP53  STK11    0.340011 1.828336e-10 4.604699e-10      49     3754      102
  TP53    ATM    0.367956 1.180830e-28 1.784365e-27     199     3604      360
  TP53    MGA    0.425992 2.562552e-17 1.290767e-16     166     3637      267
  TP53   EGFR    0.434646 2.303489e-14 8.950699e-14     140     3663      223
  TP53   AKT3    0.513608 5.447253e-04 7.334915e-04      48     3755       67
  TP53   KRAS    0.521329 1.589838e-17 8.316074e-17     344     3459      442
  KRAS   EGFR    0.552140 2.667865e-03 3.328712e-03      26      760      337
  KRAS  PBRM1    0.574941 2.307323e-03 2.905518e-03      31      755      385
PIK3CA  STK11    0.579947 2.471514e-02 2.897637e-02      20     1336      131
  TP53 CREBBP    0.593877 6.463970e-08 1.274058e-07     216     3587      254
  TP53  PTPRD    0.629957 5.388207e-07 9.642054e-07     247     3556      274
```

## Strongest co-occurrence signals (high odds ratio, small q)

```
 gene_a  gene_b  odds_ratio      p_value   q_value_bh  n_both  g1_only  g2_only
    ATM    AKT3    8.011977 1.660914e-21 1.188865e-20      47      512       68
   AKT3  CREBBP    7.462946 3.147661e-18 1.783674e-17      40       75      430
   AKT3     MGA    7.248349 7.550587e-17 3.540965e-16      37       78      396
   AKT3   PTPRD    7.168359 2.567008e-18 1.517883e-17      42       73      479
   AKT3 SMARCA4    5.735384 1.011040e-11 2.806152e-11      29       86      358
   AKT3   STK11    5.287141 1.141470e-05 1.784367e-05      12      103      139
  ERBB4   PTPRD    5.286653 7.263108e-38 1.646305e-36     124      337      397
   AKT3   ERBB4    5.164092 2.806980e-11 7.634986e-11      31       84      430
 CREBBP     MGA    5.055296 1.194185e-31 2.030115e-30     105      365      328
   AKT3     APC    4.843590 7.224769e-14 2.655591e-13      46       69      780
   EGFR    AKT3    4.751921 1.070140e-08 2.274048e-08      24      339       91
    RB1    AKT3    4.618391 1.740997e-09 4.013145e-09      28      420       87
  ERBB4 SMARCA4    4.498102 3.640193e-24 3.536187e-23      87      374      300
   AKT3  CTNNB1    4.349699 2.830537e-08 5.745568e-08      25       90      387
SMARCA4   PTPRD    4.318661 1.129397e-24 1.279983e-23      94      293      427
```

## Outputs

- `outputs/tables/step5_pancancer_pairwise_fisher.tsv`
- `outputs/tables/step5_mutation_frequency_by_gene.tsv`
- `outputs/tables/step5_mutation_frequency_by_study_gene.tsv` (if study labels available)
- `outputs/tables/step5_tp53_group_counts.tsv`
- `outputs/figures/step5_pancancer_fisher_logp_heatmap.png`
