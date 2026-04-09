import pandas as pd

def preprocess():
    data_all = pd.read_csv('data/data_mutations.txt', sep='\t') # first row is header
    # preprocessing and making mutation matrix
    # proposed gene list 
    gene_list = [
        # Main genes 
        "TP53", "RB1", "ATM",

        # RTK/RAS/PI3K pathway 
        "KRAS", "EGFR", "PIK3CA", "AKT3", "PTEN", "ERBB4",

        # Wnt pathway 
        "APC", "CTNNB1",

        # Lung Adenocarcinoma
        "STK11", "SMARCA4", "PBRM1", "CREBBP",

        # Extras
        "MGA", "PTPRD"
    ]

    # determine different types of mutations 
    unique_values = data_all['Variant_Classification'].unique()

    # list coding mutations 
    coding_mutations = [
        "Missense_Mutation", "Nonsense_Mutation",
        "Frame_Shift_Ins", "Frame_Shift_Del",
        "Splice_Site"
    ]

    # filter for only coding mutations 
    data_filtered = data_all[data_all["Variant_Classification"].isin(coding_mutations)]

    # filter for only genes in the gene list 
    data_filtered = data_filtered[data_filtered["Hugo_Symbol"].isin(gene_list)]


    # filter for only Lung Adenocarcinoma patients by joining another table with only those cancer types 
    clinical_data_all = pd.read_csv('data/lung_msk_mind_2020_clinical_data (1).tsv', sep='\t')
    clinical_data_LA = clinical_data_all[(clinical_data_all["Cancer Type"] == "Lung Adenocarcinoma") & (clinical_data_all["Cancer Type Detailed"] == "Lung Adenocarcinoma")].copy()

    # merge tables 
    df = pd.merge(
        data_filtered,
        clinical_data_LA,
        left_on='Tumor_Sample_Barcode',
        right_on='Sample ID',
        how='inner'
    )

    # pivot table: one row per sample, one column per gene
    # 1 = mutation at gene; 0 = no mutation at gene 
    mutation_matrix = (
        df.assign(mut=1)
        .pivot_table(
            index="Tumor_Sample_Barcode",
            columns="Hugo_Symbol",
            values="mut",
            aggfunc="max",
            fill_value=0
        )
    )

    # verify process of joining tables 
    print(df.shape) # (55, 122)
    print(df["Tumor_Sample_Barcode"].nunique()) # 12
    print(df["Sample ID"].nunique()) # 12 
    # there are only 12 rows/patients instead of 17 b/c patients were excluded for not having coding-mutations
    # there are 55 rows instead of 12 because each row = one mutation event, and each sample has multiple mutations


    # make sure all LUAD samples are kept, including samples with 0 retained mutations
    all_luad_samples = clinical_data_LA["Sample ID"].unique()
    mutation_matrix = mutation_matrix.reindex(all_luad_samples, fill_value=0)

    # sanity checks
    print(mutation_matrix.shape) # (17, 15)
    # 17 rows = all patients; 15 columnds = all genes retained 

    print(mutation_matrix.head())
    # shows that the table has rows of tumor samples and columns of gene names 

    # all rows' IDs are unique
    print(mutation_matrix.index.is_unique) # True

    # check for no NaNs 
    print(sorted(pd.unique(mutation_matrix.values.ravel())))

    # check for samples with at least one retained mutation
    print((mutation_matrix.sum(axis=1) > 0).sum()) # 12
    # means there are 5 tumor samples with no coding/non-silent mutations in selected genes

    # check which genes were retained after joining tables
    print(mutation_matrix.columns.tolist())
    # lost two genes: PTEN (RTK/RAS/PI3K pathway) and CTNNB1 (Wnt pathway)

    # for each sample-gene pair in the merged long table, the matrix entry should be 1
    check_pairs = df[["Tumor_Sample_Barcode", "Hugo_Symbol"]].drop_duplicates()

    all_checks_pass = True
    for _, row in check_pairs.iterrows():
        sample = row["Tumor_Sample_Barcode"]
        gene = row["Hugo_Symbol"]
        if mutation_matrix.loc[sample, gene] != 1:
            all_checks_pass = False
            print(f"Mismatch found: sample={sample}, gene={gene}")

    print(all_checks_pass) # True
    # means that for every mutation event in the raw data correctly appears as a 1 in the matrix
    # confirms no data loss/pivot errors 

    mutation_matrix["TP53_status"] = mutation_matrix["TP53"].apply(
        lambda x: "Mut" if x == 1 else "WT"
    )

    print("TP53 status counts:")
    print(mutation_matrix["TP53_status"].value_counts())

    # split by TP53 status
    tp53_mut = mutation_matrix[mutation_matrix["TP53_status"] == "Mut"].copy()
    tp53_wt  = mutation_matrix[mutation_matrix["TP53_status"] == "WT"].copy()

    return tp53_mut, tp53_wt

if __name__ == "__main__":
    tp53_mut, tp53_wt = preprocess()