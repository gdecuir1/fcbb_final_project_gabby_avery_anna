import csv

def convert_gmt_to_csv(gmt_filename, output_csv_filename):
    with open(gmt_filename, 'r') as gmt_file, open(output_csv_filename, 'w', newline='') as csv_file:
        writer = csv.writer(csv_file)
        # Write the header Step 7 expects
        writer.writerow(['pathway', 'gene'])
        
        for line in gmt_file:
            # GMT files are tab-separated
            parts = line.strip().split('\t')
            if len(parts) < 3:
                continue
                
            pathway_name = parts[0]
            # parts[1] is usually a description or URL, so we skip it.
            # The actual genes start at index 2
            genes = parts[2:]
            
            for gene in genes:
                if gene:  # Make sure the gene string isn't empty
                    writer.writerow([pathway_name, gene])

if __name__ == "__main__":
    # Replace with the exact name of the file you downloaded
    input_gmt = "config/pathways/gene_symbols.gmt" 
    output_csv = "config/pathways/pathway_to_genes.csv"
    
    convert_gmt_to_csv(input_gmt, output_csv)
    print(f"Successfully converted {input_gmt} to {output_csv}!")