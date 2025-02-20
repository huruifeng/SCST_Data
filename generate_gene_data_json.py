import pandas as pd
import json
import os

# %% =============================================
project = "ST"
# Load normalized expression data (sparse format)
print("Reading data...")
expression_data = pd.read_csv(project + "/normalized_expression_sparse.csv",index_col=None, header=0)

# %% ==========================================
# Group data by gene
grouped_by_gene = expression_data.groupby("Gene")

all_genes = grouped_by_gene.groups.keys()
all_genes = [gene_i.replace("/", "_") for gene_i in list(set(all_genes))]
with open(project + "/gene_list.json", "w") as f:
    json.dump(sorted(all_genes), f)

stop
# %% =============================================
# Create directory for genes
os.makedirs(project+ "/gene_jsons", exist_ok=True)

# Save each gene as a separate JSON file
print("Saving gene...")
i=0
for gene, df in grouped_by_gene:
    try:
        i += 1
        print(i,gene)
        gene_dict = dict(zip(df["Spot"], df["Expression"]))
        
        safe_gene_name = gene.replace("/", "_")

        # Create JSON file
        file_name = f"{project}/gene_jsons/{safe_gene_name}.json"
        if not os.path.exists(file_name):
            with open(file_name, "w") as f:
                json.dump(gene_dict, f, indent=4)
        else:
            print(f"File {file_name} already exists. Skipping...")

    except Exception as e:
        print(f"Error in processing {gene} !!! Check the error_gene.txt")
        with open(project + "/error_gene_json.txt","a") as f_err:
          f_err.write(gene + "\n")
