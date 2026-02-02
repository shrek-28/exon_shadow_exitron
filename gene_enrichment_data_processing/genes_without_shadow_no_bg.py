import pandas as pd

files = {
    "bonobo": "data/exon_shadow/gene_enrichment/bonobo_gene_enrichment/genes_without_shadow_nobg_enrichment.csv",
    "chimpanzee": "data/exon_shadow/gene_enrichment/chimpanzee_gene_enrichment/genes_without_shadow_no_bg_enrichment.csv",
    "gorilla": "data/exon_shadow/gene_enrichment/gorillla_gene_enrichment/genes_without_shadow_nobg_enrichment.csv",
    "borangutan": "data/exon_shadow/gene_enrichment/borangutan_gene_enrichment/genes_without_shadow_no_bg_enrichment.csv",
    "human": "data/exon_shadow/gene_enrichment/human_gene_enrichment/genes_without_shadow_no_bg_enrichment.csv",
    # "gibbon": "PATH/TO/gibbon.csv",
    "sorangutan": "data/exon_shadow/gene_enrichment/sorangutan_gene_enrichment/genes_without_shadow_no_bg_enrichment.csv",
}

dfs = []

for species, path in files.items():
    df = pd.read_csv(path)

    if "Pathway" not in df.columns:
        raise KeyError(f"'Pathway' column missing in {path}")

    df["Pathway"] = df["Pathway"].str.replace(r"^[^ ]+\s+", "", regex=True)

    dfs.append(df[["Pathway"]].rename(columns={"Pathway": species}))

final_df = pd.concat(dfs, axis=1)
final_df.to_csv("gene_enrichment_data_processing/results/combined_list_genes_without_shadow_no_bg.csv", index=False)

print(final_df.head())