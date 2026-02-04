import pandas as pd 
import re 
import os 
from collections import defaultdict

def gene_enrichment_processor(files, final_path):

    def clean_pathway(pathway, species_name=None):
        if pd.isna(pathway):
            return pathway

        pathway = str(pathway).strip()

        if species_name is not None:
            species_lower = species_name.lower()

            if species_lower == "gorilla":
                # do nothing at all
                return pathway.strip()

            elif species_lower == "gibbon":
                # remove KEGG prefix first
                pathway = re.sub(r"^[^ ]+\s+", "", pathway)
                # remove gibbon-specific suffix
                if "Nomascus leucogenys" in pathway:
                    pathway = re.sub(r"\s*-\s*Nomascus leucogenys.*$", "", pathway)

            else:
                # other species → remove KEGG prefix only
                pathway = re.sub(r"^[^ ]+\s+", "", pathway)

        return pathway.strip()

    def parse_genes(gene_string):
        if pd.isna(gene_string) or gene_string == "":
            return set()
        return set(re.split(r"[|\s]+", str(gene_string).strip()))

    # property → species → gene set
    property_species_genes = defaultdict(dict)

    # property → species → enrichment quantity
    property_species_metric = defaultdict(dict)

    for species, path in files.items():
        if not os.path.exists(path):
            continue

        df = pd.read_csv(path)

        if "Pathway" not in df.columns or "Genes" not in df.columns:
            raise KeyError(f"Missing required columns in {path}")

        # Check for required columns for metric
        if not {"nGenes", "Pathway Genes", "Fold Enrichment"}.issubset(df.columns):
            raise KeyError(f"{path} missing nGenes, Pathway Genes, or Fold Enrichment")

        df["Pathway"] = df["Pathway"].apply(clean_pathway, species_name=species)

        for _, row in df.iterrows():
            prop = row["Pathway"]
            genes = parse_genes(row["Genes"])
            if genes:
                property_species_genes[prop][species] = genes

            # --- compute enrichment_quantity ---
            if pd.notna(row["nGenes"]) and pd.notna(row["Pathway Genes"]) and pd.notna(row["Fold Enrichment"]):
                if row["Fold Enrichment"] > 0 and row["nGenes"] > 0 and row["Pathway Genes"] > 0:
                    metric = ((row["nGenes"] / row["Pathway Genes"]) * row["Fold Enrichment"]) ** 0.5
                else:
                    metric = 0
            else:
                metric = 0

            property_species_metric[prop][species] = metric

    rows = []

    for prop, species_genes in property_species_genes.items():
        species_list = sorted(species_genes.keys())
        gene_sets = list(species_genes.values())

        # common genes
        common_genes = set.intersection(*gene_sets) if gene_sets else set()

        # gene → species map
        gene_species_map = defaultdict(set)
        for sp, genes in species_genes.items():
            for g in genes:
                gene_species_map[g].add(sp)

        # enrichment quantity per species
        metric_dict = property_species_metric[prop]

        rows.append({
            "property": prop,
            "n_species": len(species_list),
            "species": "; ".join(species_list),

            "n_Genes_per_species": "; ".join(
                f"{sp}:{len(species_genes[sp])}"
                for sp in species_list
            ),

            "genes_common_in_all_species": "|".join(sorted(common_genes)),
            "genes_common_in_all_species_count": len(common_genes),

            "genes_unique_to_one_species": "; ".join(
                f"{sp}:{'|'.join(sorted(species_genes[sp] - common_genes))}"
                for sp in species_list
            ),

            # gene → number of species
            "gene_species_count": "|".join(
                f"{g}:{len(gene_species_map[g])}"
                for g in sorted(gene_species_map)
            ),

            # gene → species list
            "gene_species_list": "|".join(
                f"{g}:{','.join(sorted(gene_species_map[g]))}"
                for g in sorted(gene_species_map)
            ),

            # --- enrichment quantity ---
            "enrichment_quantity_per_species": "; ".join(
                f"{sp}:{metric_dict.get(sp, 0):.4f}" for sp in species_list
            )
        })

    pd.DataFrame(rows).to_csv(final_path, index=False)

files = {
        "bonobo":"exitron_gene_enrichment/gsea_raw_results/bonobo_results.csv", 
        "borangutan": "exitron_gene_enrichment/gsea_raw_results/borangutan_results.csv",
        "chimpanzee": "exitron_gene_enrichment/gsea_raw_results/chimpanzee_results.csv", 
        "gibbon": "exitron_gene_enrichment/gsea_raw_results/gibbon_results.csv",
        "gorilla": "exitron_gene_enrichment/gsea_raw_results/gorilla_results.csv",
        "human": "exitron_gene_enrichment/gsea_raw_results/human_results.csv", 
        "sorangutan": "exitron_gene_enrichment/gsea_raw_results/sorangutan_results.csv"
}


final_path = "exitron_gene_enrichment/gsea_processed_results/plot.csv"

gene_enrichment_processor(files, final_path)