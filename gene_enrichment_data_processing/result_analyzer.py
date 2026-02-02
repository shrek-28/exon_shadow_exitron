import os
import pandas as pd

folder_path = "gene_enrichment_data_processing/results"  
out_folder_path = "gene_enrichment_data_processing/results/analysed_results"

for file in os.listdir(folder_path):
    if not file.endswith(".csv"):
        continue

    df = pd.read_csv(os.path.join(folder_path, file))

    property_map = {}

    for species in df.columns:
        values = df[species].dropna().astype(str)

        for prop in values:
            if prop not in property_map:
                property_map[prop] = set()
            property_map[prop].add(species)

    result = pd.DataFrame({
        "property": property_map.keys(),
        "n_species": [len(v) for v in property_map.values()],
        "species": ["; ".join(sorted(v)) for v in property_map.values()]
    })

    out_file = file.replace(".csv", "_property_occurrence.csv")
    result.to_csv(os.path.join(out_folder_path, out_file), index=False)
