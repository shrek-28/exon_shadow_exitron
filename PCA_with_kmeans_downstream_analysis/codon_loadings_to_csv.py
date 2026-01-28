import pandas as pd
import os

species_list = [
    "bonobo", "borangutan", "chimpanzee",
    "sorangutan", "human", "gibbon", "gorilla"
]

windows = [5000, 7500, 10000]
regions = ["with_coding_region", "without_coding_region"]

base_dir = "PCA_with_kmeans"
out_dir = "DATA/pca_loadings"

os.makedirs(out_dir, exist_ok=True)

for species in species_list:
    for window in windows:
        for region in regions:

            in_path = os.path.join(
                base_dir,
                species,
                str(window),
                region,
                "PCA_loadings.tsv"
            )

            if not os.path.exists(in_path):
                continue

            df = pd.read_csv(in_path, sep="\t")

            region_fn = "coding" if region == "with_coding_region" else "noncoding"

            outname = f"{species}_{window}_{region_fn}_pca_loadings.csv"
            out_path = os.path.join(out_dir, outname)

            df.to_csv(out_path, index=False)