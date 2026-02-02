import os
import pandas as pd

folder_path = "gene_enrichment_data_processing/results"  # CHANGE THIS

for file in os.listdir(folder_path):
    if file.endswith(".csv"):
        file_path = os.path.join(folder_path, file)
        df = pd.read_csv(file_path)

        if "gibbon" in df.columns:
            df["gibbon"] = (
                df["gibbon"]
                .astype(str)
                .str.split("-", n=1)
                .str[0]
            )

            df.to_csv(file_path, index=False)
