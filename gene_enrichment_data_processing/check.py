import os

BASE_DIR = "data/exon_shadow/gene_enrichment"

TARGET_FILES = (
    "geneswshadow_no_bg_enrichment.csv",
    "geneswshadow_nobg_enrichment.csv"
)

species_to_file = {}

for folder in os.listdir(BASE_DIR):
    folder_path = os.path.join(BASE_DIR, folder)

    if not os.path.isdir(folder_path):
        continue
    if not folder.endswith("_gene_enrichment"):
        continue

    species = folder.replace("_gene_enrichment", "")

    for fname in TARGET_FILES:
        full_path = os.path.join(folder_path, fname)
        if os.path.exists(full_path):
            species_to_file[species] = full_path
            break  # stop after first match

# result
print(species_to_file)
