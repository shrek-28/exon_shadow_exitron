import pandas as pd

# --------------------------------------------------
# MANUALLY DEFINE FILES
# --------------------------------------------------
# Each tuple: (file_path, dataset_name, mode)
# mode ∈ {"count", "length", "shadow"}

files = [
    (
        r"final_strand_shadow_asymmetry/asymmetry/all_species_direction_total_shadow_count_asymmetry.tsv",
        "Shadow Count",
        "shadow"
    ),
    (
        r"final_strand_shadow_asymmetry/asymmetry/all_species_direction_total_shadow_length_asymmetry.tsv",
        "Shadow Length",
        "shadow"
    ),
    (
        r"final_strand_shadow_asymmetry/asymmetry/all_species_strand_cds_count_asymmetry.tsv",
        "CDS",
        "count"
    ),
    (
        r"final_strand_shadow_asymmetry/asymmetry/all_species_strand_cds_length_asymmetry.tsv",
        "CDS", 
        "length"
    ), 
    (
        r"final_strand_shadow_asymmetry/asymmetry/all_species_strand_down_shadow_length_asymmetry.tsv",
        "Downstream Shadow",
        "length"
    ), 
    (
        r"final_strand_shadow_asymmetry/asymmetry/all_species_strand_down_shadow_positive_count_asymmetry.tsv",
        "Downstream Shadow",
        "count"
    ),
    (
        r"final_strand_shadow_asymmetry/asymmetry/all_species_strand_scfr_count_asymmetry.tsv",
        "SCFR",
        "count"
    ),
    (
        r"final_strand_shadow_asymmetry/asymmetry/all_species_strand_scfr_length_asymmetry.tsv",
        "SCFR",
        "length"
    ),
    (
        r"final_strand_shadow_asymmetry/asymmetry/all_species_strand_total_shadow_count_asymmetry.tsv",
        "Total Shadow", 
        "count"
    ), 
    (
        r"final_strand_shadow_asymmetry/asymmetry/all_species_strand_total_shadow_length_asymmetry.tsv",
        "Total Shadow", 
        "length"
    ), 
    (
        r"final_strand_shadow_asymmetry/asymmetry/all_species_strand_up_shadow_length_asymmetry.tsv",
        "Upstream Shadow",
        "length"
    ),
    (
        r"final_strand_shadow_asymmetry/asymmetry/all_species_strand_up_shadow_positive_count_asymmetry.tsv",
        "Upstream Shadow",
        "count"
    )
    # add more entries if needed
]

# --------------------------------------------------
# SPECIES ORDER + FINAL DISPLAY NAMES
# --------------------------------------------------
species_map = {
    "human": "Human",
    "bonobo": "Bonobo",
    "chimpanzee": "Chimpanzee",
    "gorilla": "Gorilla",
    "borangutan": "Bornean Orangutan",
    "sorangutan": "Sumatran Orangutan",
    "gibbon": "Gibbon"
}

species_order = list(species_map.keys())

# --------------------------------------------------
# STORAGE
# --------------------------------------------------
count_rows = []
length_rows = []
shadow_rows = []

# --------------------------------------------------
# PROCESS FILES
# --------------------------------------------------
for file_path, dataset_name, mode in files:
    df = pd.read_csv(file_path, sep="\t", header=None)
    last_col = df.columns[-1]  # asymmetry value

    row = {"dataset": dataset_name}

    for sp in species_order:
        val = df.loc[df.iloc[:, 0] == sp, last_col].values
        row[species_map[sp]] = val[0] if len(val) > 0 else None

    if mode == "count":
        count_rows.append(row)
    elif mode == "length":
        length_rows.append(row)
    elif mode == "shadow":
        shadow_rows.append(row)
    else:
        raise ValueError(f"Invalid mode: {mode}")

# --------------------------------------------------
# FINAL DATAFRAMES
# --------------------------------------------------
df_count = pd.DataFrame(count_rows)
df_length = pd.DataFrame(length_rows)
df_shadow = pd.DataFrame(shadow_rows)

# --------------------------------------------------
# SAVE CSVs
# --------------------------------------------------
df_count.to_csv("final_strand_shadow_asymmetry/intermed_data/strand_count_asymmetry_summary.csv", index=False)
df_length.to_csv("final_strand_shadow_asymmetry/intermed_data/strand_length_asymmetry_summary.csv", index=False)
df_shadow.to_csv("final_strand_shadow_asymmetry/intermed_data/shadow_asymmetry_summary.csv", index=False)
