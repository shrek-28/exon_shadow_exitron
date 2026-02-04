import pandas as pd
import itertools

# ================================
# 1. USER INPUT (NO AUTODETECTION)
# ================================

# Each tuple: (csv_path, foreground_label, background_label)
FILES = [
   ("exitron_gene_enrichment/gsea_processed_results/plot.csv", "exitrons", "no_bg")
    # add all your files explicitly
]

FOREGROUNDS = [
    "exitrons"
]

BACKGROUNDS = [
    "no_bg"
]

dfs = []

# ================================
# 2. PROCESS FILES
# ================================
for path, fg, bg in FILES:
    df = pd.read_csv(path)

    # Clean pathway names if gibbon info is present
    df["property"] = df.apply(
    lambda row: str(row["property"]).split("-", 1)[0].strip()
    if "gibbon" in str(row["species"]).lower()
    else str(row["property"]),
    axis=1
    )

    # ================================
    # STANDARDIZE SPECIFIC PATHWAY NAMES
    # ================================

    df["property"] = df["property"].replace({
        "Pathways of neurodegeneration-multiple diseases":
            "Pathways of Neurodegeneration - Multiple Diseases",
        "Chemical carcinogenesis-reactive oxygen species":
            "Chemical Carcinogens - ROS"
    })


    required = {"property", "species", "enrichment_quantity_per_species"}
    missing_cols = required - set(df.columns)
    if missing_cols:
        raise ValueError(f"{path} missing columns: {missing_cols}")

    df["foreground"] = fg
    df["background"] = bg

    # Fix "gorillla" typo
    df["species"] = df["species"].str.replace(r"\bgorillla\b", "gorilla", regex=True)

    # Explode rows for multiple species in the 'species' column
    if df["species"].str.contains(";").any():
        df["species"] = df["species"].str.split(";")
        df = df.explode("species")

    df["species"] = df["species"].str.strip()

    # ================================
    # 2a. EXTRACT ENRICHMENT QUANTITY FOR THE SPECIES
    # ================================
    def extract_enrichment(row):
        text = row["enrichment_quantity_per_species"]
        if pd.isna(text) or text.strip() == "":
            return 0.0

        # Split by semicolon
        entries = text.split(";")
        q_dict = {}
        for e in entries:
            if ":" in e:
                sp, val = e.split(":")
                sp = sp.strip().replace("gorillla", "gorilla")
                try:
                    val = float(val.strip())
                except:
                    val = 0.0
                q_dict[sp] = val

        # Return value for this row’s species
        return q_dict.get(row["species"], 0.0)

    df["enrichment_quantity"] = df.apply(extract_enrichment, axis=1)

    dfs.append(df[[
        "property", "species", "foreground", "background", "enrichment_quantity"
    ]])

# ================================
# 3. CONCAT ALL CSVs
# ================================
df_present = pd.concat(dfs, ignore_index=True)

# ================================
# 4. BUILD FULL PERMUTATION SPACE
# ================================
all_pathways = df_present["property"].unique()
all_species = df_present["species"].unique()

full_grid = pd.DataFrame(
    itertools.product(
        all_pathways,
        all_species,
        FOREGROUNDS,
        BACKGROUNDS
    ),
    columns=["property", "species", "foreground", "background"]
)

# ================================
# 5. MERGE + FILL MISSING
# ================================
df_full = full_grid.merge(
    df_present,
    on=["property", "species", "foreground", "background"],
    how="left"
)

df_full["enrichment_quantity"] = df_full["enrichment_quantity"].fillna(0.0)

# ================================
# 6. HELPER COLUMNS
# ================================
df_full["enriched"] = (df_full["enrichment_quantity"] > 0).astype(int)
df_full["enrichment_state"] = df_full["enriched"].map({1: "enriched", 0: "not_enriched"})

df_full = df_full.drop(['foreground', 'background'], axis=1)
# ================================
# 7. EXPORT
# ================================
df_full.to_csv("exitron_gene_enrichment/gsea_processed_results/enrichment_long_format_with_quantity.csv", index=False)

print("DONE")
print("Rows:", df_full.shape[0])
print("Unique pathways:", df_full['property'].nunique())
print("Unique species:", df_full['species'].nunique())