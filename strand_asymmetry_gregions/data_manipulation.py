import pandas as pd

# ----------------------------------
# GLOBAL SPECIES ORDER (FOR PLOTTING)
# ----------------------------------

ORDER = [
    "Gibbon", 
    "Gorilla",
    "Human",
    "Bonobo",
    "Chimpanzee",
    "Bornean Orangutan",
    "Sumatran Orangutan",
]

# ----------------------------------
# SPECIES NORMALIZATION (ROBUST)
# ----------------------------------

def normalize_species(s):
    s = str(s).strip().lower()

    mapping = {
        "human": "Human",
        "bonobo": "Bonobo",
        "chimpanzee": "Chimpanzee",
        "gibbon": "Gibbon",
        "gorilla": "Gorilla",
        "borangutan": "Bornean Orangutan",
        "sorangutan": "Sumatran Orangutan",
        "bornean orangutan": "Bornean Orangutan",
        "sumatran orangutan": "Sumatran Orangutan",
    }

    return mapping.get(s, s.title())


# ----------------------------------
# GENERIC ASYMMETRY LOADER (SAFE)
# ----------------------------------

def load_asymmetry(path, value_col):
    df = pd.read_csv(path, sep="\t", header=None)

    df = df.rename(columns={
        0: "species",
        1: "pos_count",
        2: "pos",
        3: "neg_count",
        4: "neg",
        5: "difference",
        6: value_col,
    })

    # Drop strand sign columns
    df = df.drop(columns=["pos", "neg"])

    # Normalize species names
    df["species"] = df["species"].apply(normalize_species)

    # Keep only required columns
    df = df[["species", value_col]]

    # Reindex safely (NaN if species absent)
    return (
        df.set_index("species")
          .reindex(ORDER)
          .reset_index()
    )


# ----------------------------------
# LOAD ALL DATASETS
# ----------------------------------

cds_df = load_asymmetry(
    "data/strand_asymmetry/all_species_cds_strand_count_asymmetry.tsv",
    "cds_strand_asymmetry",
)

strand_count_df = load_asymmetry(
    "data/strand_asymmetry/all_species_scfr_strand_count_asymmetry.tsv",
    "scfr_strandcount_asymmetry",
)

composite_exon_df = load_asymmetry(
    "data/strand_asymmetry/all_species_composite_exon_shadow_strand_count_asymmetry.tsv",
    "comp_exon_strand_asymmetry",
)

multi_exon_df = load_asymmetry(
    "data/strand_asymmetry/all_species_multi_exon_shadow_strand_count_asymmetry.tsv",
    "multi_exon_strand_asymmetry",
)

single_exon_df = load_asymmetry(
    "data/strand_asymmetry/all_species_single_exon_shadow_strand_count_asymmetry.tsv",
    "single_exon_strand_asymmetry",
)

all_species_df = load_asymmetry(
    "data/strand_asymmetry/all_species_shadow_strand_count_asymmetry.tsv", 
    "all_species_shadow_asymmetry"
)

# ----------------------------------
# FINAL MERGE
# ----------------------------------

final_df = (
    pd.concat(
        [
            cds_df,
            strand_count_df,
            composite_exon_df,
            multi_exon_df,
            single_exon_df,
            all_species_df
        ],
        axis=1
    )
    .loc[:, ~pd.concat(
        [
            cds_df,
            strand_count_df,
            composite_exon_df,
            multi_exon_df,
            single_exon_df,
            all_species_df
        ],
        axis=1
    ).columns.duplicated()]
)

final_df.to_csv("strand_asymmetry_gregions/intermediate_asymmetry_final_data.csv")