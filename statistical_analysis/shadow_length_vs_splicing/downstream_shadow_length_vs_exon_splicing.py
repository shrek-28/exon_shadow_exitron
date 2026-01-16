import pandas as pd
import glob
import os
import re
import itertools
from scipy.stats import kruskal, mannwhitneyu
from statsmodels.stats.multitest import multipletests


# -----------------------------
# Settings
# -----------------------------
folder = "data/exon_shadow/filtered_combined_non_zero"
files = glob.glob(os.path.join(folder, "*.tsv"))

valid_categories = ['alternative', 'constitutive', 'unique']

summary_results = []
posthoc_results = []

# -----------------------------
# Loop over files
# -----------------------------
for file in files:
    df = pd.read_csv(file, sep="\t", low_memory=False)

    # Extract species name
    species = os.path.basename(file).split("_")[0]

    # Keep relevant columns
    df = df[['last_exon_splicing', 'downstream_len_in_scfr']]

    print(species)
    print("df:", len(df))

    df = df[df['downstream_len_in_scfr'] != 0]
    print("after dropping zeros:", len(df))
  
    n_total = len(df)

    # Group data
    groups = {
        cat: df[df['last_exon_splicing'] == cat]['downstream_len_in_scfr'].values
        for cat in valid_categories
    }

    # Skip if fewer than 2 non-empty groups
    non_empty = [v for v in groups.values() if len(v) > 0]
    if len(non_empty) < 2:
        continue

    # -----------------------------
    # Kruskal–Wallis
    # -----------------------------
    H, p_kw = kruskal(*non_empty)

    k = len(non_empty)
    epsilon2 = (H - k + 1) / (n_total - k)

    summary_results.append({
        "species": species,
        "n_total": n_total,
        "Kruskal_H": H,
        "Kruskal_p": p_kw,
        "epsilon2": epsilon2
    })

    # -----------------------------
    # Pairwise Mann–Whitney U
    # -----------------------------
    for g1, g2 in itertools.combinations(valid_categories, 2):
        if len(groups[g1]) == 0 or len(groups[g2]) == 0:
            continue

        U, p_u = mannwhitneyu(
            groups[g1],
            groups[g2],
            alternative='two-sided'
        )

        n1, n2 = len(groups[g1]), len(groups[g2])

        r_rb = 1 - (2*U)/(n1*n2)

        posthoc_results.append({
            "species": species,
            "group1": g1,
            "group2": g2,
            "U_stat": U,
            "p_value": p_u, 
            "rank_biserial_r": r_rb
        })

# -----------------------------
# Save results
# -----------------------------
summary_df = pd.DataFrame(summary_results)

valid_mask = summary_df["Kruskal_p"].notna()

summary_df.loc[valid_mask, "Kruskal_corrected"] = multipletests(
    summary_df.loc[valid_mask, "Kruskal_p"],
    method="fdr_bh"
)[1]

posthoc_df = pd.DataFrame(posthoc_results)

summary_df.to_csv(
    "statistical_analysis/shadow_length_vs_splicing/results/last_exon_splicing_upstream_shadow_summary.csv", index=False
)
posthoc_df.to_csv(
    "statistical_analysis/shadow_length_vs_splicing/results/last_exon_splicing_upstream_posthoc_tests.csv", index=False
)

print("summary df")
print(summary_df)

print("posthoc df")
print(posthoc_df)