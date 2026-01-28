import pandas as pd
import glob
import os
import re
from scipy.stats import kruskal, spearmanr
from statsmodels.stats.multitest import multipletests


# -----------------------------
# Settings
# -----------------------------
folder = "data/exon_shadow/individual"
files = glob.glob(os.path.join(folder, "*_multi_exon_filtered.tsv"))

results_summary = []

# -----------------------------
# Loop over species files
# -----------------------------
for file in files:
    # Read TSV
    df = pd.read_csv(file, sep="\t", low_memory=False)

    # Extract species name from filepath
    m = re.search(r"\\([^_]+)", file)
    species = m.group(1) if m else os.path.basename(file).split("_")[0]

    # Keep relevant columns and convert to numeric
    df = df[['last_exon_tx_count', 'downstream_len_in_scfr']]

    print(species)
    print("df:", len(df))

    df = df[df['downstream_len_in_scfr'] != 0]
    print("after dropping zeros:", len(df))
  
    n_total = len(df)

    # Skip species for tests only if <2 last_exon_tx_count groups
    if df['last_exon_tx_count'].nunique() < 2:
        results_summary.append({
            "species": species,
            "n_total": n_total,
            "Kruskal_H": float('nan'),
            "Kruskal_p": float('nan'),
            "Kruskal_epsilon2": float('nan'),
            "ANOVA_F": float('nan'),
            "ANOVA_p": float('nan'),
        })
        continue

    # -----------------------------
    # Group data by last_exon_tx_count
    # -----------------------------
    groups = [group['downstream_len_in_scfr'].values for _, group in df.groupby('last_exon_tx_count')]

    # -----------------------------
    # Kruskal-Wallis H test
    # -----------------------------
    kruskal_stat, kruskal_p = kruskal(*groups)

    # Epsilon-squared effect size for Kruskal-Wallis
    k = len(groups)
    kruskal_epsilon2 = (kruskal_stat - k + 1) / (n_total - k)


    # -----------------------------
    # Save results
    # -----------------------------
    results_summary.append({
        "species": species,
        "n_total": n_total,
        "Kruskal_H": kruskal_stat,
        "Kruskal_p": kruskal_p,
        "Kruskal_epsilon2": kruskal_epsilon2,
    })

# -----------------------------
# Convert results to DataFrame
# -----------------------------
results_df = pd.DataFrame(results_summary)

valid_mask = results_df["Kruskal_p"].notna()

results_df.loc[valid_mask, "Kruskal_corrected"] = multipletests(
    results_df.loc[valid_mask, "Kruskal_p"],
    method="fdr_bh"
)[1]

print(results_df)

# Optionally save to CSV
# results_df.to_csv("statistical_analysis/shadow_length_vs_transcript_count/results/downstream_exon_shadow_vs_last_exon_tx_count.csv", index=False)
