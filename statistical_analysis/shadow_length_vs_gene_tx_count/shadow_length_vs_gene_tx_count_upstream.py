import pandas as pd
import glob
import os
import re
from scipy.stats import kruskal, spearmanr
from statsmodels.stats.multitest import multipletests


# -----------------------------
# Settings
# -----------------------------
folder = "data/exon_shadow/filtered_combined_non_zero"
files = glob.glob(os.path.join(folder, "*_filtered_combined_non_zero.tsv"))

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
    df = df[['gene_tx_count', 'upstream_len_in_scfr']]

    print(species)
    print("df:", len(df))

    df = df[df['upstream_len_in_scfr'] != 0]
    print("after dropping zeros:", len(df))
  
    n_total = len(df)

    # Skip species for tests only if <2 gene_tx_count groups
    if df['gene_tx_count'].nunique() < 2:
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
    # Group data by gene_tx_count
    # -----------------------------
    groups = [group['upstream_len_in_scfr'].values for _, group in df.groupby('gene_tx_count')]

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
results_df.to_csv("statistical_analysis/shadow_length_vs_gene_tx_count/results/upstream_exon_shadow_vs_gene_tx_count.csv", index=False)
