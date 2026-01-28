import pandas as pd
import glob
import os
import re
import numpy as np
from scipy.stats import pearsonr, spearmanr, kendalltau
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
    df = pd.read_csv(file, sep="\t", low_memory=False)

    # Extract species
    m = re.search(r"\\([^_]+)", file)
    species = m.group(1) if m else os.path.basename(file).split("_")[0]

    df = df[['last_exon_sharing', 'downstream_len_in_scfr']].apply(
        pd.to_numeric, errors='coerce'
    ).dropna()

    df = df[df['downstream_len_in_scfr'] != 0]

    n_total = len(df)

    print(species, "n =", n_total)

    # Need at least 3 points for correlations to not be garbage
    if n_total < 3:
        results_summary.append({
            "species": species,
            "n_total": n_total,
            "pearson_r": np.nan,
            "pearson_p": np.nan,
            "spearman_r": np.nan,
            "spearman_p": np.nan,
            "kendall_tau": np.nan,
            "kendall_p": np.nan
        })
        continue

    x = df['last_exon_sharing'].values
    y = df['downstream_len_in_scfr'].values

    # -----------------------------
    # Correlation tests
    # -----------------------------
    pearson_r, pearson_p = pearsonr(x, y)
    spearman_r, spearman_p = spearmanr(x, y)
    kendall_tau, kendall_p = kendalltau(x, y)

    results_summary.append({
        "species": species,
        "n_total": n_total,
        "pearson_r": pearson_r,
        "pearson_p": pearson_p,
        "spearman_r": spearman_r,
        "spearman_p": spearman_p,
        "kendall_tau": kendall_tau,
        "kendall_p": kendall_p
    })

# -----------------------------
# Results table
# -----------------------------
results_df = pd.DataFrame(results_summary)

# -----------------------------
# FDR correction (separately per test)
# -----------------------------
for test in ["pearson", "spearman", "kendall"]:
    p_col = f"{test}_p"
    valid = results_df[p_col].notna()

    results_df.loc[valid, f"{test}_fdr"] = multipletests(
        results_df.loc[valid, p_col],
        method="fdr_bh"
    )[1]

print(results_df)

results_df.to_csv(
    "statistical_analysis/exon_sharing_vs_shadow_length/results/"
    "downstream_exon_shadow_vs_last_exon_sharing_correlations.csv",
    index=False
)
