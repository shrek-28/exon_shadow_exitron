import pandas as pd
import numpy as np
import glob
import os

# Folder containing the files
folder = "C:/Users/admin/Desktop/Coding/exon_shadow_exitron/exitron_length_dist/raw_exitron_data"  # replace with your folder path

files = glob.glob(os.path.join(folder, "*_exitron_candidates_filtered.tsv"))

log_rows = []
orig_stats = []

for file in files:
    species = os.path.basename(file).split("_exitron_candidates_filtered.tsv")[0]
    df = pd.read_csv(file, sep="\t")
    
    # Detect intron_length column case-insensitively
    intron_col = [col for col in df.columns if "intron_length" in col.lower()]
    if not intron_col:
        continue
    
    intron_vals = df[intron_col[0]].dropna()
    if len(intron_vals) == 0:
        continue
    
    # Log-transform intron lengths for plotting
    for val in intron_vals:
        log_rows.append({
            "species": species,
            "log_intron_length": np.log10(val)  # log-transform
        })
    
    # Original statistics for annotation
    orig_stats.append({
        "species": species,
        "orig_mean": intron_vals.mean(),
        "orig_median": intron_vals.median(),
        "orig_min": intron_vals.min(),
        "orig_max": intron_vals.max()
    })

# Create separate dataframes
log_data = pd.DataFrame(log_rows)
stats_data = pd.DataFrame(orig_stats)

# Save CSVs
log_data.to_csv("log_intron_lengths.csv", index=False)  # for boxplot
stats_data.to_csv("intron_stats.csv", index=False)      # for annotations

print("CSV files created:")
print("- log_intron_lengths.csv (log-transformed lengths per species)")
print("- intron_stats.csv (original statistics per species)")