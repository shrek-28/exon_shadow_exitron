import pandas as pd
import numpy as np
import glob
import os
import sys

# -----------------------------
# Command-line arguments
# -----------------------------
if len(sys.argv) != 3:
    print("Usage: python compile_metrics.py <input_folder> <output_folder>")
    sys.exit(1)

input_folder = sys.argv[1]
output_folder = sys.argv[2]

if not os.path.exists(output_folder):
    os.makedirs(output_folder)

# -----------------------------
# Subfolder mapping
# -----------------------------
subfolders = {
    "upstream": "_upstream_with_metric.tsv",
    "downstream": "_downstream_with_metric.tsv",
    "exitron": "_exitron_with_metric.tsv",
    "SCFR": "_scfr_with_metric.tsv"
}

# -----------------------------
# Collect metrics
# -----------------------------
records = []

for seq_type, suffix in subfolders.items():
    folder_path = os.path.join(input_folder, seq_type.lower())
    if not os.path.exists(folder_path):
        print(f"Warning: folder {folder_path} does not exist, skipping.")
        continue
    
    files = glob.glob(os.path.join(folder_path, f"*{suffix}"))
    
    for f in files:
        species = os.path.basename(f).split("_")[0]  # species name as first part
        df = pd.read_csv(f, sep="\t")
        if 'Metric' not in df.columns:
            continue
        
        metrics = pd.to_numeric(df['Metric'], errors='coerce').dropna()
        if len(metrics) == 0:
            continue
        
        records.append({
            'species': species,
            'seq_type': seq_type,
            'metrics': metrics
        })

# -----------------------------
# Compute statistics
# -----------------------------
stats_list = []
for rec in records:
    metrics = np.array(rec['metrics'])
    stats_list.append({
        'species': rec['species'],
        'seq_type': rec['seq_type'],
        'N': len(metrics),
        'min': metrics.min(),
        'max': metrics.max(),
        'mean': metrics.mean(),
        'q1': np.percentile(metrics, 25),
        'median': np.median(metrics),
        'q3': np.percentile(metrics, 75)
    })

stats_df = pd.DataFrame(stats_list)

# -----------------------------
# Save statistics
# -----------------------------
stats_file = os.path.join(output_folder, "metrics_statistics.tsv")
stats_df.to_csv(stats_file, sep="\t", index=False)
print(f"Statistics saved to: {stats_file}")
