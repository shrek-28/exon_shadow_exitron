import pandas as pd
import numpy as np
import sys

# -----------------------------
# Command-line arguments
# -----------------------------
if len(sys.argv) != 3:
    print("Usage: python add_metric.py <input_file> <output_file>")
    sys.exit(1)

input_file = sys.argv[1]
output_file = sys.argv[2]

# -----------------------------
# Load file
# -----------------------------
df = pd.read_csv(input_file, sep="\t")

# -----------------------------
# Compute row-level metric
# -----------------------------
def compute_metric(freq_str, mag_str):
    try:
        freqs = np.array([float(x) for x in str(freq_str).split(";")])
        mags  = np.array([float(x) for x in str(mag_str).split(";")])
        if len(freqs) != len(mags) or mags.sum() == 0:
            return np.nan
        return np.sum(freqs * mags) / np.sum(mags)
    except:
        return np.nan

df['Metric'] = df.apply(lambda row: compute_metric(row['Frequencies'], row['Magnitudes']), axis=1)

# -----------------------------
# Save output
# -----------------------------
df.to_csv(output_file, sep="\t", index=False)
print(f"Saved with metric column to: {output_file}")
