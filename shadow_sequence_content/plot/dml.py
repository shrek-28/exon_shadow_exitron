import pandas as pd
import os

# -----------------------------
# PATHS
# -----------------------------
BASE_INPUT_DIR = r"C:\Users\admin\Desktop\Coding\exon_shadow_exitron\shadow_sequence_content\output_files"
BASE_OUTPUT_DIR = r"C:\Users\admin\Desktop\Coding\exon_shadow_exitron\shadow_sequence_content\plot\inter_data"

os.makedirs(BASE_OUTPUT_DIR, exist_ok=True)

SPECIES_ORDER = [
    "gibbon",
    "gorilla",
    "human",
    "bonobo",
    "chimpanzee",
    "borangutan",
    "sorangutan"
]

SPECIES_LABELS = {
    "gibbon": "Gibbon",
    "gorilla": "Gorilla",
    "human": "Human",
    "bonobo": "Bonobo",
    "chimpanzee": "Chimpanzee",
    "borangutan": "Bornean Orangutan",
    "sorangutan": "Sumatran Orangutan"
}

def summarize_upstream_downstream():
    rows = []

    for sp in SPECIES_ORDER:
        for region in ["upstream", "downstream"]:
            fpath = os.path.join(BASE_INPUT_DIR, f"{sp}_{region}.csv")
            df = pd.read_csv(fpath, usecols=["GC_percent"])
            
            stats = {
                "species": SPECIES_LABELS[sp],
                "region": region.capitalize(),  # For plot fill/color
                "min": df["GC_percent"].min(),
                "q1": df["GC_percent"].quantile(0.25),
                "median": df["GC_percent"].median(),
                "q3": df["GC_percent"].quantile(0.75),
                "max": df["GC_percent"].max(),
                "mean": df["GC_percent"].mean(),
                "count": len(df['GC_percent'])
            }

            rows.append(stats)

    out = pd.DataFrame(rows)
    
    # Maintain fixed species order
    out["species"] = pd.Categorical(
        out["species"],
        categories=[SPECIES_LABELS[s] for s in SPECIES_ORDER],
        ordered=True
    )

    # Save combined CSV for side-by-side boxplot
    out.to_csv(
        os.path.join(BASE_OUTPUT_DIR, "boxplot_stats.csv"),
        index=False
    )

if __name__ == "__main__":
    summarize_upstream_downstream()