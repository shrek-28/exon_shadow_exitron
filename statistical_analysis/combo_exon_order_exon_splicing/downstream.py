import os 
import glob
import pandas as pd 
import statsmodels.formula.api as smf
import statsmodels.api as sm

folder = "data/exon_shadow/filtered_combined_non_zero"
files = glob.glob(os.path.join(folder, "*.tsv"))

anova_results = []

for file in files:
    df = pd.read_csv(file, sep="\t", low_memory=False)

    # Extract species name
    species = os.path.basename(file).split("_")[0]

    df = df[['downstream_len_in_scfr', 'last_exon_splicing', 'last_exon_order']]

    print(species)
    print("df:", len(df))
    
    df = df[df['downstream_len_in_scfr'] != 0]
    print("post-removal:", len(df))

    n_total = len(df)

    # rank transform
    df["rank_upstream"] = df["downstream_len_in_scfr"].rank(method="average")

    model = smf.ols(
        "rank_upstream ~ C(last_exon_splicing) * C(last_exon_order)",
        data=df
    ).fit()

    table = sm.stats.anova_lm(model, typ=2)

    table = table.reset_index()
    table["species"] = species
    table["n_total"] = n_total

    anova_results.append(table)

anova_df = pd.concat(anova_results, ignore_index=True)
anova_df.to_csv("statistical_analysis/combo_exon_order_exon_splicing/results/downstream_shadow_exon_order_splicing.csv")
print(anova_df)