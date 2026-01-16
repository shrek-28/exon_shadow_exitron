import pandas as pd 

bonobo = pd.read_csv("data/exon_shadow/bonobo_composite_exon_filtered.tsv", sep="\t")
borangutan = pd.read_csv("data/exon_shadow/borangutan_composite_exon_filtered.tsv", sep="\t")
chimpanzee = pd.read_csv("data/exon_shadow/chimpanzee_composite_exon_filtered.tsv", sep="\t")
gibbon = pd.read_csv("data/exon_shadow/gibbon_composite_exon_filtered.tsv", sep="\t")
gorilla = pd.read_csv("data/exon_shadow/gorilla_composite_exon_filtered.tsv", sep="\t")
human = pd.read_csv("data/exon_shadow/human_composite_exon_filtered.tsv", sep="\t")
sorangutan = pd.read_csv("data/exon_shadow/sorangutan_composite_exon_filtered.tsv", sep="\t")

bonobo = bonobo[['exon_strand']]
borangutan = borangutan[['exon_strand']]
chimpanzee = chimpanzee[['exon_strand']]
gibbon = gibbon[['exon_strand']]
gorilla = gorilla[['exon_strand']]
human = human[['exon_strand']]
sorangutan = sorangutan[['exon_strand']]

bonobo.rename(columns={'exon_strand': 'Bonobo'}, inplace=True)
borangutan.rename(columns={'exon_strand': 'Bornean Orangutan'}, inplace=True)
chimpanzee.rename(columns={'exon_strand': 'Chimpanzee'}, inplace=True)
gibbon.rename(columns={'exon_strand': 'Gibbon'}, inplace=True)
gorilla.rename(columns={'exon_strand': 'Gorilla'}, inplace=True)
human.rename(columns={'exon_strand': 'Human'}, inplace=True)
sorangutan.rename(columns={'exon_strand': 'Sumatran Orangutan'}, inplace=True)

df = pd.concat([gibbon, gorilla, human, bonobo, chimpanzee, borangutan, sorangutan], axis=1)

df_long = df.melt(var_name='species', value_name='strand')
strand_counts = df_long.groupby(['species', 'strand']).size().reset_index(name='count')

species_order = [
    "Gibbon", 
    "Gorilla", 
    "Human",
    "Bonobo",
    "Chimpanzee",
    "Bornean Orangutan",
    "Sumatran Orangutan"
]

strand_counts["species"] = pd.Categorical(
    strand_counts["species"],
    categories=species_order,
    ordered=True
)

strand_counts = strand_counts.sort_values("species")

print(strand_counts)
strand_counts.to_csv("strand_asymmetry_distribution/summary_table_composite_exon.csv")