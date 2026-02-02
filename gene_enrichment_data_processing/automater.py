import pandas as pd 

def gene_enrichment_processor(files, final_path):
    dfs = []

    for species, path in files.items():
        df = pd.read_csv(path)

        if "Pathway" not in df.columns:
            raise KeyError(f"'Pathway' column missing in {path}")

        # remove everything before and including first space
        df["Pathway"] = df["Pathway"].str.replace(
            r"^[^ ]+\s+", "", regex=True
        )

        dfs.append(df[["Pathway"]].rename(columns={"Pathway": species}))

    final_df = pd.concat(dfs, axis=1)
    final_df.to_csv(final_path, index=False)

## FOR GENES WITH NO SHADOW - ALL GENES AS BACKGROUND 
genes_no_shadow_allgenes_bg = {
    "bonobo": "data/exon_shadow/gene_enrichment/bonobo_gene_enrichment/genes_without_shadow_nobg_enrichment.csv",
    "chimpanzee": "data/exon_shadow/gene_enrichment/chimpanzee_gene_enrichment/genes_without_shadow_no_bg_enrichment.csv",
    "gorilla": "data/exon_shadow/gene_enrichment/gorillla_gene_enrichment/genes_without_shadow_nobg_enrichment.csv",
    "borangutan": "data/exon_shadow/gene_enrichment/borangutan_gene_enrichment/genes_without_shadow_no_bg_enrichment.csv",
    "human": "data/exon_shadow/gene_enrichment/human_gene_enrichment/genes_without_shadow_no_bg_enrichment.csv",
    # "gibbon": "PATH/TO/gibbon.csv",
    "sorangutan": "data/exon_shadow/gene_enrichment/sorangutan_gene_enrichment/genes_without_shadow_no_bg_enrichment.csv",
}

genes_no_shadow_allgenes_bg_path = r"gene_enrichment_data_processing/results/combined_list_genes_no_shadow_allgenes_bg.csv"

## GENES WITH NO SHADOW - NO GENES AS BACKGROUND 
genes_no_shadow_nobg_files = {
    "bonobo": "data/exon_shadow/gene_enrichment/bonobo_gene_enrichment/genes_without_shadow_nobg_enrichment.csv",
    "chimpanzee": "data/exon_shadow/gene_enrichment/chimpanzee_gene_enrichment/genes_without_shadow_no_bg_enrichment.csv",
    "gorilla": "data/exon_shadow/gene_enrichment/gorillla_gene_enrichment/genes_without_shadow_nobg_enrichment.csv",
    "borangutan": "data/exon_shadow/gene_enrichment/borangutan_gene_enrichment/genes_without_shadow_no_bg_enrichment.csv",
    "human": "data/exon_shadow/gene_enrichment/human_gene_enrichment/genes_without_shadow_no_bg_enrichment.csv",
    # "gibbon": "PATH/TO/gibbon.csv",
    "sorangutan": "data/exon_shadow/gene_enrichment/sorangutan_gene_enrichment/genes_without_shadow_no_bg_enrichment.csv",
}
genes_no_shadow_nobg_path = "gene_enrichment_data_processing/results/combined_list_genes_without_shadow_no_bg.csv"

## NO DOWNSTREAM SHADOW - NO BACKGROUND GENES
no_ds_shadow_no_bg_files = {'bonobo': 'data/exon_shadow/gene_enrichment/bonobo_gene_enrichment/no_downstream_no_bg_enrichment.csv', 
                            'borangutan': 'data/exon_shadow/gene_enrichment/borangutan_gene_enrichment/no_downstream_no_bg_enrichment.csv', 
                            'chimpanzee': 'data/exon_shadow/gene_enrichment/chimpanzee_gene_enrichment/no_downstream_shadow_no_bg_enrichment.csv', 
                            'gibbon': 'data/exon_shadow/gene_enrichment/gibbon_gene_enrichment/no_downstream_no_bg_enrichment.csv', 
                            'gorillla': 'data/exon_shadow/gene_enrichment/gorillla_gene_enrichment/no_downstream_no_bg_enrichment.csv', 
                            'human': 'data/exon_shadow/gene_enrichment/human_gene_enrichment/no_downstream_no_bg_enrichment.csv', 
                            'sorangutan': 'data/exon_shadow/gene_enrichment/sorangutan_gene_enrichment/no_downstream_no_bg_enrichment.csv'}
no_ds_shadow_no_bg_path = "gene_enrichment_data_processing/results/no_ds_shadow_no_bg.csv"

## NO UPSTREAM SHADOW ALL GENES BACKGROUND ENRICHMENT 
no_us_shadow_allgenes_bg_files = {
    'bonobo': 'data/exon_shadow/gene_enrichment/bonobo_gene_enrichment/no_upstream_allgenes_bg_enrichment.csv', 
    'borangutan': 'data/exon_shadow/gene_enrichment/borangutan_gene_enrichment/no_upstream_allgenes_bg_enrichment.csv', 
    'chimpanzee': 'data/exon_shadow/gene_enrichment/chimpanzee_gene_enrichment/no_upstream_allgenes_bg_enrichment.csv', 
    'gibbon': 'data/exon_shadow/gene_enrichment/gibbon_gene_enrichment/no_upstream_allgenes_bg_enrichment.csv', 
    'gorillla': 'data/exon_shadow/gene_enrichment/gorillla_gene_enrichment/no_upstream_allgenes_bg_enrichment.csv', 
    'human': 'data/exon_shadow/gene_enrichment/human_gene_enrichment/no_upstream_allgenes_bg_enrichment.csv', 
    'sorangutan': 'data/exon_shadow/gene_enrichment/sorangutan_gene_enrichment/no_upstream_allgenes_bg_enrichment.csv'
    } 
no_us_shadow_allgenes_bg_path = "gene_enrichment_data_processing/results/no_us_shadow_allgenes_bg.csv"

## NO UPSTREAM SHADOW NO BACKGROUND 
no_us_shadow_nobg_files = {
    'bonobo': 'data/exon_shadow/gene_enrichment/bonobo_gene_enrichment/no_upstream_no_bg_enrichment.csv', 
    'borangutan': 'data/exon_shadow/gene_enrichment/borangutan_gene_enrichment/no_upstream_no_bg_enrichment.csv', 
    'chimpanzee': 'data/exon_shadow/gene_enrichment/chimpanzee_gene_enrichment/no_upstream_no_bg_enrichment.csv', 
    'gibbon': 'data/exon_shadow/gene_enrichment/gibbon_gene_enrichment/no_upstream_no_bg_enrichment.csv', 
    'gorillla': 'data/exon_shadow/gene_enrichment/gorillla_gene_enrichment/no_upstream_no_bg_enrichment.csv', 
    'human': 'data/exon_shadow/gene_enrichment/human_gene_enrichment/no_upstream_no_bg_enrichment.csv', 
    'sorangutan': 'data/exon_shadow/gene_enrichment/sorangutan_gene_enrichment/no_upstream_no_bg_enrichment.csv'
}
no_us_shadow_nobg_path = "gene_enrichment_data_processing/results/no_us_shadow_no_bg.csv"

## NO UPSTREAM SHADOW - ALL GENES WITH SHADOW BACKGROUND
no_us_shadow_geneswshadow_bg_files = {
    'bonobo': 'data/exon_shadow/gene_enrichment/bonobo_gene_enrichment/no_upstream_geneswshadow_enrichment.csv', 
    'borangutan': 'data/exon_shadow/gene_enrichment/borangutan_gene_enrichment/no_upstream_geneswshadow_enrichment.csv', 
    'chimpanzee': 'data/exon_shadow/gene_enrichment/chimpanzee_gene_enrichment/no_upstream_geneswshadow_enrichment.csv', 
    'gibbon': 'data/exon_shadow/gene_enrichment/gibbon_gene_enrichment/no_upstream_geneswshadow_enrichment.csv', 
    'gorilla': 'data/exon_shadow/gene_enrichment/gorillla_gene_enrichment/no_upstream_geneswshadow_enrichment.csv', 
    'human': 'data/exon_shadow/gene_enrichment/human_gene_enrichment/no_upstream_geneswshadow_enrichment.csv', 
    'sorangutan': 'data/exon_shadow/gene_enrichment/sorangutan_gene_enrichment/no_upstream_geneswshadow_enrichment.csv'
}
no_us_shadow_geneswshadow_bg_path = "gene_enrichment_data_processing/results/no_us_shadow_geneswshadow_bg.csv"

## SYMMETRIC SHADOW ALL GENES BG ENRICHMENT 
sym_shadow_allgenes_bg_files = {
    'bonobo': 'data/exon_shadow/gene_enrichment/bonobo_gene_enrichment/symmetric_allgenes_bg_enrichment.csv', 
    'borangutan': 'data/exon_shadow/gene_enrichment/borangutan_gene_enrichment/symmetric_allgenes_bg_enrichment.csv', 
    'chimpanzee': 'data/exon_shadow/gene_enrichment/chimpanzee_gene_enrichment/symmetric_allgenes_bg_enrichment.csv', 
    'gibbon': 'data/exon_shadow/gene_enrichment/gibbon_gene_enrichment/symmetric_allgenes_bg_enrichment.csv', 
    'human': 'data/exon_shadow/gene_enrichment/human_gene_enrichment/symmetric_allgenes_bg_enrichment.csv', 
    'sorangutan': 'data/exon_shadow/gene_enrichment/sorangutan_gene_enrichment/symmetric_allgenes_bg_enrichment.csv'
}
sym_shadow_allgenes_path = "gene_enrichment_data_processing/results/sym_shadow_allgenes_bg.csv"

## SYMMETRIC SHADOWS NO BACKGROUND ENRICHMENT 
sym_shadow_nobg_files = {
    'bonobo': 'data/exon_shadow/gene_enrichment/bonobo_gene_enrichment/symmetric_no_bg_enrichment.csv', 
    'borangutan': 'data/exon_shadow/gene_enrichment/borangutan_gene_enrichment/symmetric_no_bg_enrichment.csv', 
    'chimpanzee': 'data/exon_shadow/gene_enrichment/chimpanzee_gene_enrichment/symmetric_no_bg_enrichment.csv', 
    'gibbon': 'data/exon_shadow/gene_enrichment/gibbon_gene_enrichment/symmetric_no_bg_enrichment.csv', 
    'human': 'data/exon_shadow/gene_enrichment/human_gene_enrichment/symmetric_no_bg_enrichment.csv', 
    'sorangutan': 'data/exon_shadow/gene_enrichment/sorangutan_gene_enrichment/symmetric_no_bg_enrichment.csv',
    'gorilla': 'data\exon_shadow/gene_enrichment/gorillla_gene_enrichment/shadow_symmetric_no_bg_enrichment.csv'
}
sym_shadow_nobg_path = "gene_enrichment_data_processing/results/sym_shadow_nobg.csv"

## SYMMETRIC SHADOW ALL GENES WITH SHADOW ENRICHMENT 
sym_shadow_geneswshadow_files = {
    'bonobo': 'data/exon_shadow/gene_enrichment/bonobo_gene_enrichment/symmetric_geneswshadow_bg_enrichment.csv', 
    'borangutan': 'data/exon_shadow/gene_enrichment/borangutan_gene_enrichment/symmetric_geneswshadow_bg_enrichment.csv', 
    'chimpanzee': 'data/exon_shadow/gene_enrichment/chimpanzee_gene_enrichment/symmetric_geneswshadow_bg_enrichment.csv', 
    'gibbon': 'data/exon_shadow/gene_enrichment/gibbon_gene_enrichment/symmetric_geneswshadow_bg_enrichment.csv', 
    'human': 'data/exon_shadow/gene_enrichment/human_gene_enrichment/symmetric_geneswshadow_bg_enrichment.csv'
}
## not there in gorilla and sumatran orangutan 
sym_shadow_geneswshadow_path = "gene_enrichment_data_processing/results/sym_shadow_geneswshadow_bg.csv"

## GENES WITH SHADOW/ NO BACKGROUND SPECIFIED
geneswshadow_no_bg_files = {
    'bonobo': 'data/exon_shadow/gene_enrichment/bonobo_gene_enrichment/geneswshadow_nobg_enrichment.csv', 
    'borangutan': 'data/exon_shadow/gene_enrichment/borangutan_gene_enrichment/geneswshadow_nobg_enrichment.csv', 
    'chimpanzee': 'data/exon_shadow/gene_enrichment/chimpanzee_gene_enrichment/geneswshadow_no_bg_enrichment.csv', 
    'gibbon': 'data/exon_shadow/gene_enrichment/gibbon_gene_enrichment/geneswshadow_no_bg_enrichment.csv', 
    'gorilla': 'data/exon_shadow/gene_enrichment/gorillla_gene_enrichment/geneswshadow_no_bg_enrichment.csv', 
    'human': 'data/exon_shadow/gene_enrichment/human_gene_enrichment/geneswshadow_nobg_enrichment.csv', 
    'sorangutan': 'data/exon_shadow/gene_enrichment/sorangutan_gene_enrichment/geneswshadow_no_bg_enrichment.csv'
}
geneswshadow_no_bg_path = "gene_enrichment_data_processing/results/geneswsahdow_no_bg.csv"

#----------------------------------------------------------------------------------------
#----------------------------------------------------------------------------------------
## FUNCTION CALLS

## genes with no shadow: all genes as background 
gene_enrichment_processor(files=genes_no_shadow_allgenes_bg, final_path=genes_no_shadow_allgenes_bg_path)

## genes with no shadow: no background (default background used)
gene_enrichment_processor(files=genes_no_shadow_nobg_files, final_path=genes_no_shadow_nobg_path)

## no downstream shadow: no background (default background)
gene_enrichment_processor(files=no_ds_shadow_no_bg_files, final_path=no_ds_shadow_no_bg_path)

## no upstream shadow: all genes as background 
gene_enrichment_processor(files=no_us_shadow_allgenes_bg_files, final_path=no_us_shadow_allgenes_bg_path)

## no upstream shadow: no background (default background)
gene_enrichment_processor(files=no_us_shadow_nobg_files, final_path=no_us_shadow_nobg_path)

## no upstream shadow: all genes with shadow as background 
gene_enrichment_processor(files=no_us_shadow_geneswshadow_bg_files, final_path=no_us_shadow_geneswshadow_bg_path)

## symmetric shadow: all genes as background 
gene_enrichment_processor(files=sym_shadow_allgenes_bg_files, final_path=sym_shadow_allgenes_path)

## symmetric shadow: no background 
gene_enrichment_processor(files=sym_shadow_nobg_files, final_path=sym_shadow_nobg_path)

## symmetric shadow: all genes with shadow as background
gene_enrichment_processor(files=sym_shadow_geneswshadow_files, final_path=sym_shadow_geneswshadow_path)

## genes with shadow enrichment with no background 
gene_enrichment_processor(files=geneswshadow_no_bg_files, final_path=geneswshadow_no_bg_path)