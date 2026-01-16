# Clear workspace
rm(list = ls())

# Load required libraries
library(ggplot2)
library(tidyr)
library(dplyr)
library(stringr)

# ---- Read data ----
df <- read.csv(
  "C:\\Users\\admin\\Documents\\internship\\strand_asymmetry_gregions\\intermediate_asymmetry_final_data.csv",
  stringsAsFactors = FALSE, 
  row.names = 1
)

# Wrap species labels
df$species <- str_wrap(df$species, width = 10)

# ---- Convert to long format ----
df_long <- df %>%
  pivot_longer(
    cols = -species,
    names_to = "region",
    values_to = "strand_asymmetry"
  )

# Ensure numeric
df_long$strand_asymmetry <- as.numeric(df_long$strand_asymmetry)

# Preserve species order
df_long$species <- factor(df_long$species, levels = unique(df$species))

regions_to_plot <- c(
  "cds_strand_asymmetry",
  "scfr_strandcount_asymmetry",
  "all_species_shadow_asymmetry",
  "single_exon_strand_asymmetry",
  "multi_exon_strand_asymmetry",
  "comp_exon_strand_asymmetry"
)

df_long_sub <- df_long %>%
  filter(region %in% regions_to_plot)

# ---- Plot ----
plot <- ggplot(
  df_long_sub,
  aes(
    x = species,
    y = strand_asymmetry,
    color = region,
    group = region
  )
) +
  geom_line(size = 0.75) +
  geom_point(size = 2) +
  scale_y_continuous(limits = c(-0.25, 0.25)) +
  scale_color_manual(
    values = c(
      cds_strand_asymmetry = "#1b9e77",
      scfr_strandcount_asymmetry = "#d95f02",
      comp_exon_strand_asymmetry = "#7570b3",
      multi_exon_strand_asymmetry = "#e7298a",
      single_exon_strand_asymmetry = "#66a61e",
      all_species_shadow_asymmetry = "lightgreen"
    ),
    labels = c(
      cds_strand_asymmetry = "CDS",
      scfr_strandcount_asymmetry = "SCFR",
      comp_exon_strand_asymmetry = "Composite Exon Shadow",
      multi_exon_strand_asymmetry = "Multi-exon Shadow",
      single_exon_strand_asymmetry = "Single-exon Shadow",
      all_species_shadow_asymmetry = "Complete Shadow"
    )
  ) +
  theme_minimal() +
  theme(
    legend.title = element_blank(),
    axis.text.x = element_text(hjust = 0.5)
  ) +
  labs(
    title = "Strand Asymmetry Across Species and Genomic Regions",
    x = "Species",
    y = "Strand Asymmetry"
  )

plot

ggsave("C:\\Users\\admin\\Documents\\internship\\strand_asymmetry_gregions\\strand_asymmetry_plot.pdf", width=10, height=6, dpi=300)