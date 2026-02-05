# =========================================
# Grouped + Stacked Barplot (Single Panel)
# =========================================

suppressPackageStartupMessages({
  library(tidyverse)
  library(stringr)
})

# -------------------------
# 1. Handle arguments
# -------------------------
if (interactive()) {
  # ---- RStudio usage ----
  INPUT_FILE  <- "C://Users//admin//Desktop//Coding//exon_shadow_exitron//shadow_distributions_regions//no_of_shadows.tsv"
  OUTPUT_FILE <- "C://Users//admin//Desktop//Coding//exon_shadow_exitron//shadow_distributions_regions//output_plot.pdf"
} else {
  # ---- Command-line usage ----
  args <- commandArgs(trailingOnly = TRUE)
  
  if (length(args) != 2) {
    stop("Usage: Rscript plot_code.R <input_file> <output_file>")
  }
  
  INPUT_FILE  <- args[1]
  OUTPUT_FILE <- args[2]
}

# -------------------------
# 2. Load data
# -------------------------
df <- read.delim(INPUT_FILE, stringsAsFactors = FALSE)

# -------------------------
# 3. Species naming + order
# -------------------------
df$Species <- df$Species |>
  str_to_lower() |>
  recode(
    "borangutan" = "Bornean Orangutan",
    "sorangutan" = "Sumatran Orangutan"
  ) |>
  str_to_title()

SPECIES_ORDER <- c(
  "Gibbon",
  "Gorilla",
  "Human",
  "Bonobo",
  "Chimpanzee",
  "Bornean Orangutan",
  "Sumatran Orangutan"
)

df$Species <- factor(df$Species, levels = SPECIES_ORDER)

totals_df <- df |>
  select(
    Species,
    total_genes,
    total_tx,
    total_exons
  ) |>
  pivot_longer(
    -Species,
    names_to = "Feature",
    values_to = "Total"
  ) |>
  mutate(
    Feature = recode(
      Feature,
      "total_genes"        = "Genes",
      "total_tx"  = "Transcripts",
      "total_exons"        = "Exons"
    ),
    Feature = factor(Feature, levels = c("Genes", "Transcripts", "Exons")),
    y = 115   # fixed position ABOVE with-shadow %
  )

# -------------------------
# 4. Reshape data
# -------------------------
# -------------------------
# 5. Plot
# -------------------------
p <- ggplot(
  plot_df,
  aes(x = Feature, y = Percentage, fill = Shadow)
) +
  geom_col(width = 0.7) +
  geom_text(
    data = totals_df,
    aes(
      x = Feature,
      y = y,
      label = scales::comma(Total)
    ),
    inherit.aes = FALSE,
    size = 3,
    fontface = "bold"
  )+ 
  
  # ---- WITHOUT shadow: below bar
  geom_text(
    data = subset(plot_df, Shadow == "Without Shadow"),
    aes(
      y = -3,
      label = sprintf("%.1f", Percentage)
    ),
    size = 3.5
  ) +
  
  # ---- WITH shadow: above bar
  geom_text(
    data = subset(plot_df, Shadow == "With Shadow"),
    aes(
      y = 103,
      label = sprintf("%.1f", Percentage)
    ),
    size = 3.5
  ) +
  
  facet_wrap(
    ~ Species,
    nrow = 1,
    strip.position = "bottom",
    labeller = labeller(Species = label_wrap_gen(width = 0))
  ) +
  
  scale_y_continuous(
    limits = c(-8, 120),
    breaks = c(0, 25, 50, 75, 100),
    expand = c(0, 0)
  ) +
  
  labs(
    title="Percentage of Sequence Types with Shadow",
    x = NULL,
    y = "Percentage (%)",
    fill = "Shadow Status"
  ) +
  
  theme_bw(base_size = 12) +
  theme(
    panel.border = element_blank(),
    strip.background = element_blank(),
    strip.placement = "outside",
    strip.text = element_text(face = "bold", size = 13),
    panel.grid.major.x = element_blank(),
    panel.grid.minor = element_blank(),
    legend.position = "top",
    axis.text.x = element_text(color="black", size=10, angle=45, hjust=1),
    plot.title = element_text(hjust=0, face="bold")
  )

print(p)
ggsave(OUTPUT_FILE, plot=p, width=14, height=6, dpi=300)