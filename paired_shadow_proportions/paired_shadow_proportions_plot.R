# =========================================
# Grouped + Stacked Barplot: Shadow Types
# =========================================

suppressPackageStartupMessages({
  library(tidyverse)
  library(stringr)
})

# -------------------------
# 1. Handle arguments
# -------------------------
if (interactive()) {
  INPUT_FILE  <- "C://Users//admin//Desktop//Coding//exon_shadow_exitron//paired_shadow_proportions//paired_shadow_proportions.tsv"
  OUTPUT_FILE <- "C://Users//admin//Desktop//Coding//exon_shadow_exitron//paired_shadow_proportions//paired_shadow_proportions_plot.pdf"
} else {
  args <- commandArgs(trailingOnly = TRUE)
  if (length(args) != 2) stop("Usage: Rscript plot_code.R <input_file> <output_file>")
  INPUT_FILE  <- args[1]
  OUTPUT_FILE <- args[2]
}

# -------------------------
# 2. Load data
# -------------------------
df <- read.delim(INPUT_FILE, stringsAsFactors = FALSE)

# -------------------------
# 3. Species naming (single controlled path)
# -------------------------
species_map <- c(
  "gibbon"      = "Gibbon",
  "gorilla"     = "Gorilla",
  "human"       = "Human",
  "bonobo"      = "Bonobo",
  "chimpanzee"  = "Chimpanzee",
  "borangutan"  = "Bornean Orangutan",
  "sorangutan"  = "Sumatran Orangutan"
)

df$Species <- df$Species |>
  str_to_lower() |>
  recode(!!!species_map)

SPECIES_ORDER <- unname(species_map)
df$Species <- factor(df$Species, levels = SPECIES_ORDER)

# -------------------------
# 4. Shadow type renaming (authoritative)
# -------------------------
shadow_map <- c(
  pc_only_up              = "Only Upstream",
  pc_only_down            = "Only Downstream",
  pc_both_shadow          = "Both Shadows",
  pc_shadow_with_exitron  = "With Exitron",
  pc_shadow_without_exitron = "Without Exitron"
)

shadow_levels <- c(
  "Only Upstream",
  "Only Downstream",
  "Both Shadows",
  "With Exitron",
  "Without Exitron"
)

# -------------------------
# 5. Build plot dataframe
# -------------------------
plot_df <- df %>%
  select(Species, all_of(names(shadow_map))) %>%
  pivot_longer(
    cols = -Species,
    names_to = "Raw",
    values_to = "Percentage"
  ) %>%
  mutate(
    ShadowType = recode(Raw, !!!shadow_map),
    ShadowType = factor(ShadowType, levels = shadow_levels),
    Bar = if_else(
      ShadowType %in% c("Only Upstream","Only Downstream","Both Shadows"),
      "Total",
      "Paired"
    ),
    Bar = factor(Bar, levels = c("Total","Paired"))
  )

# -------------------------
# 6. Totals above bars
# -------------------------
totals_df <- df %>%
  transmute(
    Species,
    Total  = total_shadows,
    Paired = both_shadow
  ) %>%
  pivot_longer(
    -Species,
    names_to = "Bar",
    values_to = "Total"
  ) %>%
  mutate(
    Bar = factor(Bar, levels = c("Total","Paired")),
    y = 110
  )


# -------------------------
# 7. Manual color scale (EXACT MATCH)
# -------------------------
shadow_colors <- c(
  "Only Upstream"   = "#1b9e77",
  "Only Downstream" = "#d95f02",
  "Both Shadows"    = "#7570b3",
  "With Exitron"    = "#e7298a",
  "Without Exitron" = "#66a61e"
)

# -------------------------
# 8. Plot
# -------------------------
p <- ggplot(plot_df, aes(x = Bar, y = Percentage, fill = ShadowType)) +
  geom_col(width = 0.55) +
  
  geom_text(
    aes(label = sprintf("%.1f", Percentage)),
    position = position_stack(vjust = 0.5),
    size = 3,
    fontface = "bold",
    color = "black"
  ) +
  
  geom_text(
    data = totals_df,
    aes(x = Bar, y = y, label = scales::comma(Total)),
    inherit.aes = FALSE,
    size = 3,
    fontface = "bold"
  ) +
  
  facet_wrap(
    ~ Species,
    nrow = 1,
    strip.position = "bottom",
    labeller = labeller(Species = label_wrap_gen(width = 10))
  ) +
  
  scale_fill_manual(values = shadow_colors) +
  
  scale_y_continuous(
    limits = c(0,120),
    breaks = seq(0,100,25),
    expand = c(0,0)
  ) +
  
  labs(
    title = "Proportion of Shadow Types Across Species",
    x = NULL,
    y = "Percentage (%)",
    fill = "Shadow Type"
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
    plot.title = element_text(hjust = 0, face = "bold"),
    axis.text.x = element_text(angle = 45, hjust = 1, vjust = 1)
  )

print(p)

ggsave(OUTPUT_FILE, plot=p, width=12, height=5, dpi=300)
