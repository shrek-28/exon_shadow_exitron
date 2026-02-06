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
  INPUT_FILE  <- "C://Users//admin//Desktop//Coding//exon_shadow_exitron//shadow_types//shadow_types.tsv"
  OUTPUT_FILE <- "C://Users//admin//Desktop//Coding//exon_shadow_exitron//shadow_types//shadow_types_plot.pdf"
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
# Species naming + order
# -------------------------
df$Species <- df$Species |>
  str_to_lower() |>
  recode(
    "borangutan" = "Bornean Orangutan",
    "sorangutan" = "Sumatran Orangutan"
  ) |>
  str_to_title()

SPECIES_ORDER <- c(
  "Gibbon","Gorilla","Human","Bonobo",
  "Chimpanzee","Bornean Orangutan","Sumatran Orangutan"
)
df$Species <- factor(df$Species, levels = SPECIES_ORDER)


# -------------------------
# Build plot dataframe
# -------------------------
plot_df <- df %>%
  select(
    Species,
    pc_positive_shadows,
    pc_zero_shadows,
    pc_negative_shadows,
    pc_up_shadows,
    pc_down_shadows
  ) %>%
  pivot_longer(-Species, names_to = "Raw", values_to = "Percentage") %>%
  mutate(
    ShadowType = recode(
      Raw,
      "pc_positive_shadows" = "Positive",
      "pc_zero_shadows"     = "Zero",
      "pc_negative_shadows" = "Negative",
      "pc_up_shadows"       = "Upstream",
      "pc_down_shadows"     = "Downstream"
    ),
    Bar = case_when(
      ShadowType %in% c("Positive","Zero","Negative") ~ "Total",
      TRUE                                            ~ "Positive"
    ),
    Bar = factor(Bar, levels = c("Total","Positive"))
  )

# -------------------------
# Global stacking order
# -------------------------
plot_df$ShadowType <- factor(
  plot_df$ShadowType,
  levels = c("Positive","Zero","Negative","Upstream","Downstream")
)

# -------------------------
# Totals above bars
# -------------------------
totals_df <- df %>%
  transmute(
    Species,
    Total    = total_shadows,
    Positive = positive_shadows
  ) %>%
  pivot_longer(-Species, names_to = "Bar", values_to = "Total") %>%
  mutate(
    Bar = factor(Bar, levels = c("Total","Positive")),
    y = 110
  )

# -------------------------
# Plot
# -------------------------
p <- ggplot(plot_df, aes(x = Bar, y = Percentage, fill = ShadowType)) +
  geom_col(width = 0.55) +
  
  # Stack-centered percentage labels
  geom_text(
    aes(label = sprintf("%.1f", Percentage)),
    position = position_stack(vjust = 0.5),
    fontface = "bold",
    size = 3,
    color = "black"
  ) +
  
  # Totals above bars
  geom_text(
    data = totals_df,
    aes(x = Bar, y = y, label = scales::comma(Total)),
    inherit.aes = FALSE,
    fontface = "bold",
    size = 3
  ) +
  
  facet_wrap(
    ~ Species,
    nrow = 1,
    strip.position = "bottom",
    labeller =  labeller(Species = label_wrap_gen(width = 10))
  ) +
  
  scale_fill_manual(values = c(
    "Positive"   = "#1b9e77",
    "Zero"       = "#7570b3",
    "Negative"   = "#d95f02",
    "Upstream"   = "#e7298a",
    "Downstream" = "#66a61e"
  )) +
  
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
